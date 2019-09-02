/******************************************************************************
 * exact_parallel_minimum_cut.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <functional>
#include <memory>
#include <unordered_map>
#include <vector>

#include "algorithms/global_mincut/minimum_cut_helpers.h"
#include "algorithms/global_mincut/noi_minimum_cut.h"
#include "algorithms/global_mincut/viecut.h"
#include "coarsening/test_wrapper.h"
#include "common/configuration.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/fifo_node_bucket_pq.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "data_structure/priority_queues/node_bucket_pq.h"
#include "tools/random_functions.h"
#include "tools/timer.h"

#ifdef PARALLEL

#include "parallel/coarsening/contract_graph.h"
#include "parallel/coarsening/contraction_tests.h"
#include "parallel/coarsening/sparsify.h"
#include "parallel/data_structure/union_find.h"

#else
#include "coarsening/contract_graph.h"
#include "coarsening/contraction_tests.h"
#include "coarsening/sparsify.h"
#include "data_structure/union_find.h"
#endif

class exact_parallel_minimum_cut : public minimum_cut {
 public:
    exact_parallel_minimum_cut() { }

    ~exact_parallel_minimum_cut() { }

    static constexpr bool debug = false;
    static constexpr bool timing = true;

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G) {
        return perform_minimum_cut(G, false);
    }

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G,
                                   bool indirect) {
        if (!minimum_cut_helpers::graphValid(G))
            return -1;
        std::vector<std::shared_ptr<graph_access> > graphs;
        timer t;
        EdgeWeight mincut = G->getMinDegree();
#ifdef PARALLEL
        viecut heuristic_mc;
        mincut = heuristic_mc.perform_minimum_cut(G, true);
        LOGC(timing) << "VieCut found cut " << mincut
                     << " [Time: " << t.elapsed() << "s]";
#endif

        graphs.push_back(G);

        // if PARALLEL is set, NodeInCut are already set to the result of viecut
        // This is what we want.
#ifndef PARALLEL
        minimum_cut_helpers::setInitialCutValues(graphs);
#endif

        while (graphs.back()->number_of_nodes() > 2 && mincut > 0) {
            std::shared_ptr<graph_access> curr_g = graphs.back();
            timer ts;
#ifdef PARALLEL

            noi_minimum_cut noi;

            auto uf = parallel_modified_capforest(curr_g, mincut);
            if (uf.n() == curr_g->number_of_nodes()) {
                uf = noi.modified_capforest(curr_g, mincut);
                LOG1 << "seq capforest needed";
            }

#else
            LOG1 << "Error: Running exact_parallel_minimum_cut without PARALLEL"
                 << " Using normal noi_minimum_cut instead!";

            noi_minimum_cut noi;
            auto uf = noi.modified_capforest(curr_g, mincut);
#endif

            if (uf.n() > 1) {
                std::vector<NodeID> mapping(curr_g->number_of_nodes());
                std::vector<NodeID> part(curr_g->number_of_nodes(),
                                         UNDEFINED_NODE);
                std::vector<std::vector<NodeID> > reverse_mapping;
                NodeID current_pid = 0;
                for (NodeID n : curr_g->nodes()) {
                    NodeID part_id = uf.Find(n);

                    if (part[part_id] == UNDEFINED_NODE) {
                        part[part_id] = current_pid++;
                    }

                    mapping[n] = part[part_id];
                    curr_g->setPartitionIndex(n, part[part_id]);
                }

                graphs.push_back(
                    contraction::contractGraph(curr_g, mapping,
                                               current_pid, reverse_mapping));

                mincut = minimum_cut_helpers::updateCut(graphs, mincut);
            } else {
                break;
            }
        }

        if (!indirect && configuration::getConfig()->save_cut)
            minimum_cut_helpers::retrieveMinimumCut(graphs);

        return mincut;
    }

    std::vector<NodeID> randomStartNodes(std::shared_ptr<graph_access> G) {
        std::vector<NodeID> start_nodes;
        for (int i = 0; i < omp_get_max_threads(); ++i)
            start_nodes.push_back(
                random_functions::next() % G->number_of_nodes());

        return start_nodes;
    }

    std::vector<NodeID> bfsStartNodes(std::shared_ptr<graph_access> G) {
        NodeID starting_node = random_functions::next() % G->number_of_nodes();
        std::vector<NodeID> start_nodes;
        start_nodes.push_back(starting_node);

        for (int i = 1; i < omp_get_max_threads(); ++i) {
            std::deque<NodeID> bfs;
            std::vector<bool> nodes(G->number_of_nodes(), false);
            size_t found = i;

            for (auto el : start_nodes) {
                bfs.push_back(el);
                nodes[el] = true;
            }

            while (!bfs.empty() && found < G->number_of_nodes()) {
                NodeID no = bfs.front();
                bfs.pop_front();
                for (EdgeID e : G->edges_of(no)) {
                    NodeID tgt = G->getEdgeTarget(e);
                    if (!nodes[tgt]) {
                        found++;
                        nodes[tgt] = true;
                        bfs.push_back(tgt);
                        if (found == G->number_of_nodes()) {
                            start_nodes.push_back(tgt);
                            break;
                        }
                    }
                }
            }
        }
        return start_nodes;
    }

    union_find parallel_modified_capforest(
        std::shared_ptr<graph_access> G, const EdgeWeight mincut) {
        union_find uf(G->number_of_nodes());
        LOG << "Contract all edges with value at least " << mincut;
        timer t;

        timer timer2;
        std::vector<NodeID> start_nodes = randomStartNodes(G);

        // std::vector<bool> would be bad for thread-safety
        std::vector<uint8_t> visited(G->number_of_nodes(), false);
        std::vector<size_t> times(G->number_of_nodes(), 0);

#pragma omp parallel for
        for (int i = 0; i < omp_get_num_threads(); ++i) {
            fifo_node_bucket_pq pq(G->number_of_nodes(), mincut + 1);
            std::vector<bool> blacklisted(G->number_of_nodes(), false);
            std::vector<NodeID> r_v(G->number_of_nodes(), 0);

            NodeID starting_node = start_nodes[i];
            NodeID current_node = starting_node;

            pq.insert(current_node, 0);

            timer t;
            size_t elements = 0;

            while (!pq.empty()) {
                current_node = pq.deleteMax();

                blacklisted[current_node] = true;
                if (visited[current_node]) {
                    continue;
                } else {
                    visited[current_node] = true;
                }

                elements++;

                for (EdgeID e : G->edges_of(current_node)) {
                    NodeID tgt = G->getEdgeTarget(e);

                    if (r_v[tgt] < mincut) {
                        if ((r_v[tgt] + G->getEdgeWeight(e)) >= mincut) {
                            if (!blacklisted[tgt]) {
                                uf.Union(current_node, tgt);
                            }
                        }

                        if (!visited[tgt]) {
                            size_t new_rv =
                                std::min(r_v[tgt] + G->getEdgeWeight(e),
                                         mincut);
                            r_v[tgt] = new_rv;
                            if (!visited[tgt] && !blacklisted[tgt]) {
                                if (pq.contains(tgt)) {
                                    pq.increaseKey(tgt, new_rv);
                                } else {
                                    pq.insert(tgt, new_rv);
                                }
                            }
                        }
                    }
                }
            }
        }
        return uf;
    }
};
