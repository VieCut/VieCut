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

#include "algorithms/global_mincut/noi_minimum_cut.h"
#include "algorithms/global_mincut/viecut.h"
#include "coarsening/test_wrapper.h"
#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/fifo_node_bucket_pq.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "data_structure/priority_queues/node_bucket_pq.h"
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

#include "definitions.h"
#include "tools/random_functions.h"

#include "algorithms/global_mincut/minimum_cut_helpers.h"
#include "data_structure/priority_queues/bucket_pq.h"
#include "data_structure/priority_queues/node_bucket_pq.h"
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <memory>
#include <unordered_map>

class exact_parallel_minimum_cut : public minimum_cut
{
public:
    exact_parallel_minimum_cut() { }

    ~exact_parallel_minimum_cut() { }

    static constexpr bool debug = false;
    static constexpr bool timing = true;

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G, bool save_cut) {
        return perform_minimum_cut(G, save_cut, false);
    }

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G,
                                   std::string pq,
                                   bool save_cut,
                                   bool disable_limiting) {
        return perform_minimum_cut(G, save_cut, false, pq, disable_limiting);
    }

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G,
                                   bool save_cut,
                                   bool __attribute__ ((unused)) indirect,
                                   std::string __attribute__ ((unused)) pq = "default",
                                   bool __attribute__ ((unused)) disable_limiting = false) {

        if (!minimum_cut_helpers::graphValid(G))
            return -1;
        std::vector<std::shared_ptr<graph_access> > graphs;
        timer t;
        EdgeWeight global_mincut = G->getMinDegree();
#ifdef PARALLEL
        viecut heuristic_mc;
        global_mincut = heuristic_mc.perform_minimum_cut(G, true);
        LOGC(timing) << "VieCut found cut " << global_mincut << " [Time: " << t.elapsed() << "s]";
#endif

        graphs.push_back(G);

        // if PARALLEL is set, NodeInCut are already set to the result of viecut. This is what we want.
#ifndef PARALLEL
        minimum_cut_helpers::setInitialCutValues(graphs, save_cut);
#endif

        while (graphs.back()->number_of_nodes() > 2 && global_mincut > 0) {
            std::shared_ptr<graph_access> curr_g = graphs.back();
            timer ts;
            union_find uf(curr_g->number_of_nodes());
#ifdef PARALLEL

            noi_minimum_cut noi;

            if (graphs.size() == 1 || curr_g->number_of_nodes() < graphs[graphs.size() - 2]->number_of_nodes()) {
                global_mincut = std::min(global_mincut,
                                         parallel_modified_capforest(curr_g, global_mincut,
                                                                     uf, graphs, save_cut, pq, !disable_limiting));
                if (uf.n() == curr_g->number_of_nodes()) {
                    global_mincut = std::min(global_mincut,
                                             noi.modified_capforest(curr_g, global_mincut,
                                                                    uf, graphs, save_cut, pq));
                    LOG1 << "seq capforest needed";
                }
            }
            else {
                global_mincut = std::min(global_mincut,
                                         noi.modified_capforest(curr_g, global_mincut,
                                                                uf, graphs, save_cut, pq));
            }

#else
            std::cerr << "Error: Running exact_parallel_minimum_cut without PARALLEL defined."
                " Using normal noi_minimum_cut instead!" << std::endl;

            noi_minimum_cut noi;
            global_mincut = std::min(global_mincut,
                                     noi.modified_capforest(curr_g, global_mincut,
                                                            uf, graphs, save_cut, pq));
#endif

            if (uf.n() > 1) {
                std::vector<NodeID> mapping(curr_g->number_of_nodes());
                std::vector<NodeID> part(curr_g->number_of_nodes(), UNDEFINED_NODE);
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

                graphs.push_back(contraction::contractGraph(curr_g, mapping, current_pid, reverse_mapping));

                global_mincut = minimum_cut_helpers::updateCutValueAfterContraction(graphs, global_mincut, save_cut);
            }
            else {
                break;
            }
        }

        if (!indirect && save_cut)
            minimum_cut_helpers::retrieveMinimumCut(graphs);

        return global_mincut;
    }

    std::vector<NodeID> randomStartNodes(std::shared_ptr<graph_access> G) {
        std::vector<NodeID> start_nodes;
        for (int i = 0; i < omp_get_max_threads(); ++i)
            start_nodes.push_back(random_functions::next() % G->number_of_nodes());

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

    EdgeWeight
    parallel_modified_capforest(
        std::shared_ptr<graph_access> G, EdgeWeight mincut,
        union_find& uf,
        std::vector<std::shared_ptr<graph_access> >& all_graphs,
        bool save_cut,
        std::string pq_str = "default",
        bool limit = true) {

        size_t alphamin = mincut;
        all_graphs.size();
        timer timer2;
        std::vector<NodeID> start_nodes = randomStartNodes(G);

        std::vector<size_t> times(G->number_of_nodes(), 0);
#pragma omp parallel for
        for (int i = 0; i < omp_get_max_threads(); ++i) {

            priority_queue_interface* pq = noi_minimum_cut::selectPq(pq_str, mincut, G, limit);
            std::vector<bool> visited(G->number_of_nodes(), false);
            std::vector<bool> seen(G->number_of_nodes(), false);
            std::vector<bool> blacklisted(G->number_of_nodes(), false);
            std::vector<EdgeWeight> r_v(G->number_of_nodes(), 0);

            NodeID starting_node = start_nodes[i];
            NodeID current_node = starting_node;

            pq->insert(current_node, 0);
            EdgeWeight alpha = 0;

            timer t;
            size_t elements = 0;

            while (!pq->empty()) {
                current_node = pq->deleteMax();

                visited[current_node] = true;
                if (times[current_node]++ >= 1) {
                    blacklisted[current_node] = true;
                    continue;
                }

                elements++;

                alpha += G->getWeightedNodeDegree(current_node) - (2 * r_v[current_node]);
                if (alpha < alphamin && (alpha > 0 || elements < G->number_of_nodes())) {

                    alphamin = alpha;
                    if (save_cut) {
                        for (NodeID idx : all_graphs[0]->nodes()) {
                            NodeID coarseID = idx;
                            for (size_t lv = 0; lv < all_graphs.size() - 1; ++lv) {
                                coarseID = all_graphs[lv]->getPartitionIndex(coarseID);
                            }

                            if (!visited[coarseID]) {
                                all_graphs[0]->setNodeInCut(idx, false);
                            }
                            else {
                                all_graphs[0]->setNodeInCut(idx, true);
                            }
                        }
                    }
                }

                if (limit) {
                    processVertexLimited(G, current_node, visited, r_v, blacklisted, seen, mincut, uf, pq);
                }
                else {
                    processVertexUnlimited(G, current_node, visited, r_v, blacklisted, seen, mincut, uf, pq);
                }
            }
            delete pq;
        }

        return alphamin;
    }

    void processVertexLimited(std::shared_ptr<graph_access> G,
                              NodeID current_node,
                              std::vector<bool>& visited,
                              std::vector<EdgeWeight>& r_v,
                              std::vector<bool>& blacklisted,
                              std::vector<bool>& seen,
                              EdgeWeight mincut,
                              union_find& uf,
                              priority_queue_interface* pq) {
        for (EdgeID e : G->edges_of(current_node)) {
            NodeID tgt = G->getEdgeTarget(e);

            if (!visited[tgt]) {

                bool increase = false;

                if (r_v[tgt] < mincut) {
                    increase = true;
                    if ((r_v[tgt] + G->getEdgeWeight(e)) >= mincut) {
                        if (!blacklisted[tgt])
                            uf.Union(current_node, tgt);
                    }
                }

                r_v[tgt] += G->getEdgeWeight(e);

                size_t new_rv = std::min(r_v[tgt], mincut);

                if (seen[tgt]) {
                    if (increase && !visited[tgt] && !blacklisted[tgt]) {
                        pq->increaseKey(tgt, new_rv);
                    }
                }
                else {
                    seen[tgt] = true;
                    if (!blacklisted[tgt]) {
                        pq->insert(tgt, new_rv);
                    }
                }
            }
        }
    }

    void processVertexUnlimited(std::shared_ptr<graph_access> G,
                                NodeID current_node,
                                std::vector<bool>& visited,
                                std::vector<EdgeWeight>& r_v,
                                std::vector<bool>& blacklisted,
                                std::vector<bool>& seen,
                                EdgeWeight mincut,
                                union_find& uf,
                                priority_queue_interface* pq) {

        for (EdgeID e : G->edges_of(current_node)) {
            NodeID tgt = G->getEdgeTarget(e);

            if (!visited[tgt]) {
                if (r_v[tgt] < mincut) {
                    if ((r_v[tgt] + G->getEdgeWeight(e)) >= mincut) {
                        if (!blacklisted[tgt])
                            uf.Union(current_node, tgt);
                    }
                }

                r_v[tgt] += G->getEdgeWeight(e);

                if (seen[tgt]) {
                    if (!visited[tgt] && !blacklisted[tgt]) {
                        pq->increaseKey(tgt, r_v[tgt]);
                    }
                }
                else {
                    seen[tgt] = true;
                    if (!blacklisted[tgt]) {
                        pq->insert(tgt, r_v[tgt]);
                    }
                }
            }
        }
    }
};
