/******************************************************************************
 * noi_minimum_cut.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2017-2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "algorithms/global_mincut/minimum_cut.h"
#include "algorithms/global_mincut/minimum_cut_helpers.h"
#include "algorithms/multicut/multicut_problem.h"
#include "common/configuration.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "data_structure/priority_queues/bucket_pq.h"
#include "data_structure/priority_queues/fifo_node_bucket_pq.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "data_structure/priority_queues/node_bucket_pq.h"
#include "data_structure/priority_queues/vecMaxNodeHeap.h"
#include "tools/random_functions.h"
#include "tools/timer.h"

#ifdef PARALLEL
#include "parallel/coarsening/contract_graph.h"
#include "parallel/coarsening/contraction_tests.h"
#include "parallel/data_structure/union_find.h"
#else
#include "coarsening/contract_graph.h"
#include "coarsening/contraction_tests.h"
#include "data_structure/union_find.h"
#endif

class noi_minimum_cut : public minimum_cut {
 public:
    noi_minimum_cut() { }
    ~noi_minimum_cut() { }
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
        graphs.push_back(G);
        minimum_cut_helpers::setInitialCutValues(graphs);

        while (graphs.back()->number_of_nodes() > 2 && mincut > 0) {
            auto uf = modified_capforest(graphs.back(), mincut);
            graphs.emplace_back(contraction::fromUnionFind(graphs.back(), &uf));
            mincut = minimum_cut_helpers::updateCut(graphs, mincut);
        }

        if (!indirect && configuration::getConfig()->save_cut)
            minimum_cut_helpers::retrieveMinimumCut(graphs);

        return mincut;
    }

    static priority_queue_interface * selectPq(std::shared_ptr<graph_access> G,
                                               EdgeWeight mincut) {
        priority_queue_interface* pq;

        size_t max_deg = G->getMaxDegree();

        if (configuration::getConfig()->disable_limiting) {
            mincut = max_deg;
        }

        if (configuration::getConfig()->pq == "default") {
            if (mincut > 10000 && mincut > G->number_of_nodes()) {
                pq = new vecMaxNodeHeap(G->number_of_nodes());
            } else {
                pq = new fifo_node_bucket_pq(G->number_of_nodes(), mincut);
            }
        } else {
            if (configuration::getConfig()->pq == "bqueue") {
                pq = new fifo_node_bucket_pq(G->number_of_nodes(), mincut);
            } else {
                if (configuration::getConfig()->pq == "heap") {
                    pq = new vecMaxNodeHeap(G->number_of_nodes());
                } else {
                    if (configuration::getConfig()->pq == "bstack") {
                        pq = new node_bucket_pq(G->number_of_nodes(), mincut);
                    } else {
                        std::cerr << "unknown pq type "
                                  << configuration::getConfig()->pq
                                  << std::endl;
                        exit(1);
                        pq = new maxNodeHeap();
                    }
                }
            }
        }
        return pq;
    }

    union_find modified_capforest(std::shared_ptr<graph_access> G,
                                  EdgeWeight mincut) {
        union_find uf(G->number_of_nodes());

        priority_queue_interface* pq = selectPq(G, mincut);

        std::vector<bool> visited(G->number_of_nodes(), false);
        std::vector<bool> seen(G->number_of_nodes(), false);
        std::vector<EdgeWeight> r_v(G->number_of_nodes(), 0);

        NodeID starting_node = random_functions::next() % G->number_of_nodes();

        NodeID current_node = starting_node;

        pq->insert(current_node, 0);

        while (!pq->empty()) {
            current_node = pq->deleteMax();
            visited[current_node] = true;
            if (!configuration::getConfig()->disable_limiting) {
                for (EdgeID e : G->edges_of(current_node)) {
                    NodeID tgt = G->getEdgeTarget(e);
                    if (!visited[tgt]) {
                        bool increase = false;

                        if (r_v[tgt] < mincut) {
                            increase = true;
                            if ((r_v[tgt] + G->getEdgeWeight(e)) >= mincut) {
                                uf.Union(current_node, tgt);
                            }
                        }

                        r_v[tgt] += G->getEdgeWeight(e);

                        size_t new_rv = std::min(r_v[tgt], mincut);

                        if (seen[tgt]) {
                            if (increase && !visited[tgt]) {
                                pq->increaseKey(tgt, new_rv);
                            }
                        } else {
                            seen[tgt] = true;
                            pq->insert(tgt, new_rv);
                        }
                    }
                }
            } else {
                for (EdgeID e : G->edges_of(current_node)) {
                    NodeID tgt = G->getEdgeTarget(e);

                    if (!visited[tgt]) {
                        if (r_v[tgt] < mincut) {
                            if ((r_v[tgt] + G->getEdgeWeight(e)) >= mincut) {
                                uf.Union(current_node, tgt);
                            }
                        }

                        r_v[tgt] += G->getEdgeWeight(e);

                        if (seen[tgt]) {
                            if (!visited[tgt]) {
                                pq->increaseKey(tgt, r_v[tgt]);
                            }
                        } else {
                            seen[tgt] = true;
                            pq->insert(tgt, r_v[tgt]);
                        }
                    }
                }
            }
        }
        delete pq;
        return uf;
    }

    // Terrible code duplication, but differences are too big for reuse
    union_find modified_capforest(problemPointer mcp,
                                  EdgeWeight mincut) {
        mutableGraphPtr G = mcp->graph;
        std::unordered_set<NodeID> terminal_set;

        for (const auto& t : mcp->terminals) {
            terminal_set.insert(t.position);
        }

        union_find uf(G->n());

        vecMaxNodeHeap pq(G->n());
        std::vector<bool> visited(G->number_of_nodes(), false);
        std::vector<bool> seen(G->number_of_nodes(), false);
        std::vector<EdgeWeight> r_v(G->number_of_nodes(), 0);

        NodeID starting_node = random_functions::next() % G->number_of_nodes();
        NodeID current_node = starting_node;
        pq.insert(current_node, 0);
        timer t;

        while (!pq.empty()) {
            current_node = pq.deleteMax();
            visited[current_node] = true;
            for (EdgeID e : G->edges_of(current_node)) {
                NodeID tgt = G->getEdgeTarget(current_node, e);
                if (!visited[tgt]) {
                    bool increase = false;
                    if (r_v[tgt] < mincut) {
                        increase = true;
                        if ((r_v[tgt] + G->getEdgeWeight(current_node, e)) >=
                            mincut) {
                            if (terminal_set.count(tgt) == 0
                                && terminal_set.count(current_node) == 0) {
                                uf.Union(current_node, tgt);
                            }
                        }
                    }

                    r_v[tgt] += G->getEdgeWeight(current_node, e);

                    size_t new_rv = std::min(r_v[tgt], mincut);

                    if (seen[tgt]) {
                        if (increase && !visited[tgt]) {
                            pq.increaseKey(tgt, new_rv);
                        }
                    } else {
                        seen[tgt] = true;
                        pq.insert(tgt, new_rv);
                    }
                }
            }
        }
        return uf;
    }
};
