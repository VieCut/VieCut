/******************************************************************************
 * contraction_tests.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <atomic>
#include <memory>
#include <utility>
#include <vector>

#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "tlx/logger.hpp"

#ifdef PARALLEL
#include "parallel/data_structure/union_find.h"
#else
#include "data_structure/union_find.h"
#endif

class tests {
 private:
    static void sort3(EdgeWeight w1, EdgeWeight w2, EdgeWeight w3) {
        if (w1 < w2) {
            if (w2 < w3) {
                return;
            } else if (w1 < w3) {
                std::swap(w2, w3);
            } else {
                EdgeWeight tmp = std::move(w1);
                w1 = std::move(w3);
                w3 = std::move(w2);
                w2 = std::move(tmp);
            }
        } else {
            if (w1 < w3) {
                std::swap(w1, w2);
            } else if (w3 < w2) {
                std::swap(w1, w3);
            } else {
                EdgeWeight tmp = std::move(w1);
                w1 = std::move(w2);
                w2 = std::move(w3);
                w3 = std::move(tmp);
            }
        }
    }

 public:
    static union_find prTests12(std::shared_ptr<graph_access> G,
                                EdgeWeight weight_limit,
                                bool find_all_cuts = false) {
        union_find uf(G->number_of_nodes());
        G->computeDegrees();

        // workaround for std::vector<bool> not being usable in parallel
        std::vector<uint8_t> contracted(G->number_of_nodes(), false);

        NodeID end = G->number_of_nodes();
#pragma omp parallel for schedule(dynamic, 100)
        for (NodeID n = 0; n < end; ++n) {
            NodeWeight n_wgt = G->getWeightedNodeDegree(n);

            for (EdgeID e : G->edges_of(n)) {
                EdgeWeight wgt = G->getEdgeWeight(e);
                NodeID tgt = G->getEdgeTarget(e);
                NodeWeight tgt_wgt = G->getWeightedNodeDegree(tgt);
                if (wgt >= weight_limit) {
                    contracted[n] = true;
                    contracted[tgt] = true;
                    uf.Union(n, tgt);
                }

                // if we want to find all cuts
                // we are not allowed to contract an edge
                // when an incident vertex has degree mincut
                // (as the singleton cut might be important)
                if (!find_all_cuts ||
                    (n_wgt >= weight_limit && tgt_wgt >= weight_limit)) {
                    // node degrees change when we contract edges.
                    // thus, we only use PR 2 or 3
                    // when the incident vertices haven't been contracted yet
                    // keeping a data structure with current degrees
                    // would be too expensive in parallel
                    if ((2 * wgt) > n_wgt) {
                        if (__sync_bool_compare_and_swap(&contracted[n],
                                                         false, true)) {
                            contracted[tgt] = true;
                            uf.Union(n, tgt);
                        }
                    } else if ((2 * wgt) > tgt_wgt) {
                        if (__sync_bool_compare_and_swap(&contracted[tgt],
                                                         false, true)) {
                            contracted[n] = true;
                            uf.Union(n, tgt);
                        }
                    }
                }
            }
        }
        return uf;
    }

    static union_find prTests34(std::shared_ptr<graph_access> G,
                                EdgeWeight weight_limit,
                                bool find_all_cuts = false) {
        union_find uf(G->number_of_nodes());
        G->computeDegrees();
        std::vector<uint8_t> finished(G->number_of_nodes(), false);
        std::vector<uint8_t> contracted(G->number_of_nodes(), 0);
#pragma omp parallel
        {
            std::vector<EdgeID> marked(G->number_of_nodes(), UNDEFINED_EDGE);
#pragma omp for schedule(dynamic, 100)
            for (NodeID n = 0; n < G->number_of_nodes(); ++n) {
                if (finished[n])
                    continue;

                finished[n] = true;

                for (EdgeID e : G->edges_of(n)) {
                    NodeID tgt = G->getEdgeTarget(e);
                    if (tgt > n) {
                        marked[tgt] = e;
                    }
                }

                NodeID deg_n = G->getWeightedNodeDegree(n);
                for (EdgeID e1 : G->edges_of(n)) {
                    NodeID tgt = G->getEdgeTarget(e1);
                    NodeID deg_tgt = G->getWeightedNodeDegree(tgt);
                    EdgeWeight w1 = G->getEdgeWeight(e1);
                    if (finished[tgt]) {
                        marked[tgt] = UNDEFINED_EDGE;
                        continue;
                    }

                    finished[tgt] = true;
                    EdgeWeight wgt_sum = G->getEdgeWeight(e1);
                    if (tgt > n) {
                        for (EdgeID e2 : G->edges_of(tgt)) {
                            NodeID tgt2 = G->getEdgeTarget(e2);
                            if (marked[tgt2] == UNDEFINED_EDGE)
                                continue;

                            if (marked[tgt2] < G->get_first_edge(n)
                                || marked[tgt2] >= G->get_first_invalid_edge(n))
                                continue;

                            EdgeWeight w2 = G->getEdgeWeight(e2);
                            EdgeWeight w3 = G->getEdgeWeight(marked[tgt2]);

                            wgt_sum += std::min(w2, w3);

                            bool contractible_one_cut =
                                !find_all_cuts
                                && 2 * (w1 + w3) >= deg_n
                                && 2 * (w1 + w2) >= deg_tgt;

                            bool contractible_all_cuts =
                                find_all_cuts
                                && 2 * (w1 + w3) > deg_n
                                && 2 * (w1 + w2) > deg_tgt
                                && deg_n >= weight_limit
                                && deg_tgt >= weight_limit;

                            if (contractible_one_cut || contractible_all_cuts) {
                                // node degrees change when we contract edges.
                                // thus, we only use PR 2 or 3 when the
                                // incident vertices haven't been contracted yet
                                // keeping a data structure with current
                                // degrees would be too expensive in parallel
                                if (__sync_bool_compare_and_swap(
                                        &contracted[n], false, true)) {
                                    if (__sync_bool_compare_and_swap(
                                            &contracted[tgt], false, true)) {
                                        uf.Union(n, tgt);
                                        break;
                                    }
                                }
                            }
                        }

                        if (wgt_sum >= weight_limit) {
                            contracted[n] = true;
                            contracted[tgt] = true;
                            uf.Union(n, tgt);
                        }
                        marked[tgt] = UNDEFINED_EDGE;
                    }
                }
            }
        }
        return uf;
    }

    static void findHeavyEdges(std::shared_ptr<graph_access> G,
                               union_find* uf,
                               EdgeWeight weight_limit) {
#pragma omp parallel for schedule(guided)
        for (NodeID n = 0; n < G->number_of_nodes(); ++n) {
            for (EdgeID e : G->edges_of(n)) {
                if (G->getEdgeWeight(e) > weight_limit) {
                    if (!uf->SameSet(n, G->getEdgeTarget(e)))
                        uf->Union(n, G->getEdgeTarget(e));
                }
            }
        }
    }

    static void findHeavyTriangles(std::shared_ptr<graph_access> G,
                                   union_find* uf,
                                   EdgeWeight weight_limit) {
#pragma omp parallel
        {
            std::vector<bool> marked(G->number_of_nodes(), false);

#pragma omp for schedule(guided)
            for (NodeID n = 0; n < G->number_of_nodes(); ++n) {
                for (EdgeID e : G->edges_of(n)) {
                    NodeID tgt = G->getEdgeTarget(e);
                    if (tgt > n) {
                        marked[tgt] = true;
                    }
                }

                for (EdgeID e1 : G->edges_of(n)) {
                    NodeID tgt = G->getEdgeTarget(e1);
                    if (tgt > n) {
                        for (EdgeID e2 : G->edges_of(tgt)) {
                            NodeID tgt2 = G->getEdgeTarget(e2);
                            if (marked[tgt2]) {
                                for (EdgeID e3 : G->edges_of(n)) {
                                    if (G->getEdgeTarget(e3) == tgt2) {
                                        EdgeWeight w1 = G->getEdgeWeight(e1);
                                        EdgeWeight w2 = G->getEdgeWeight(e2);
                                        EdgeWeight w3 = G->getEdgeWeight(e3);

                                        sort3(w1, w2, w3);

                                        if (w1 + w2 > weight_limit) {
#pragma omp critical
                                            {
                                                uf->Union(n, tgt);
                                                uf->Union(n, tgt2);
                                            }
                                        }
                                        break;
                                    }
                                }
                            }
                        }
                        marked[tgt] = false;
                    }
                }
            }
        }
    }
};
