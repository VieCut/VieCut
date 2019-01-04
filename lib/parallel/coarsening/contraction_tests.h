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

#include "data_structure/graph_access.h"
#include "definitions.h"
#include "tlx/logger.hpp"

#ifdef PARALLEL
#include "parallel/data_structure/union_find.h"
#else
#include "data_structure/union_find.h"
#endif

#include <cstdint>
#include <cstdlib>

class tests
{

private:
    static void sort3(EdgeWeight& w1, EdgeWeight& w2, EdgeWeight& w3) {
        if (w1 < w2) {
            if (w2 < w3) {
                return;
            }
            else if (w1 < w3) {
                std::swap(w2, w3);
            }
            else {
                EdgeWeight tmp = std::move(w1);
                w1 = std::move(w3);
                w3 = std::move(w2);
                w2 = std::move(tmp);
            }
        }
        else {
            if (w1 < w3) {
                std::swap(w1, w2);
            }
            else if (w3 < w2) {
                std::swap(w1, w3);
            }
            else {
                EdgeWeight tmp = std::move(w1);
                w1 = std::move(w2);
                w2 = std::move(w3);
                w3 = std::move(tmp);
            }
        }
    }

public:
    static union_find prTests12(std::shared_ptr<graph_access> G, EdgeWeight weight_limit) {

        union_find uf(G->number_of_nodes());

        NodeID end = G->number_of_nodes();
#pragma omp parallel for schedule(dynamic, 100)
        for (NodeID n = 0; n < end; ++n) {
            for (EdgeID e : G->edges_of(n)) {
                EdgeWeight wgt = G->getEdgeWeight(e);
                if (wgt >= weight_limit) {
                    uf.Union(n, G->getEdgeTarget(e));
                }
                // PR Test 2 cannot be performed in parallel without too much overhead and thus will not be performed in shared-memory parallel.
            }
        }

        return uf;
    }

    static union_find prTests34(std::shared_ptr<graph_access> G, EdgeWeight weight_limit) {

        union_find uf(G->number_of_nodes());

        std::vector<bool> donezo(G->number_of_nodes(), false);

#pragma omp parallel
        {
            std::vector<EdgeID> marked(G->number_of_nodes(), false);
#pragma omp for schedule(dynamic, 100)
            for (NodeID n = 0; n < G->number_of_nodes(); ++n) {

                if (donezo[n])
                    continue;

                donezo[n] = true;

                for (EdgeID e : G->edges_of(n)) {
                    NodeID tgt = G->getEdgeTarget(e);
                    if (tgt > n) {
                        marked[tgt] = e;
                    }
                }

                EdgeWeight deg_n = G->getWeightedNodeDegree(n);
                for (EdgeID e1 : G->edges_of(n)) {
                    NodeID tgt = G->getEdgeTarget(e1);
                    EdgeWeight deg_tgt = G->getWeightedNodeDegree(tgt);
                    if (donezo[tgt])
                        continue;

                    donezo[tgt] = true;
                    EdgeWeight wgt_sum = G->getEdgeWeight(e1);
                    if (tgt > n) {
                        for (EdgeID e2 : G->edges_of(tgt)) {
                            NodeID tgt2 = G->getEdgeTarget(e2);
                            if (!marked[tgt2])
                                continue;

                            EdgeWeight w1 = G->getEdgeWeight(e1);
                            EdgeWeight w2 = G->getEdgeWeight(e2);
                            EdgeWeight w3 = G->getEdgeWeight(marked[tgt2]);

                            wgt_sum += std::min(w2, w3);

                            if (2 * (w1 + w3) > deg_n &&
                                2 * (w1 + w2) > deg_tgt) {
                                if (!uf.SameSet(n, tgt))
                                    uf.Union(n, tgt);
                            }
                        }

                        if (wgt_sum >= weight_limit) {
                            if (!uf.SameSet(n, tgt))
                                uf.Union(n, tgt);
                        }
                        marked[tgt] = 0;
                    }
                }
            }
        }
        return uf;
    }

    static void findHeavyEdges(graph_access& G,
                               union_find& uf,
                               EdgeWeight weight_limit) {

#pragma omp parallel for schedule(guided)
        for (NodeID n = 0; n < G.number_of_nodes(); ++n) {
            for (EdgeID e : G.edges_of(n)) {
                if (G.getEdgeWeight(e) > weight_limit) {
                    if (!uf.SameSet(n, G.getEdgeTarget(e)))
                        uf.Union(n, G.getEdgeTarget(e));
                }
            }
        }
    }

    static void findHeavyTriangles(graph_access& G,
                                   union_find& uf,
                                   EdgeWeight weight_limit) {

#pragma omp parallel
        {
            std::vector<bool> marked(G.number_of_nodes(), false);

#pragma omp for schedule(guided)
            for (NodeID n = 0; n < G.number_of_nodes(); ++n) {

                for (EdgeID e : G.edges_of(n)) {
                    NodeID tgt = G.getEdgeTarget(e);
                    if (tgt > n) {
                        marked[tgt] = true;
                    }
                }

                for (EdgeID e1 : G.edges_of(n)) {
                    NodeID tgt = G.getEdgeTarget(e1);
                    if (tgt > n) {
                        for (EdgeID e2 : G.edges_of(tgt)) {
                            NodeID tgt2 = G.getEdgeTarget(e2);
                            if (marked[tgt2]) {
                                for (EdgeID e3 : G.edges_of(n)) {
                                    if (G.getEdgeTarget(e3) == tgt2) {
                                        EdgeWeight w1 = G.getEdgeWeight(e1);
                                        EdgeWeight w2 = G.getEdgeWeight(e2);
                                        EdgeWeight w3 = G.getEdgeWeight(e3);

                                        sort3(w1, w2, w3);

                                        if (w1 + w2 > weight_limit) {
#pragma omp critical
                                            {
                                                uf.Union(n, tgt);
                                                uf.Union(n, tgt2);
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
