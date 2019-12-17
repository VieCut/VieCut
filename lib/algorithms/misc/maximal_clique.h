/******************************************************************************
 * maximal_clique.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/
#pragma once

#include <algorithm>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "algorithms/misc/equal_neighborhood.h"
#include "algorithms/multicut/multicut_problem.h"
#include "data_structure/mutable_graph.h"

class maximal_clique {
 public:
    maximal_clique() { }

    std::vector<std::vector<NodeID> > findCliques(
        std::shared_ptr<mutable_graph> graph) {
        for (NodeID n : graph->nodes()) {
            std::vector<std::pair<NodeID, EdgeWeight> > P;
            std::vector<NodeID> R = { n };
            std::vector<std::pair<NodeID, EdgeWeight> > X;

            size_t external = 0;
            for (EdgeID e : graph->edges_of(n)) {
                auto [t, w] = graph->getEdge(n, e);
                if (t > n) {
                    P.emplace_back(graph->getEdge(n, e));
                } else {
                    X.emplace_back(graph->getEdge(n, e));
                    external += w;
                }
            }

            std::sort(P.begin(), P.end(), pairs);
            std::sort(X.begin(), X.end(), pairs);

            bronKerbosch(graph, &P, &R, &X, 0, external);
        }
        return cliques;
    }

    std::tuple<std::vector<std::pair<NodeID, EdgeWeight> >,
               std::vector<std::pair<NodeID, EdgeWeight> >, size_t>
    neighborhoodIntersection(std::vector<std::pair<NodeID, EdgeWeight> >* N,
                             std::vector<std::pair<NodeID, EdgeWeight> >* P,
                             std::vector<NodeID>* R,
                             std::vector<std::pair<NodeID, EdgeWeight> >* X,
                             size_t external_weight) {
        size_t n_id = 0;
        size_t p_id = 0;
        size_t r_id = 0;
        size_t x_id = 0;
        std::vector<std::pair<NodeID, EdgeWeight> > intersect_p;
        std::vector<std::pair<NodeID, EdgeWeight> > intersect_x;
        while (n_id < N->size()) {
            NodeID n = (*N)[n_id].first;
            NodeID p;
            if (p_id == P->size()) {
                p = UNDEFINED_NODE;
            } else {
                p = (*P)[p_id].first;
            }
            NodeID x;
            if (x_id == X->size()) {
                x = UNDEFINED_NODE;
            } else {
                x = (*X)[x_id].first;
            }

            if (n <= p && n <= x) {
                ++n_id;
            }

            if (n < p && n < x) {
                while ((*R)[r_id] < n) {
                    ++r_id;
                }

                if ((*R)[r_id] != n) {
                    external_weight += (*N)[n_id - 1].second;
                }
                continue;
            }

            if (n == x) {
                auto new_w = (*N)[n_id].second + (*X)[x_id].second;
                intersect_x.emplace_back(n, new_w);
                ++x_id;
            }

            if (n > x) {
                ++x_id;
            }

            if (n == p) {
                auto new_w = (*N)[n_id].second + (*P)[p_id].second;
                intersect_p.emplace_back(n, new_w);
                ++p_id;
            }

            if (n > p) {
                ++p_id;
            }
        }

        return std::make_tuple(intersect_p, intersect_x, external_weight);
    }

    bool bronKerbosch(std::shared_ptr<mutable_graph> graph,
                      std::vector<std::pair<NodeID, EdgeWeight> >* P,
                      std::vector<NodeID>* R,
                      std::vector<std::pair<NodeID, EdgeWeight> >* X,
                      size_t depth,
                      size_t ex_weight) {
        if (ex_weight > (P->size() + R->size() - 1)) {
            return false;
        }

        if (P->size() == 0) {
            if (X->size() == 0) {
                cliques.emplace_back(*R);
                LOG1 << "new clique " << cliques.back()
                     << " with ex-weight " << ex_weight;
            }
            return true;
        }

        NodeID n = (*P)[0].first;
        std::vector<std::pair<NodeID, EdgeWeight> > neighborhood;
        for (EdgeID e : graph->edges_of(n)) {
            neighborhood.emplace_back(graph->getEdge(n, e));
        }
        std::sort(neighborhood.begin(), neighborhood.end(), pairs);
        size_t n_id = 0;
        size_t p_id = 0;
        std::vector<std::pair<NodeID, EdgeWeight> > intersect;
        while (p_id < P->size()) {
            if (n_id == neighborhood.size() ||
                neighborhood[n_id].first > (*P)[p_id].first) {
                intersect.emplace_back((*P)[p_id]);
                ++p_id;
                continue;
            }

            if (neighborhood[n_id].first < (*P)[p_id].first) {
                ++n_id;
                continue;
            }

            n_id++;
            p_id++;
        }

        for (size_t i = 0; i < intersect.size(); ++i) {
            n = intersect[i].first;

            neighborhood.clear();
            for (EdgeID e : graph->edges_of(n)) {
                neighborhood.emplace_back(graph->getEdge(n, e));
            }
            std::sort(neighborhood.begin(), neighborhood.end(), pairs);

            auto [nextP, nextX, current_ex_weight] =
                neighborhoodIntersection(&neighborhood, P, R, X, ex_weight);

            if (current_ex_weight >= 2 * (P->size() + R->size() - 1)) {
                return false;
            }

            std::vector<NodeID> nextR = *R;
            nextR.emplace_back(n);

            if (!bronKerbosch(graph, &nextP, &nextR,
                              &nextX, depth + 1, current_ex_weight)) {
                return false;
            }

            // remove n from P
            auto f = std::find_if(P->begin(), P->end(),
                                  [=](const auto el1) {
                                      return (el1.first == n);
                                  });

            if (f != P->end()) {
                P->erase(f);
            }

            auto f2 = std::lower_bound(X->begin(), X->end(), n,
                                       [](const auto& p1, const auto& p2) {
                                           return p1.first < p2;
                                       });

            if (f2 == X->end() || (*f2).first != n) {
                X->insert(f2, std::make_pair(n, 0));
            }
        }
        return true;
    }

    std::vector<std::vector<NodeID> > cliques;
};
