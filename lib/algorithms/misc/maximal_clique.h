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

            for (EdgeID e : graph->edges_of(n)) {
                NodeID t = graph->getEdgeTarget(n, e);
                if (t > n) {
                    P.emplace_back(graph->getEdge(n, e));
                } else {
                    X.emplace_back(graph->getEdge(n, e));
                }
            }

            std::sort(P.begin(), P.end(), pairs);
            std::sort(X.begin(), X.end(), pairs);

            bronKerbosch(graph, &P, &R, &X, 0);
        }
        return cliques;
    }

    std::vector<std::pair<NodeID, EdgeWeight> > neighborhoodIntersection(
        std::vector<std::pair<NodeID, EdgeWeight> >* neighborhood,
        std::vector<std::pair<NodeID, EdgeWeight> >* current) {
        size_t n_id = 0;
        size_t p_id = 0;
        std::vector<std::pair<NodeID, EdgeWeight> > intersect;
        while (n_id < neighborhood->size() && p_id < current->size()) {
            if ((*neighborhood)[n_id].first < (*current)[p_id].first) {
                ++n_id;
                continue;
            }
            if ((*neighborhood)[n_id].first > (*current)[p_id].first) {
                ++p_id;
                continue;
            }
            NodeID vtx = (*neighborhood)[n_id].first;
            auto new_w = (*neighborhood)[n_id].second + (*current)[p_id].second;
            intersect.emplace_back(vtx, new_w);
            n_id++;
            p_id++;
        }

        return intersect;
    }

    void bronKerbosch(std::shared_ptr<mutable_graph> graph,
                      std::vector<std::pair<NodeID, EdgeWeight> >* P,
                      std::vector<NodeID>* R,
                      std::vector<std::pair<NodeID, EdgeWeight> >* X,
                      size_t depth) {
        if (P->size() == 0) {
            if (X->size() == 0) {
                cliques.emplace_back(*R);
                LOG1 << "new clique " << cliques.back();
            }
            return;
        }

        NodeID n = (*P)[0].first;
        std::vector<std::pair<NodeID, EdgeWeight> > neighborhood;
        for (EdgeID e : graph->edges_of(n)) {
            neighborhood.emplace_back(graph->getEdge(n, e));
        }
        std::sort(neighborhood.begin(), neighborhood.end(), pairs);

        std::vector<NodeID> non_neighbors;

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

            auto nextP = neighborhoodIntersection(&neighborhood, P);
            auto nextX = neighborhoodIntersection(&neighborhood, X);

            std::vector<NodeID> nextR = *R;
            nextR.emplace_back(n);

            bronKerbosch(graph, &nextP, &nextR, &nextX, depth + 1);
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
    }

    std::vector<std::vector<NodeID> > cliques;
};
