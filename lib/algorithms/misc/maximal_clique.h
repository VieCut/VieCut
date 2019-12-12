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

#include <memory>

#include "algorithms/multicut/multicut_problem.h"
#include "data_structure/mutable_graph.h"



class maximal_clique {
 public:
    maximal_clique() { };

    void findCliques(std::shared_ptr<multicut_problem> problem) {
        for (NodeID n : problem->graph->nodes()) {
            std::vector<std::pair<NodeID, EdgeWeight> > P;
            std::vector<NodeID> R = {n};
            std::vector<std::pair<NodeID, EdgeWeight> > X;

            for (EdgeID e : problem->graph->edges_of(n)) {
                P.emplace_back(problem->graph->getEdge(n, e));
            }

            std::sort(P.begin(), P.end(), pairs);

            X = P;

           // LOG1 << "starting" << P << R << X;



            bronKerbosch(problem->graph, &P, &R, &X);
        }
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
                      std::vector<std::pair<NodeID, EdgeWeight> >* X) {
        if (P->size() == 0 && X->size() == 0) {
            LOG1 << "clique" << *R;
        }

        std::vector<std::pair<NodeID, EdgeWeight> > neighborhood;
        while (P->size() > 0) {
            NodeID n = (*P)[0].first;
            for (EdgeID e : graph->edges_of(n)) {
                neighborhood.emplace_back(graph->getEdge(n, e));
            }

            std::sort(neighborhood.begin(), neighborhood.end(), pairs);
            auto nextP = neighborhoodIntersection(&neighborhood, P);
            auto nextX = neighborhoodIntersection(&neighborhood, X);
            
            std::vector<NodeID> nextR = *R;
            nextR.emplace_back(n);
            
            bronKerbosch(graph, &nextP, &nextR, &nextX);
            // remove n from P
            auto f = std::find_if(nextP.begin(), nextP.end(), 
                [=](const auto el1) {
                    return (el1.first == n);
                });

            if (f != nextP.end()) {
                nextP.erase(f);
            }

            *P = nextP;

            auto f2 = std::lower_bound(nextX.begin(), nextX.end(), n,
                [](const auto& p1, const auto& p2) {
                    return p1.first < p2;
                });

            if (f2 == nextX.end() || (*f2).first != n) {
                nextX.insert(f2, std::make_pair(n, 0));
            }

            *X = nextX;
        }


    }

};