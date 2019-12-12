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
                auto [tgt, wgt] = problem->graph->getEdge(n, e);
                P.emplace_back(tgt, wgt);
            }

            std::sort(P.begin(), P.end(), [](const auto& p1, const auto& p2) {
                return p1.first < p2.first;
            });

            X = P;



            bronKerbosch(problem->graph, &P, &R, &X);
        }
    }

    void bronKerbosch(std::shared_ptr<mutable_graph> graph, 
                      std::vector<std::pair<NodeID, EdgeWeight> >* P,
                      std::vector<NodeID>* R,
                      std::vector<std::pair<NodeID, EdgeWeight> >* X) {
        
        
    }

};