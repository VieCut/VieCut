/******************************************************************************
 * graph_algorithms.h
 * 
 * Source of VieCut
 * 
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 * 
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <data_structure/graph_access.h>
#include <definitions.h>
#include <memory>


class graph_algorithms
{
public:
    static std::vector<NodeID> top_k_degrees(std::shared_ptr<graph_access> G, size_t k) {

        std::vector<std::pair<NodeID, EdgeWeight> > all_degrees;
        for (NodeID n : G->nodes()) {
            all_degrees.emplace_back(n, G->getNodeDegree(n));
        }

        return find_top_k(all_degrees, k);
    }

    static std::vector<NodeID> weighted_top_k_degrees(std::shared_ptr<graph_access> G, size_t k) {

        std::vector<std::pair<NodeID, EdgeWeight> > all_degrees;
        for (NodeID n : G->nodes()) {
            all_degrees.emplace_back(n, G->getWeightedNodeDegree(n));
        }

        return find_top_k(all_degrees, k);
    }

private:
    static std::vector<NodeID> find_top_k(std::vector<std::pair<NodeID, EdgeWeight> >& in, size_t k) {
        std::nth_element(in.begin(), in.end() - k, in.end(), [](auto a1, auto a2) {
                             return a1.second < a2.second;
                         });

        std::vector<NodeID> out;
        for (size_t i = in.size() - k; i < in.size(); ++i) {
            out.emplace_back(in[i].first);
        }
        return out;
    }
};
