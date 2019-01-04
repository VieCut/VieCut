/******************************************************************************
 * label_propagation.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@univie.ac.at>
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include "data_structure/graph_access.h"
#include "definitions.h"
#include "tlx/logger.hpp"
#include "tools/random_functions.h"

#include <cstdint>
#include <cstdlib>

#include <random>

class label_propagation
{
    static constexpr bool debug = false;
    static constexpr bool timing = true;

public:
    label_propagation() { }
    virtual ~label_propagation() { }

    std::vector<NodeID> propagate_labels(std::shared_ptr<graph_access> G) {

        std::vector<NodeID> cluster_id(G->number_of_nodes());
        std::vector<NodeID> permutation(G->number_of_nodes());

        // use this vector to store connection strengths to blocks.
        // in order to only reset elements that are actually written in we store (in second)
        // which vertex last used this element.
        std::vector<std::pair<EdgeWeight, NodeID> > hash_vec(G->number_of_nodes(),
                                                             std::make_pair<EdgeWeight, NodeID>(0, 0));

        for (size_t i = 0; i < cluster_id.size(); ++i) {
            cluster_id[i] = i;
        }

        random_functions::permutate_vector_local(permutation, true);

        size_t iterations = 2;

        LOG << "Number of iterations: " << iterations;

        timer t;

        // size_t change = G.number_of_nodes();
        for (size_t j = 0; j < iterations; j++) {
            for (NodeID node : G->nodes()) {
                NodeID n = permutation[node];
                PartitionID max_block = cluster_id[n];
                EdgeWeight max_value = 0;

                for (EdgeID e : G->edges_of(n)) {

                    NodeID target = G->getEdgeTarget(e);
                    PartitionID block = cluster_id[target];

                    // we want to store connection of vertex n to block 'block'
                    // thus reset now that we look at outgoing edges from n
                    if (hash_vec[block].second != n) {
                        hash_vec[block].first = 0;
                        hash_vec[block].second = n;
                    }

                    hash_vec[block].first += G->getEdgeWeight(e);

                    // set strongest connected block if higher connection strength, random if equal
                    if (hash_vec[block].first > max_value ||
                        (hash_vec[block].first == max_value &&
                         random_functions::next() % 2)) {
                        max_value = hash_vec[block].first;
                        max_block = block;
                    }
                }

                cluster_id[n] = max_block;
            }
            LOGC(timing) << "LP: Iteration " << j << " -  Timer: " << t.elapsedToZero();
        }

        return cluster_id;
    }

    void certify_clusters(
        std::shared_ptr<graph_access> G,
        std::vector<NodeID>& cluster_id,
        std::vector<std::vector<NodeID> >& reverse_mapping) {

        for (size_t p = 0; p < reverse_mapping.size(); ++p) {
            std::vector<EdgeWeight> in (reverse_mapping[p].size());
            std::vector<EdgeWeight> out(reverse_mapping[p].size());

            for (size_t i = 0; i < reverse_mapping[p].size(); ++i) {
                for (EdgeID e : G->edges_of(reverse_mapping[p][i])) {
                    if (cluster_id[G->getEdgeTarget(e)] == p) {
                        in[i] += G->getEdgeWeight(e);
                    }
                    else {
                        out[i] += G->getEdgeWeight(e);
                    }
                }
            }
        }
    }
};
