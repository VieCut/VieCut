/******************************************************************************
 * measurements.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2020 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/
#pragma once

#include <vector>

#include "common/definitions.h"
#include "data_structure/mutable_graph.h"

class measurements {
 public:
    static FlowType flowValue(const mutable_graph& original_graph,
                              const std::vector<NodeID>& original_terminals,
                              bool verbose, const std::vector<NodeID>& sol) {
        std::vector<size_t> block_sizes(original_terminals.size(), 0);

        EdgeWeight total_weight = 0;
        for (NodeID n : original_graph.nodes()) {
            block_sizes[sol[n]]++;
            for (EdgeID e : original_graph.edges_of(n)) {
                NodeID tgt = original_graph.getEdgeTarget(n, e);
                EdgeWeight wgt = original_graph.getEdgeWeight(n, e);
                if (sol[n] != sol[tgt]) {
                    total_weight += wgt;
                }
            }
        }
        total_weight /= 2;

        for (size_t i = 0; i < block_sizes.size(); ++i) {
            size_t terminals = 0;
            LOGC(verbose) << "Block " << i << " has size " << block_sizes[i];
            if (block_sizes[i] == 0) {
                LOG1 << "ERROR: Empty block " << i;
                exit(1);
            }

            for (NodeID ot : original_terminals) {
                if (sol[ot] == i) {
                    terminals++;
                }
            }

            if (terminals != 1) {
                LOG1 << terminals << " terminals in block " << i;
                exit(1);
            }
        }
        return total_weight;
    }
};
