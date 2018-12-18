/******************************************************************************
 * graph_extractor.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once

#include "data_structure/graph_access.h"
#include "definitions.h"
#include <memory>

class graph_extractor
{
public:
    graph_extractor() { }
    virtual ~graph_extractor() { }

    std::pair<std::shared_ptr<graph_access>, std::vector<NodeID> > extract_block(std::shared_ptr<graph_access> G,
                                                                                 PartitionID block) {

        std::shared_ptr<graph_access> extracted_block = std::make_shared<graph_access>();
        // build reverse mapping
        std::vector<NodeID> reverse_mapping;
        NodeID nodes = 0;
        NodeID dummy_node = G->number_of_nodes() + 1;
        for (NodeID node : G->nodes()) {
            if (G->getPartitionIndex(node) == block) {
                reverse_mapping.push_back(nodes++);
            }
            else {
                reverse_mapping.push_back(dummy_node);
            }
        }

        extracted_block->start_construction(nodes, G->number_of_edges());

        for (NodeID node : G->nodes()) {
            if (G->getPartitionIndex(node) == block) {
                NodeID new_node = extracted_block->new_node();
                // extracted_block.setNodeWeight( new_node, G.getNodeWeight(node));
                for (EdgeID e : G->edges_of(node)) {
                    NodeID target = G->getEdgeTarget(e);
                    if (G->getPartitionIndex(target) == block) {
                        EdgeID new_edge = extracted_block->new_edge(new_node, reverse_mapping[target]);
                        extracted_block->setEdgeWeight(new_edge, G->getEdgeWeight(e));
                    }
                }
            }
        }

        extracted_block->finish_construction();
        return std::make_pair(extracted_block, reverse_mapping);
    }
};
