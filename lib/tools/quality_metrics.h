/******************************************************************************
 * quality_metrics.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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

#ifndef QUALITY_METRICS_10HC2I5M
#define QUALITY_METRICS_10HC2I5M

#include "data_structure/graph_access.h"

class quality_metrics
{
public:
    quality_metrics() { }

    ~quality_metrics() { }

    EdgeWeight edge_cut(graph_access& G) {
        EdgeWeight edgeCut = 0;
        for (NodeID n : G.nodes()) {
            PartitionID partitionIDSource = G.getPartitionIndex(n);
            for (EdgeID e : G.edges_of(n)) {
                NodeID targetNode = G.getEdgeTarget(e);
                PartitionID partitionIDTarget = G.getPartitionIndex(targetNode);

                if (partitionIDSource != partitionIDTarget) {
                    edgeCut += G.getEdgeWeight(e);
                }
            }
        }
        return edgeCut / 2;
    }
};

#endif /* end of include guard: QUALITY_METRICS_10HC2I5M */
