/******************************************************************************
 * minimum_cut_helpers.h
 *
 * Source of VieCut
 *
 ******************************************************************************
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

#include <data_structure/graph_access.h>
#include <memory>

class minimum_cut_helpers
{
private:
    // Get index of minimum degree vertex
    static size_t minimumIndex(std::shared_ptr<graph_access> G) {
        size_t minimum_index = 0;
        for (NodeID n : G->nodes()) {
            if (G->getWeightedNodeDegree(n) == G->getMinDegree()) {
                minimum_index = n;
                break;
            }
        }
        return minimum_index;
    }

public:
    static bool graphValid(std::shared_ptr<graph_access> G) {

        // graph does not exist
        if (!G.use_count())
            return false;

        // graph has no nodes
        if (!G->number_of_nodes())
            return false;

        // in process of building graph
        if (G->currentlyBuildingGraph())
            return false;

        return true;
    }

    // set minimum cut to initial value (one minimum degree vertex) - this cut will be updated later in the global_mincut
    static void setInitialCutValues(__attribute__ ((unused)) std::vector<std::shared_ptr<graph_access> > graphs) {
#ifdef SAVECUT
        size_t minimum_index = minimumIndex(graphs.back());

        for (NodeID idx : graphs[0]->nodes()) {
            if (idx == minimum_index) {
                graphs[0]->setNodeInCut(idx, true);
            }
            else {
                graphs[0]->setNodeInCut(idx, false);
            }
        }
#endif
    }

    static EdgeWeight updateCutValueAfterContraction(__attribute__ ((unused)) std::vector<std::shared_ptr<graph_access> > graphs,
                                                     EdgeWeight previous_mincut) {
#ifdef SAVECUT
        std::shared_ptr<graph_access> new_graph = graphs.back();
        if (new_graph->number_of_nodes() > 1) {
            if (new_graph->getMinDegree() < previous_mincut) {
                size_t minimum_index = minimumIndex(graphs.back());

                for (NodeID idx : graphs[0]->nodes()) {
                    NodeID coarseID = idx;
                    for (size_t lv = 0; lv < graphs.size() - 1; ++lv) {
                        coarseID = graphs[lv]->getPartitionIndex(coarseID);
                    }
                    if (coarseID == minimum_index) {
                        graphs[0]->setNodeInCut(idx, true);
                    }
                    else {
                        graphs[0]->setNodeInCut(idx, false);
                    }
                }
            }
        }
#endif
        if (graphs.back()->number_of_nodes() > 1) {
            return std::min<EdgeWeight>(previous_mincut, graphs.back()->getMinDegree());
        }
        else {
            return previous_mincut;
        }
    }

    static void retrieveMinimumCut(__attribute__ ((unused)) std::vector<std::shared_ptr<graph_access> > graphs) {

#ifdef SAVECUT
        std::shared_ptr<graph_access> G = graphs[0];

        size_t inside = 0, outside = 0;
        for (NodeID n : G->nodes()) {
            if (G->getNodeInCut(n)) {
                inside++;
                G->setPartitionIndex(n, 0);
            }
            else {
                outside++;
                G->setPartitionIndex(n, 1);
            }
        }

        size_t __attribute__ ((unused)) smaller = 0;
        if (inside > outside) {
            smaller = 1;
        }
#ifndef NDEBUG
        for (NodeID n : G->nodes()) {
            if (G->getPartitionIndex(n) == smaller) {
                std::cout << "n " << n << std::endl;
            }
        }
#endif

        LOG1 << "smaller side of cut has " << std::min(inside, outside) << " nodes.";
#endif // SAVECUT
    }

    static std::vector<std::vector<NodeID> > remap_cluster(
        std::shared_ptr<graph_access> G, std::vector<NodeID>& cluster_id) {

        std::vector<std::vector<NodeID> > reverse_mapping;

        PartitionID cur_no_clusters = 0;
        std::unordered_map<PartitionID, PartitionID> remap;

        std::vector<NodeID> part(G->number_of_nodes(), UNDEFINED_NODE);

        for (NodeID node : G->nodes()) {
            PartitionID cur_cluster = cluster_id[node];
            // check whether we already had that
            if (part[cur_cluster] == UNDEFINED_NODE) {
                part[cur_cluster] = cur_no_clusters++;
                reverse_mapping.emplace_back();
            }

            cluster_id[node] = part[cur_cluster];
#ifdef SAVECUT
            G->setPartitionIndex(node, part[cur_cluster]);
#endif
            reverse_mapping[part[cur_cluster]].push_back(node);
        }

        return reverse_mapping;
    }
};
