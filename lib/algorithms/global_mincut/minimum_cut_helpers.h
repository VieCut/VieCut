/******************************************************************************
 * minimum_cut_helpers.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include "common/configuration.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"

class minimum_cut_helpers {
 private:
    static constexpr bool debug = false;
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

    // set minimum cut to initial value (one minimum degree vertex)
    // - this cut will be updated later in the global_mincut
    static void setInitialCutValues(
        const std::vector<std::shared_ptr<graph_access> >& graphs) {
        if (configuration::getConfig()->save_cut) {
            size_t minimum_index = minimumIndex(graphs.back());

            for (NodeID idx : graphs[0]->nodes()) {
                if (idx == minimum_index) {
                    graphs[0]->setNodeInCut(idx, true);
                } else {
                    graphs[0]->setNodeInCut(idx, false);
                }
            }
        }
    }

    static EdgeWeight updateCut(
        const std::vector<std::shared_ptr<graph_access> >& graphs,
        EdgeWeight previous_mincut) {
        if (configuration::getConfig()->save_cut) {
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
                        } else {
                            graphs[0]->setNodeInCut(idx, false);
                        }
                    }
                }
            }
        }
        if (graphs.back()->number_of_nodes() > 1) {
            return std::min(previous_mincut, graphs.back()->getMinDegree());
        } else {
            return previous_mincut;
        }
    }

    static void retrieveMinimumCut(
        std::vector<std::shared_ptr<graph_access> > graphs) {
        std::shared_ptr<graph_access> G = graphs[0];
        size_t inside = 0, outside = 0;
        for (NodeID n : G->nodes()) {
            if (G->getNodeInCut(n)) {
                inside++;
                G->setPartitionIndex(n, 0);
            } else {
                outside++;
                G->setPartitionIndex(n, 1);
            }
        }

        [[maybe_unused]] size_t smaller = 0;
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

        LOG << "smaller side has " << std::min(inside, outside) << " nodes.";
    }

    static std::pair<std::vector<NodeID>,
                     std::vector<std::vector<NodeID> > > remap_cluster(
        std::shared_ptr<graph_access> G,
        const std::vector<NodeID>& cluster_id) {

        std::vector<NodeID> mapping;
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

            mapping.emplace_back(part[cur_cluster]);

            if (configuration::getConfig()->save_cut) {
                G->setPartitionIndex(node, part[cur_cluster]);
            }
            reverse_mapping[part[cur_cluster]].push_back(node);
        }

        return std::make_pair(mapping, reverse_mapping);
    }

    static std::pair<std::vector<NodeID>, std::unordered_map<NodeID, NodeID> >
    reInsertVertices(std::shared_ptr<mutable_graph> out_graph,
                     const std::vector<std::shared_ptr<graph_access> >& graphs,
                     const std::vector<size_t>& ge_ids,
                     const std::vector<std::vector<
                                           std::pair<NodeID, NodeID> > >&
                     guaranteed_edges,
                     const EdgeWeight global_mincut) {
        NodeID n0 = graphs[0]->number_of_nodes();
        std::vector<NodeID> out_graph_mapping;

        for (size_t i = 0; i < out_graph->getOriginalNodes(); ++i) {
            out_graph_mapping.emplace_back(out_graph->getCurrentPosition(i));
        }

        std::unordered_map<NodeID, NodeID> deleted_vertex_mappings;

        for (size_t i = ge_ids.size(); i--; ) {
            size_t current_id = ge_ids[i];
            for (auto e : guaranteed_edges[i]) {
                graphs[current_id]->setPartitionIndex(e.first, UNDEFINED_NODE);
                NodeID vertex_id = n0 * current_id + e.first;
                NodeID new_node = out_graph->new_empty_node();
                deleted_vertex_mappings[vertex_id] = new_node;

                NodeID neighbour = e.second;

                for (size_t j = current_id; j < graphs.size() - 1; ++j) {
                    neighbour = graphs[j]->getPartitionIndex(neighbour);
                    if (neighbour == UNDEFINED_NODE) {
                        NodeID deleted_id = n0 * j + neighbour;
                        neighbour = deleted_vertex_mappings[deleted_id];
                        break;
                    }
                }

                out_graph->new_edge_order(new_node, neighbour, global_mincut);
                for (NodeID n : out_graph->containedVertices(neighbour)) {
                    out_graph->setCurrentPosition(n, UNDEFINED_NODE);
                }
            }
        }
        return std::make_pair(out_graph_mapping, deleted_vertex_mappings);
    }

    static void setVertexLocations(
        std::shared_ptr<mutable_graph> out_graph,
        const std::vector<std::shared_ptr<graph_access> >& graphs,
        const std::vector<NodeID>& out_graph_mapping,
        const std::unordered_map<NodeID, NodeID>& deleted_vertex_mappings) {
        if (configuration::getConfig()->save_cut) {
            NodeID n0 = graphs[0]->number_of_nodes();
            out_graph->setOriginalNodes(graphs[0]->number_of_nodes());
            std::vector<NodeID> final;

            for (NodeID n : out_graph->nodes()) {
                std::vector<NodeID> empty;
                out_graph->setContainedVertices(n, empty);
            }

            std::vector<NodeID> b(out_graph->n(), 0);

            bool broken = false;
            for (NodeID idx : graphs[0]->nodes()) {
                broken = false;
                NodeID coarseID = idx;
                for (size_t lv = 0; lv < graphs.size() - 1; ++lv) {
                    if (graphs[lv]->getPartitionIndex(coarseID)
                        == UNDEFINED_NODE) {
                        NodeID deleted_id = (n0 * lv) + coarseID;
                        coarseID = deleted_vertex_mappings.at(deleted_id);
                        broken = true;
                        break;
                    }
                    coarseID = graphs[lv]->getPartitionIndex(coarseID);
                }

                // if the vertex was deleted previously, it was readded
                // and deleted_vertex_mappings holds the correct id
                if (!broken) {
                    coarseID = out_graph_mapping[coarseID];
                }

                out_graph->setCurrentPosition(idx, coarseID);
                out_graph->addContainedVertex(coarseID, idx);

                ++b[coarseID];
            }

            static constexpr bool cut_logs = true;

            if (cut_logs) {
                std::sort(b.begin(), b.end());
                printLogs(b);
            }
        }
    }

    static void printLogs(const std::vector<NodeID>& b) {
        NodeID empty = std::lower_bound(b.begin(), b.end(), 1) - b.begin();
        NodeID id1 = b.end() - std::upper_bound(b.begin(), b.end(), 1);
        NodeID id10 = b.end() - std::upper_bound(b.begin(), b.end(), 9);
        NodeID id100 = b.end() - std::upper_bound(b.begin(), b.end(), 99);
        NodeID id1000 = b.end() - std::upper_bound(b.begin(), b.end(), 999);
        NodeID id10000 = b.end() - std::upper_bound(b.begin(), b.end(), 9999);

        LOG1 << "--------------------------------------------------";
        LOG1 << "Cut stats:";
        LOG1 << "Largest block: " << b.back();
        LOG1 << "Number of empty blocks: " << empty;
        LOG1 << "Number of nontrivial blocks: " << id1;
        LOG1 << "Number of blocks of size >=10: " << id10;
        LOG1 << "Number of blocks of size >=100: " << id100;
        LOG1 << "Number of blocks of size >=1000: " << id1000;
        LOG1 << "Number of blocks of size >=10000: " << id10000;
        LOG1 << "--------------------------------------------------";
    }
};
