/******************************************************************************
 * matula_approx.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 * 
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "minimum_cut.h"
#include "noi_minimum_cut.h"

#include "tools/timer.h"

#ifdef PARALLEL
#include "parallel/coarsening/contract_graph.h"
#else
#include "coarsening/contract_graph.h"
#endif

#include "definitions.h"
#include "minimum_cut.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <memory>
#include <unordered_map>

#include <cstdint>
#include <cstdlib>
class matula_approx : public minimum_cut
{
public:
    matula_approx() { }
    virtual ~matula_approx() { }
    static constexpr bool debug = false;

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G) {

        if (!minimum_cut_helpers::graphValid(G))
            return -1;

        std::vector<std::shared_ptr<graph_access> > graphs;
        graphs.push_back(G);
        EdgeWeight global_mincut = G->getMinDegree();

#ifdef SAVECUT
        size_t minindex = 0;
        for (NodeID n : G->nodes()) {
            if (G->getWeightedNodeDegree(n) == G->getMinDegree()) {
                minindex = n;
                break;
            }
        }

        for (NodeID idx : graphs[0]->nodes()) {
            if (idx == minindex) {
                graphs[0]->setNodeInCut(idx, true);
            }
            else {
                graphs[0]->setNodeInCut(idx, false);
            }
        }

#endif
        while (graphs.back()->number_of_nodes() > 2 && global_mincut > 0) {

#ifdef SAVECUT

            if (graphs.back()->getMinDegree() < global_mincut) {
                size_t minindex = 0;
                for (NodeID n : graphs.back()->nodes()) {
                    if (graphs.back()->getWeightedNodeDegree(n) == graphs.back()->getMinDegree()) {
                        minindex = n;
                        break;
                    }
                }

                for (NodeID idx : graphs[0]->nodes()) {
                    NodeID coarseID = idx;
                    for (size_t lv = 0; lv < graphs.size() - 1; ++lv) {
                        coarseID = graphs[lv]->getPartitionIndex(coarseID);
                    }
                    if (coarseID == minindex) {
                        graphs[0]->setNodeInCut(idx, true);
                    }
                    else {
                        graphs[0]->setNodeInCut(idx, false);
                    }
                }
            }

#endif
            std::shared_ptr<graph_access> curr_g = graphs.back();
            noi_minimum_cut noi;
            std::vector<std::pair<NodeID, NodeID> > contractable;
            size_t j = global_mincut / 2;

            if (!j)
                return global_mincut;

            union_find uf(curr_g->number_of_nodes());
            noi.modified_capforest(curr_g, j, uf, graphs);
            global_mincut = std::min(global_mincut, graphs.back()->getMinDegree());

            LOG << "Global mincut now: " << global_mincut;

            size_t num_contractable = 0;

            for (auto edge : contractable) {

                NodeID first = uf.Find(edge.first);
                NodeID second = uf.Find(edge.second);
                if (first != second) {
                    uf.Union(first, second);
                    ++num_contractable;
                }
            }

            num_contractable = curr_g->number_of_nodes() - uf.n();
            std::vector<std::vector<NodeID> > reverse_mapping(
                curr_g->number_of_nodes() - num_contractable);
            std::vector<NodeID> mapping(curr_g->number_of_nodes());
            std::vector<NodeID> part(curr_g->number_of_nodes(), UNDEFINED_NODE);
            NodeID current_pid = 0;

            for (NodeID n : curr_g->nodes()) {
                NodeID part_id = uf.Find(n);
                if (part[part_id] == UNDEFINED_NODE) {
                    part[part_id] = current_pid++;
                }
                mapping[n] = part[part_id];
                curr_g->setPartitionIndex(n, part[part_id]);
                reverse_mapping[part[part_id]].push_back(n);
            }

            assert(current_pid == reverse_mapping.size());
            std::shared_ptr<graph_access> new_graph = contraction::contractGraph(
                curr_g, mapping,
                reverse_mapping.size(),
                reverse_mapping);

            graphs.push_back(new_graph);

            if (new_graph->number_of_nodes() > 1) {
#ifdef SAVECUT
                if (new_graph->getMinDegree() < global_mincut) {
                    size_t minindex = 0;

                    for (NodeID n : G->nodes()) {
                        if (new_graph->getWeightedNodeDegree(n) ==
                            new_graph->getMinDegree()) {
                            minindex = n;
                            break;
                        }
                    }

                    for (NodeID idx : graphs[0]->nodes()) {
                        NodeID coarseID = idx;
                        for (size_t lv = 0; lv < graphs.size() - 1; ++lv) {
                            coarseID = graphs[lv]->getPartitionIndex(coarseID);
                        }
                        if (coarseID == minindex) {
                            graphs[0]->setNodeInCut(idx, true);
                        }
                        else {
                            graphs[0]->setNodeInCut(idx, false);
                        }
                    }
                }
#endif
                global_mincut = std::min(graphs.back()->getMinDegree(), global_mincut);
            }
        }

        minimum_cut_helpers::retrieveMinimumCut(graphs);

        return global_mincut;
    }
};
