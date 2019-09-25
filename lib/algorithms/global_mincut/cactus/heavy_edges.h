/******************************************************************************
 * heavy_edges.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common/definitions.h"
#include "data_structure/mutable_graph.h"

class heavy_edges {
 public:
    static void removeHeavyEdges(
        std::shared_ptr<mutable_graph> G,
        std::vector<std::tuple<NodeID, std::vector<NodeID> > >* cactusEdge,
        EdgeWeight mincut) {
        std::unordered_map<NodeID, std::vector<NodeID> > contract;
        std::vector<NodeID> markForCactus;
        for (NodeID n : G->nodes()) {
            if (G->isEmpty(n))
                continue;

            for (EdgeID e : G->edges_of(n)) {
                EdgeWeight wgt = G->getEdgeWeight(n, e);
                NodeID target = G->getEdgeTarget(n, e);

                if (G->isEmpty(target))
                    continue;

                if (wgt > mincut) {
                    NodeID v1 = G->containedVertices(n)[0];
                    NodeID v2 = G->containedVertices(target)[0];
                    NodeID min = std::min(v1, v2);
                    NodeID max = std::max(v1, v2);

                    if (contract.find(min) == contract.end()) {
                        contract[min] = { max };
                    } else {
                        contract[min].emplace_back(max);
                    }
                }

                if (wgt == mincut) {
                    if (G->get_first_invalid_edge(n) == 1) {
                        // each edge is seen from both adjacent nodes
                        // so we get all edges
                        markForCactus.emplace_back(G->containedVertices(n)[0]);
                    }
                }
            }
        }

        for (const auto& [lowest, others] : contract) {
            std::unordered_set<NodeID> vtxset;
            vtxset.insert(G->getCurrentPosition(lowest));
            for (const auto& v : others) {
                vtxset.insert(G->getCurrentPosition(v));
            }
            if (vtxset.size() > 1) {
                G->contractVertexSet(vtxset);
            }
        }

        for (const NodeID& e : markForCactus) {
            if (G->n() > 2) {
                NodeID n = G->getCurrentPosition(e);
                VIECUT_ASSERT_EQ(G->get_first_invalid_edge(n), 1);
                NodeID t = G->getEdgeTarget(n, 0);
                if (G->isEmpty(t)) {
                    continue;
                }
                NodeID vtx_in_t = G->containedVertices(t)[0];
                cactusEdge->emplace_back(vtx_in_t, G->containedVertices(n));
                G->deleteVertex(n);
            }
        }
    }

    static void reInsertVertices(
        std::shared_ptr<mutable_graph> G,
        std::vector<std::tuple<NodeID, std::vector<NodeID> > > toInsert,
        EdgeWeight mincut) {
        for (size_t i = toInsert.size(); i-- > 0; ) {
            const auto& [t, cont] = toInsert[i];
            NodeID curr = G->getCurrentPosition(t);
            NodeID vtx = G->new_empty_node();
            G->new_edge_order(curr, vtx, mincut);
            G->setContainedVertices(vtx, cont);
            for (const auto& e : G->containedVertices(vtx)) {
                G->setCurrentPosition(e, vtx);
            }
        }
    }
};
