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

    using NeighboursAndContents =
        std::vector<std::tuple<std::vector<std::pair<NodeID, EdgeWeight> >,
                               std::vector<NodeID>, NodeID> >;
    heavy_edges(EdgeWeight mincut) : mincut(mincut) { }

    std::vector<std::tuple<NodeID, std::vector<NodeID> > > removeHeavyEdges(
        std::shared_ptr<mutable_graph> G) {
        std::vector<std::tuple<NodeID, std::vector<NodeID> > > cactusEdge;
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
                cactusEdge.emplace_back(vtx_in_t, G->containedVertices(n));
                G->deleteVertex(n);
            }
        }
        return cactusEdge;
    }

    std::vector<std::tuple<std::pair<NodeID, NodeID>, std::vector<NodeID> > >
    contractCycleEdges(std::shared_ptr<mutable_graph> G) {
        std::vector<std::tuple<std::pair<NodeID, NodeID>,
                               std::vector<NodeID> > > cycleEdges;
        if (mincut % 2 != 0)
            return cycleEdges;

        // as we contract edges, we use basic for loop so G->n() can update
        for (NodeID n = 0; n < G->n(); ++n) {
            if (G->get_first_invalid_edge(n) == 2
                && G->getWeightedNodeDegree(n) == mincut) {
                NodeID n0 = G->getEdgeTarget(n, 0);
                NodeID n1 = G->getEdgeTarget(n, 1);
                if (G->isEmpty(n0) || G->isEmpty(n1))
                    continue;
                // if the edges have different weights, the heavier of them will
                // be contracted in local routines
                if (G->getEdgeWeight(n, 0) != mincut / 2
                    || G->getEdgeWeight(n, 1) != mincut / 2)
                    continue;

                NodeID p0 = G->containedVertices(n0)[0];
                NodeID p1 = G->containedVertices(n1)[0];
                auto contained = G->containedVertices(n);
                G->setContainedVertices(n, { });
                for (const auto& c : contained) {
                    G->setCurrentPosition(c, UNDEFINED_NODE);
                }

                G->contractEdgeSparseTarget(n0, G->getReverseEdge(n, 0));
                cycleEdges.emplace_back(std::make_pair(p0, p1), contained);
            }
        }

        return cycleEdges;
    }

    // Return type: std::vector<std::pair<std::vector<std::pair<NodeID,
    // EdgeWeight> >, std::vector<NodeID> > >;
    NeighboursAndContents contractMinimumDegreeVertices(
        std::shared_ptr<mutable_graph> G) {
        NeighboursAndContents nbc;
        if (mincut % 2 == 0)
            return nbc;
        for (NodeID n = 0; n < G->n(); ++n) {
            if (G->getWeightedNodeDegree(n) == mincut) {
                bool empty_neighbor = false;
                EdgeID non_minimum_neighbour = UNDEFINED_EDGE;
                for (EdgeID e : G->edges_of(n)) {
                    NodeID ngbr = G->getEdgeTarget(n, e);
                    if (G->getWeightedNodeDegree(ngbr) != mincut
                        && non_minimum_neighbour == UNDEFINED_EDGE) {
                        non_minimum_neighbour = e;
                    }

                    if (G->isEmpty(ngbr)) {
                        empty_neighbor = true;
                        LOG1 << "empty neighbour!";
                        break;
                    }
                }
                if (empty_neighbor || non_minimum_neighbour == UNDEFINED_EDGE)
                    continue;

                LOG1 << "k";
                nbc.emplace_back();
                auto& neighbours = std::get<0>(nbc.back());
                for (EdgeID e : G->edges_of(n)) {
                    NodeID t = G->containedVertices(G->getEdgeTarget(n, e))[0];
                    EdgeWeight w = G->getEdgeWeight(n, e);
                    neighbours.emplace_back(t, w);
                }
                auto& contents = std::get<1>(nbc.back());
                contents = G->containedVertices(n);
                G->setContainedVertices(n, { });
                for (const auto& c : contents) {
                    G->setCurrentPosition(c, UNDEFINED_NODE);
                }
                auto& ctr = std::get<2>(nbc.back());
                ctr = non_minimum_neighbour;
                // TODO(anoe): use reverse edge when checked that it works!
                G->contractEdge(n, non_minimum_neighbour);
            }
        }

        return nbc;
    }

    void reInsertMinimumDegree(std::shared_ptr<mutable_graph> G,
                               NeighboursAndContents nbc) {
        for (size_t i = nbc.size(); i-- > 0; ) {
            const auto& [neighbours, contents, ctr_ngbr] = nbc[i];
            NodeID p = neighbours[ctr_ngbr].first;
            NodeID n = G->getCurrentPosition(p);
            NodeID other_ngbr = UNDEFINED_NODE;
            EdgeWeight wgt = 0;

            for (size_t j = 0; j < neighbours.size(); ++j) {
                if (j == ctr_ngbr)
                    continue;
                NodeID pj = neighbours[j].first;
                NodeID nj = G->getCurrentPosition(pj);
                if (nj != n) {
                    wgt += neighbours[j].second;
                    other_ngbr = nj;
                }
            }
            NodeID reIns = G->new_empty_node();
            /*if (other_ngbr != UNDEFINED_NODE) {
                VIECUT_ASSERT_EQ(wgt, mincut / 2);
                EdgeID e = UNDEFINED_EDGE;
                for (EdgeID arc : G->edges_of(n)) {
                    if (G->getEdgeTarget(n, arc) == other_ngbr) {
                        e = arc;
                        break;
                    }
                }
                G->new_edge_order(n, reIns, mincut / 2);
                G->new_edge_order(other_ngbr, reIns, mincut / 2);

                EdgeWeight w01 = G->getEdgeWeight(n, e);
                if (w01 == (mincut / 2)) {
                    G->deleteEdge(n, e);
                } else {
                    G->setEdgeWeight(n, e, w01 - (mincut / 2));
                }
            } else {*/
            G->new_edge_order(n, reIns, mincut);
            // }

            G->setContainedVertices(reIns, contents);
            for (NodeID v : contents) {
                G->setCurrentPosition(v, reIns);
            }
        }
    }

    void reInsertCycles(
        std::shared_ptr<mutable_graph> G,
        std::vector<std::tuple<std::pair<NodeID, NodeID>,
                               std::vector<NodeID> > > toInsert) {
        for (size_t i = toInsert.size(); i-- > 0; ) {
            const auto& [p, cont] = toInsert[i];
            NodeID n0 = G->getCurrentPosition(p.first);
            NodeID n1 = G->getCurrentPosition(p.second);
            NodeID reIns = G->new_empty_node();
            if (n0 == n1) {
                G->new_edge_order(n0, reIns, mincut);
                G->setContainedVertices(reIns, cont);
                for (NodeID v : cont) {
                    G->setCurrentPosition(v, reIns);
                }
            } else {
                EdgeID e = UNDEFINED_EDGE;
                for (EdgeID arc : G->edges_of(n0)) {
                    if (G->getEdgeTarget(n0, arc) == n1) {
                        e = arc;
                        break;
                    }
                }
                G->new_edge_order(n0, reIns, mincut / 2);
                G->new_edge_order(n1, reIns, mincut / 2);
                EdgeWeight w01 = G->getEdgeWeight(n0, e);
                if (w01 == (mincut / 2)) {
                    G->deleteEdge(n0, e);
                } else {
                    G->setEdgeWeight(n0, e, w01 - (mincut / 2));
                }
            }
            G->setContainedVertices(reIns, cont);
            for (NodeID v : cont) {
                G->setCurrentPosition(v, reIns);
            }
        }
    }

    void reInsertVertices(
        std::shared_ptr<mutable_graph> G,
        std::vector<std::tuple<NodeID, std::vector<NodeID> > > toInsert) {
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

 private:
    EdgeWeight mincut;
};
