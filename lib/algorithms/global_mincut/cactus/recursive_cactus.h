/******************************************************************************
 * recursive_cactus.h
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
#include <map>
#include <memory>
#include <queue>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/flow/push_relabel.h"
#include "algorithms/global_mincut/noi_minimum_cut.h"
#include "algorithms/misc/graph_algorithms.h"
#include "algorithms/misc/strongly_connected_components.h"
#include "algorithms/multicut/multicut_problem.h"
#include "common/configuration.h"
#include "common/definitions.h"
#include "data_structure/mutable_graph.h"
#include "data_structure/priority_queues/node_bucket_pq.h"
#include "tlx/logger.hpp"
#include "tools/string.h"

#ifdef PARALLEL
#include "parallel/coarsening/contract_graph.h"
#include "parallel/data_structure/union_find.h"
#else
#include "coarsening/contract_graph.h"
#include "data_structure/union_find.h"
#endif

class recursive_cactus {
 public:
    recursive_cactus() { }
    explicit recursive_cactus(EdgeWeight mincut) : mincut(mincut) { }
    ~recursive_cactus() { }

    static constexpr bool debug = false;
    bool timing = configuration::getConfig()->verbose;

    void setMincut(EdgeWeight mc) {
        mincut = mc;
    }

    std::shared_ptr<mutable_graph> flowMincut(
        const std::vector<std::shared_ptr<graph_access> >& graphs) {
        std::vector<std::shared_ptr<mutable_graph> > flow_graphs;
        flow_graphs.emplace_back(
            mutable_graph::from_graph_access(graphs.back()));
        std::vector<std::vector<std::pair<NodeID, NodeID> > > packed;
        std::vector<std::vector<std::pair<NodeID, NodeID> > > deleted;
        bool continuing = true;
        size_t num_vertices = graphs.back()->number_of_nodes();

        while (continuing) {
            auto flow_graph =
                std::make_shared<mutable_graph>(*flow_graphs.back());
            packed.emplace_back(flow_graph->number_of_nodes(),
                                std::make_pair(UNDEFINED_NODE, UNDEFINED_NODE));
            deleted.emplace_back(flow_graph->number_of_nodes(),
                                 std::make_pair(UNDEFINED_NODE,
                                                UNDEFINED_NODE));
            continuing = false;

            for (NodeID n = flow_graph->n(); n-- != 0; ) {
                if (flow_graph->get_first_invalid_edge(n) == 1
                    && num_vertices > 2) {
                    NodeID tgt = flow_graph->containedVertices(
                        flow_graph->getEdgeTarget(n, 0))[0];
                    NodeID original = flow_graph->containedVertices(n)[0];
                    continuing = true;
                    num_vertices--;

                    if (flow_graph->getEdgeWeight(n, 0) == mincut) {
                        packed.back()[n] = std::make_pair(tgt, original);
                    } else {
                        deleted.back()[n] = std::make_pair(tgt, original);
                    }
                }
            }

            for (NodeID n = flow_graph->number_of_nodes(); n-- != 0; ) {
                if (packed.back()[n].first != UNDEFINED_NODE
                    || deleted.back()[n].first != UNDEFINED_NODE) {
                    flow_graph->deleteVertex(n);
                }
            }
            flow_graphs.push_back(flow_graph);
            LOGC(timing) << "t " << t.elapsed() << " contracted from "
                         << packed.back().size() << " to "
                         << flow_graph->number_of_nodes() << " vertices!";
        }

        auto in_graph = flow_graphs.back()->simplify();

        auto out_graph = recursiveCactus(in_graph, 0);

        VIECUT_ASSERT_TRUE(isCNCR(out_graph));
        LOGC(timing) << "t " << t.elapsed() << " cactus n "
                     << out_graph->n() << " m " << out_graph->n();

        out_graph->setOriginalNodes(flow_graphs[0]->number_of_nodes());

        for (NodeID n = 0; n < out_graph->number_of_nodes(); ++n) {
            std::vector<NodeID> cv;
            for (NodeID v : out_graph->containedVertices(n)) {
                NodeID p = flow_graphs.back()->containedVertices(v)[0];
                cv.emplace_back(p);
                out_graph->setCurrentPosition(p, n);
            }
            out_graph->setContainedVertices(n, cv);
        }

        // unpack
        for (size_t lv = packed.size(); lv-- > 0; ) {
            for (size_t n = 0; n < packed[lv].size(); ++n) {
                if (packed[lv][n].first != UNDEFINED_NODE) {
                    auto [t, orig] = packed[lv][n];
                    NodeID newnode = out_graph->new_empty_node();
                    NodeID t_now = out_graph->getCurrentPosition(t);

                    out_graph->addContainedVertex(newnode, orig);
                    out_graph->setCurrentPosition(orig, newnode);
                    out_graph->new_edge_order(t_now, newnode, mincut);
                }

                if (deleted[lv][n].first != UNDEFINED_NODE) {
                    auto [t, orig] = deleted[lv][n];
                    NodeID t_now = out_graph->getCurrentPosition(t);

                    out_graph->addContainedVertex(t_now, orig);
                    out_graph->setCurrentPosition(orig, t_now);
                }
            }
        }

        LOGC(timing) << "t " << t.elapsed() << " cactus_unpack n "
                     << out_graph->n() << " m " << out_graph->m();
        return out_graph;
    }

 private:
    std::tuple<NodeID, EdgeID, NodeID> centralFlowEdge(
        std::shared_ptr<mutable_graph> G) {
        NodeID random_vtx = random_functions::nextInt(0, G->n() - 1);
        NodeID v1 = std::get<2>(graph_algorithms::bfsDistances(G, random_vtx));
        auto [parent, distance, v2] = graph_algorithms::bfsDistances(G, v1);

        uint32_t max_distance = distance[v2];
        for (uint32_t d = max_distance; d > (max_distance + 1) / 2; --d) {
            if (distance[v2] != d) {
                LOG1 << distance[v2] << " of " << v2 << " is not " << d;
                exit(1);
            }
            v2 = parent[v2];
        }

        for (EdgeID e : G->edges_of(v2)) {
            if (G->getEdgeTarget(v2, e) == parent[v2]) {
                return std::make_tuple(v2, e, parent[v2]);
            }
        }

        LOG1 << "Central flow edge didn't find an edge!";
        exit(1);
    }

    std::tuple<NodeID, EdgeID, NodeID> maximumFlowEdge(
        std::shared_ptr<mutable_graph> G) {
        NodeWeight max_degree = 0;
        NodeID s = UNDEFINED_NODE;

        for (NodeID n : G->nodes()) {
            if (G->getNodeDegree(n) > max_degree &&
                !G->isEmpty(n)) {
                max_degree = G->getNodeDegree(n);
                s = n;
            }
        }

        NodeID t = UNDEFINED_NODE;
        EdgeID e = UNDEFINED_EDGE;
        NodeWeight max_ngbr = 0;

        for (EdgeID edge : G->edges_of(s)) {
            NodeID ngbr = G->getEdgeTarget(s, edge);
            if (G->getNodeDegree(ngbr) > max_ngbr &&
                !G->isEmpty(ngbr)) {
                max_ngbr = G->getNodeDegree(ngbr);
                t = ngbr;
                e = edge;
            }
        }

        if (t == UNDEFINED_NODE) {
            LOG1 << "Heaviest vertex has only empty neighbours!";
            return findFlowEdge(G);
        } else {
            return std::make_tuple(s, e, t);
        }
    }

    std::tuple<NodeID, EdgeID, NodeID> findFlowEdge(
        std::shared_ptr<mutable_graph> G) {
        NodeID s = random_functions::nextInt(0, G->n() - 1);
        NodeID tgt = 0;
        NodeID max_edge = G->get_first_invalid_edge(s) - 1;
        EdgeID e = random_functions::nextInt(0, max_edge);
        bool edge_found = false;
        while (!edge_found) {
            while (G->isEmpty(s)) {
                if (s + 1 >= G->n()) {
                    s = 0;
                } else {
                    s++;
                }
            }
            e = 0;
            while (e < G->get_first_invalid_edge(s)
                   && (G->isEmpty(G->getEdgeTarget(s, e)))) {
                e++;
            }

            if (e < G->get_first_invalid_edge(s)) {
                edge_found = true;
            } else {
                if (s + 1 >= G->n()) {
                    s = 0;
                } else {
                    s++;
                }
            }
        }

        tgt = G->getEdgeTarget(s, e);
        return std::make_tuple(s, e, tgt);
    }

    std::vector<std::tuple<NodeID, std::vector<NodeID> > >
    removeHeavyEdges(std::shared_ptr<mutable_graph> G) {
        std::unordered_map<NodeID, std::vector<NodeID> > contract;
        std::vector<std::tuple<NodeID, std::vector<NodeID> > > cactusEdge;
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

    std::shared_ptr<mutable_graph> recursiveCactus(
        std::shared_ptr<mutable_graph> G, size_t depth) {
        if (depth % 100 == 0) {
            LOGC(configuration::getConfig()->verbose)
                << "G n " << G->n() << " m " << G->m() << " depth " << depth;
        }

        std::vector<std::tuple<NodeID, std::vector<NodeID> > > cactusEdge;

        if (depth % 10 == 0) {
            size_t previous = UNDEFINED_NODE;
            cactusEdge = removeHeavyEdges(G);

            // implicit do-while loop
            while (previous > G->n()) {
                previous = G->n();
                // create empty multicut problem to be able to run mod_capforest

                multicut_problem mcp(G);
                auto problem = std::make_shared<multicut_problem>(mcp);
                noi_minimum_cut noi;
                auto uf = noi.modified_capforest(problem, mincut + 1);

                if (uf.n() < G->n()) {
                    G = contraction::fromUnionFind(G, &uf);
                }

                auto uf12 = allCutsPrTests12(G, mincut);
                if (uf12.n() < G->n()) {
                    G = contraction::fromUnionFind(G, &uf12);
                }

                auto uf34 = allCutsPrTests34(G, mincut);
                if (uf34.n() < G->n()) {
                    G = contraction::fromUnionFind(G, &uf34);
                }
            }
        }

        if (G->number_of_nodes() == 1 || G->number_of_edges() == 0) {
            reInsertVertices(G, cactusEdge);
            return G;
        }

        VIECUT_ASSERT_TRUE(isCNCR(G));
        FlowType max_flow;

        // auto [s, e, tgt] = findFlowEdge(G);
        auto [s, e, tgt] = maximumFlowEdge(G);
        // auto [s, e, tgt] = centralFlowEdge(G);

        {
            std::vector<NodeID> vtcs = { s, tgt };
            push_relabel pr;
            max_flow = pr.solve_max_flow_min_cut(
                G, vtcs, 0, false, false).first;
        }

        if (max_flow > (FlowType)mincut) {
            LOG << "max flow is larger " << max_flow;

            VIECUT_ASSERT_EQ(G->getEdgeTarget(s, e), tgt);
            G->contractEdge(s, e);
            G = recursiveCactus(G, depth + 1);
            reInsertVertices(G, cactusEdge);
            VIECUT_ASSERT_TRUE(isCNCR(G));

            return G;
        } else {
            if (G->number_of_nodes() == 2) {
                reInsertVertices(G, cactusEdge);
                VIECUT_ASSERT_TRUE(isCNCR(G));

                return G;
            }
            // contract
            strongly_connected_components scc;
            auto [v, num_comp] = scc.strong_components(G);

            if (num_comp == 2
                && (G->getWeightedNodeDegree(s) == mincut
                    || G->getWeightedNodeDegree(tgt) == mincut)) {
                std::vector<int> empty;
                v.swap(empty);
                NodeID ctr = (G->getWeightedNodeDegree(s) == mincut) ? s : tgt;
                VIECUT_ASSERT_EQ(G->getWeightedNodeDegree(ctr), mincut);
                NodeID other = (ctr == s) ? tgt : s;

                auto elementsInCtr = G->containedVertices(ctr);
                auto elementsInOther = G->containedVertices(other);

                G->contractEdge(s, e);
                NodeID contracted_v = G->getCurrentPosition(elementsInCtr[0]);
                G->setContainedVertices(contracted_v, elementsInOther);
                for (NodeID n : elementsInOther) {
                    G->setCurrentPosition(n, contracted_v);
                }

                VIECUT_ASSERT_TRUE(isCNCR(G));

                auto ret = recursiveCactus(G, depth + 1);

                NodeID other_now = ret->getCurrentPosition(elementsInOther[0]);
                NodeID new_node = ret->new_empty_node();

                ret->new_edge(other_now, new_node, mincut);
                ret->setContainedVertices(new_node, elementsInCtr);

                for (NodeID n : elementsInCtr) {
                    ret->setCurrentPosition(n, new_node);
                }

                reInsertVertices(ret, cactusEdge);

                return ret;
            }

            auto STCactus = findSTCactus(v, G, s, num_comp);

            for (int c = 0; c < static_cast<int>(num_comp); ++c) {
                STCactus = mergeCactusWithComponent(STCactus, G, depth, c, v);
            }

            reInsertVertices(STCactus, cactusEdge);
            VIECUT_ASSERT_TRUE(isCNCR(STCactus));
            return STCactus;
        }
    }

    std::shared_ptr<mutable_graph> mergeCactusWithComponent(
        std::shared_ptr<mutable_graph> STCactus,
        std::shared_ptr<mutable_graph> G,
        size_t depth, int component,
        const std::vector<int>& scc_result) {
        // find a node in G that is contracted
        // and one that is not contracted,
        // use their location in contracted graphs
        // to re-find nodes as IDs swap around
        NodeID uncontracted_base_vertex = UNDEFINED_NODE;
        NodeID contracted_base_vertex = UNDEFINED_NODE;

        std::unordered_set<NodeID> all_ctr;
        for (size_t i = 0; i < scc_result.size(); ++i) {
            if (scc_result[i] != component) {
                all_ctr.insert(i);
                if (contracted_base_vertex == UNDEFINED_NODE) {
                    if (!G->containedVertices(i).empty()) {
                        contracted_base_vertex = G->containedVertices(i)[0];
                    }
                }
            } else {
                if (uncontracted_base_vertex == UNDEFINED_NODE) {
                    if (!G->containedVertices(i).empty()) {
                        uncontracted_base_vertex = G->containedVertices(i)[0];
                    }
                }
            }
        }

        auto graph = std::make_shared<mutable_graph>(*G);
        graph->contractVertexSet(all_ctr);
        auto n_i = recursiveCactus(graph, depth + 1);
        NodeID merge_vtx_in_cactus = STCactus->getCurrentPosition(
            uncontracted_base_vertex);
        NodeID nibar = n_i->getCurrentPosition(contracted_base_vertex);
        STCactus = mergeGraphs(STCactus, merge_vtx_in_cactus, n_i, nibar);
        VIECUT_ASSERT_TRUE(isCNCR(STCactus));
        return STCactus;
    }

    void reInsertVertices(std::shared_ptr<mutable_graph> G,
                          std::vector<
                              std::tuple<NodeID,
                                         std::vector<NodeID> > > toInsert) {
        for (const auto& [t, cont] : toInsert) {
            NodeID curr = G->getCurrentPosition(t);
            NodeID vtx = G->new_empty_node();
            G->new_edge_order(curr, vtx, mincut);
            G->setContainedVertices(vtx, cont);
            for (const auto& e : G->containedVertices(vtx)) {
                G->setCurrentPosition(e, vtx);
            }
        }
    }

    bool isCNCR(std::shared_ptr<mutable_graph> G) {
        for (NodeID n : G->nodes()) {
            if (G->isEmpty(n)) {
                EdgeWeight deg = 0;
                NodeID num = 0;
                for (EdgeID e : G->edges_of(n)) {
                    deg += G->getEdgeWeight(n, e);
                    num++;
                }

                if (deg % mincut != 0) {
                    LOG1 << "bad degree in vertex " << n;
                    return false;
                }

                if (deg / mincut == 3) {
                    LOG1 << "empty 3junction vertex " << n;
                    return false;
                }

                if (deg / mincut == 2) {
                    if (num < 4) {
                        LOG1 << "empty 2junction node with 2cycle " << n;
                        return false;
                    }
                }
            }
        }
        return true;
    }

    std::shared_ptr<mutable_graph>
    mergeGraphs(std::shared_ptr<mutable_graph> G1,
                NodeID v1,
                std::shared_ptr<mutable_graph> G2,
                NodeID v2) {
        NodeID v1_n = G1->number_of_nodes();

        if (G2->number_of_nodes() == 1) {
            return G1;
        }

        if (G1->number_of_nodes() == 1) {
            return G2;
        }

        std::vector<NodeID> empty;
        G1->setContainedVertices(v1, empty);

        for (size_t i = 0; i < G2->number_of_nodes(); ++i) {
            if (i == v2) {
                for (NodeID c : G2->containedVertices(i)) {
                    if (G1->getCurrentPosition(c) == v1) {
                        G1->setCurrentPosition(c, v1);
                        G1->addContainedVertex(v1, c);
                    }
                }
            } else {
                NodeID v = G1->new_empty_node();

                for (NodeID c : G2->containedVertices(i)) {
                    G1->setCurrentPosition(c, v);
                }

                G1->setContainedVertices(v, G2->containedVertices(i));
            }
        }

        for (NodeID n : G2->nodes()) {
            NodeID new_id = n + v1_n - (n > v2);
            if (n == v2) {
                new_id = v1;
            }

            for (EdgeID e : G2->edges_of(n)) {
                NodeID t = G2->getEdgeTarget(n, e);
                if (t != v2) {
                    NodeID t_id = t + v1_n - (t > v2);
                    G1->new_edge(new_id, t_id, G2->getEdgeWeight(n, e));
                }
            }
        }

        if (G1->isEmpty(v1)) {
            canonizeCactus(G1, v1);
        }

        return G1;
    }

    bool canonizeCactus(std::shared_ptr<mutable_graph> cactus, NodeID vertex) {
        size_t junctions = 0;
        std::vector<std::tuple<NodeID, EdgeID, EdgeID> > neighbors_tree;
        std::vector<std::pair<NodeID, EdgeID> > neighbors_cycle;
        for (EdgeID e : cactus->edges_of(vertex)) {
            EdgeWeight wgt = cactus->getEdgeWeight(vertex, e);
            if (wgt == mincut) {
                neighbors_tree.emplace_back(cactus->getEdgeTarget(vertex, e),
                                            cactus->getReverseEdge(vertex, e),
                                            e);
            } else {
                VIECUT_ASSERT_EQ(mincut % 2, 0);
                VIECUT_ASSERT_EQ(mincut / 2, wgt);
                neighbors_cycle.emplace_back(cactus->getEdgeTarget(vertex, e),
                                             cactus->getReverseEdge(vertex, e));
            }
            junctions += wgt;
        }

        if (junctions % mincut != 0) {
            LOG1 << mincut << " " << junctions;
            LOG1 << cactus;
            LOG1 << "ERROR! NOT EVEN CYCLED. EXITING!";
            exit(1);
        }

        junctions /= mincut;

        if (junctions == 2 && !neighbors_tree.empty()) {
            EdgeID ed = std::get<2>(neighbors_tree[0]);
            cactus->contractEdge(vertex, ed);
            return true;
        }

        if (junctions == 3) {
            std::vector<std::pair<std::pair<NodeID, EdgeID>,
                                  std::pair<NodeID, EdgeID> > > pairs;
            while (neighbors_cycle.size() > 2) {
                // find matching vertices in cycle, remove from array
                // (only up to 6 elements in vector, so O(n) delete is fine)
                std::map<NodeID, size_t> neighbors;
                for (size_t i = 1; i < neighbors_cycle.size(); ++i) {
                    neighbors.emplace(neighbors_cycle[i].first, i);
                }

                std::queue<NodeID> q;
                std::vector<bool> seen(cactus->number_of_nodes(), false);
                // if cactus is actually a cactus,
                // the cycles are only connected through 'vertex'
                seen[vertex] = true;
                q.push(neighbors_cycle[0].first);

                while (!q.empty()) {
                    NodeID n = q.front();
                    q.pop();
                    for (EdgeID e : cactus->edges_of(n)) {
                        NodeID tgt = cactus->getEdgeTarget(n, e);
                        if (neighbors.count(tgt)) {
                            pairs.emplace_back(neighbors_cycle[0],
                                               neighbors_cycle[neighbors[tgt]]);
                            std::queue<NodeID> empty;
                            std::swap(q, empty);
                            neighbors_cycle.erase(neighbors_cycle.begin()
                                                  + neighbors[tgt]);
                            neighbors_cycle.erase(neighbors_cycle.begin());
                            break;
                        }

                        if (!seen[tgt]) {
                            q.push(tgt);
                        }
                    }
                }
            }
            if (neighbors_cycle.size() == 2) {
                pairs.emplace_back(neighbors_cycle[0], neighbors_cycle[1]);
            }

            VIECUT_ASSERT_EQ(pairs.size() + neighbors_tree.size(), 3);

            // delete vertex and replace with 3 new empty 2-junction vertices
            // merge the ones that are in a 2-cycle

            std::vector<NodeID> new_nodes;

            for (auto p : pairs) {
                new_nodes.emplace_back(cactus->new_empty_node());
                cactus->new_edge_order(new_nodes.back(),
                                       p.first.first, mincut / 2);
                cactus->new_edge_order(new_nodes.back(),
                                       p.second.first, mincut / 2);
            }

            for (auto n : neighbors_tree) {
                new_nodes.emplace_back(std::get<0>(n));
            }

            cactus->new_edge_order(new_nodes[0], new_nodes[1], mincut / 2);
            cactus->new_edge_order(new_nodes[0], new_nodes[2], mincut / 2);
            cactus->new_edge_order(new_nodes[1], new_nodes[2], mincut / 2);

            cactus->deleteVertex(vertex);
            return true;
        }

        return false;
    }

    std::shared_ptr<mutable_graph> findSTCactus(
        const std::vector<int>& v,
        std::shared_ptr<mutable_graph> G,
        NodeID s,
        int num_comp) {
        auto contract = std::make_shared<mutable_graph>();
        contract->start_construction(num_comp);
        NodeID contained = G->containedVertices(s)[0];
        contract->setOriginalNodes(G->getOriginalNodes());

        for (NodeID n : contract->nodes()) {
            contract->setContainedVertices(n, { });
        }

        for (NodeID n = 0; n < G->getOriginalNodes(); ++n) {
            NodeID n_in_contract = v[G->getCurrentPosition(n)];
            contract->addContainedVertex(n_in_contract, n);
            contract->setCurrentPosition(n, n_in_contract);
        }

        for (NodeID n : contract->nodes()) {
            for (NodeID m : contract->nodes()) {
                if (n < m) {
                    contract->new_edge(n, m, 0);
                }
            }
        }

        for (NodeID n : G->nodes()) {
            for (EdgeID e : G->edges_of(n)) {
                NodeID t = G->getEdgeTarget(n, e);
                EdgeWeight wgt = G->getEdgeWeight(n, e);
                int ctr = v[t];
                if (v[n] > ctr) {
                    EdgeID e_ctr = ctr - (ctr > v[n]);
                    auto wgt_ctr = wgt + contract->getEdgeWeight(v[n], e_ctr);
                    contract->setEdgeWeight(v[n], e_ctr, wgt_ctr);
                }
            }
        }

        for (NodeID n : contract->nodes()) {
            for (EdgeID e = contract->getNodeDegree(n); e-- > 0; ) {
                if (contract->getEdgeWeight(n, e) == 0) {
                    contract->deleteEdge(n, e);
                }
            }
        }

        auto stcactus = std::make_shared<mutable_graph>();
        NodeID num_vertices = contract->n();
        stcactus->start_construction(num_vertices);
        stcactus->resizePositions(contract->getOriginalNodes());

        s = contract->getCurrentPosition(contained);

        node_bucket_pq pq(num_vertices, mincut + 1);

        for (NodeID n = 0; n < num_vertices; ++n) {
            stcactus->new_node();

            if (n != s)
                pq.insert(n, 0);
        }

        for (EdgeID e : contract->edges_of(s)) {
            NodeID tgt = contract->getEdgeTarget(s, e);
            pq.increaseKey(tgt, contract->getEdgeWeight(s, e));
        }

        std::vector<NodeID> node_mapping(contract->number_of_nodes());
        std::vector<NodeID> rev_node_mapping(contract->number_of_nodes());

        stcactus->setContainedVertices(0, contract->containedVertices(s));
        node_mapping[s] = 0;

        for (NodeID v : stcactus->containedVertices(0)) {
            stcactus->setCurrentPosition(v, 0);
        }

        for (NodeID n = 1; n < num_vertices; ++n) {
            NodeID next = pq.deleteMax();
            for (EdgeID e : contract->edges_of(next)) {
                NodeID tgt = contract->getEdgeTarget(next, e);
                EdgeWeight wgt = pq.getKey(tgt);
                if (pq.contains(tgt)) {
                    pq.increaseKey(
                        tgt,
                        std::min(wgt + contract->getEdgeWeight(next, e),
                                 mincut));
                }
            }

            node_mapping[next] = n;
            rev_node_mapping[n] = next;

            stcactus->setContainedVertices(
                n, contract->containedVertices(next));

            for (NodeID v : stcactus->containedVertices(n)) {
                stcactus->setCurrentPosition(v, n);
            }
        }

        std::vector<std::vector<NodeID> > A;
        std::vector<NodeID> B;
        std::vector<bool> order;

        size_t i = 1;
        B.emplace_back(0);
        order.emplace_back(false);

        while (i < (contract->number_of_nodes() - 1)) {
            EdgeWeight cycle_degree = 0;
            std::unordered_set<NodeID> curr_cycle;
            NodeID n = rev_node_mapping[i];

            while ((cycle_degree == 0 || cycle_degree == mincut)
                   && (i + 1 < contract->number_of_nodes())) {
                n = rev_node_mapping[i];
                for (EdgeID e : contract->edges_of(n)) {
                    NodeID tgt = contract->getEdgeTarget(n, e);
                    EdgeWeight wgt = contract->getEdgeWeight(n, e);
                    if (curr_cycle.count(tgt)) {
                        cycle_degree -= wgt;
                    } else {
                        cycle_degree += wgt;
                    }
                }

                if (cycle_degree == mincut) {
                    i++;
                    curr_cycle.insert(n);
                }
            }

            if (curr_cycle.size() > 0) {
                A.emplace_back();
                order.emplace_back(true);
                for (size_t v = i - curr_cycle.size(); v < i; ++v) {
                    A.back().emplace_back(v);
                }
            } else {
                i++;
                B.emplace_back(node_mapping[n]);
                order.emplace_back(false);
            }
        }

        order.emplace_back(false);
        B.emplace_back(num_vertices - 1);

        NodeID previous = 0;
        size_t a_index = 0, b_index = 0;
        VIECUT_ASSERT_EQ(order.size(), A.size() + B.size());
        for (size_t i = 0; i < (A.size() + B.size() - 1); ++i) {
            if (order[i]) {
                // make cycle
                for (size_t j = 0; j < A[a_index].size(); ++j) {
                    if (j > 0) {
                        stcactus->new_edge_order(A[a_index][j - 1],
                                                 A[a_index][j], mincut / 2);
                    } else {
                        stcactus->new_edge_order(previous,
                                                 A[a_index][0], mincut / 2);
                    }
                    if (j == A[a_index].size() - 1) {
                        // last vertex, connect with next cycle or ordered vtx
                        NodeID next;

                        if (order[i + 1] == true) {
                            next = stcactus->new_empty_node();
                        } else {
                            next = B[b_index];
                        }

                        stcactus->new_edge_order(A[a_index][j], next,
                                                 mincut / 2);
                        stcactus->new_edge_order(previous, next, mincut / 2);
                        previous = next;
                    }
                }
                a_index++;
            } else {
                // make ordered vtx
                if (!order[i + 1]) {
                    stcactus->new_edge_order(B[b_index],
                                             B[b_index + 1], mincut);
                }
                previous = B[b_index];
                b_index++;
            }
        }
        stcactus->finish_construction();

        return stcactus;
    }

    union_find allCutsPrTests12(std::shared_ptr<mutable_graph> G,
                                EdgeWeight limit) {
        union_find uf(G->number_of_nodes());
        std::vector<EdgeWeight> degrees;
        std::vector<bool> contracted(G->number_of_nodes(), false);

        for (NodeID n : G->nodes()) {
            degrees.push_back(G->getWeightedNodeDegree(n));
        }

        for (NodeID n : G->nodes()) {
            for (EdgeID e : G->edges_of(n)) {
                NodeID target = uf.Find(G->getEdgeTarget(n, e));
                NodeID source = uf.Find(n);
                EdgeWeight wgt = G->getEdgeWeight(n, e);

                if (target != source
                    && (wgt > limit
                        || 2 * wgt > degrees[source]
                        || 2 * wgt > degrees[target])
                    && degrees[source] > limit
                    && degrees[target] > limit) {
                    EdgeWeight new_w =
                        degrees[source] + degrees[target] - (2 * wgt);

                    degrees[source] = new_w;
                    degrees[target] = new_w;
                    uf.Union(source, target);
                }
            }
        }
        return uf;
    }

    union_find allCutsPrTests34(std::shared_ptr<mutable_graph> G,
                                EdgeWeight weight_limit) {
        union_find uf(G->number_of_nodes());
        std::vector<EdgeID> marked(G->number_of_nodes(), UNDEFINED_EDGE);
        std::vector<bool> finished(G->number_of_nodes(), false);
        std::vector<bool> contracted(G->number_of_nodes(), false);

        for (NodeID n : G->nodes()) {
            if (finished[n])
                continue;

            finished[n] = true;
            for (EdgeID e : G->edges_of(n)) {
                NodeID tgt = G->getEdgeTarget(n, e);
                if (tgt > n) {
                    marked[tgt] = e;
                }
            }

            EdgeWeight deg_n = G->getWeightedNodeDegree(n);
            for (EdgeID e1 : G->edges_of(n)) {
                NodeID tgt = G->getEdgeTarget(n, e1);
                EdgeWeight deg_tgt = G->getWeightedNodeDegree(tgt);
                if (finished[tgt]) {
                    marked[tgt] = UNDEFINED_EDGE;
                    continue;
                }

                finished[tgt] = true;
                EdgeWeight wgt_sum = G->getEdgeWeight(n, e1);
                if (tgt > n) {
                    for (EdgeID e2 : G->edges_of(tgt)) {
                        NodeID tgt2 = G->getEdgeTarget(tgt, e2);
                        if (marked[tgt2] == UNDEFINED_EDGE)
                            continue;

                        if (marked[tgt2] >= G->get_first_invalid_edge(n)
                            || marked[tgt2] < G->get_first_edge(n)) {
                            continue;
                        }

                        EdgeWeight w1 = G->getEdgeWeight(n, e1);
                        EdgeWeight w2 = G->getEdgeWeight(tgt, e2);
                        EdgeWeight w3 = G->getEdgeWeight(n, marked[tgt2]);

                        wgt_sum += std::min(w2, w3);

                        bool contractible =
                            2 * (w1 + w3) > deg_n && 2 * (w1 + w2) > deg_tgt
                            && deg_n > weight_limit && deg_tgt > weight_limit;

                        if (contractible
                            && !contracted[n] && !contracted[tgt]) {
                            uf.Union(n, tgt);
                            contracted[n] = true;
                            contracted[tgt] = true;
                        }
                    }

                    if (wgt_sum > weight_limit) {
                        uf.Union(n, tgt);
                        contracted[n] = true;
                        contracted[tgt] = true;
                    }
                    marked[tgt] = UNDEFINED_EDGE;
                }
            }
        }
        return uf;
    }

    timer t;
    EdgeWeight mincut;
};
