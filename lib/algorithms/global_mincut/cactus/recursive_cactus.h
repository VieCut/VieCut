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
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/flow/push_relabel.h"
#include "algorithms/global_mincut/cactus/all_cut_local_red.h"
#include "algorithms/global_mincut/cactus/graph_modification.h"
#include "algorithms/global_mincut/cactus/heavy_edges.h"
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
        auto in_graph = mutable_graph::from_graph_access(graphs.back());
        auto out_graph = recursiveCactus(in_graph, 0);
        VIECUT_ASSERT_TRUE(graph_modification::isCNCR(out_graph, mincut));
        LOGC(timing) << "t " << t.elapsed() << " cactus n "
                     << out_graph->n() << " m " << out_graph->m();
        return out_graph;
    }

 private:
    std::shared_ptr<mutable_graph> recursiveCactus(
        std::shared_ptr<mutable_graph> G, size_t depth) {
        switch (configuration::getConfig()->optimization) {
        case 5:
        {
            heavy_edges he(mincut);
            auto cactusEdges = he.removeHeavyEdges(G);
            G = internalRecursiveCactus(G, depth);
            he.reInsertVertices(G, cactusEdges);
            return G;
        }
        case 6:
        {
            heavy_edges he(mincut);
            auto cactusEdges = he.removeHeavyEdges(G);
            auto cycleEdges = he.contractCycleEdges(G);
            G = internalRecursiveCactus(G, depth);
            he.reInsertCycles(G, cycleEdges);
            he.reInsertVertices(G, cactusEdges);
            return G;
        }
        default:
            return internalRecursiveCactus(G, depth);
        }
    }

    std::shared_ptr<mutable_graph> internalRecursiveCactus(
        std::shared_ptr<mutable_graph> G, size_t depth) {
        if (t.elapsed() > 1800) {
            // TODO(anoe): this is only temporary!
            LOG1 << "TIMEOUT!";
            return G;
        }

        if (depth % 100 == 0) {
            LOGC(configuration::getConfig()->verbose)
                << "G n " << G->n() << " m " << G->m() << " depth " << depth;
        }
        if (configuration::getConfig()->optimization >= 4) {
            if (depth % 10 == 0) {
                size_t previous = UNDEFINED_NODE;
                // implicit do-while loop
                while (previous > G->n()) {
                    previous = G->n();
                    // create empty multicut problem to run modified_capforest
                    multicut_problem mcp(G);
                    auto problem = std::make_shared<multicut_problem>(mcp);
                    noi_minimum_cut noi;
                    auto uf = noi.modified_capforest(problem, mincut + 1);
                    G = contraction::fromUnionFind(G, &uf);
                    auto uf12 = all_cut_local_red::allCutsPrTests12(G, mincut);
                    G = contraction::fromUnionFind(G, &uf12);
                    auto uf34 = all_cut_local_red::allCutsPrTests34(G, mincut);
                    G = contraction::fromUnionFind(G, &uf34);
                }
            }
        }
        if (G->number_of_nodes() == 1 || G->number_of_edges() == 0) {
            return G;
        }
        VIECUT_ASSERT_TRUE(graph_modification::isCNCR(G, mincut));
        FlowType max_flow;
        NodeID s = 0;
        NodeID tgt = 0;
        EdgeID e = 0;

        std::string es = configuration::getConfig()->edge_selection;

        if (es == "heavy")
            std::tie(s, e, tgt) = maximumFlowEdge(G);
        if (es == "heavy_weighted" || es == "" || es == "heavy_vertex")
            std::tie(s, e, tgt) = maximumWeightedFlowEdge(G);
        if (es == "central")
            std::tie(s, e, tgt) = centralFlowEdge(G);
        if (es == "random")
            std::tie(s, e, tgt) = findFlowEdge(G);

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
            return G;
        } else {
            if (G->number_of_nodes() == 2) {
                return G;
            }
            strongly_connected_components scc;
            auto [v, num_comp, blocksizes] = scc.strong_components(G);
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

                VIECUT_ASSERT_TRUE(graph_modification::isCNCR(G, mincut));
                auto ret = recursiveCactus(G, depth + 1);
                NodeID other_now = ret->getCurrentPosition(elementsInOther[0]);
                NodeID new_node = ret->new_empty_node();
                ret->new_edge(other_now, new_node, mincut);
                ret->setContainedVertices(new_node, elementsInCtr);
                for (NodeID n : elementsInCtr) {
                    ret->setCurrentPosition(n, new_node);
                }
                return ret;
            }

            auto STCactus = findSTCactus(v, G, s, num_comp);
            double g_n = static_cast<double>(G->n());
            // first the small blocks, last the big one. then we don't need to
            // copy graphs as the small ones are newly generated
            for (int c = 0; c < static_cast<int>(num_comp); ++c) {
                if (static_cast<double>(blocksizes[c]) <= (g_n / 2.0)) {
                    STCactus = mergeCactusWithComponent(
                        STCactus, G, depth, c, v, blocksizes[c]);
                }
            }
            for (int c = 0; c < static_cast<int>(num_comp); ++c) {
                if (static_cast<double>(blocksizes[c]) > (g_n / 2.0)) {
                    STCactus = mergeCactusWithComponent(
                        STCactus, G, depth, c, v, blocksizes[c]);
                }
            }
            VIECUT_ASSERT_TRUE(graph_modification::isCNCR(STCactus, mincut));
            return STCactus;
        }
    }

    std::shared_ptr<mutable_graph> mergeCactusWithComponent(
        std::shared_ptr<mutable_graph> STCactus,
        std::shared_ptr<mutable_graph> G,
        size_t depth, int component,
        const std::vector<int>& scc_result, size_t blocksize) {
        NodeID uncontracted_base_vertex = UNDEFINED_NODE;
        NodeID contracted_base_vertex = UNDEFINED_NODE;
        std::shared_ptr<mutable_graph> graph;
        if (static_cast<double>(blocksize) <=
            (static_cast<double>(G->n()) / 2.0)) {
            graph = std::make_shared<mutable_graph>();
            graph->start_construction(blocksize + 1);
            graph->setOriginalNodes(G->getOriginalNodes());
            graph->new_empty_node();
            for (NodeID n : graph->nodes()) {
                graph->setContainedVertices(n, { });
            }

            std::vector<NodeID> contained;
            NodeID vtx = 0;

            for (NodeID n : G->nodes()) {
                if (scc_result[n] == component) {
                    graph->new_empty_node();
                    contained.emplace_back(vtx++);
                    if (uncontracted_base_vertex == UNDEFINED_NODE &&
                        !G->containedVertices(n).empty()) {
                        uncontracted_base_vertex = G->containedVertices(n)[0];
                    }
                    for (NodeID con : G->containedVertices(n)) {
                        graph->addContainedVertex(contained.back(), con);
                        graph->setCurrentPosition(con, contained.back());
                    }
                } else {
                    contained.emplace_back(blocksize);
                    if (contracted_base_vertex == UNDEFINED_NODE &&
                        !G->containedVertices(n).empty()) {
                        contracted_base_vertex = G->containedVertices(n)[0];
                    }
                    for (NodeID con : G->containedVertices(n)) {
                        graph->addContainedVertex(blocksize, con);
                        graph->setCurrentPosition(con, blocksize);
                    }
                }
            }

            for (NodeID n : G->nodes()) {
                if (contained[n] != blocksize) {
                    EdgeWeight to_contracted = 0;
                    for (EdgeID e : G->edges_of(n)) {
                        NodeID t = G->getEdgeTarget(n, e);
                        EdgeWeight wgt = G->getEdgeWeight(n, e);
                        if (contained[t] == blocksize) {
                            to_contracted += wgt;
                        } else if (contained[n] < contained[t]) {
                            graph->new_edge(contained[n], contained[t], wgt);
                        }
                    }

                    if (to_contracted > 0) {
                        graph->new_edge(contained[n], blocksize, to_contracted);
                    }
                }
            }
        } else {
            // find a node in G that is contracted
            // and one that is not contracted,
            // use their location in contracted graphs
            // to re-find nodes as IDs swap around
            std::unordered_set<NodeID> all_ctr;
            for (size_t i = 0; i < scc_result.size(); ++i) {
                if (scc_result[i] != component) {
                    all_ctr.insert(i);
                    if (contracted_base_vertex == UNDEFINED_NODE &&
                        !G->containedVertices(i).empty()) {
                        contracted_base_vertex = G->containedVertices(i)[0];
                    }
                } else {
                    if (uncontracted_base_vertex == UNDEFINED_NODE &&
                        !G->containedVertices(i).empty()) {
                        uncontracted_base_vertex = G->containedVertices(i)[0];
                    }
                }
            }
            graph = G;
            graph->contractVertexSet(all_ctr);
        }
        auto n_i = recursiveCactus(graph, depth + 1);
        NodeID merge_vtx_in_cactus = STCactus->getCurrentPosition(
            uncontracted_base_vertex);
        NodeID nibar = n_i->getCurrentPosition(contracted_base_vertex);
        STCactus = graph_modification::mergeGraphs(
            STCactus, merge_vtx_in_cactus, n_i, nibar, mincut);
        VIECUT_ASSERT_TRUE(graph_modification::isCNCR(STCactus, mincut));
        return STCactus;
    }

    std::shared_ptr<mutable_graph> findSTCactus(
        const std::vector<int>& v, std::shared_ptr<mutable_graph> G,
        NodeID s, int num_comp) {
        auto contract = std::make_shared<mutable_graph>();
        contract->start_construction(num_comp);
        NodeID contained = G->containedVertices(s)[0];
        contract->setOriginalNodes(G->getOriginalNodes());
        for (NodeID n : contract->nodes()) {
            contract->setContainedVertices(n, { });
        }
        for (NodeID n = 0; n < G->getOriginalNodes(); ++n) {
            NodeID pos = G->getCurrentPosition(n);
            if (pos < G->n()) {
                NodeID n_in_contract = v[G->getCurrentPosition(n)];
                contract->addContainedVertex(n_in_contract, n);
                contract->setCurrentPosition(n, n_in_contract);
            }
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
                    EdgeWeight new_wgt = wgt + contract->getEdgeWeight(next, e);
                    pq.increaseKey(tgt, std::min(new_wgt, mincut));
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

    std::tuple<NodeID, EdgeID, NodeID> maximumWeightedFlowEdge(
        std::shared_ptr<mutable_graph> G) {
        NodeWeight max_degree = 0;
        NodeID s = UNDEFINED_NODE;

        for (NodeID n : G->nodes()) {
            if (G->getWeightedNodeDegree(n) > max_degree &&
                !G->isEmpty(n)) {
                max_degree = G->getWeightedNodeDegree(n);
                s = n;
            }
        }

        NodeID t = UNDEFINED_NODE;
        EdgeID e = UNDEFINED_EDGE;
        NodeWeight max_ngbr = 0;
        for (EdgeID edge : G->edges_of(s)) {
            NodeID ngbr = G->getEdgeTarget(s, edge);
            if (G->getWeightedNodeDegree(ngbr) > max_ngbr &&
                !G->isEmpty(ngbr)) {
                max_ngbr = G->getWeightedNodeDegree(ngbr);
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

    timer t;
    EdgeWeight mincut;
};
