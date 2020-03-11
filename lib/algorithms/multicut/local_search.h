/******************************************************************************
 * local_search.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2018-2020 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include "common/configuration.h"
#include "common/definitions.h"

class local_search {
 public:

    static EdgeWeight flowBetweenBlocks(const mutable_graph& original_graph,
                                        const std::vector<bool>& fixed_vertex,
                                        std::vector<NodeID>* sol,
                                        NodeID terminal1,
                                        NodeID terminal2) {
        std::vector<NodeID>& solution = *sol;

        std::vector<NodeID> mapping(original_graph.n(), UNDEFINED_NODE);
        auto G = std::make_shared<mutable_graph>();
        EdgeWeight sol_weight = 0;

        NodeID id = 2;
        for (NodeID n : original_graph.nodes()) {
            if (solution[n] != terminal1 && solution[n] != terminal2)
                continue;
            if (fixed_vertex[n]) {
                mapping[n] = (solution[n] == terminal1 ? 0 : 1);
            } else {
                mapping[n] = id;
                id++;
            }
        }
        G->start_construction(id);
        std::unordered_map<NodeID, EdgeWeight> edgesToFixed0;
        std::unordered_map<NodeID, EdgeWeight> edgesToFixed1;
        for (NodeID n : original_graph.nodes()) {
            if (solution[n] != terminal1 && solution[n] != terminal2)
                continue;
            NodeID m_n = mapping[n];
            for (EdgeID e : original_graph.edges_of(n)) {
                auto [t, w] = original_graph.getEdge(n, e);
                NodeID m_t = mapping[t];
                if ((solution[t] != terminal1 && solution[t] != terminal2)
                    || m_n >= m_t || m_t < 2)
                    continue;

                if (solution[t] != solution[n]) {
                    sol_weight += w;
                }

                if (m_n < 2) {
                    if (m_n == 0) {
                        if (edgesToFixed0.count(m_t) > 0) {
                            edgesToFixed0[m_t] += w;
                        } else {
                            edgesToFixed0[m_t] = w;
                        }
                    } else {
                        if (edgesToFixed1.count(m_t) > 0) {
                            edgesToFixed1[m_t] += w;
                        } else {
                            edgesToFixed1[m_t] = w;
                        }
                    }
                } else {
                    G->new_edge_order(m_n, m_t, w);
                }
            }
        }
        for (auto [n, w] : edgesToFixed0) {
            G->new_edge_order(n, 0, w);
        }
        for (auto [n, w] : edgesToFixed1) {
            G->new_edge_order(n, 1, w);
        }

        std::vector<NodeID> terminals = { 0, 1 };
        push_relabel pr;
        auto [f, s] = pr.solve_max_flow_min_cut(G, terminals, 0, true);
        std::unordered_set<NodeID> zero;

        for (NodeID v : s) {
            zero.insert(v);
        }

        size_t improvement = (sol_weight - f);
        for (NodeID n : original_graph.nodes()) {
            if (solution[n] == terminal1 || solution[n] == terminal2) {
                if (fixed_vertex[n]) {
                    if (zero.count(mapping[n]) != (solution[n] == terminal1)) {
                        LOG1 << "DIFFERENT";
                        exit(1);
                    }
                }
                NodeID map = mapping[n];
                if (zero.count(map) > 0) {
                    solution[n] = terminal1;
                } else {
                    solution[n] = terminal2;
                }
            }
        }
        return improvement;
    }

    static EdgeWeight flowLocalSearch(const mutable_graph& original_graph,
                                      const std::vector<NodeID>& original_terms,
                                      const std::vector<bool>& fixed_vertex,
                                      std::vector<NodeID>* sol) {
        std::vector<NodeID>& solution = *sol;
        std::vector<std::vector<EdgeWeight> > blockConnectivity(
            original_terms.size());

        EdgeWeight improvement = 0;

        std::vector<std::pair<NodeID, NodeID> > neighboringBlocks;

        for (auto& b : blockConnectivity) {
            b.resize(original_terms.size(), 0);
        }

        for (NodeID n : original_graph.nodes()) {
            NodeID blockn = solution[n];
            for (EdgeID e : original_graph.edges_of(n)) {
                auto [t, w] = original_graph.getEdge(n, e);
                if (solution[t] > blockn) {
                    blockConnectivity[blockn][solution[t]] += w;
                }
            }
        }

        for (size_t i = 0; i < blockConnectivity.size(); ++i) {
            for (size_t j = 0; j < blockConnectivity[i].size(); ++j) {
                // TODO (anoe): find lower limit for when flow is useful.
                if (blockConnectivity[i][j] > 0) {
                    neighboringBlocks.emplace_back(i, j);
                }
            }
        }

        random_functions::permutate_vector_good(&neighboringBlocks);

        for (auto [a, b] : neighboringBlocks) {
            improvement += flowBetweenBlocks(
                original_graph, fixed_vertex, sol, a, b);
        }

        return improvement;
    }

    static EdgeWeight improveSolution(const mutable_graph& original_graph,
                                      const std::vector<NodeID>& original_terms,
                                      const std::vector<bool>& fixed_vertex,
                                      std::vector<NodeID>* sol,
                                      EdgeWeight valueBefore) {
        std::vector<NodeID>& current_solution = *sol;
        EdgeWeight ls_bound = valueBefore;

        bool change_found = true;

        while (change_found) {
            change_found = false;

            FlowType imp = flowLocalSearch(original_graph, original_terms,
                                           fixed_vertex, sol);
            ls_bound -= imp;
            if (imp > 0) {
                change_found = true;
            }
            std::vector<NodeID> permute(original_graph.n(), 0);
            std::vector<bool> inBoundary(original_graph.n(), true);
            std::vector<std::pair<NodeID, int64_t> > nextBest(
                original_graph.n(), { UNDEFINED_NODE, 0 });

            random_functions::permutate_vector_good(&permute, true);

            for (NodeID v : original_graph.nodes()) {
                NodeID n = permute[v];
                if (fixed_vertex[n] || !inBoundary[n])
                    continue;

                std::vector<EdgeWeight> blockwgt(
                    configuration::getConfig()->num_terminals, 0);
                NodeID ownBlockID = current_solution[n];
                for (EdgeID e : original_graph.edges_of(n)) {
                    auto [t, w] = original_graph.getEdge(n, e);
                    NodeID block = current_solution[t];
                    blockwgt[block] += w;
                }

                EdgeWeight ownBlockWgt = blockwgt[ownBlockID];
                NodeID maxBlockID = 0;
                EdgeWeight maxBlockWgt = 0;
                for (size_t i = 0; i < blockwgt.size(); ++i) {
                    if (i != ownBlockID) {
                        if (blockwgt[i] > maxBlockWgt) {
                            maxBlockID = i;
                            maxBlockWgt = blockwgt[i];
                        }
                    }
                }

                if (maxBlockWgt) {
                    inBoundary[n] = false;
                }

                int64_t gain = static_cast<int64_t>(maxBlockWgt)
                               - static_cast<int64_t>(ownBlockWgt);

                bool notDoublemoved = true;
                for (EdgeID e : original_graph.edges_of(n)) {
                    auto [t, w] = original_graph.getEdge(n, e);
                    auto [nbrBlockID, nbrGain] = nextBest[t];
                    int64_t movegain = nbrGain + gain + 2 * w;
                    if (current_solution[t] == current_solution[n] &&
                        nbrBlockID == maxBlockID && movegain > 0
                        && movegain > gain) {
                        current_solution[n] = maxBlockID;
                        current_solution[t] = maxBlockID;
                        ls_bound -= movegain;
                        if (movegain > 0) {
                            change_found = true;
                        }
                        notDoublemoved = false;

                        for (EdgeID e : original_graph.edges_of(n)) {
                            NodeID b = original_graph.getEdgeTarget(n, e);
                            nextBest[b] = std::make_pair(UNDEFINED_NODE, 0);
                            inBoundary[b] = true;
                        }
                        nextBest[t] = std::make_pair(UNDEFINED_NODE, 0);
                        for (EdgeID e : original_graph.edges_of(t)) {
                            NodeID b = original_graph.getEdgeTarget(t, e);
                            nextBest[b] = std::make_pair(UNDEFINED_NODE, 0);
                            inBoundary[b] = true;
                        }
                    }
                }

                if (!notDoublemoved) {
                    continue;
                }

                if (gain >= 0) {
                    current_solution[n] = maxBlockID;
                    ls_bound -= gain;
                    for (EdgeID e : original_graph.edges_of(n)) {
                        NodeID t = original_graph.getEdgeTarget(n, e);
                        nextBest[t] = std::make_pair(UNDEFINED_NODE, 0);
                        inBoundary[t] = true;
                    }
                } else {
                    nextBest[n] = std::make_pair(maxBlockID, gain);
                }
            }
        }
        return ls_bound;
    }

 private:
};
