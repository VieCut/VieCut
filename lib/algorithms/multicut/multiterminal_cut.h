/******************************************************************************
 * multiterminal_cut.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <mpi.h>

#include <algorithm>
#include <memory>
#include <queue>
#include <unordered_set>
#include <vector>

#include "algorithms/misc/strongly_connected_components.h"
#include "algorithms/multicut/branch_multicut.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "tlx/string/split.hpp"

class multiterminal_cut {
 public:
    static constexpr bool debug = false;
    multiterminal_cut() { }

    std::vector<NodeID> setOriginalTerminals(std::shared_ptr<mutable_graph> G) {
        auto config = configuration::getConfig();
        std::vector<NodeID> terminals;
        if (config->partition_file.empty()) {
            if (config->top_k > 0)
                terminals = topKTerminals(G);
            if (config->random_k > 0)
                terminals = randomTerminals(G);
            if (terminals.empty())
                terminals = terminalsByID(G);

            if (config->preset_percentage > 0) {
                NodeID blocksize = G->n() / terminals.size();
                config->bfs_size =
                    std::max(
                        static_cast<size_t>(
                            blocksize * config->preset_percentage / 100),
                        static_cast<size_t>(1));
            }
            return terminals;
        } else {
            auto s = tlx::split(".", config->partition_file);
            config->bfs_size = 1;
            if (s.size() < 2) {
                LOG1 << "partition filename invalid";
                exit(1);
            }

            if (s.back() == "pre") {
                config->total_terminals = std::stoi(s[s.size() - 2]);
                return presetFileTerminals(G);
            } else {
                config->total_terminals = std::stoi(s.back());
                return orderFileTerminals(G);
            }
        }
    }

    size_t multicut(std::shared_ptr<mutable_graph> G,
                    std::vector<NodeID> terminals) {
        strongly_connected_components cc;
        auto cfg = configuration::getConfig();
        auto problems = splitConnectedComponents(G, terminals);
        FlowType flow_sum = 0;
        for (auto& problem : problems) {
            if (debug) {
                graph_algorithms::checkGraphValidity(problem.graph);
            }

            std::vector<NodeID> terminals;
            for (size_t i = 0; i < problem.terminals.size(); ++i) {
                terminals.emplace_back(problem.terminals[i].position);
            }

            branch_multicut bmc(problem.graph, terminals);
            auto problem_pointer = std::make_shared<multicut_problem>(problem);
            addSurroundingAreaToTerminals(problem_pointer, terminals);
            flow_sum += bmc.find_multiterminal_cut(problem_pointer);
        }
        return flow_sum;
    }

 private:
    static void addSurroundingAreaToTerminals(
        std::shared_ptr<multicut_problem> mcp,
        std::vector<NodeID> terminals) {
        auto config = configuration::getConfig();
        if (config->bfs_size) {
            std::vector<bool> already_in_block(
                mcp->graph->number_of_nodes(), false);
            std::vector<std::vector<NodeID> > blocks;

            for (NodeID t : terminals) {
                already_in_block[t] = true;
            }

            for (NodeID t : terminals) {
                blocks.emplace_back();
                size_t current_block_size = 1;
                std::queue<NodeID> bfs_q;
                bfs_q.emplace(t);
                blocks.back().emplace_back(t);

                while (!bfs_q.empty()
                       && current_block_size < config->bfs_size) {
                    NodeID n = bfs_q.front();
                    bfs_q.pop();
                    for (EdgeID e : mcp->graph->edges_of(n)) {
                        NodeID tgt = mcp->graph->getEdgeTarget(n, e);
                        if (!already_in_block[tgt]) {
                            already_in_block[tgt] = true;
                            bfs_q.push(tgt);
                            blocks.back().emplace_back(tgt);
                            if (++current_block_size >= config->bfs_size) {
                                break;
                            }
                        }
                    }
                }
            }

            graph_contraction::contractIsolatingBlocks(mcp, blocks);

            for (NodeID n : mcp->graph->nodes()) {
                mcp->graph->setPartitionIndex(n, 0);
            }
            for (size_t l = 0; l < terminals.size(); ++l) {
                NodeID v = mcp->graph->getCurrentPosition(terminals[l]);
                mcp->graph->setPartitionIndex(v, l);
            }
        }
    }

    static std::vector<multicut_problem> splitConnectedComponents(
        std::shared_ptr<mutable_graph> G, std::vector<NodeID> all_terminals) {
        std::vector<multicut_problem> problems;
        std::vector<int> t_comp;
        strongly_connected_components cc;

        auto [components, num_comp, unused] = cc.strong_components(G);
        (void)unused;
        std::vector<std::vector<terminal> > terminals(num_comp);
        std::vector<NodeID> ctr(num_comp, 0);
        std::vector<NodeID> num_terminals(num_comp, 0);

        // space to create all connected components.
        // we only actually create the components
        // that contain at least one terminal
        std::vector<std::shared_ptr<mutable_graph> >
        component_subgraph(num_comp);

        for (NodeID t : all_terminals) {
            num_terminals[components[t]]++;
        }

        for (NodeID t : all_terminals) {
            int terminal_component = components[t];

            if (num_terminals[terminal_component] == 1)
                continue;

            if (!vector::contains(t_comp, terminal_component)) {
                graph_extractor ge;

                auto [G_out, reverse_mapping] =
                    ge.extract_block(G, terminal_component, components);

                component_subgraph[terminal_component] = G_out;

                for (NodeID i : all_terminals) {
                    if (components[i] == terminal_component) {
                        terminals[terminal_component].emplace_back(
                            reverse_mapping[i], ctr[terminal_component]++);
                    }
                }
            }
            t_comp.emplace_back(terminal_component);
        }

        for (size_t i = 0; i < num_comp; ++i) {
            if (terminals[i].size() > 1) {
                problems.emplace_back(component_subgraph[i], terminals[i]);
            }
        }

        return problems;
    }

    std::vector<NodeID> topKTerminals(std::shared_ptr<mutable_graph> G) {
        auto v = graph_algorithms::top_k_degrees(
            G, configuration::getConfig()->top_k);
        std::vector<NodeID> terminals;
        for (auto vtx : v) {
            terminals.push_back(vtx);
        }
        return terminals;
    }

    std::vector<NodeID> randomTerminals(std::shared_ptr<mutable_graph> G) {
        std::vector<NodeID> terminals;
        for (int i = 0; i < configuration::getConfig()->random_k; ++i) {
            terminals.emplace_back(random_functions::nextInt(0, G->n() - 1));
            LOG << "Set random terminal " << terminals.back();
        }
        return terminals;
    }

    std::vector<NodeID> terminalsByID(std::shared_ptr<mutable_graph> G) {
        std::vector<NodeID> terminals;
        for (auto term : configuration::getConfig()->term_strings) {
            try {
                NodeID terminal = std::stoi(term);
                if (terminal < G->number_of_nodes()) {
                    terminals.emplace_back(std::stoi(term));
                } else {
                    LOG1 << term << " >= " << G->n();
                }
            } catch (...) {
                LOG1 << term << " is not a valid terminal! Continuing without.";
            }
        }
        return terminals;
    }

    std::vector<NodeID> presetFileTerminals(std::shared_ptr<mutable_graph> G) {
        strongly_connected_components cc;
        auto config = configuration::getConfig();
        auto [components, num_comp, unused] = cc.strong_components(G);
        (void)unused;
        auto v = graph_io::readVector<NodeID>(config->partition_file);
        std::vector<NodeID> term;
        std::vector<NodeID> terminals;
        size_t n_start = G->n();

        for (int c = 0; c < static_cast<int>(num_comp); ++c) {
            for (size_t i = 0; i < config->total_terminals; ++i) {
                std::unordered_set<NodeID> contractSet;
                for (size_t n = 0; n < n_start; ++n) {
                    if (v[n] == i && components[n] == c) {
                        if (contractSet.size() == 0) {
                            term.emplace_back(n);
                        }
                        contractSet.emplace(G->getCurrentPosition(n));
                    }
                }
                if (contractSet.size() > 1) {
                    G->contractVertexSet(contractSet);
                }
            }
        }

        for (size_t i = 0; i < term.size(); ++i) {
            terminals.emplace_back(G->getCurrentPosition(term[i]));
        }

        config->bfs_size = 1;
        return terminals;
    }

    std::vector<NodeID> orderFileTerminals(std::shared_ptr<mutable_graph> G) {
        auto config = configuration::getConfig();
        auto v = graph_io::readVector<NodeID>(config->partition_file);
        auto o = graph_io::readVector<NodeID>(config->partition_file + ".pos");
        std::vector<NodeID> term;
        std::vector<NodeID> terminals;

        NodeID terminal_size = 1;
        NodeID start_n = G->n();
        if (config->preset_percentage > 0) {
            NodeID blocksize = G->n() / config->total_terminals;
            terminal_size =
                std::max(static_cast<size_t>(blocksize
                                             * config->preset_percentage / 100),
                         static_cast<size_t>(1));
        }

        for (size_t i = 0; i < config->total_terminals; ++i) {
            std::unordered_set<NodeID> contractSet;
            for (size_t n = 0; n < start_n; ++n) {
                if (v[n] == i && o[n] <= terminal_size) {
                    if (term.size() == i) {
                        term.emplace_back(n);
                    }
                    contractSet.emplace(G->getCurrentPosition(n));
                }
            }

            if (contractSet.size() > 1) {
                G->contractVertexSet(contractSet);
            }
        }

        for (size_t i = 0; i < config->total_terminals; ++i) {
            terminals.emplace_back(G->getCurrentPosition(term[i]));
        }
        return terminals;
    }
};
