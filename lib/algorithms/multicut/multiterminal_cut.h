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

#include <memory>
#include <queue>
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
            if (config->top_k > 0) {
                auto v = graph_algorithms::top_k_degrees(G, config->top_k);
                for (auto vtx : v) {
                    terminals.push_back(vtx);
                }
            } else {
                for (auto term : config->term_strings) {
                    try {
                        NodeID terminal = std::stoi(term);
                        if (terminal < G->number_of_nodes()) {
                            terminals.emplace_back(std::stoi(term));
                        } else {
                            LOG1 << term
                                 << " is not a valid terminal! Continuing without.";
                        }
                    } catch (...) {
                        LOG1 << term
                             << " is not a valid terminal! Continuing without.";
                    }
                }

                if (config->random_k > 0) {
                    for (int i = 0; i < config->random_k; ++i) {
                        terminals.emplace_back(
                            random_functions::nextInt(0,
                                                      G->number_of_nodes() - 1));
                        LOG << "Set random terminal " << terminals.back();
                    }
                }
            }

            config->total_terminals = terminals.size();

            if (config->preset_percentage > 0) {
                NodeID blocksize = G->n() / terminals.size();
                config->bfs_size = blocksize * config->preset_percentage / 100;
            }
        } else {
            auto s = tlx::split(".", config->partition_file);
            if (s.size() < 2) {
                LOG1 << "partition filename invalid";
                exit(1);
            }

            size_t maxn = G->n();
            if (s.back() == "pre") {
                NodeID num_partitions = std::stoi(s[s.size() - 2]);
                config->total_terminals = num_partitions;

                strongly_connected_components cc;
                auto [components, num_comp, unused] = cc.strong_components(G);
                (void)unused;
                std::vector<NodeID> v =
                    graph_io::readVector<NodeID>(config->partition_file);
                std::vector<NodeID> term;

                for (int c = 0; c < static_cast<int>(num_comp); ++c) {
                    for (size_t i = 0; i < num_partitions; ++i) {
                        bool terminal_is_set = false;
                        size_t size = 0;
                        std::unordered_set<NodeID> vset;
                        for (size_t n = 0; n < maxn; ++n) {
                            if (v[n] == i && components[n] == c) {
                                size++;
                                if (!terminal_is_set) {
                                    terminal_is_set = true;
                                    term.emplace_back(n);
                                }
                                vset.emplace(G->getCurrentPosition(n));
                            }
                        }

                        if (vset.size() > 1) {
                            G->contractVertexSet(vset);
                        }
                    }
                }

                for (size_t i = 0; i < term.size(); ++i) {
                    terminals.emplace_back(G->getCurrentPosition(term[i]));
                }

                config->bfs_size = 1;
            } else {
                std::vector<NodeID> v =
                    graph_io::readVector<NodeID>(config->partition_file);
                std::vector<NodeID> order =
                    graph_io::readVector<NodeID>(config->partition_file + ".pos");
                size_t num_partitions = std::stoi(s.back());
                config->total_terminals = num_partitions;
                std::vector<NodeID> term;
                size_t maxn = G->n();

                if (config->preset_percentage > 0) {
                    NodeID blocksize = G->n() / num_partitions;
                    config->bfs_size = blocksize * config->preset_percentage / 100;
                }

                for (size_t i = 0; i < num_partitions; ++i) {
                    size_t size = 0;
                    std::unordered_set<NodeID> vset;
                    for (size_t n = 0; n < maxn; ++n) {
                        if (v[n] == i && order[n] <= config->bfs_size) {
                            size++;
                            if (term.size() == i) {
                                term.emplace_back(n);
                            }
                            vset.emplace(G->getCurrentPosition(n));
                        }
                    }

                    if (vset.size() > 1) {
                        G->contractVertexSet(vset);
                    }
                }

                for (size_t i = 0; i < num_partitions; ++i) {
                    terminals.emplace_back(G->getCurrentPosition(term[i]));
                }
                config->bfs_size = 1;
            }
        }
        return terminals;
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
};
