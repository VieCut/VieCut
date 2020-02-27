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
#include <tuple>
#include <unordered_set>
#include <vector>

#include "algorithms/misc/strongly_connected_components.h"
#include "algorithms/multicut/branch_multicut.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "io/graph_io.h"
#include "tlx/string/split.hpp"

class multiterminal_cut {
 public:
    static constexpr bool debug = false;
    multiterminal_cut() { }

    std::vector<NodeID> setOriginalTerminals(std::shared_ptr<mutable_graph> G) {
        auto config = configuration::getConfig();
        std::vector<NodeID> terminalMapping;
        if (config->partition_file.empty()) {
            std::vector<NodeID> terminals;
            if (config->top_k > 0)
                terminals = topKTerminals(G);
            if (config->random_k > 0)
                terminals = randomTerminals(G);
            if (config->distant_terminals > 0)
                terminals = distantTerminals(G);
            if (terminals.empty())
                terminals = terminalsByID(G);

            config->num_terminals = terminals.size();

            if (config->preset_percentage > 0) {
                NodeID blocksize = G->n() / terminals.size();
                size_t desired = blocksize * config->preset_percentage / 100;
                size_t minimum_size = 1;
                config->bfs_size = std::max(desired, minimum_size);
            }

            terminalMapping = addSurroundingAreaToTerminals(G, terminals);
            return terminalMapping;
        } else {
            auto s = tlx::split(".", config->partition_file);
            config->bfs_size = 1;
            if (s.size() < 2) {
                LOG1 << "partition filename invalid";
                exit(1);
            }

            if (s.back() == "pre") {
                try {
                    config->num_terminals = std::stoi(s[s.size() - 2]);
                } catch (const std::exception&) {
                    LOG1 << "ERROR: " << s[s.size() - 2]
                         << " is not a valid number of terminals. Exiting!";
                    exit(1);
                }

                return presetFileTerminals(G);
            } else {
                try {
                    config->num_terminals = std::stoi(s.back());
                } catch (const std::exception&) {
                    LOG1 << "ERROR: " << s.back()
                         << " is not a valid number of terminals. Exiting!";
                    exit(1);
                }
                return orderFileTerminals(G);
            }
        }
    }

    size_t multicut(std::shared_ptr<mutable_graph> G,
                    std::vector<NodeID> terminals, NodeID num_terminals) {
        auto cfg = configuration::getConfig();
        cfg->num_terminals = num_terminals;

        // todo: use orig_graph here so we can actually have correct map
        auto [problems, originalGraphs, nodeProblemMapping, positionInProblem,
            globalTerminalIndex] = splitConnectedComponents(G, terminals);

        std::vector<std::vector<NodeID> > solutions;
        std::vector<NodeID> globalSolution;

        FlowType flow_sum = 0;
        for (size_t p = 0; p < problems.size(); ++p) {
            auto& problem = problems[p];
            if (debug) {
                graph_algorithms::checkGraphValidity(problem.graph);
            }

            std::vector<NodeID> terminals;
            for (size_t i = 0; i < problem.terminals.size(); ++i) {
                terminals.emplace_back(problem.terminals[i].position);
            }

            branch_multicut bmc(originalGraphs[p], terminals);
            auto p_pointer = std::make_shared<multicut_problem>(problem);
            auto [sol, flow] = bmc.find_multiterminal_cut(p_pointer);
            flow_sum += flow;

            if (cfg->write_solution || cfg->inexact) {
                solutions.emplace_back(sol);
            }
        }

        if (cfg->write_solution || cfg->inexact) {
            std::vector<NodeID> blocksize(cfg->num_terminals, 0);
            for (NodeID i = 0; i < G->n(); ++i) {
                if (nodeProblemMapping[i] == UNDEFINED_NODE) {
                    globalSolution.emplace_back(0);
                    blocksize[0]++;
                } else {
                    NodeID pr = nodeProblemMapping[i];
                    NodeID localId = positionInProblem[i];
                    NodeID sol = solutions[pr][localId];
                    globalSolution.emplace_back(globalTerminalIndex[pr][sol]);

                    blocksize[globalSolution[i]]++;
                }
            }
            LOG1 << "size: " << blocksize;

            if (cfg->write_solution) {
                graph_io gio;
                gio.writeVector(globalSolution, "solution");
            }

            if (cfg->inexact) {
                // local search for better solutions
                for (NodeID n : G->nodes()) { }
            }
        }

        return flow_sum;
    }

 private:
    static std::vector<NodeID> addSurroundingAreaToTerminals(
        std::shared_ptr<mutable_graph> graph,
        std::vector<NodeID> terminals) {
        auto config = configuration::getConfig();
        std::vector<NodeID> terminalMapping(graph->n(), UNDEFINED_NODE);

        for (size_t i = 0; i < terminals.size(); ++i) {
            terminalMapping[terminals[i]] = i;
        }

        if (config->bfs_size > 1) {
            for (size_t i = 0; i < terminals.size(); ++i) {
                NodeID t = terminals[i];
                size_t curr_blocksize = 1;
                std::queue<NodeID> bfs_q;
                bfs_q.emplace(t);
                terminalMapping[t] = i;

                while (!bfs_q.empty() && curr_blocksize < config->bfs_size) {
                    NodeID n = bfs_q.front();
                    bfs_q.pop();
                    for (EdgeID e : graph->edges_of(n)) {
                        NodeID tgt = graph->getEdgeTarget(n, e);
                        if (terminalMapping[tgt] == UNDEFINED_NODE) {
                            terminalMapping[tgt] = i;
                            bfs_q.push(tgt);
                            if (++curr_blocksize >= config->bfs_size) {
                                break;
                            }
                        }
                    }
                }
            }
        }
        return terminalMapping;
    }

    static std::tuple<std::vector<multicut_problem>,
                      std::vector<mutable_graph>,
                      std::vector<NodeID>,
                      std::vector<NodeID>,
                      std::vector<std::vector<NodeID> > >
    splitConnectedComponents(std::shared_ptr<mutable_graph> G,
                             const std::vector<NodeID>& terminalMapping) {
        std::vector<multicut_problem> problems;
        std::vector<mutable_graph> originalGraphs;
        std::vector<int> t_comp;
        strongly_connected_components cc;

        auto config = configuration::getConfig();

        auto [components, num_comp, blocksizes] = cc.strong_components(G);
        (void)blocksizes;

        std::vector<NodeID> nodeProblemMapping(G->n(), UNDEFINED_NODE);
        std::vector<NodeID> positionInProblem(G->n(), UNDEFINED_NODE);
        std::vector<std::vector<bool> > terminalInProblem(num_comp);
        std::vector<size_t> numTerminalsInProblem(num_comp);

        for (size_t i = 0; i < terminalInProblem.size(); ++i) {
            terminalInProblem[i].resize(config->num_terminals, false);
        }

        for (size_t i = 0; i < components.size(); ++i) {
            size_t c = static_cast<size_t>(components[i]);
            nodeProblemMapping[i] = c;
            if (terminalMapping[i] != UNDEFINED_NODE) {
                if (!terminalInProblem[c][terminalMapping[i]]) {
                    numTerminalsInProblem[c]++;
                    terminalInProblem[c][terminalMapping[i]] = true;
                }
            }
        }

        std::vector<NodeID> ctr(num_comp, 0);
        std::vector<NodeID> num_terminals(num_comp, 0);

        std::vector<std::vector<NodeID> > globalTerminalIndex;

        for (size_t problem = 0; problem < num_comp; problem++) {
            if (numTerminalsInProblem[problem] >= 1) {
                globalTerminalIndex.emplace_back();
                graph_extractor ge;
                auto [G_out, mapping, reverse_mapping] =
                    ge.extract_block(G, problem, components);

                originalGraphs.emplace_back(*G_out);

                std::vector<std::vector<NodeID> > termBlocks(
                    config->num_terminals);

                for (size_t i = 0; i < mapping.size(); ++i) {
                    NodeID n = mapping[i];
                    if (positionInProblem[n] != UNDEFINED_NODE) {
                        LOG1 << "ERROR: NODE IN MULTIPLE CCs";
                        exit(1);
                    }
                    positionInProblem[n] = i;
                    if (terminalMapping[n] != UNDEFINED_NODE) {
                        termBlocks[terminalMapping[n]].emplace_back(i);
                    }
                }

                std::vector<terminal> problemTerms;
                for (size_t t = 0; t < termBlocks.size(); ++t) {
                    if (termBlocks[t].size() > 0) {
                        problemTerms.emplace_back(termBlocks[t][0],
                                                  problemTerms.size());
                        globalTerminalIndex[problem].emplace_back(t);
                    }
                }

                problems.emplace_back(G_out, problemTerms);
                auto p = std::make_shared<multicut_problem>(problems.back());
                graph_contraction::contractIsolatingBlocks(p, termBlocks);
            }
        }

        return std::make_tuple(problems,
                               originalGraphs,
                               nodeProblemMapping,
                               positionInProblem,
                               globalTerminalIndex);
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

    std::vector<NodeID> distantTerminals(std::shared_ptr<mutable_graph> G) {
        std::vector<NodeID> terminals;
        size_t r = random_functions::nextInt(0, G->n() - 1);
        int dt = configuration::getConfig()->distant_terminals;
        for (int i = 0; i < dt; ++i) {
            std::queue<NodeID> Q;
            std::vector<bool> found(G->n(), false);

            Q.push(r);
            found[r] = true;
            for (auto t : terminals) {
                Q.push(t);
                found[t] = true;
            }

            while (!Q.empty()) {
                NodeID t = Q.front();
                Q.pop();
                for (EdgeID e : G->edges_of(t)) {
                    NodeID n = G->getEdgeTarget(t, e);
                    if (!found[n]) {
                        Q.push(n);
                        found[n] = true;
                    }
                }

                if (Q.empty()) {
                    terminals.emplace_back(t);
                    break;
                }
            }
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

    std::vector<NodeID> presetFileTerminals(std::shared_ptr<mutable_graph>) {
        auto config = configuration::getConfig();
        auto v = graph_io::readVector<NodeID>(config->partition_file);
        std::vector<NodeID> terminalMapping;

        for (NodeID n : v) {
            if (n == config->num_terminals) {
                terminalMapping.emplace_back(UNDEFINED_NODE);
            } else {
                terminalMapping.emplace_back(n);
            }
        }

        return terminalMapping;
    }

    std::vector<NodeID> orderFileTerminals(std::shared_ptr<mutable_graph> G) {
        auto config = configuration::getConfig();
        auto v = graph_io::readVector<NodeID>(config->partition_file);
        auto o = graph_io::readVector<NodeID>(config->partition_file + ".pos");
        std::vector<NodeID> terminalMapping(G->n(), UNDEFINED_NODE);
        NodeID terminal_size = 1;
        if (config->preset_percentage > 0) {
            NodeID blocksize = G->n() / config->num_terminals;
            size_t desired_size = blocksize * config->preset_percentage / 100;
            size_t minimum_size = 1;
            terminal_size = std::max(desired_size, minimum_size);
        }

        for (size_t i = 0; i < config->num_terminals; ++i) {
            for (size_t n = 0; n < G->n(); ++n) {
                if (v[n] == i && o[n] <= terminal_size) {
                    terminalMapping[n] = i;
                }
            }
        }

        return terminalMapping;
    }
};
