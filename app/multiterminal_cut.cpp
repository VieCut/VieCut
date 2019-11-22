/******************************************************************************
 * multiterminal_cut.cpp
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2018-2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include <omp.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>

#include "algorithms/multicut/multicut_problem.h"
#include "algorithms/multicut/multiterminal_cut.h"
#include "data_structure/graph_access.h"
#include "gperftools/malloc_extension.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tlx/string/split.hpp"
#include "tools/random_functions.h"
#include "tools/string.h"
#include "tools/timer.h"

int main(int argn, char** argv) {
    static constexpr bool debug = false;
    tlx::CmdlineParser cmdl;
    auto config = configuration::getConfig();
    size_t num_terminals;

    cmdl.add_param_string("graph", config->graph_filename,
                          "path to graph file");
    cmdl.add_stringlist('t', "terminal", config->term_strings,
                        "add terminal vertex");
    cmdl.add_int('k', "top_k", config->top_k,
                 "multiterminal cut between top k vertices (invalidates t)");
    cmdl.add_int('r', "random_k", config->random_k,
                 "multiterminal cut between k random vertices");
    cmdl.add_size_t('p', "proc", config->threads, "number of threads");
    cmdl.add_size_t('s', "seed", config->seed, "random seed");
    cmdl.add_size_t('b', "bfs_size", config->bfs_size,
                    "use bfs to find and start on terminal block");
    cmdl.add_size_t('n', "preset_percentage", config->preset_percentage,
                    "percentag of vertices that are preset");
    cmdl.add_string('e', "edge_selection", config->edge_selection,
                    "edge selection rule");
    cmdl.add_string('q', "queue_type", config->queue_type,
                    "Type of priority queue used");
    cmdl.add_size_t('c', "contraction_type", config->contraction_type,
                    "Disable contraction mechanisms");
    cmdl.add_string('f', "partition_file", config->partition_file,
                    "Partition file");

    if (!cmdl.process(argn, argv))
        return -1;

    // MallocExtension::instance()->SetMemoryReleaseRate(0.0);

    random_functions::setSeed(config->seed);

    std::vector<NodeID> terminals;
    std::shared_ptr<mutable_graph> G = mutable_graph::from_graph_access(
        graph_io::readGraphWeighted(config->graph_filename));

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

        num_terminals = terminals.size();
        config->total_terminals = num_terminals;

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
            num_terminals = num_partitions;
            config->total_terminals = num_terminals;

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
            num_terminals = num_partitions;
            config->total_terminals = num_terminals;
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

    if (terminals.size() < 2) {
        std::cerr << "ERROR: Number of terminals ("
                  << terminals.size()
                  << ") too small! Exiting..." << std::endl;
        exit(-1);
    }

    multiterminal_cut mc;
    timer t;
    FlowType flow = mc.multicut(G, terminals);
    std::cout << "RESULT selection_rule=" << config->edge_selection
              << " pq=" << config->queue_type
              << " contraction_type=" << config->contraction_type
              << " graph=" << config->graph_filename
              << " time=" << t.elapsed()
              << " terminals=" << num_terminals
              << " cut=" << flow
              << " n=" << G->number_of_nodes()
              << " m=" << G->number_of_edges() / 2
              << " processes=" << config->threads
              << " seed=" << config->seed << std::endl;
}
