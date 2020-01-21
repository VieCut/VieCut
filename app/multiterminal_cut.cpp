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

#include <mpi.h>
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
#include "tools/random_functions.h"
#include "tools/string.h"
#include "tools/timer.h"

int main(int argn, char** argv) {
    MPI_Init(&argn, &argv);
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    tlx::CmdlineParser cmdl;
    auto config = configuration::getConfig();

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
    cmdl.add_string('f', "partition_file", config->partition_file,
                    "Partition file");
    cmdl.add_flag('i', "use_ilp", config->use_ilp, "Use ILP");

    if (!cmdl.process(argn, argv))
        return -1;

    random_functions::setSeed(config->seed);

    std::vector<NodeID> terminals;

    std::shared_ptr<mutable_graph> G;
    multiterminal_cut mc;
    if (mpi_rank == 0) {
        G = mutable_graph::from_graph_access(
            graph_io::readGraphWeighted(config->graph_filename));
        terminals = mc.setOriginalTerminals(G);
    }

    size_t termsize = terminals.size();
    MPI_Bcast(&termsize, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    terminals.resize(termsize);
    MPI_Bcast(&terminals.front(), termsize, MPI_INT, 0, MPI_COMM_WORLD);

    if (config->total_terminals == 0)
        config->total_terminals = terminals.size();

    if (terminals.size() < 2) {
        std::cerr << "ERROR: Number of terminals (" << terminals.size()
                  << ") too small! Exiting..." << std::endl;
        exit(-1);
    }

    timer t;
    FlowType flow = mc.multicut(G, terminals);
    std::cout << "RESULT selection_rule=" << config->edge_selection
              << " pq=" << config->queue_type
              << " graph=" << config->graph_filename
              << " time=" << t.elapsed()
              << " terminals=" << config->total_terminals
              << " cut=" << flow
              << " n=" << G->number_of_nodes()
              << " m=" << G->number_of_edges() / 2
              << " use_ilp=" << config->use_ilp
              << " difference=" << config->bound_difference
              << " contract_n=" << config->n
              << " contract_m=" << config->m
              << " processes=" << config->threads
              << " seed=" << config->seed << std::endl;
}
