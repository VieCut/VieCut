/******************************************************************************
 * multiterminal_cut.cpp
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2018-2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
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
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tools/random_functions.h"
#include "tools/string.h"
#include "tools/timer.h"

#ifdef USE_TCMALLOC
#include "gperftools/malloc_extension.h"
#endif

int main(int argn, char** argv) {
    MPI_Init(&argn, &argv);
    int mpi_rank, mpi_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    tlx::CmdlineParser cmdl;
    auto config = configuration::getConfig();

    cmdl.add_param_string("graph", config->graph_filename,
                          "path to graph file");

    cmdl.add_size_t('b', "bfs_size", config->bfs_size,
                    "use bfs to find and start on terminal block");
    cmdl.add_size_t('c', "contraction_depth_around_terminals",
                    config->contractionDepthAroundTerminal,
                    "Contract vertices close to heaviest terminal [only -X]");
    cmdl.add_flag('d', "disable_cpu_affinity", config->disable_cpu_affinity,
                  "Default CPU affinity (i.e. let the OS decide where to run)");
    cmdl.add_int('D', "distant_terminals", config->distant_terminals,
                 "Use terminals with high distance from each other");
    cmdl.add_string('e', "edge_selection", config->edge_selection,
                    "edge selection rule");
    cmdl.add_string('f', "partition_file", config->partition_file,
                    "Partition file");
    cmdl.add_flag('i', "use_ilp", config->use_ilp, "Use ILP");
    cmdl.add_double('l', "removeTerminalsBeforeBranch",
                    config->removeTerminalsBeforeBranch,
                    "Remove low degree terminals before branch [only -X]");
    cmdl.add_int('k', "top_k", config->top_k,
                 "multiterminal cut between top k vertices (invalidates t)");
    cmdl.add_double('n', "preset_percentage", config->preset_percentage,
                    "percentag of vertices that are preset");
    cmdl.add_string('o', "first_branch_path", config->first_branch_path,
                    "Print graph at time of first branching, then terminate.");
    cmdl.add_size_t('p', "proc", config->threads, "number of threads");
    cmdl.add_string('q', "queue_type", config->queue_type,
                    "Type of priority queue used");
    cmdl.add_int('r', "random_k", config->random_k,
                 "multiterminal cut between k random vertices");
    cmdl.add_size_t('s', "seed", config->seed, "random seed");
    cmdl.add_stringlist('t', "terminal", config->term_strings,
                        "add terminal vertex");
    cmdl.add_flag('w', "write_solution", config->write_solution,
                  "Print best solution");
    cmdl.add_flag('X', "inexact", config->inexact, "Apply inexact heuristics");

    if (!cmdl.process(argn, argv))
        return -1;

    config->random_flows = config->threads / 2;
    config->high_distance_flows = config->threads / 2;

    random_functions::setSeed(config->seed);

    std::vector<NodeID> terminals;

    std::shared_ptr<mutable_graph> G;
    FlowType flow;
    timer t;
    try {
        multiterminal_cut mc;
        G = mutable_graph::from_graph_access(
            graph_io::readGraphWeighted(config->graph_filename));
        terminals = mc.setOriginalTerminals(G);

        if (config->num_terminals < 2) {
            std::cerr << "ERROR: Number of terminals (" << terminals.size()
                      << ") too small! Exiting..." << std::endl;
            exit(-1);
        }

        LOG1 << config->graph_filename << " terminals:"
             << config->num_terminals << " contraction:"
             << config->preset_percentage;
        t.restart();
        flow = mc.multicut(G, terminals, config->num_terminals);

#ifdef USE_GUROBI
    } catch (GRBException e) {
        LOG1 << e.getMessage();
#endif
    } catch (const std::exception& e) {
        LOG1 << e.what();
    }
    std::cout << "RESULT selection_rule=" << config->edge_selection
              << " pq=" << config->queue_type
              << " graph=" << config->graph_filename
              << " time=" << t.elapsed()
              << " terminals=" << config->num_terminals
              << " cut=" << flow
              << " n=" << G->number_of_nodes()
              << " m=" << G->number_of_edges() / 2
              << " use_ilp=" << config->use_ilp
              << " difference=" << config->bound_difference
              << " contract_n=" << config->n
              << " contract_m=" << config->m
              << " processes=" << config->threads
              << " inexact=" << config->inexact
              << " seed=" << config->seed << std::endl;
    MPI_Finalize();
}
