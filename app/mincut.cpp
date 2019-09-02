/******************************************************************************
 * mincut.cpp
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2017-2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <sstream>

#include "algorithms/global_mincut/algorithms.h"
#include "algorithms/global_mincut/minimum_cut.h"
#include "algorithms/global_mincut/minimum_cut_helpers.h"
#include "data_structure/graph_access.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tools/random_functions.h"
#include "tools/string.h"
#include "tools/timer.h"

int main(int argn, char** argv) {
    static constexpr bool debug = false;

    tlx::CmdlineParser cmdl;
    size_t num_iterations = 1;

    auto cfg = configuration::getConfig();

    cmdl.add_param_string("graph", cfg->graph_filename, "path to graph file");
#ifdef PARALLEL
    std::vector<std::string> procs;
    cmdl.add_stringlist('p', "proc", procs, "number of processes");
    cmdl.add_param_string("algo", cfg->algorithm, "algorithm name");
#else
    cmdl.add_param_string("algo", cfg->algorithm, "algorithm name");
#endif

    cmdl.add_string('q', "pq", cfg->queue_type,
                    "name of priority queue implementation");
    cmdl.add_size_t('i', "iter", num_iterations, "number of iterations");
    cmdl.add_bool('l', "disable_limiting", cfg->disable_limiting,
                  "disable limiting of PQ values");
    cmdl.add_bool('s', "save_cut", cfg->save_cut,
                  "compute and store minimum cut");
    cmdl.add_double('c', "contraction_factor", cfg->contraction_factor,
                    "contraction factor for pre-run of viecut");
    cmdl.add_string('k', "sampling_type", cfg->sampling_type,
                    "sampling variant for pre-run of viecut");
    cmdl.add_flag('b', "balanced", cfg->find_most_balanced_cut,
                  "find most balanced minimum cut");

    if (!cmdl.process(argn, argv))
        return -1;

    std::vector<int> numthreads;

    std::shared_ptr<graph_access> G =
        graph_io::readGraphWeighted(
            configuration::getConfig()->graph_filename);
    timer t;

    if (G->getMinDegree() == 0) {
        LOG1 << "empty nodes are bad, exiting";
        exit(1);
    }

    LOG1 << "io time: " << t.elapsed();
    // ***************************** perform cut *****************************
#ifdef PARALLEL
    LOG1 << "PARALLEL DEFINED!";

    size_t i;
    try {
        for (i = 0; i < procs.size(); ++i) {
            numthreads.emplace_back(std::stoi(procs[i]));
        }
    } catch (...) {
        LOG1 << procs[i]
             << " is not a valid number of workers! Continuing without.";
    }
#else
    LOG1 << "PARALLEL NOT DEFINED";
#endif
    if (numthreads.empty())
        numthreads.emplace_back(1);
    timer tdegs;

    for (size_t i = 0; i < num_iterations; ++i) {
        for (int numthread : numthreads) {
            auto seed = i;
            LOG << "random seed " << seed;
            random_functions::setSeed(seed);
            configuration::getConfig()->seed = seed;

            NodeID n = G->number_of_nodes();
            EdgeID m = G->number_of_edges();

            minimum_cut* mc = selectMincutAlgorithm(cfg->algorithm);
            omp_set_num_threads(numthread);
            cfg->threads = numthread;

            t.restart();
            EdgeWeight cut;
            cut = mc->perform_minimum_cut(G);

            if (cfg->save_cut && false) {
                graph_io::writeCut(G, cfg->graph_filename + ".cut");
            }

            std::string graphname = string::basename(cfg->graph_filename);

            std::string algprint = cfg->algorithm;
#ifdef PARALLEL
            algprint += "par";
#endif
            algprint += cfg->pq;

            if (cfg->disable_limiting) {
                algprint += "unlimited";
            }

            std::cout << "RESULT source=taa algo=" << algprint << " graph="
                      << graphname << " time=" << t.elapsed()
                      << " cut=" << cut << " n=" << n
                      << " m=" << m / 2 << " processes="
                      << numthread << std::endl;
        }
    }
}
