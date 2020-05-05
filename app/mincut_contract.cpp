/******************************************************************************
 * mincut_contract.cpp
 *
 * Source of VieCut.
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
#include "data_structure/graph_access.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tlx/logger.hpp"
#include "tools/string.h"

#ifdef PARALLEL
#include "parallel/coarsening/sparsify.h"
#else
#include "coarsening/sparsify.h"
#endif

int main(int argn, char** argv) {
    static constexpr bool debug = false;
    static const bool timing = true;

    tlx::CmdlineParser cmdl;

    int iter = 1;
    double contraction_factor;

    cmdl.add_param_string("graph", configuration::getConfig()->graph_filename,
                          "path to graph file");

#ifdef PARALLEL
    std::vector<std::string> procs;
    cmdl.add_stringlist('p', "proc", procs, "number of processes");
    cmdl.add_param_string("algo", configuration::getConfig()->algorithm,
                          "algorithm name ('vc', 'exact')");
#else
    cmdl.add_param_string("algo", configuration::getConfig()->algorithm,
                          "algorithm name ('vc', 'noi', 'matula', 'pr', 'ks')");
#endif

    cmdl.add_param_string("graph", configuration::getConfig()->graph_filename,
                          "path to graph file");
    cmdl.add_string('q', "pq", configuration::getConfig()->pq,
                    "name of priority queue implementation");
    cmdl.add_int('i', "iter", iter, "number of iterations");
    cmdl.add_double('c', "contraction factor", contraction_factor,
                    "contract until only 1 - factor of vertices are left");
    cmdl.add_bool('s', "save_cut", configuration::getConfig()->save_cut,
                  "compute and store minimum cut");
    cmdl.add_string('k', "sampling_type",
                    configuration::getConfig()->sampling_type,
                    "sampling variant for pre-run of viecut");

    if (!cmdl.process(argn, argv))
        return -1;

    timer t;
    graphAccessPtr G = graph_io::readGraphWeighted(
        configuration::getConfig()->graph_filename);

    if (G->getMinDegree() == 0) {
        LOG1 << "empty nodes are bad, exiting";
        exit(1);
    }

    LOG << "io time: " << t.elapsed();
    t.restart();

    sparsify sf;

    std::vector<size_t> numthreads;

#ifdef PARALLEL
    LOG1 << "PARALLEL DEFINED";
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

    for (int i = 0; i < iter; ++i) {
        for (size_t numthread : numthreads) {
            random_functions::setSeed(i);
            configuration::getConfig()->seed = i;
            minimum_cut* mc = selectMincutAlgorithm(
                configuration::getConfig()->algorithm);
            t.restart();
            omp_set_num_threads(numthread);

            auto G2 = sf.one_ks(G);

            LOGC(timing) << "Contraction: " << t.elapsed() << "s";
            EdgeWeight cut = mc->perform_minimum_cut(G2);

            if (configuration::getConfig()->save_cut) {
                graph_io::writeCut(G,
                                   configuration::getConfig()->graph_filename
                                   + ".cut");
            }

            double time = t.elapsed();
            std::cout << "RESULT source=taa algo=contract"
                      << static_cast<int>(100) * contraction_factor
                      << " graph="
                      << configuration::getConfig()->graph_filename
                      << " time=" << time
                      << " cut=" << cut
                      << " n=" << G->number_of_nodes()
                      << " m=" << G->number_of_edges() << " processes="
                      << omp_get_max_threads() << std::endl;
        }
    }
}
