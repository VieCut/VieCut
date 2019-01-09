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

#include <algorithm>
#include <algorithms/global_mincut/minimum_cut_helpers.h>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <omp.h>
#include <sstream>

#include "algorithms/global_mincut/algorithms.h"
#include "algorithms/global_mincut/minimum_cut.h"
#include "data_structure/graph_access.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tools/random_functions.h"
#include "tools/string.h"
#include "tools/timer.h"

int main(int argn, char** argv) {

    static constexpr bool debug = false;

    tlx::CmdlineParser cmdl;

    std::string graph_filename;
    std::string algo;
    std::string pq = "default";
    int iter = 1;
    bool disable_limiting = false;

    cmdl.add_param_string("graph", graph_filename, "path to graph file");
#ifdef PARALLEL
    std::vector<std::string> procs;
    cmdl.add_stringlist('p', "proc", procs, "number of processes");
    cmdl.add_param_string("algo", algo, "algorithm name ('vc', 'exact')");
#else
    cmdl.add_param_string("algo", algo, "algorithm name ('vc', 'noi', 'matula', 'pr', 'ks')");
#endif

    cmdl.add_string('q', "pq", pq, "name of priority queue implementation");
    cmdl.add_int('i', "iter", iter, "number of iterations");
    cmdl.add_bool('l', "disable_limiting", disable_limiting, "disable limiting of PQ values");

    if (!cmdl.process(argn, argv))
        return -1;

    std::vector<int> numthreads;

    std::shared_ptr<graph_access> G = graph_io::readGraphWeighted(graph_filename);
    timer t;

    if (G->getMinDegree() == 0) {
        LOG1 << "empty nodes are bad, exiting";
        exit(1);
    }

    LOG1 << "io time: " << t.elapsed();
    // ***************************** perform cut ***************************************
#ifdef PARALLEL
    LOG1 << "PARALLEL DEFINED!";

    size_t i;
    try {
        for (i = 0; i < procs.size(); ++i) {
            numthreads.emplace_back(std::stoi(procs[i]));
        }
    } catch (...) {
        LOG1 << procs[i] << " is not a valid number of workers! Continuing without.";
    }
#else
    LOG1 << "PARALLEL NOT DEFINED";
#endif
    if (numthreads.empty())
        numthreads.emplace_back(1);
    timer tdegs;

    for (int i = 0; i < iter; ++i) {
        for (int numthread : numthreads) {
            auto seed = i;// std::chrono::system_clock::now().time_since_epoch().count();
            LOG << "random seed " << seed;
            random_functions::setSeed(seed);

            minimum_cut* mc = selectMincutAlgorithm(algo);
            omp_set_num_threads(numthread);

            t.restart();
            EdgeWeight cut;
            if (algo == "noi") {
                cut = mc->perform_minimum_cut(G, pq, disable_limiting);
            }
            else {
                cut = mc->perform_minimum_cut(G);
            }

#ifdef SAVECUT
            graph_io::writeCut(G, graph_filename + ".cut");
#endif

            std::string graphname = string::basename(graph_filename);

            std::string algprint = algo;
#ifdef PARALLEL
            algprint += "par";
#endif
            algprint += pq;

            if (disable_limiting) {
                algprint += "unlimited";
            }

            std::cout << "RESULT source=taa algo=" << algprint << " graph="
                      << graphname << " time=" << t.elapsed()
                      << " cut=" << cut << " mindeg=" << G->getMinDegree() << " n=" << G->number_of_nodes()
                      << " m=" << G->number_of_edges() / 2 << " processes="
                      << numthread << std::endl;
        }
    }
}
