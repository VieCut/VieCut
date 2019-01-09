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

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <omp.h>
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

    std::string graph_filename;
    std::string algo;
    std::string pq = "default";
    int iter = 1;
    double contraction_factor;

    cmdl.add_param_string("graph", graph_filename, "path to graph file");

#ifdef PARALLEL
    std::vector<std::string> procs;
    cmdl.add_stringlist('p', "proc", procs, "number of processes");
    cmdl.add_param_string("algo", algo, "algorithm name ('vc', 'exact')");
#else
    cmdl.add_param_string("algo", algo, "algorithm name ('vc', 'noi', 'matula', 'pr', 'ks')");
#endif

    cmdl.add_param_string("graph", graph_filename, "path to graph file");
    cmdl.add_param_string("algo", algo, "global_mincut name");
    cmdl.add_string('q', "pq", pq, "name of priority queue implementation");
    cmdl.add_int('i', "iter", iter, "number of iterations");
    cmdl.add_double('c', "contraction factor", contraction_factor, "contract until only n*(1-contraction_factor) vertices are left");

    if (!cmdl.process(argn, argv))
        return -1;

    std::shared_ptr<graph_access> G = graph_io::readGraphWeighted(graph_filename);
    timer t;

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
        LOG1 << procs[i] << " is not a valid number of workers! Continuing without.";
    }
#else
    LOG1 << "PARALLEL NOT DEFINED";
#endif
    if (numthreads.empty())
        numthreads.emplace_back(1);

    for (int i = 0; i < iter; ++i) {
        for (size_t numthread : numthreads) {
            random_functions::setSeed(i);
            minimum_cut* mc = selectMincutAlgorithm(algo);
            t.restart();
            omp_set_num_threads(numthread);

            auto G2 = sf.one_ks(G, contraction_factor, i);

            LOGC(timing) << "Contraction: " << t.elapsed() << "s";
            EdgeWeight cut = mc->perform_minimum_cut(G2);

#ifdef SAVECUT
            graph_io::writeCut(G, graph_filename + ".cut");
#endif

            std::string graphname = string::basename(graph_filename);

            double time = t.elapsed();
            std::cout << "RESULT source=taa algo=contract" << (int)100 * contraction_factor << " graph="
                      << graphname << " time=" << time
                      << " cut=" << cut
                      << " n=" << G->number_of_nodes()
                      << " m=" << G->number_of_edges() << " processes="
                      << omp_get_max_threads() << std::endl;
        }
    }
}
