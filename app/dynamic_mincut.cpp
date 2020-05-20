/******************************************************************************
 * dynamic_mincut.cpp
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2020 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <stddef.h>

#include <ext/alloc_traits.h>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#ifdef PARALLEL
#include "parallel/algorithm/parallel_cactus.h"
#else
#include "algorithms/global_mincut/cactus/cactus_mincut.h"
#endif

#include "algorithms/global_mincut/dynamic/dynamic_mincut.h"
#include "common/configuration.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tlx/logger.hpp"
#include "tools/random_functions.h"

int main(int argn, char** argv) {
    tlx::CmdlineParser cmdl;
    size_t num_iterations = 1;
    auto cfg = configuration::getConfig();
    cmdl.add_param_string("graph", cfg->graph_filename, "path to graph file");
#ifdef PARALLEL
    size_t procs = 1;
    cmdl.add_size_t('p', "proc", procs, "number of processes");
#endif
    cmdl.add_size_t('i', "iter", num_iterations, "number of iterations");
    cmdl.add_flag('v', "verbose", cfg->verbose, "more verbose logs");
    cmdl.add_size_t('r', "seed", cfg->seed, "random seed");

    if (!cmdl.process(argn, argv))
        return -1;

#ifdef PARALLEL
    LOGC(cfg->verbose) << "PARALLEL DEFINED, USING " << procs << " THREADS";
    omp_set_num_threads(procs);
    parallel_cactus<mutableGraphPtr> cactus;
#else
    LOGC(cfg->verbose) << "PARALLEL NOT DEFINED";
    cactus_mincut<mutableGraphPtr> cactus;
#endif

    random_functions::setSeed(cfg->seed);

    cfg->save_cut = true;

    auto [numV, tempEdges] = graph_io::readTemporalGraph(cfg->graph_filename);
    mutableGraphPtr G = std::make_shared<mutable_graph>();
    G->start_construction(numV);
    G->finish_construction();

    // auto G = graph_io::readGraphWeighted<mutable_graph>(cfg->graph_filename);

    dynamic_mincut dynmc;
    dynmc.initialize(G);

    size_t ctr = 0;
    timer t;
    for (auto [s, t, w, timestamp] : tempEdges) {
        LOG1 << ctr++ << " of " << tempEdges.size();
        if (w > 0) {
            dynmc.addEdge(s, t, w);
        } else {
            dynmc.removeEdge(s, t);
        }
    }
    LOG1 << t.elapsed();
    mutableGraphPtr mgp = dynmc.getOriginalGraph();
    graph_io::writeGraphWeighted(mgp->to_graph_access(), "output.graph");
}
