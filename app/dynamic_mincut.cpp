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
    auto cfg = configuration::getConfig();
    std::string initial_graph = "";
    std::string dynamic_edges = "";
    cmdl.add_string('i', "initial_graph", initial_graph, "path to graph file");
    cmdl.add_string('d', "dynamic_edges", dynamic_edges, "path to edge list");
#ifdef PARALLEL
    size_t procs = 1;
    cmdl.add_size_t('p', "proc", procs, "number of processes");
#endif
    cmdl.add_flag('v', "verbose", cfg->verbose, "more verbose logs");
    cmdl.add_size_t('r', "seed", cfg->seed, "random seed");

    if (!cmdl.process(argn, argv))
        return -1;

    random_functions::setSeed(cfg->seed);

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

    mutableGraphPtr G;

    if (dynamic_edges == "") {
        LOG1 << "ERROR: No list of dynamic edges given! Use parameter -d!";
        exit(1);
    }

    auto [numV, tempEdges] = graph_io::readTemporalGraph(dynamic_edges);

    if (initial_graph == "") {
        G = std::make_shared<mutable_graph>();
        G->start_construction(numV);
        G->finish_construction();
    } else {
        G = graph_io::readGraphWeighted<mutable_graph>(initial_graph);
    }

    timer t;

    dynamic_mincut dynmc;
    dynmc.initialize(G);

    size_t ctr = 0;
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
