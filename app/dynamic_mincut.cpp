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
#include "parallel/algorithm/exact_parallel_minimum_cut.h"
#include "parallel/algorithm/parallel_cactus.h"
#else
#include "algorithms/global_mincut/cactus/cactus_mincut.h"
#include "algorithms/global_mincut/noi_minimum_cut.h"
#endif

#include "algorithms/global_mincut/dynamic/dynamic_mincut.h"
#include "common/configuration.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tlx/logger.hpp"
#include "tlx/string.hpp"
#include "tools/random_functions.h"

int main(int argn, char** argv) {
    tlx::CmdlineParser cmdl;
    auto cfg = configuration::getConfig();
    std::string initial_graph = "";
    std::string dynamic_edges = "";
    bool run_static = false;
    cmdl.add_string('i', "initial_graph", initial_graph, "path to graph file");
    cmdl.add_string('d', "dynamic_edges", dynamic_edges, "path to edge list");
#ifdef PARALLEL
    size_t procs = 1;
    cmdl.add_size_t('p', "proc", procs, "number of processes");
#endif
    cmdl.add_size_t('r', "seed", cfg->seed, "random seed");
    cmdl.add_bool('s', "static", run_static, "run static algorithm");
    cmdl.add_flag('v', "verbose", cfg->verbose, "more verbose logs");

    if (!cmdl.process(argn, argv))
        return -1;

    auto edgenames = tlx::split('.', dynamic_edges);
    std::string seedstr = "";
    std::string delstr = "";
    std::string insstr = "";
    if (edgenames.size() >= 4) {
        seedstr = edgenames[edgenames.size() - 1];
        delstr = edgenames[edgenames.size() - 2];
        insstr = edgenames[edgenames.size() - 3];
    }

    random_functions::setSeed(cfg->seed);
    bool verbose = configuration::getConfig()->verbose;

#ifdef PARALLEL
    LOGC(verbose) << "PARALLEL DEFINED, USING " << procs << " THREADS";
    omp_set_num_threads(procs);
    parallel_cactus<mutableGraphPtr> cactus;
#else
    LOGC(verbose) << "PARALLEL NOT DEFINED";
    cactus_mincut<mutableGraphPtr> cactus;
#endif

    cfg->save_cut = true;

    mutableGraphPtr G;

    if (dynamic_edges == "") {
        LOG1 << "ERROR: No list of dynamic edges given! Use parameter -d!";
        exit(1);
    }

    auto [numV, tempEdges] = graph_io::readTemporalGraph(dynamic_edges);

    LOG1 << "graph with " << numV << " vertices!";
    if (initial_graph == "") {
        G = std::make_shared<mutable_graph>();
        G->start_construction(numV);
        G->finish_construction();
    } else {
        G = graph_io::readGraphWeighted<mutable_graph>(initial_graph);
    }

    timer t;
    size_t ctr = 0;
    size_t cutchange = 0;
    size_t staticruns = 0;

    if (run_static) {
#ifdef PARALLEL
        exact_parallel_minimum_cut<mutableGraphPtr> static_alg;
#else
        noi_minimum_cut<mutableGraphPtr> static_alg;
#endif
        EdgeID previous_timestamp = std::get<3>(tempEdges[0]);
        EdgeWeight previous_cut = static_alg.perform_minimum_cut(G);
        size_t edgesInBatch = 0;
        for (auto [s, t, w, timestamp] : tempEdges) {
            if (timestamp != previous_timestamp) {
                edgesInBatch = 0;
                previous_timestamp = timestamp;
                staticruns++;
                EdgeWeight current_cut = static_alg.perform_minimum_cut(G);
                if (current_cut != previous_cut) {
                    previous_cut = current_cut;
                    cutchange++;
                }
            }
            edgesInBatch++;
            ctr++;
            if (ctr % 1000 == 0) {
                LOG1 << ctr << " out of " << tempEdges.size()
                     << " -> " << staticruns << " batches so far";
            }
            if (s == t) continue;
            if (w > 0) {
                G->new_edge_order(s, t, w);
            } else {
                EdgeID eToT = UNDEFINED_EDGE;
                for (EdgeID e : G->edges_of(s)) {
                    if (G->getEdgeTarget(s, e) == t) {
                        eToT = e;
                        break;
                    }
                }
                if (eToT == UNDEFINED_EDGE) {
                    LOG1 << "Warning: delete edge that doesn't exist!";
                } else {
                    G->deleteEdge(s, eToT);
                }
            }
        }
    } else {
        dynamic_mincut dynmc;
        EdgeWeight previous_cut = dynmc.initialize(G);
        EdgeWeight current_cut = 0;
        for (auto [s, t, w, timestamp] : tempEdges) {
            ctr++;
            if (s == t) continue;
            if (w > 0) {
                current_cut = dynmc.addEdge(s, t, w);
            } else {
                current_cut = dynmc.removeEdge(s, t);
            }
            if (current_cut != previous_cut) {
                previous_cut = current_cut;
                cutchange++;
            }
        }
        staticruns = dynmc.getCallsOfStaticAlgorithm();
    }

    std::string graph = initial_graph;
    if (graph == "") {
        graph = dynamic_edges;
    }

    LOG1 << "RESULT"
         << " graph=" << graph
         << " time=" << t.elapsed()
         << " cutchange=" << cutchange
         << " insert=" << insstr
         << " delete=" << delstr
         << " seed=" << seedstr
         << " static=" << run_static
         << " runs_static_alg=" << staticruns;
}
