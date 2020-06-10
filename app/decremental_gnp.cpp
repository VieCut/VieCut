/******************************************************************************
 * decremental_gnp.cpp
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

bool existsEdge(mutableGraphPtr G, NodeID s, NodeID t) {
    if (s == t) {
        return true;
    }

    for (EdgeID e : G->edges_of(s)) {
        NodeID tgt = G->getEdgeTarget(s, e);
        if (tgt == t) {
            return true;
        }
    }
    return false;
}

int main(int argn, char** argv) {
    tlx::CmdlineParser cmdl;
    auto cfg = configuration::getConfig();
    std::string initial_graph = "";
    std::string dynamic_edges = "";
    bool run_static = false;
    size_t edges_per_vertex = 0;
    cmdl.add_size_t('e', "edges_per_vertex", edges_per_vertex,
                    "num edges to add to each vertex to overlay gnp");
    cmdl.add_string('i', "initial_graph", initial_graph, "path to graph file");
#ifdef PARALLEL
    size_t procs = 1;
    cmdl.add_size_t('p', "proc", procs, "number of processes");
#endif
    cmdl.add_size_t('r', "seed", cfg->seed, "random seed");
    cmdl.add_bool('s', "static", run_static, "run static algorithm");
    cmdl.add_flag('v', "verbose", cfg->verbose, "more verbose logs");

    if (!cmdl.process(argn, argv))
        return -1;

    auto G = graph_io::readGraphWeighted<mutable_graph>(initial_graph);
    std::vector<std::pair<NodeID, NodeID> > decrementalEdges;

    LOG1 << "Creating edges...";
    for (auto v : G->nodes()) {
        size_t i = 0;
        while (i < edges_per_vertex) {
            NodeID r = random_functions::nextInt(0, G->n() - 1);
            if (!existsEdge(G, v, r)) {
                G->new_edge_order(v, r, 1);
                decrementalEdges.emplace_back(v, r);
                ++i;
            }
        }
    }

    LOG1 << "Permuting...";

    random_functions::permutate_vector_good(&decrementalEdges);

    LOG1 << "Starting decremental...";

    timer time;
    bool verbose = cfg->verbose;

    size_t ctr = 0;
    if (run_static) {
#ifdef PARALLEL
        exact_parallel_minimum_cut<mutableGraphPtr> static_alg;
#else
        noi_minimum_cut<mutableGraphPtr> static_alg;
#endif
        EdgeWeight previous_cut = static_alg.perform_minimum_cut(G);
        for (auto [s, t] : decrementalEdges) {
            LOGC(verbose) << ctr << " after " << time.elapsed();
            ctr++;
            EdgeID eToT = UNDEFINED_EDGE;
            for (EdgeID e : G->edges_of(s)) {
                if (G->getEdgeTarget(s, e) == t) {
                    eToT = e;
                    break;
                }
            }
            if (eToT == UNDEFINED_EDGE) {
                LOG1 << "Warning: delete edge that doesn't exist!";
            }
            G->deleteEdge(s, eToT);
            EdgeWeight current_cut = static_alg.perform_minimum_cut(G);
            if (current_cut != previous_cut) {
                previous_cut = current_cut;
                LOG1 << "cut changed to " << current_cut << " after "
                     << time.elapsed() << " and " << ctr;
            }
        }
    } else {
        dynamic_mincut dynmc;
        EdgeWeight previous_cut = dynmc.initialize(G);
        EdgeWeight current_cut = 0;
        for (auto [s, t] : decrementalEdges) {
            LOGC(verbose) << ctr << " after " << time.elapsed();
            ctr++;
            current_cut = dynmc.removeEdge(s, t);

            if (current_cut != previous_cut) {
                previous_cut = current_cut;
                LOG1 << "cut changed to " << current_cut << " after "
                     << time.elapsed() << " and " << ctr;
            }
        }
    }

    LOG1 << time.elapsed();
}
