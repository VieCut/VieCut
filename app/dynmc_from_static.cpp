/******************************************************************************
 * dynmc_from_static.cpp
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
    std::string init_graph = "";
    bool run_static = false;
    double insert_factor = 0.0;
    double remove_factor = 0.0;
    size_t batch_size = 1;
    cmdl.add_string('g', "initial_graph", init_graph, "path to graph file");
    cmdl.add_double('i', "insert_factor", insert_factor,
                    "factor of edges that are added dynamically");
    cmdl.add_double('r', "remove_factor", remove_factor,
                    "factor of edges that are removed dynamically");
    cmdl.add_size_t('b', "batch_size", batch_size, "batch size");
    cmdl.add_size_t('r', "seed", configuration::getConfig()->seed, "rnd seed");
    cmdl.add_bool('s', "run_static", run_static, "Run static algorithm");
    cmdl.add_bool('v', "vbs", configuration::getConfig()->verbose, "verbose");
    cmdl.add_size_t('d', "depth",
                    configuration::getConfig()->depthOfPartialRelabeling,
                    "depth of partial relabeling");

    configuration::getConfig()->save_cut = true;

    if (!cmdl.process(argn, argv))
        return -1;

    random_functions::setSeed(configuration::getConfig()->seed);
    graphAccessPtr G = graph_io::readGraphWeighted<graph_access>(init_graph);

    size_t insert_edges =
        std::floor(insert_factor * static_cast<double>(G->m()) * 0.5);
    size_t remove_edges =
        std::floor(remove_factor * static_cast<double>(G->m()) * 0.5);

    if (insert_factor + remove_factor >= 0.999) {
        LOG1 << "ERROR: Insert factor + remove factor >= 1. Exiting!";
        exit(1);
    }

    size_t ie = 0;
    size_t re = 0;
    std::vector<bool> insertMarked(G->m(), false);
    std::vector<bool> removeMarked(G->m(), false);

    EdgeID nextEdge = random_functions::nextInt(0, G->m() - 1);
    while (ie < insert_edges) {
        while (insertMarked[nextEdge]
               || G->getEdgeTarget(nextEdge) < G->getEdgeSource(nextEdge)) {
            nextEdge = random_functions::nextInt(0, G->m() - 1);
        }
        insertMarked[nextEdge] = true;
        ie++;
    }
    while (re < remove_edges) {
        while (insertMarked[nextEdge] || removeMarked[nextEdge]
               || G->getEdgeTarget(nextEdge) < G->getEdgeSource(nextEdge)) {
            nextEdge = random_functions::nextInt(0, G->m() - 1);
        }
        removeMarked[nextEdge] = true;
        re++;
    }

    std::vector<std::tuple<NodeID, NodeID, int64_t, uint64_t> > edges;
    LOG1 << "insert " << insert_edges << " remove " << remove_edges;
    mutableGraphPtr dynG = std::make_shared<mutable_graph>();
    std::vector<std::tuple<NodeID, NodeID, int64_t, bool> > dynEdges;

    dynG->start_construction(G->n());
    for (EdgeID e : G->edges()) {
        NodeID src = G->getEdgeSource(e);
        NodeID tgt = G->getEdgeTarget(e);
        EdgeWeight wgt = G->getEdgeWeight(e);
        if (tgt < src) continue;
        if (!insertMarked[e]) {
            dynG->new_edge_order(src, tgt, wgt);
        } else {
            dynEdges.emplace_back(src, tgt, wgt, true);
        }
        if (removeMarked[e]) {
            dynEdges.emplace_back(src, tgt, wgt, false);
        }
    }
    dynG->finish_construction();
    std::shuffle(dynEdges.begin(), dynEdges.end(), random_functions::getRand());

    timer run_timer;
    size_t ctr = 0;
    size_t cutchange = 0;
    size_t staticruns = 1;
    LOG1 << dynG->m();
    size_t insert = 0;
    size_t remove = 0;
    if (run_static) {
        size_t currentBatchSize = 0;
#ifdef PARALLEL
        exact_parallel_minimum_cut<mutableGraphPtr> static_alg;
#else
        noi_minimum_cut<mutableGraphPtr> static_alg;
#endif
        EdgeWeight previous_cut = static_alg.perform_minimum_cut(dynG);
        for (auto [s, t, w, isInsert] : dynEdges) {
            if (currentBatchSize >= batch_size) {
                currentBatchSize = 0;
                staticruns++;
                EdgeWeight current_cut = static_alg.perform_minimum_cut(dynG);
                if (current_cut != previous_cut) {
                    previous_cut = current_cut;
                    cutchange++;
                    LOG1 << "after " << ctr << " " << current_cut;
                }
            }
            currentBatchSize++;
            ctr++;
            if (s == t) continue;
            if (isInsert) {
                dynG->new_edge_order(s, t, w);
                insert++;
            } else {
                EdgeID eToT = UNDEFINED_EDGE;
                for (EdgeID e : dynG->edges_of(s)) {
                    if (dynG->getEdgeTarget(s, e) == t) {
                        eToT = e;
                        break;
                    }
                }
                if (eToT == UNDEFINED_EDGE) {
                    LOG1 << "Warning: delete edge that doesn't exist!";
                } else {
                    dynG->deleteEdge(s, eToT);
                    remove++;
                }
            }
        }
        staticruns++;
        EdgeWeight current_cut = static_alg.perform_minimum_cut(dynG);
        if (current_cut != previous_cut) {
            LOG1 << "at end, cut " << current_cut;
            cutchange++;
        }
    } else {
        dynamic_mincut dynmc;
        EdgeWeight previous_cut = dynmc.initialize(dynG);
        EdgeWeight current_cut = 0;
        for (auto [s, t, w, isInsert] : dynEdges) {
            ctr++;
            if (s == t) continue;
            if (isInsert) {
                insert++;
                current_cut = dynmc.addEdge(s, t, w);
            } else {
                remove++;
                current_cut = dynmc.removeEdge(s, t);
            }
            if (current_cut != previous_cut) {
                previous_cut = current_cut;
                cutchange++;
                LOG1 << "after " << ctr << " " << current_cut;
            }
        }
        staticruns = dynmc.getCallsOfStaticAlgorithm();
    }
    LOG1 << insert << " " << remove;
    LOG1 << dynG->m();

    LOG1 << "RESULT"
         << " graph=" << init_graph
         << " time=" << run_timer.elapsed()
         << " cutchange=" << cutchange
         << " insert=" << insert_edges
         << " delete=" << remove_edges
         << " seed=" << configuration::getConfig()->seed
         << " static=" << run_static
         << " runs_static_alg=" << staticruns;
}
