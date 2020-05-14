/******************************************************************************
 * make_graph_weighted.cpp
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <stdlib.h>
#include <time.h>

#include <cstdint>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "common/configuration.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "io/graph_io.h"
#include "parallel/coarsening/contract_graph.h"
#include "tlx/cmdline_parser.hpp"
#include "tlx/logger.hpp"
#include "tools/random_functions.h"

int main(int argn, char** argv) {
    static const bool debug = true;

    tlx::CmdlineParser cmdl;

    random_functions::setSeed(time(NULL));

    auto cfg = configuration::getConfig();

    std::string output_path = "";
    size_t min = 0;
    size_t max = 0;
    size_t multi = 0;

    cfg->seed = UNDEFINED_NODE;

    cmdl.add_param_string("graph", cfg->graph_filename, "path to graph file");
    cmdl.add_param_string("output_path", output_path, "output path");
    cmdl.add_param_size_t("min", min, "minimum edge weight");
    cmdl.add_param_size_t("max", max, "maximum edge weight");
    cmdl.add_param_size_t(
        "multi", multi,
        "multiply weight of one (outgoing) edge per vertex by this factor");
    cmdl.add_size_t('s', "seed", cfg->seed,
                    "random seed (using time otherwise)");

    if (!cmdl.process(argn, argv)) {
        LOG << "Error in command line processing!";
        exit(1);
    }

    if (cfg->seed != UNDEFINED_NODE) {
        random_functions::setSeed(cfg->seed);
    }

    LOG1 << "Using seed " << random_functions::getSeed();

    auto G = graph_io::readGraphWeighted(cfg->graph_filename);

    std::vector<bool> covered(G->number_of_nodes());
    std::unordered_set<uint64_t> heavy_edges;
    std::unordered_map<uint64_t, EdgeWeight> edge_weights;

    size_t edges = 0;
    for (NodeID n : G->nodes()) {
        // if (!covered[n]) {
        EdgeID e_highweight = random_functions::nextInt(
            G->get_first_edge(n), G->get_first_invalid_edge(n) - 1);
        NodeID t = G->getEdgeTarget(e_highweight);

        edges++;
        uint64_t pair = contraction::get_uint64_from_pair(n, t);
        heavy_edges.insert(pair);
        covered[t] = true;
        covered[n] = true;

//        }
    }

    for (NodeID n : G->nodes()) {
        if (n % 100000 == 0) {
            LOG1 << n;
        }
        if (n % 10000 == 0) {
            std::cout << "." << std::flush;
        }

        for (EdgeID e : G->edges_of(n)) {
            NodeID tgt = G->getEdgeTarget(e);
            uint64_t pair = contraction::get_uint64_from_pair(n, tgt);
            if (tgt > n) {
                EdgeWeight wgt = random_functions::nextInt(min, max);
                if (heavy_edges.count(pair) > 0) {
                    wgt *= multi;
                }
                G->setEdgeWeight(e, wgt);
                edge_weights[pair] = wgt;
            } else {
                EdgeWeight wgt = edge_weights[pair];
                G->setEdgeWeight(e, wgt);
            }
        }
    }

    std::cout << std::endl;
    LOG1 << "Writing graph to " << output_path;
    graph_io::writeGraphWeighted(G, output_path);
}
