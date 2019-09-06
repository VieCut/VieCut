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

#include "common/configuration.h"
#include "io/graph_io.h"
#include "parallel/coarsening/contract_graph.h"
#include "tlx/cmdline_parser.hpp"
#include "tools/random_functions.h"

int main(int argn, char** argv) {
    static const bool debug = true;

    tlx::CmdlineParser cmdl;

    auto cfg = configuration::getConfig();

    std::string output_path = "";
    size_t min = 0;
    size_t max = 0;
    size_t multi = 0;

    cmdl.add_param_string("graph", cfg->graph_filename, "path to graph file");
    cmdl.add_param_string("output_path", output_path, "output path");
    cmdl.add_param_size_t("min", min, "minimum edge weight");
    cmdl.add_param_size_t("max", max, "maximum edge weight");
    cmdl.add_param_size_t(
        "multi", multi,
        "multiply weight of one (outgoing) edge per vertex by this factor");

    if (!cmdl.process(argn, argv)) {
        LOG << "Error in command line processing!";
        exit(1);
    }

    std::unordered_set<NodeID> heavy_edges;
    std::unordered_map<NodeID, EdgeWeight> edge_weights;

    auto G = graph_io::readGraphWeighted(cfg->graph_filename);

    size_t edges = 0;
    for (NodeID n : G->nodes()) {
        EdgeID e_highweight = random_functions::nextInt(
            G->get_first_edge(n), G->get_first_invalid_edge(n) - 1);

        edges++;
        NodeID t = G->getEdgeTarget(e_highweight);
        heavy_edges.insert(contraction::get_uint64_from_pair(n, t));
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
