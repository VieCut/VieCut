/******************************************************************************
 * kcore.cpp
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
#include <cstdint>
#include <iostream>
#include <memory>
#include <sstream>

#include "algorithms/global_mincut/noi_minimum_cut.h"
#include "algorithms/misc/core_decomposition.h"
#include "algorithms/misc/strongly_connected_components.h"
#include "data_structure/graph_access.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tools/graph_extractor.h"

int main(int argn, char** argv) {
    static const bool debug = false;

    tlx::CmdlineParser cmdl;

    std::string graph_filename;
    bool no_cut = false;
    bool lowest_core = false;

    std::vector<std::string> cores;
    cmdl.add_stringlist('k', "cores", cores, "kCores");
    cmdl.add_bool('l', "lowest_core", lowest_core,
                  "Search for lowest core where cut is not min degree");
    cmdl.add_bool('c', "no_cut", no_cut, "Disable minimum cut testing.");
    cmdl.add_param_string("graph", graph_filename, "path to graph file");

    if (!cmdl.process(argn, argv))
        return -1;

    timer t;
    std::shared_ptr<graph_access> G =
        graph_io::readGraphWeighted(graph_filename);

    LOG << "io time: " << t.elapsed();
    t.restart();
    k_cores kCores = core_decomposition::batagelj_zaversnik(G);
    std::vector<size_t> tgts;
    uint32_t max_core =
        kCores.degrees[kCores.vertices[G->number_of_nodes() - 1]];

    if (lowest_core) {
        size_t min_core = 2;
        for (size_t core = min_core; core < max_core; core++) {
            if (kCores.buckets[core] != kCores.buckets[core - 1]
                || (core == min_core)) {
                tgts.push_back(core);
            }
        }
    } else {
        size_t i;
        try {
            for (i = 0; i < cores.size(); ++i) {
                tgts.emplace_back(std::stoi(cores[i]));
            }
        } catch (...) {
            LOG1 << cores[i] << " is not a valid kCore. Continuing without.";
        }
    }

    if (tgts.empty()) {
        LOG1 << "ERROR: No target kCore. Set with -l or -k";
        exit(1);
    }

    for (size_t target_core : tgts) {
        if (max_core >= target_core) {
            auto connected_graph = core_decomposition::createCoreGraph(
                kCores, target_core, G);
            std::string outputname = graph_filename
                                     + "_core_" + std::to_string(target_core);

            LOG1 << "output graph: strongly connected component with core "
                 << target_core << " nodes: "
                 << connected_graph->number_of_nodes()
                 << " edges: " << connected_graph->number_of_edges();

            size_t result = 0;
            if (!no_cut) {
                minimum_cut* mc = new noi_minimum_cut();
                result = mc->perform_minimum_cut(connected_graph);
            }

            if (result < connected_graph->getMinDegree()) {
                std::string out_path = graph_filename + "_core_" +
                                       std::to_string(target_core);
                LOG1 << "small min cut " << result << " -> saving";
                graph_io::writeGraph(connected_graph, out_path);
                graph_io::writeGraphDimacsKS(connected_graph,
                                             out_path + "_ks", "cut");
                if (lowest_core)
                    exit(0);
            } else {
                LOG1 << "minimum degree equals minimum cut "
                     << result << " " << connected_graph->getMinDegree();
            }
        }
    }
}
