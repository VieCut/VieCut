/******************************************************************************
 * balanced_small_cut.cpp
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#ifdef PARALLEL
#include "parallel/algorithm/parallel_cactus.h"
#else
#include "algorithms/global_mincut/cactus/cactus_mincut.h"
#endif
#include "algorithms/global_mincut/viecut.h"
#include "algorithms/misc/strongly_connected_components.h"
#include "data_structure/graph_access.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tlx/logger.hpp"
#include "tools/graph_extractor.h"
#include "tools/macros_assertions.h"
#include "tools/timer.h"

int main(int argn, char** argv) {
    timer t;
    tlx::CmdlineParser cmdl;
    auto cfg = configuration::getConfig();
    bool output = false;
    cmdl.add_param_string("graph", cfg->graph_filename, "path to graph file");
    cmdl.add_bool('o', "output", output, "Write intermediate graphs to disk");
    cmdl.add_flag('v', "verbose", cfg->verbose, "Verbose output.");
    cmdl.add_size_t('p', "processes", cfg->threads, "Number of processes!");

    cfg->find_most_balanced_cut = true;
    cfg->save_cut = true;

    if (!cmdl.process(argn, argv))
        return -1;

    std::shared_ptr<graph_access> original_graph =
        graph_io::readGraphWeighted(cfg->graph_filename);

#ifdef PARALLEL
    parallel_cactus mc;
#else
    cactus_mincut mc;
#endif
    LOG1 << "io time: " << t.elapsed();

    EdgeWeight cut = 0;
    timer t_this;
    size_t ct = 0;

    auto G = original_graph;
    std::vector<std::shared_ptr<graph_access> > graph_vec = { original_graph };

    union_find uf(original_graph->number_of_nodes());

    while (G->number_of_nodes() > 1) {
        t_this.restart();
        auto [current_cut, mg] = mc.findAllMincuts(graph_vec);

        NodeWeight largest_block = 0;
        NodeID largest_id = 0;
        for (NodeID n : mg->nodes()) {
            if (mg->containedVertices(n).size() > largest_block) {
                largest_block = mg->containedVertices(n).size();
                largest_id = n;
            }
        }

        for (NodeID n : mg->nodes()) {
            if (n != largest_id) {
                if (mg->containedVertices(n).size() > 0) {
                    NodeID v0 = mg->containedVertices(n)[0];
                    for (NodeID v : mg->containedVertices(n)) {
                        uf.Union(v0, v);
                        for (EdgeID e : original_graph->edges_of(v)) {
                            NodeID t = original_graph->getEdgeTarget(e);
                            uf.Union(v, t);
                        }
                    }
                }
            }
        }

        G = contraction::fromUnionFind(original_graph, &uf);
        graph_vec = { original_graph, G };

        cut = current_cut;

        if (output) {
            std::string name = cfg->graph_filename + "_" + std::to_string(ct++);
            graph_io::writeGraphWeighted(G, name);
        }

        LOG1 << "Minimum cut " << cut;
        LOG1 << "-------------------------------";
        LOG1 << "Next graph: G=(" << G->number_of_nodes()
             << "," << G->number_of_edges() << ")";
    }
}
