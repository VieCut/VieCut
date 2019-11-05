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
    // cfg->find_lowest_conductance = true;
    cfg->save_cut = true;
    cfg->optimization = 6;

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
        auto [current_cut, mg, mb_edges] = mc.findAllMincuts(graph_vec);

        NodeWeight largest_block = 0;
        NodeID largest_id = 0;
        for (NodeID n : mg->nodes()) {
            if (mg->containedVertices(n).size() > largest_block) {
                largest_block = mg->containedVertices(n).size();
                largest_id = n;
            }
        }

        std::vector<std::tuple<NodeID, NodeID, NodeWeight> > heaviest_neighbors(
            mg->n(), std::make_tuple(0, 0, 0));

        for (NodeID n : original_graph->nodes()) {
            NodeID pos = mg->getCurrentPosition(n);
            if (pos != largest_id) {
                for (EdgeID e : original_graph->edges_of(n)) {
                    NodeID t = original_graph->getEdgeTarget(e);
                    NodeID pos_t = mg->getCurrentPosition(t);
                    if (pos == pos_t) {
                        uf.Union(n, t);
                    } else {
                        if (mb_edges.count(e) == 0) {
                            uf.Union(n, t);
                        }
                    }
                }
            }
        }

        auto random_it =
            std::next(std::begin(mb_edges), random_functions::nextInt(
                          0, mb_edges.size() - 1));

        EdgeID e = *random_it;
        original_graph->setEdgeWeight(e, original_graph->getEdgeWeight(e) + 1);
        NodeID in_n = original_graph->getEdgeSource(e);
        NodeID in_t = original_graph->getEdgeTarget(e);

        for (EdgeID e : original_graph->edges_of(in_t)) {
            if (original_graph->getEdgeTarget(e) == in_n) {
                original_graph->setEdgeWeight(
                    e, original_graph->getEdgeWeight(e) + 1);
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
