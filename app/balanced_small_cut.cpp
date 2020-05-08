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

static void augmentMostBalancedCut(
    graphAccessPtr original_graph,
    const std::vector<std::pair<NodeID, EdgeID> >& mb_edges) {
    auto random_it = std::next(std::begin(mb_edges), random_functions::nextInt(
                                   0, mb_edges.size() - 1));

    auto [n, e] = *random_it;
    original_graph->setEdgeWeight(e, original_graph->getEdgeWeight(e) + 1);
    NodeID in_n = original_graph->getEdgeSource(e);
    NodeID in_t = original_graph->getEdgeTarget(e);

    for (EdgeID e : original_graph->edges_of(in_t)) {
        if (original_graph->getEdgeTarget(e) == in_n) {
            original_graph->setEdgeWeight(
                e, original_graph->getEdgeWeight(e) + 1);
        }
    }
}

static void augmentCutEdges(graphAccessPtr original_graph,
                            mutableGraphPtr mg,
                            NodeID largest_id) {
    std::vector<std::vector<EdgeID> > interblock_edges;
    for (NodeID n : mg->nodes()) {
        interblock_edges.emplace_back();
        if (n == largest_id)
            continue;

        for (NodeID v : mg->containedVertices(n)) {
            size_t edges = 0;
            size_t nonedges = 0;
            for (EdgeID e : original_graph->edges_of(v)) {
                edges++;
                NodeID t = original_graph->getEdgeTarget(e);
                if (mg->getCurrentPosition(t) != n) {
                    nonedges++;
                    interblock_edges.back().emplace_back(e);
                }
            }
        }
    }

    std::vector<bool> increased(mg->n(), false);
    for (NodeID n : original_graph->nodes()) {
        NodeID block = mg->getCurrentPosition(n);
        if (!increased[block] && block != largest_id) {
            auto random_it = std::next(
                std::begin(interblock_edges[block]),
                random_functions::nextInt(
                    0, interblock_edges[block].size() - 1));

            increased[block] = true;
            EdgeID e = *random_it;

            original_graph->setEdgeWeight(
                e, original_graph->getEdgeWeight(e) + 1);

            NodeID n = original_graph->getEdgeSource(e);
            NodeID in_t = original_graph->getEdgeTarget(e);
            increased[mg->getCurrentPosition(in_t)] = true;
            for (EdgeID e : original_graph->edges_of(in_t)) {
                if (original_graph->getEdgeTarget(e) == n) {
                    original_graph->setEdgeWeight(
                        e, original_graph->getEdgeWeight(e) + 1);
                    break;
                }
            }
        }
    }
}

int main(int argn, char** argv) {
    timer t;
    tlx::CmdlineParser cmdl;
    auto cfg = configuration::getConfig();
    bool output = false;
    bool augment_all = false;
    bool augment_most_balanced = false;
    cmdl.add_param_string("graph", cfg->graph_filename, "path to graph file");
    cmdl.add_bool('o', "output", output, "Write intermediate graphs to disk");
    cmdl.add_flag('v', "verbose", cfg->verbose, "Verbose output.");
    cmdl.add_size_t('p', "processes", cfg->threads, "Number of processes!");
    cmdl.add_flag('c', "find lowest conductance", cfg->find_lowest_conductance,
                  "Find cut with lowest conductance");
    cmdl.add_flag('a', "augment", augment_all, "augment all minimum cuts");
    cmdl.add_flag('b', "augment_most_balanced", augment_most_balanced,
                  "augment most balanced minimum cut");

    cfg->find_most_balanced_cut = true;
    cfg->save_cut = true;
    cfg->set_node_in_cut = true;

    if (!cmdl.process(argn, argv))
        return -1;

    graphAccessPtr original_graph =
        graph_io::readGraphWeighted(cfg->graph_filename);

#ifdef PARALLEL
    parallel_cactus<graphAccessPtr> mc;
#else
    cactus_mincut<graphAccessPtr> mc;
#endif
    LOG1 << "io time: " << t.elapsed();

    EdgeWeight cut = 0;
    timer t_this;
    size_t ct = 0;

    auto G = original_graph;
    std::vector<graphAccessPtr> graph_vec = { original_graph };
    union_find uf(original_graph->number_of_nodes());

    while (G->number_of_nodes() > 1) {
        t_this.restart();
        auto [current_cut, mg, mb_edges] = mc.findAllMincuts(graph_vec);
        bool augment_mb = ((augment_most_balanced || augment_all)
                           && (current_cut > 1));
        std::vector<NodeID> largest_block;
        NodeID largest_id = 0;
        for (NodeID n : mg->nodes()) {
            if (mg->containedVertices(n).size() > largest_block.size()) {
                largest_block = mg->containedVertices(n);
                largest_id = n;
            }
        }

        for (NodeID n : original_graph->nodes()) {
            NodeID pos = mg->getCurrentPosition(n);
            if (pos != largest_id) {
                for (EdgeID e : original_graph->edges_of(n)) {
                    NodeID t = original_graph->getEdgeTarget(e);
                    NodeID pos_t = mg->getCurrentPosition(t);
                    if (pos == pos_t) {
                        uf.Union(n, t);
                    } else {
                        if (!augment_mb) {
                            uf.Union(n, t);
                        }
                    }
                }
            }
        }

        if (augment_mb) {
            if (augment_all) {
                augmentCutEdges(original_graph, mg, largest_id);
            } else {
                augmentMostBalancedCut(original_graph, mb_edges);
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
