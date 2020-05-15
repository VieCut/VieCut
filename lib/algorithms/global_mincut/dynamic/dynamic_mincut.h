/******************************************************************************
 * dynamic_minimum.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2020 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#ifdef PARALLEL
#include "parallel/algorithm/parallel_cactus.h"
#else
#include "algorithms/global_mincut/cactus/cactus_mincut.h"
#endif

#include "algorithms/global_mincut/dynamic/cactus_path.h"
#include "common/definitions.h"
#include "data_structure/mutable_graph.h"
#include "tlx/logger.hpp"
#include "tools/timer.h"

class dynamic_mincut {
 private:
    mutableGraphPtr original_graph;
    mutableGraphPtr out_cactus;
    EdgeWeight current_cut;
#ifdef PARALLEL
    parallel_cactus<mutableGraphPtr> cactus;
#else
    cactus_mincut<mutableGraphPtr> cactus;
#endif

 public:
    dynamic_mincut() { }
    ~dynamic_mincut() { }

    EdgeWeight initialize(mutableGraphPtr graph) {
        auto [cut, outgraph, balanced] = cactus.findAllMincuts(graph);
        original_graph = graph;
        out_cactus = outgraph;
        current_cut = cut;
        return cut;
    }

    EdgeWeight addEdge(NodeID s, NodeID t, EdgeWeight w) {
        timer timer;
        NodeID sCactusPos = out_cactus->getCurrentPosition(s);
        NodeID tCactusPos = out_cactus->getCurrentPosition(t);
        original_graph->new_edge_order(s, t, w);

        if (sCactusPos != tCactusPos) {
            if (current_cut == 0) {
                if (out_cactus->n() == 2) {
                    LOG1 << "full recompute from empty";
                    auto [cut, o, bal] = cactus.findAllMincuts(original_graph);
                    out_cactus = o;
                    current_cut = cut;
                } else {
                    LOG1 << "contract two empty vtcs";
                    out_cactus->contractVertexSet({ sCactusPos, tCactusPos });
                }
            } else {
                auto vtxset = cactus_path::findPath(out_cactus, sCactusPos,
                                                    tCactusPos, current_cut);
                if (vtxset.size() == out_cactus->n()) {
                    LOG1 << "full recompute";
                    auto [cut, outg, bal] = cactus.findAllMincuts(original_graph);
                    out_cactus = outg;
                    current_cut = cut;
                } else {
                    LOG1 << "contract set of size " << vtxset.size();
                    out_cactus->contractVertexSet(vtxset);
                }
            }
        }
        LOG1 << "t " << timer.elapsed() << " cut " << current_cut
             << " vtcs_in_cactus " << out_cactus->n();
        return current_cut;
    }
};
