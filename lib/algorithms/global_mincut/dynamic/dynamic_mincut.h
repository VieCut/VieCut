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
    bool verbose;
    mutableGraphPtr original_graph;
    mutableGraphPtr out_cactus;
    EdgeWeight current_cut;
#ifdef PARALLEL
    parallel_cactus<mutableGraphPtr> cactus;
#else
    cactus_mincut<mutableGraphPtr> cactus;
#endif

 public:
    dynamic_mincut() {
        verbose = configuration::getConfig()->verbose;
    }

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
                    LOGC(verbose) << "full recompute from empty";
                    auto [cut, o, bal] = cactus.findAllMincuts(original_graph);
                    out_cactus = o;
                    current_cut = cut;
                } else {
                    LOGC(verbose) << "contract two empty vtcs";
                    out_cactus->contractVertexSet({ sCactusPos, tCactusPos });
                }
            } else {
                auto vtxset = cactus_path::findPath(out_cactus, sCactusPos,
                                                    tCactusPos, current_cut);
                if (vtxset.size() == out_cactus->n()) {
                    LOGC(verbose) << "full recompute";
                    auto [cut, outg, b] = cactus.findAllMincuts(original_graph);
                    out_cactus = outg;
                    current_cut = cut;
                } else {
                    LOGC(verbose) << "contract set of size " << vtxset.size();
                    out_cactus->contractVertexSet(vtxset);
                }
            }
        }
        LOGC(verbose) << "t " << timer.elapsed() << " cut " << current_cut
                      << " vtcs_in_cactus " << out_cactus->n();
        return current_cut;
    }

    EdgeWeight removeEdge(NodeID s, NodeID t) {
        timer timer;
        EdgeID eToT = UNDEFINED_EDGE;
        for (EdgeID e : original_graph->edges_of(s)) {
            if (original_graph->getEdgeTarget(s, e) == t) {
                eToT = e;
                break;
            }
        }

        if (eToT == UNDEFINED_EDGE) {
            LOG1 << "Warning: Deleting edge that does not exist! Doing nothing";
            return current_cut;
        }

        EdgeWeight wgt = original_graph->getEdgeWeight(s, eToT);

        original_graph->deleteEdge(s, eToT);

        if (wgt == 0) {
            LOGC(verbose) << "edge has zero weight, current cut remains same";
            return current_cut;
        }

        NodeID sCactusPos = out_cactus->getCurrentPosition(s);
        NodeID tCactusPos = out_cactus->getCurrentPosition(t);

        if (sCactusPos != tCactusPos) {
            LOGC(verbose) << "previously mincut between vertices, recompute";
            auto [cut, outg, b] = cactus.findAllMincuts(original_graph);
            out_cactus = outg;
            current_cut = cut;
        } else {
            push_relabel<true> pr;
            auto [flow, sourceset] = pr.solve_max_flow_min_cut(
                original_graph, { 0, 1 }, 0, false, false, current_cut);
            if (static_cast<EdgeWeight>(flow) >= current_cut) {
                LOGC(verbose) << "cut not changed!";
            } else {
                auto [cut, outg, b] = cactus.findAllMincuts(original_graph);
                out_cactus = outg;
                current_cut = cut;
                LOGC(verbose) << "recomputing, minimum cut changed to " << cut;
            }
        }
        LOGC(verbose) << "t " << timer.elapsed() << " cut " << current_cut;
        return current_cut;
    }

    mutableGraphPtr getOriginalGraph() {
        return original_graph;
    }
};
