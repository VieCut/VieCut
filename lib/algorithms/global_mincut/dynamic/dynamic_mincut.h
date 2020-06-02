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

#include <unordered_set>
#include <vector>

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
    size_t flow_problem_id;
    size_t max_cache_size = 100;

    std::vector<bool> currentlyCaching;
    std::vector<std::pair<mutableGraphPtr, std::vector<std::tuple<
                                                           NodeID, NodeID, EdgeWeight> > > > problemCache;

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
        timer t;
        auto [cut, outgraph, balanced] = cactus.findAllMincuts(graph);
        original_graph = graph;
        out_cactus = outgraph;
        current_cut = cut;
        flow_problem_id = random_functions::next();
        LOGC(verbose) << "initialize t " << t.elapsed() << " cut " << cut
                      << " cactus_vtcs " << outgraph->n();
        return cut;
    }

    EdgeWeight addEdge(NodeID s, NodeID t, EdgeWeight w) {
        timer timer;
        NodeID sCactusPos = out_cactus->getCurrentPosition(s);
        NodeID tCactusPos = out_cactus->getCurrentPosition(t);
        original_graph->new_edge_order(s, t, w);
        cacheEdge(s, t, w);
        if (sCactusPos != tCactusPos) {
            if (current_cut == 0) {
                if (out_cactus->n() == 2) {
                    LOGC(verbose) << "full recompute from empty";
                    auto [cut, o, bal] = cactus.findAllMincuts(original_graph);
                    out_cactus = o;
                    current_cut = cut;
                    // RECREATE FROM CACHE - SOMETHING IS SUUUUPER SLOW RIGHT NOW
                    LOG1 << currentlyCaching.size() << " " << cut;
                    if (currentlyCaching.size() > cut && currentlyCaching[cut]) {
                        for (auto [s, t, w] : problemCache[cut].second) {
                            sLOG1 << s << t << w;
                        }
                    }
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
                    LOG1 << currentlyCaching.size() << " " << cut;
                    if (currentlyCaching.size() > cut && currentlyCaching[cut]) {
                        for (auto [s, t, w] : problemCache[cut].second) {
                            sLOG1 << s << t << w;
                        }
                    }
                    LOG1 << "ins2 " << timer.elapsed();
                } else {
                    LOGC(verbose) << "contract set of size " << vtxset.size();
                    contractVertexSet(vtxset);
                }
            }
        }
        LOGC(verbose) << "t " << timer.elapsed() << " cut " << current_cut
                      << " vtcs_in_cactus " << out_cactus->n();
        return current_cut;
    }

    void contractVertexSet(const std::unordered_set<NodeID>& vtxset) {
        // if one vertex has high degree and all others don't, it is faster
        // to explicitly contract others into this high degree vertex instead of
        // standard set contraction (also check that no vertex is empty so we
        // have handles on the vertices)
        bool alternativeContract = true;

        NodeID high_degree = UNDEFINED_NODE;
        size_t numNonlow = 0;

        for (auto v : vtxset) {
            if (out_cactus->numContainedVertices(v) == 0) {
                alternativeContract = false;
                break;
            }
            if (out_cactus->getUnweightedNodeDegree(v) > 100) {
                if (high_degree != UNDEFINED_NODE) {
                    alternativeContract = false;
                    break;
                } else {
                    high_degree = v;
                }
            }

            if (out_cactus->getUnweightedNodeDegree(v) > 10) {
                numNonlow++;
            }
            if (numNonlow > 1) {
                alternativeContract = false;
                break;
            }
        }

        if (alternativeContract && high_degree != UNDEFINED_NODE) {
            NodeID high_origid = out_cactus->containedVertices(high_degree)[0];
            std::vector<NodeID> orig_ids;
            for (auto v : vtxset) {
                if (v != high_degree) {
                    orig_ids.emplace_back(out_cactus->containedVertices(v)[0]);
                }
            }

            for (auto v : orig_ids) {
                NodeID s = out_cactus->getCurrentPosition(high_origid);
                NodeID t = out_cactus->getCurrentPosition(v);
                EdgeID conn_edge = UNDEFINED_EDGE;
                for (EdgeID e : out_cactus->edges_of(t)) {
                    if (out_cactus->getEdgeTarget(t, e) == s) {
                        conn_edge = out_cactus->getReverseEdge(t, e);
                        break;
                    }
                }
                if (conn_edge == UNDEFINED_EDGE) {
                    out_cactus->contractSparseTargetNoEdge(s, t);
                } else {
                    out_cactus->contractEdgeSparseTarget(s, conn_edge);
                }
            }
        } else {
            out_cactus->contractVertexSet(vtxset);
        }
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
        NodeID sCactusPos = out_cactus->getCurrentPosition(s);
        NodeID tCactusPos = out_cactus->getCurrentPosition(t);

        if (wgt == 0) {
            LOGC(verbose) << "edge has zero weight, current cut remains same";
            return current_cut;
        }

        if (current_cut == 0) {
            LOGC(verbose) << "previously multiple CCs already, cut remains 0";
            return current_cut;
        }

        if (sCactusPos != tCactusPos) {
            LOGC(verbose) << "previously mincut between vertices, recompute";
            /* auto [cut, outg, b] = cactus.findAllMincuts(original_graph);
            out_cactus = outg;
            current_cut = cut;   */
            putIntoCache(out_cactus, current_cut);
            size_t fpid = flow_problem_id++;
            recursive_cactus<mutableGraphPtr> rc;
            push_relabel<false> pr;
            size_t flow = pr.solve_max_flow_min_cut(
                original_graph, { s, t }, 0, false, false, 0, fpid).first;
            auto new_g = rc.decrementalRebuild(original_graph, s, flow, fpid);
            current_cut = flow;
            out_cactus = new_g;
            LOG1 << "del1 t " << timer.elapsed();
        } else {
            push_relabel<true> pr;
            size_t fp = flow_problem_id++;
            auto [flow, sourceset] = pr.solve_max_flow_min_cut(
                original_graph, { s, t }, 0, false, false,
                current_cut, fp);
            if (static_cast<EdgeWeight>(flow) >= current_cut) {
                LOGC(verbose) << "cut not changed!";
            } else {
                putIntoCache(out_cactus, current_cut);
                recursive_cactus<mutableGraphPtr> rc;
                auto new_g = rc.decrementalRebuild(original_graph, s, flow, fp);
                current_cut = flow;
                out_cactus = new_g;
                LOG1 << "del2 t " << timer.elapsed();
                LOGC(verbose) << "recomputing, minimum cut changed to " << flow;
            }
        }
        LOGC(verbose) << "t " << timer.elapsed() << " cut " << current_cut;
        return current_cut;
    }

    mutableGraphPtr getOriginalGraph() {
        return original_graph;
    }

    mutableGraphPtr getCurrentCactus() {
        return out_cactus;
    }

    void putIntoCache(mutableGraphPtr cactusToCache, EdgeWeight cactusCut) {
        if (problemCache.size() < cactusCut + 1) {
            problemCache.resize(cactusCut + 1);
        }
        problemCache[cactusCut].first = cactusToCache;
        problemCache[cactusCut].second.clear();
        if (currentlyCaching.size() + 1 < cactusCut) {
            currentlyCaching.resize(cactusCut + 1, false);
        }
        currentlyCaching[cactusCut] = true;
    }

    void cacheEdge(NodeID s, NodeID t, EdgeWeight wgt) {
        EdgeWeight cut = current_cut;
        if (currentlyCaching.size() >= cut + 2 && currentlyCaching[cut + 1]) {
            if (problemCache[cut + 1].second.size() < max_cache_size) {
                problemCache[cut + 1].second.emplace_back(s, t, wgt);
            } else {
                currentlyCaching[cut + 1] = false;
            }
        }
    }
};
