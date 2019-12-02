/******************************************************************************
 * edge_selection.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 */

#pragma once

#include <algorithm>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/multicut/multicut_problem.h"
#include "data_structure/mutable_graph.h"
#include "tools/random_functions.h"

[[maybe_unused]] static std::tuple<NodeID, EdgeID> findHeavyEdge(
    std::shared_ptr<multicut_problem> problem) {
    EdgeWeight max_wgt = 0;
    NodeID max_id_t = 0;
    EdgeID max_id_e = 0;
    for (const auto& term : problem->terminals) {
        NodeID t = term.position;

        for (EdgeID e : problem->graph->edges_of(t)) {
            EdgeWeight curr_wgt = problem->graph->getEdgeWeight(t, e);
            if (curr_wgt > max_wgt) {
                max_wgt = curr_wgt;
                max_id_t = t;
                max_id_e = e;
            }
        }
    }
    return std::make_tuple(max_id_t, max_id_e);
}

[[maybe_unused]] static std::tuple<NodeID, EdgeID> findVertexWithHighDistance(
    std::shared_ptr<multicut_problem> problem) {
    EdgeWeight max_wgt = 0;
    NodeID max_id_t = 0;
    EdgeID max_id_e = 0;
    for (const auto& term : problem->terminals) {
        NodeID t = term.position;
        for (EdgeID e : problem->graph->edges_of(t)) {
            NodeID tgt = problem->graph->getEdgeTarget(t, e);
            EdgeWeight curr_wgt = problem->graph->getWeightedNodeDegree(tgt)
                                  - (problem->graph->getEdgeWeight(t, e));
            if (curr_wgt > max_wgt) {
                max_wgt = curr_wgt;
                max_id_t = t;
                max_id_e = e;
            }
        }
    }
    return std::make_tuple(max_id_t, max_id_e);
}

[[maybe_unused]] static std::tuple<NodeID, EdgeID> findHeavyVertex(
    std::shared_ptr<multicut_problem> problem) {
    EdgeWeight max_wgt = 0;
    NodeID max_id_t = 0;
    EdgeID max_id_e = 0;
    for (const auto& term : problem->terminals) {
        NodeID t = term.position;

        if (problem->graph->get_first_invalid_edge(t) == 1) {
            // return std::make_tuple(t, 0);
        }

        for (EdgeID e : problem->graph->edges_of(t)) {
            NodeID tgt = problem->graph->getEdgeTarget(t, e);
            EdgeWeight curr_wgt = problem->graph->getWeightedNodeDegree(tgt);
            if (curr_wgt > max_wgt) {
                max_wgt = curr_wgt;
                max_id_t = t;
                max_id_e = e;
            }
        }
    }
    return std::make_tuple(max_id_t, max_id_e);
}

[[maybe_unused]] static std::tuple<NodeID, EdgeID> findGlobalHeaviestEdge(
    std::shared_ptr<multicut_problem> problem) {
    EdgeWeight max_wgt = 0;
    NodeID max_id_t = 0;
    EdgeID max_id_e = 0;
    for (NodeID n : problem->graph->nodes()) {
        for (EdgeID e : problem->graph->edges_of(n)) {
            EdgeWeight curr_wgt = problem->graph->getEdgeWeight(n, e);
            if (curr_wgt > max_wgt) {
                max_wgt = curr_wgt;
                max_id_t = n;
                max_id_e = e;
            }
        }
    }
    return std::make_tuple(max_id_t, max_id_e);
}

[[maybe_unused]] static std::tuple<NodeID, EdgeID> findEdgeWithManyTerminals(
    const std::shared_ptr<multicut_problem> problem) {
    std::unordered_map<NodeID, uint16_t> terms;
    NodeID max_num = 0;
    NodeID max_ngbr = 0;
    std::unordered_set<NodeID> termset;
    for (const auto& term : problem->terminals) {
        NodeID t = term.position;
        termset.emplace(t);
        for (EdgeID e : problem->graph->edges_of(t)) {
            EdgeWeight wgt = problem->graph->getEdgeWeight(t, e);
            NodeID ngbr = problem->graph->getEdgeTarget(t, e);
            if (terms.find(ngbr) == terms.end()) {
                terms[ngbr] = wgt;
            } else {
                terms[ngbr] += wgt;
            }
            if (terms[ngbr] > max_num) {
                max_ngbr = ngbr;
                max_num = terms[ngbr];
            }
        }
    }

    EdgeWeight max_wgt = 0;
    NodeID max_id_t = 0;
    EdgeID max_id_e = 0;

    for (EdgeID e : problem->graph->edges_of(max_ngbr)) {
        EdgeWeight curr_wgt = problem->graph->getEdgeWeight(max_ngbr, e);
        NodeID ngbr = problem->graph->getEdgeTarget(max_ngbr, e);
        if (curr_wgt > max_wgt && termset.count(ngbr) > 0) {
            max_wgt = curr_wgt;
            max_id_t = ngbr;
            max_id_e = problem->graph->getReverseEdge(max_ngbr, e);
        }
    }
    return std::make_tuple(max_id_t, max_id_e);
}

[[maybe_unused]] static std::tuple<NodeID, EdgeID> findRandomEdge(
    std::shared_ptr<multicut_problem> problem) {
    std::vector<EdgeWeight> prefixsum(problem->terminals.size() + 1, 0);
    EdgeWeight s = 0;
    for (size_t i = 0; i < prefixsum.size() - 1; ++i) {
        s += problem->graph->getNodeDegree(problem->terminals[i].position);
        prefixsum[i + 1] = s;
    }

    EdgeID random_edge = random_functions::nextInt(0, s - 1);
    NodeID t = 0;

    while (random_edge >= prefixsum[t + 1]) {
        t++;
    }

    EdgeID e = random_edge - prefixsum[t];

    return std::make_tuple(problem->terminals[t].position, e);
}

static std::tuple<NodeID, EdgeID> findEdge(
    std::shared_ptr<multicut_problem> problem,
    const std::string& edge_selection) {
    std::pair<NodeID, EdgeID> undefined = { UNDEFINED_NODE, UNDEFINED_EDGE };
    if (problem->priority_edge != undefined) {
        auto p = problem->priority_edge;
        problem->priority_edge = undefined;
        return p;
    }

    if (edge_selection == "random")
        return findRandomEdge(problem);

    if (edge_selection == "heavy")
        return findHeavyEdge(problem);

    if (edge_selection == "heavy_global")
        return findGlobalHeaviestEdge(problem);

    if (edge_selection == "connection")
        return findEdgeWithManyTerminals(problem);

    if (edge_selection == "distance")
        return findVertexWithHighDistance(problem);

    if (edge_selection == "heavy_vertex")
        return findHeavyVertex(problem);

    return findHeavyVertex(problem);
}
