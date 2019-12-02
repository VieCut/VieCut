/******************************************************************************
 * multiterminal_cut.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "data_structure/mutable_graph.h"

struct terminal {
    terminal() { }

    terminal(NodeID position, NodeID original_id)
        : position(position), original_id(original_id), invalid_flow(true) { }

    terminal(NodeID position, NodeID original_id, bool invalid_flow)
        : position(position),
          original_id(original_id),
          invalid_flow(invalid_flow) { }

    NodeID position;
    NodeID original_id;
    bool   invalid_flow;
};

struct multicut_problem {
    multicut_problem() { }

    explicit multicut_problem(std::shared_ptr<mutable_graph> G)
        : multicut_problem(G, std::vector<terminal>()) { }

    multicut_problem(std::shared_ptr<mutable_graph> G,
                     std::vector<terminal> term)
        : multicut_problem(G,
                           term,
                           std::vector<
                               std::shared_ptr<std::vector<NodeID> > >(),
                           -1,
                           std::numeric_limits<FlowType>::max(),
                           0,
                           { UNDEFINED_NODE, UNDEFINED_EDGE },
                           "") { }

    multicut_problem(std::shared_ptr<mutable_graph> G,
                     std::vector<terminal> term,
                     std::vector<std::shared_ptr<std::vector<NodeID> > >
                     mappings,
                     FlowType lower,
                     FlowType upper,
                     EdgeWeight deleted,
                     std::pair<NodeID, EdgeID> prio,
                     std::string path) : graph(G),
                                         terminals(term),
                                         mappings(mappings),
                                         lower_bound(lower),
                                         upper_bound(upper),
                                         deleted_weight(deleted),
                                         priority_edge(prio),
                                         path(path) { }

    NodeID mapped(NodeID n) const {
        NodeID n_coarse = n;
        for (const auto& map : mappings) {
            n_coarse = (*map)[n_coarse];
        }
        return n_coarse;
    }

    std::shared_ptr<mutable_graph>                      graph;
    std::vector<terminal>                               terminals;
    std::vector<std::shared_ptr<std::vector<NodeID> > > mappings;
    FlowType                                            lower_bound;
    FlowType                                            upper_bound;
    EdgeWeight                                          deleted_weight;
    std::pair<NodeID, EdgeID>                           priority_edge;
    std::string                                         path;
};
