/******************************************************************************
 * multiterminal_cut.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once

#include <limits>
#include <memory>
#include <string>
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
                           0, "") { }

    multicut_problem(std::shared_ptr<mutable_graph> G,
                     std::vector<terminal> term,
                     std::vector<std::shared_ptr<std::vector<NodeID> > >
                     mappings,
                     FlowType lower,
                     FlowType upper,
                     EdgeWeight deleted,
                     std::string path) : graph(G),
                                         terminals(term),
                                         mappings(mappings),
                                         lower_bound(lower),
                                         upper_bound(upper),
                                         deleted_weight(deleted),
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
    std::string                                         path;
};
