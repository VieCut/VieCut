/******************************************************************************
 * parallel_cactus.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include "algorithms/global_mincut/cactus/most_balanced_minimum_cut.h"
#include "algorithms/global_mincut/cactus/recursive_cactus.h"
#include "algorithms/global_mincut/minimum_cut_helpers.h"
#include "algorithms/global_mincut/noi_minimum_cut.h"
#include "algorithms/global_mincut/viecut.h"
#include "coarsening/test_wrapper.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/bucket_pq.h"
#include "data_structure/priority_queues/fifo_node_bucket_pq.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "data_structure/priority_queues/node_bucket_pq.h"
#include "io/graph_io.h"
#include "tools/random_functions.h"
#include "tools/timer.h"

#ifdef PARALLEL
#include "parallel/algorithm/exact_parallel_minimum_cut.h"
#include "parallel/coarsening/contract_graph.h"
#include "parallel/coarsening/contraction_tests.h"
#include "parallel/coarsening/sparsify.h"
#include "parallel/data_structure/union_find.h"
#else
#include "coarsening/contract_graph.h"
#include "coarsening/contraction_tests.h"
#include "coarsening/sparsify.h"
#include "data_structure/union_find.h"
#endif

class parallel_cactus : public minimum_cut {
 public:
    parallel_cactus() { }

    ~parallel_cactus() { }

    static constexpr bool debug = false;
    bool timing = configuration::getConfig()->verbose;

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G) {
        if (!minimum_cut_helpers::graphValid(G))
            return -1;
        // compatibility with min cut interface
        return findAllMincuts(G).first;
    }

    std::pair<EdgeWeight, std::shared_ptr<mutable_graph> > findAllMincuts(
        std::shared_ptr<graph_access> G) {
        std::vector<std::shared_ptr<graph_access> > v = { G };
        return findAllMincuts(v);
    }

    std::pair<EdgeWeight, std::shared_ptr<mutable_graph> > findAllMincuts(
        std::vector<std::shared_ptr<graph_access> > graphs) {
        timer t;
        EdgeWeight mincut = graphs.back()->getMinDegree();
        recursive_cactus rc;
        exact_parallel_minimum_cut mc;
#ifdef PARALLEL
        viecut heuristic_mc;
        sparsify sf;
        auto G2 = graphs.back();
        if (configuration::getConfig()->contraction_factor > 0.0) {
            G2 = sf.one_ks(G2);
        }
        mincut = heuristic_mc.perform_minimum_cut(G2, true);
        LOGC(timing) << "VieCut found cut " << mincut
                     << " [Time: " << t.elapsed() << "s]";
#endif

        std::vector<std::vector<std::pair<NodeID, NodeID> > > guaranteed_edges;
        std::vector<size_t> ge_ids;

        // if PARALLEL is set, NodeInCut are already set to the result of viecut
        // This is what we want.
#ifndef PARALLEL
        minimum_cut_helpers::setInitialCutValues(graphs);
#endif

        NodeID previous_size = UNDEFINED_NODE;
        bool disable_blacklist = false;

        while (graphs.back()->number_of_nodes() < previous_size) {
            previous_size = graphs.back()->number_of_nodes();
            LOGC(timing) << "t " << t.elapsed() << " n "
                         << graphs.back()->number_of_nodes()
                         << " m " << graphs.back()->number_of_edges();
            std::shared_ptr<graph_access> curr_g = graphs.back();

#ifdef PARALLEL
            // all runs after first disable blacklist so that every thread
            // runs capforest on the whole graph
            auto uf = mc.parallel_modified_capforest(curr_g, mincut + 1,
                                                     disable_blacklist);
            disable_blacklist = true;
#else
            LOG1 << "Error: Running exact_parallel_minimum_cut without PARALLEL"
                 << "Using normal noi_minimum_cut instead!";

            noi_minimum_cut noi;
            auto uf = noi.modified_capforest(curr_g, mincut + 1);
#endif

            ge_ids.emplace_back(graphs.size() - 1);
            guaranteed_edges.emplace_back();
            EdgeWeight current_mincut = mincut;

            for (NodeID n : curr_g->nodes()) {
                if (curr_g->getNodeDegree(n) == 1) {
                    EdgeID e = curr_g->get_first_edge(n);
                    if ((curr_g->getEdgeWeight(e) == mincut)
                        && (guaranteed_edges.back().size() + 1 < n)) {
                        NodeID t = curr_g->getEdgeTarget(e);
                        uf.Union(n, t);
                        guaranteed_edges.back().emplace_back(n, t);
                    }
                }
            }
            LOGC(timing) << "t " << t.elapsed() << " contract "
                         << graphs.back()->number_of_nodes()
                         << " to " << uf.n();
            if (uf.n() < curr_g->number_of_nodes()) {
                graphs.push_back(contraction::fromUnionFind(curr_g, &uf));
                mincut = minimum_cut_helpers::updateCut(graphs, mincut);
            }

            auto uf12 = tests::prTests12(graphs.back(), mincut + 1, true);
            LOGC(timing) << "t12 " << t.elapsed() << " contract "
                         << graphs.back()->number_of_nodes()
                         << " to " << uf12.n();
            if (uf12.n() < graphs.back()->number_of_nodes()) {
                auto g12 = contraction::fromUnionFind(graphs.back(), &uf12);
                graphs.push_back(g12);
                mincut = minimum_cut_helpers::updateCut(graphs, mincut);
            }

            auto uf34 = tests::prTests34(graphs.back(), mincut + 1, true);
            LOGC(timing) << "t34 " << t.elapsed() << " contract "
                         << graphs.back()->number_of_nodes()
                         << " to " << uf34.n();
            if (uf34.n() < graphs.back()->number_of_nodes()) {
                auto g34 = contraction::fromUnionFind(graphs.back(), &uf34);
                graphs.push_back(g34);
                mincut = minimum_cut_helpers::updateCut(graphs, mincut);
            }

            if (current_mincut > mincut) {
                // mincut has improved, so all the edges that
                // were in guaranteed edges
                // before are now _not_ cactus edges as there is a lighter cut
                guaranteed_edges.clear();
                ge_ids.clear();
            }
        }

        // check whether there is a small cut hidden in graphs.back()
        if (graphs.back()->number_of_nodes() > 1) {
            LOGC(timing) << "Searching for minimum cut in remaining graph...";
            mincut = std::min(mincut,
                              mc.perform_minimum_cut(graphs.back(), true));
        }

        rc.setMincut(mincut);
        auto out_graph = rc.flowMincut(graphs);

        minimum_cut_helpers::setVertexLocations(out_graph, graphs, ge_ids,
                                                guaranteed_edges, mincut);

        LOGC(timing) << "t " << t.elapsed() << " unpacked - n "
                     << out_graph->n() << " m " << out_graph->m();

        if (configuration::getConfig()->find_most_balanced_cut) {
            most_balanced_minimum_cut mbmc;
            mbmc.findCutFromCactus(out_graph, mincut, graphs[0]);
        }

        return std::make_pair(mincut, out_graph);
    }
};
