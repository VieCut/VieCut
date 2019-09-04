/******************************************************************************
 * cactus_mincut.h
 *
 * Source of VieCut
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
#include "algorithms/global_mincut/minimum_cut.h"
#include "algorithms/global_mincut/noi_minimum_cut.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "io/graph_io.h"
#include "tools/string.h"
#include "tools/timer.h"

#ifdef PARALLEL
#include "parallel/coarsening/contract_graph.h"
#else
#include "coarsening/contract_graph.h"
#endif

class cactus_mincut : public minimum_cut {
 public:
    cactus_mincut() { }
    virtual ~cactus_mincut() { }
    static constexpr bool debug = false;

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G) {
        if (!minimum_cut_helpers::graphValid(G))
            return -1;
        // compatibility with min cut interface
        return findAllMincuts(G).first;
    }

    std::pair<EdgeWeight, std::shared_ptr<mutable_graph> > findAllMincuts(
        std::shared_ptr<graph_access> G) {
        std::vector<std::shared_ptr<graph_access> > graphs;
        recursive_cactus rc;
        EdgeWeight mincut = G->getMinDegree();
        graphs.push_back(G);
        noi_minimum_cut noi;
        timer t;

        std::vector<std::vector<std::pair<NodeID, NodeID> > > guaranteed_edges;
        std::vector<size_t> ge_ids;

        minimum_cut_helpers::setInitialCutValues(graphs);

        NodeID previous_size = UNDEFINED_NODE;
        while (graphs.back()->number_of_nodes() < previous_size) {
            previous_size = graphs.back()->number_of_nodes();
            auto current_graph = graphs.back();

            EdgeWeight current_mincut = mincut;

            ge_ids.emplace_back(graphs.size() - 1);
            guaranteed_edges.emplace_back();

            LOG1 << "t " << t.elapsed()
                 << " n " << graphs.back()->number_of_nodes()
                 << " m " << current_graph->number_of_edges()
                 << " cut " << mincut;

            std::vector<std::pair<NodeID, NodeID> > contractable;

            timer time;
            auto uf = noi.modified_capforest(current_graph, mincut + 1);

            for (NodeID n : current_graph->nodes()) {
                EdgeID e = current_graph->get_first_edge(n);
                if (current_graph->get_first_invalid_edge(n) - e == 1) {
                    if ((current_graph->getEdgeWeight(e) == mincut)
                        && uf.n() > 1) {
                        NodeID t = current_graph->getEdgeTarget(e);
                        uf.Union(n, t);
                        guaranteed_edges.back().emplace_back(n, t);
                    }
                }
            }

            auto newg = contraction::fromUnionFind(current_graph, &uf);
            if (newg->number_of_nodes() < current_graph->number_of_nodes())
                graphs.emplace_back(newg);
            mincut = minimum_cut_helpers::updateCut(graphs, mincut);

            union_find uf12 = tests::prTests12(graphs.back(), mincut + 1, true);
            auto g12 = contraction::fromUnionFind(graphs.back(), &uf12);
            if (g12->number_of_nodes() < newg->number_of_nodes())
                graphs.push_back(g12);
            mincut = minimum_cut_helpers::updateCut(graphs, mincut);

            union_find uf34 = tests::prTests34(graphs.back(), mincut + 1, true);
            auto g34 = contraction::fromUnionFind(graphs.back(), &uf34);
            if (g34->number_of_nodes() < g12->number_of_nodes())
                graphs.push_back(g34);
            mincut = minimum_cut_helpers::updateCut(graphs, mincut);

            if (current_mincut > mincut) {
                // mincut has improved
                // so all the edges that were in guaranteed edges
                // before are now _not_ cactus edges as there is a lighter cut
                guaranteed_edges.clear();
                ge_ids.clear();
            }
        }

        if (graphs.back()->number_of_nodes() > 1)
            mincut = noi.perform_minimum_cut(graphs.back());

        rc.setMincut(mincut);
        auto out_graph = rc.flowMincut(graphs);

        auto [out_graph_mapping, deleted_vertex_mappings] =
            minimum_cut_helpers::reInsertVertices(
                out_graph, graphs, ge_ids, guaranteed_edges, mincut);

        if (configuration::getConfig()->save_cut) {
            minimum_cut_helpers::setVertexLocations(
                out_graph, graphs, out_graph_mapping, deleted_vertex_mappings);
        }

        LOG1 << "unpacked - n " << out_graph->n() << " m " << out_graph->m();

        if (configuration::getConfig()->find_most_balanced_cut) {
            most_balanced_minimum_cut mbmc;
            mbmc.findCutFromCactus(out_graph, mincut, G);
        }

        return std::make_pair(mincut, out_graph);
    }
};
