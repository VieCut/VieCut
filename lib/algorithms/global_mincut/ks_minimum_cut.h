/******************************************************************************
 * ks_minimum_cut.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 * 
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include "data_structure/graph_access.h"
#include "minimum_cut.h"
#include "noi_minimum_cut.h"
#include "tools/random_functions.h"
#include "tools/timer.h"

#ifdef PARALLEL
#include "parallel/coarsening/contract_graph.h"
#include "parallel/data_structure/union_find.h"
#else
#include "coarsening/contract_graph.h"
#include "data_structure/union_find.h"
#endif
#include "definitions.h"
#include "minimum_cut_helpers.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <memory>
#include <unordered_map>

class ks_minimum_cut : public minimum_cut
{
public:
    static const bool debug = false;
    static const bool timing = true;

    ks_minimum_cut() { }

    ~ks_minimum_cut() { }

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G) {
        return perform_minimum_cut(G, 0);
    }

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G, size_t optimal) {

        if (!minimum_cut_helpers::graphValid(G))
            return -1;

        EdgeWeight mincut = std::numeric_limits<EdgeWeight>::max();

        timer t;
        for (size_t i = 0; i < std::log2(G->number_of_nodes()) && mincut > optimal; ++i) {
            std::vector<std::shared_ptr<graph_access> > graphs;
            graphs.push_back(G);
            EdgeWeight curr_cut = recurse(graphs, 0, i);
            mincut = std::min(mincut, curr_cut);
            LOGC(timing) << "iter=" << i << " mincut=" << mincut << " curr_cut=" << curr_cut
                         << " time=" << t.elapsed();
        }
        return mincut;
    }

    bool unionNodes(union_find& uf, NodeID src, NodeID tgt) {
        NodeID first = uf.Find(src);
        NodeID second = uf.Find(tgt);

        if (first != second) {
            uf.Union(first, second);
            return true;
        }
        else {
            return false;
        }
    }

    std::pair<std::vector<NodeID>, std::vector<std::vector<NodeID> > >
    createMappings(graph_access& G,
                   union_find& uf) {

        std::vector<std::vector<NodeID> > reverse_mapping;
        std::vector<NodeID> mapping(G.number_of_nodes());
        std::vector<NodeID> part(G.number_of_nodes(), UNDEFINED_NODE);
        NodeID current_pid = 0;

        for (NodeID n : G.nodes()) {
            NodeID part_id = uf.Find(n);

            if (part[part_id] == UNDEFINED_NODE) {
                part[part_id] = current_pid++;
                reverse_mapping.emplace_back();
            }

            mapping[n] = part[part_id];
            //  G.setPartitionIndex(n, part[part_id]);
            reverse_mapping[part[part_id]].push_back(n);
        }

        return std::make_pair(mapping, reverse_mapping);
    }

    NodeID sample_contractible(graph_access& G,
                               NodeID currentN,
                               union_find& uf,
                               double reduction,
                               size_t iteration = 0) {
        NodeID n_reduce = std::min((NodeID)((double)currentN * reduction), currentN - 2);

        EdgeID num_edges = G.number_of_edges();

        std::mt19937_64 m_mt(iteration);

        size_t reduced = 0;

        while (reduced < n_reduce) {
            EdgeWeight e_rand = m_mt() % num_edges;
            NodeID src = G.getEdgeSource(e_rand);
            NodeID tgt = G.getEdgeTarget(e_rand);
            if (uf.Union(src, tgt))
                ++reduced;
        }

        return G.number_of_nodes() - reduced;
    }

    NodeID sample_contractible_weighted(graph_access& G,
                                        NodeID currentN,
                                        union_find& uf,
                                        double reduction,
                                        size_t iteration = 0) {

        NodeID n_reduce = std::min((NodeID)((double)currentN * reduction), currentN - 2);

        timer t;
        std::vector<EdgeWeight> prefixsum;
        prefixsum.reserve(G.number_of_edges());
        size_t wgt = 0;

        for (NodeID n : G.nodes()) {
            for (EdgeID e : G.edges_of(n)) {
                wgt += G.getEdgeWeight(e);
                prefixsum.emplace_back(wgt);
            }
        }

        EdgeWeight num_edges = wgt;

        size_t contracted = 0;
        std::mt19937_64 m_mt(iteration);

        while (contracted < n_reduce) {
            EdgeWeight e_rand = m_mt() % num_edges;

            auto edge = std::lower_bound(prefixsum.begin(), prefixsum.end(), e_rand,
                                         [](const auto& in1, const EdgeWeight& e) {
                                             return in1 < e;
                                         });

            EdgeID e = edge - prefixsum.begin();

            NodeID src = G.getEdgeSource(e);
            NodeID tgt = G.getEdgeTarget(e);

            if (unionNodes(uf, src, tgt)) {
                ++contracted;
            }
        }
        return currentN - contracted;
    }

    EdgeWeight recurse(std::vector<std::shared_ptr<graph_access> >& graphs,
                       size_t current, size_t iteration) {

        std::shared_ptr<graph_access> G = graphs[current];

        if (G->number_of_nodes() > 100) {

            NodeID currentN = G->number_of_nodes();

            union_find uf(G->number_of_nodes());

            if (current > 0) {
                for (NodeID n : G->nodes()) {
                    for (EdgeID e : G->edges_of(n)) {
                        if (G->getEdgeTarget(e) > n && currentN > 2) {
                            if (G->getEdgeWeight(e) > G->getMinDegree() ||
                                G->getEdgeWeight(e) * 2 > G->getWeightedNodeDegree(n) ||
                                G->getEdgeWeight(e) * 2 > G->getWeightedNodeDegree(G->getEdgeTarget(e))) {
                                NodeID first = uf.Find(n);
                                NodeID second = uf.Find(G->getEdgeTarget(e));

                                if (first != second) {
                                    currentN--;
                                    uf.Union(first, second);
                                }
                            }
                        }
                    }
                }
            }

            sample_contractible(*G, currentN, uf, 0.4, iteration);
            // sample_contractible_weighted(G, currentN, uf, 0.4, iteration);
            auto maps = createMappings(*G, uf);
            LOG << "Contracted to " << uf.n();

            std::shared_ptr<graph_access> G2 = contraction::contractGraph(G, maps.first,
                                                                          maps.second.size(),
                                                                          maps.second);

            graphs.push_back(G2);

            size_t currvecsize = graphs.size() - 1;

            return std::min(recurse(graphs, currvecsize, iteration),
                            recurse(graphs, currvecsize, iteration + 9273 /* pseudo-random seed */));
        }
        else {
            noi_minimum_cut mc;
            return mc.perform_minimum_cut(graphs[current]);
        }
    }
};
