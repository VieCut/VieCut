/******************************************************************************
 * sparsify.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
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
#include "omp.h"

#include "data_structure/graph_access.h"
#include "parallel/coarsening/contract_graph.h"
#include "parallel/data_structure/union_find.h"
#include "tlx/logger.hpp"
#include "tools/random_functions.h"
#include "tools/timer.h"

class sparsify
{

private:
    auto createMappings(graph_access& G,
                        union_find& uf) {

        timer t;
        // std::vector<std::vector<NodeID>> reverse_mapping;
        std::vector<NodeID> mapping(G.number_of_nodes());
        std::vector<NodeID> part(G.number_of_nodes(), UNDEFINED_NODE);
        NodeID current_pid = 0;
// #pragma omp parallel for schedule(dynamic)
        for (size_t n = 0; n < G.number_of_nodes(); ++n) {
            NodeID part_id = uf.Find(n);

            if (part[part_id] == UNDEFINED_NODE) {
                part[part_id] = current_pid++;
            }

            mapping[n] = part[part_id];
        }

        // reverse_mapping.resize(current_pid);

        return std::make_tuple(mapping, current_pid);
    }

public:
    NodeID sample_contractible_weighted(graph_access& G,
                                        NodeID currentN,
                                        union_find& uf,
                                        double reduction,
                                        size_t iteration = 0) {
        NodeID n_reduce;
        timer t;
        std::vector<std::vector<EdgeWeight> > prefixsum;

        std::vector<EdgeWeight> processor_prefixes;

#pragma omp parallel
        {
            if (!omp_get_thread_num()) {
                prefixsum.resize(omp_get_num_threads());
                processor_prefixes.resize(omp_get_num_threads());
            }
#pragma omp barrier

            EdgeID m = G.number_of_edges();
            size_t wgt = 0;
            EdgeID per_thread = std::ceil((double)m / (double)omp_get_num_threads());

            auto id = omp_get_thread_num();
            prefixsum[id].reserve(per_thread);

            EdgeID my_start = id * per_thread;
            EdgeID my_end = std::min((id + 1) * per_thread, G.number_of_edges());

            for (size_t m = my_start; m < my_end; ++m) {
                wgt += G.getEdgeWeight(m);
                prefixsum[id].push_back(wgt);
            }

            processor_prefixes[id] = wgt;
#pragma omp barrier

            EdgeWeight my_prefix = 0;
            for (int i = 0; i < id; ++i) {
                my_prefix += processor_prefixes[i];
            }

            for (size_t m = my_start; m < my_end; ++m) {
                prefixsum[m / per_thread][m % per_thread] += my_prefix;
            }

            n_reduce = std::min((NodeID)((double)currentN * reduction), currentN - 2);
            size_t contracted = 0;
            std::mt19937_64 m_mt(iteration);
            n_reduce /= omp_get_num_threads();

#pragma omp barrier
            while (contracted < n_reduce) {
                auto random = m_mt();
                EdgeWeight e_rand = (random % processor_prefixes[id]) + my_prefix;

                auto edge = std::lower_bound(prefixsum[id].begin(), prefixsum[id].end(), e_rand,
                                             [](const auto& in1, const EdgeWeight& e) {
                                                 return in1 < e;
                                             });

                EdgeID e = edge - prefixsum[id].begin() + (per_thread * id);

                NodeID src = G.getEdgeSource(e);
                NodeID tgt = G.getEdgeTarget(e);

                if (!uf.SameSet(src, tgt)) {
                    if (uf.Union(src, tgt)) {
                        ++contracted;
                    }
                }
            }
        }

        return G.number_of_nodes() - n_reduce * omp_get_num_threads();
    }

    NodeID sample_contractible(graph_access& G,
                               NodeID currentN,
                               union_find& uf,
                               double reduction,
                               size_t iteration = 0) {

        NodeID n_reduce;

#pragma omp parallel
        {
            n_reduce = (NodeID)((double)currentN * reduction);
            std::mt19937_64 m_mt(iteration + omp_get_thread_num());
            EdgeID num_edges = G.number_of_edges();

            size_t my_n = 0;

            n_reduce = n_reduce / omp_get_num_threads();

            while (my_n < n_reduce) {
                // while (my_tries > num_tries++) {
                EdgeWeight e_rand = m_mt() % num_edges; // my_range + my_start;

                NodeID src = G.getEdgeSource(e_rand);
                NodeID tgt = G.getEdgeTarget(e_rand);

                if (!uf.SameSet(src, tgt)) {
                    if (uf.Union(src, tgt)) {
                        ++my_n;
                    }
                }
            }
        }

        return G.number_of_nodes() - n_reduce * omp_get_num_threads();
    }

    std::shared_ptr<graph_access> one_ks(
        std::shared_ptr<graph_access> G_in, double contraction, size_t iteration) {

        union_find uf(G_in->number_of_nodes());

        timer t;

        constexpr bool weighted = false;

        if (weighted) {
            sample_contractible_weighted(
                (*G_in), G_in->number_of_nodes(), uf, contraction, iteration);
        }
        else {
            sample_contractible(
                (*G_in), G_in->number_of_nodes(), uf, contraction, iteration);
        }

        auto [map, rev_map] = createMappings((*G_in), uf);

        std::shared_ptr<graph_access> G2 = contraction::contractGraph(G_in, map, rev_map);

        return G2;
    }

    std::shared_ptr<graph_access> random_matching(
        std::shared_ptr<graph_access> G_in) {

        std::shared_ptr<graph_access> G_out = std::make_shared<graph_access>();

        const size_t MAX_TRIES = 50;
        size_t no_of_coarse_vertices = 0;
        std::vector<NodeID> edge_matching(G_in->number_of_nodes());

        for (NodeID n : G_in->nodes()) {
            edge_matching[n] = n;
        }

        std::vector<NodeID> coarse_mapping(G_in->number_of_nodes());

        for (NodeID n : G_in->nodes()) {
            if (edge_matching[n] == n) {

                size_t no_try = 0;
                NodeID matchingPartner = n;
                while (no_try < MAX_TRIES) {

                    // match with a random neighbor
                    EdgeID s = G_in->get_first_edge(n);
                    EdgeID modulo = G_in->getNodeDegree(n);
                    EdgeID r = s + (random_functions::next() % modulo);

                    NodeID tgt = G_in->getEdgeTarget(r);

                    if (edge_matching[tgt] == tgt) {
                        matchingPartner = tgt;
                        break;
                    }

                    no_try++;
                }

                if (no_try == MAX_TRIES) {
                    coarse_mapping[n] = no_of_coarse_vertices;
                    edge_matching[n] = n;
                    no_of_coarse_vertices++;
                }
                else {
                    coarse_mapping[matchingPartner] = no_of_coarse_vertices;
                    coarse_mapping[n] = no_of_coarse_vertices;

                    edge_matching[matchingPartner] = n;
                    edge_matching[n] = matchingPartner;

                    no_of_coarse_vertices++;
                }
            }
        }

        G_out->start_construction(no_of_coarse_vertices,
                                  G_in->number_of_edges());

        std::vector<std::pair<NodeID, NodeID> > edge_positions(
            no_of_coarse_vertices,
            std::make_pair<NodeID, NodeID>(-1, -1));

        NodeID cur_no_vertices = 0;

        size_t edgerino = 0;

        for (NodeID n : G_in->nodes()) {
            // we look only at the coarser nodes
            if (coarse_mapping[n] < cur_no_vertices)
                continue;

            NodeID coarseNode = G_out->new_node();

            // do something with all outgoing edges (in auxillary graph)
            for (EdgeID e : G_in->edges_of(n)) {
                edgerino++;
                EdgeID new_coarse_edge_target = coarse_mapping[
                    G_in->getEdgeTarget(e)];

                if (new_coarse_edge_target == coarseNode) continue;

                auto pos = edge_positions[new_coarse_edge_target];

                if (pos.first != coarseNode) {
                    EdgeID coarseEdge = G_out->new_edge(coarseNode, new_coarse_edge_target);
                    G_out->setEdgeWeight(coarseEdge, G_in->getEdgeWeight(e));
                    edge_positions[new_coarse_edge_target] =
                        std::make_pair(coarseNode, coarseEdge);
                }
                else {
                    EdgeWeight new_weight = G_out->getEdgeWeight(pos.second)
                                            + G_in->getEdgeWeight(e);

                    G_out->setEdgeWeight(pos.second, new_weight);
                }
            }

            // this node was really matched
            NodeID matched_neighbor = edge_matching[n];
            if (n != matched_neighbor) {

                for (EdgeID e : G_in->edges_of(matched_neighbor)) {
                    EdgeID new_coarse_edge_target = coarse_mapping[
                        G_in->getEdgeTarget(e)];

                    if (new_coarse_edge_target == coarseNode) continue;

                    auto pos = edge_positions[new_coarse_edge_target];

                    if (pos.first != coarseNode) {
                        EdgeID coarseEdge = G_out->new_edge(
                            coarseNode, new_coarse_edge_target);
                        G_out->setEdgeWeight(coarseEdge, G_in->getEdgeWeight(e));
                        edge_positions[new_coarse_edge_target] =
                            std::make_pair(coarseNode, coarseEdge);
                    }
                    else {
                        EdgeWeight new_weight = G_out->getEdgeWeight(pos.second)
                                                + G_in->getEdgeWeight(e);

                        G_out->setEdgeWeight(pos.second, new_weight);
                    }
                }
            }

            cur_no_vertices++;
        }

        // this also resizes the edge fields ...
        G_out->finish_construction();

        return G_out;
    }

    std::shared_ptr<graph_access> remove_heavy_vertices(
        std::shared_ptr<graph_access> G_in, double percentile,
        std::vector<NodeID>& prefixsum) {

        std::shared_ptr<graph_access> G_out = std::make_shared<graph_access>();

        EdgeWeight bound_deg = G_in->getPercentile(percentile);

        // new vertex ids in sparser graph
        prefixsum.resize(G_in->number_of_nodes());

        size_t ctr = 0;
        // this is an upper bound, as edges to dense vectors are also counted.
        // we can just resize the edge vector afterwards though, no need to check
        // each edge (todo: if slow, this needs to be changed for parallelism)
        size_t existing_edges = 0;

        for (NodeID n : G_in->nodes()) {
            if (G_in->getWeightedNodeDegree(n) < bound_deg) {
                prefixsum[n] = ctr++;
                existing_edges += G_in->getNodeDegree(n);
            }
            else {
                prefixsum[n] = G_in->number_of_nodes();
            }
        }

        G_out->start_construction(ctr, existing_edges);

        size_t actually_edges = 0;

        for (NodeID n : G_in->nodes()) {
            size_t pre = prefixsum[n];
            if (pre < G_in->number_of_nodes()) {
                G_out->new_node();
                for (EdgeID e : G_in->edges_of(n)) {
                    NodeID target = G_in->getEdgeTarget(e);
                    size_t pretgt = prefixsum[target];
                    if (prefixsum[target] < G_in->number_of_nodes()) {
                        G_out->new_edge(pre, pretgt);
                        actually_edges++;
                    }
                }
            }
        }

        G_out->finish_construction();

        return G_out;
    }
};
