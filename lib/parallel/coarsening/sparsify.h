/******************************************************************************
 * sparsify.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include "data_structure/graph_access.h"
#include "parallel/coarsening/contract_graph.h"
#include "parallel/data_structure/union_find.h"
#include "tlx/logger.hpp"
#include "tools/random_functions.h"
#include "tools/timer.h"

class sparsify {
 private:
    NodeID elementsToReduce(std::shared_ptr<graph_access> G) {
        return static_cast<NodeID>(
            static_cast<double>(G->number_of_nodes())
            * configuration::getConfig()->contraction_factor);
    }

    auto createMappings(std::shared_ptr<graph_access> G,
                        union_find* uf) {
        timer t;
        std::vector<NodeID> mapping(G->number_of_nodes());
        std::vector<NodeID> part(G->number_of_nodes(), UNDEFINED_NODE);
        NodeID current_pid = 0;
        for (size_t n = 0; n < G->number_of_nodes(); ++n) {
            NodeID part_id = uf->Find(n);

            if (part[part_id] == UNDEFINED_NODE) {
                part[part_id] = current_pid++;
            }

            mapping[n] = part[part_id];
        }
        return std::make_tuple(mapping, current_pid);
    }

 public:
    NodeID sample_contractible_weighted(std::shared_ptr<graph_access> G,
                                        union_find* uf) {
        NodeID n_reduce;
        timer t;
        size_t num_threads = omp_get_num_threads();
        std::vector<std::vector<EdgeWeight> > prefixsum;

        std::vector<EdgeWeight> processor_prefixes;

#pragma omp parallel
        {
            if (!omp_get_thread_num()) {
                prefixsum.resize(num_threads);
                processor_prefixes.resize(num_threads);
            }
#pragma omp barrier

            EdgeID m = G->number_of_edges();
            size_t wgt = 0;
            EdgeID per_thread =
                std::ceil(static_cast<double>(m)
                          / static_cast<double>(num_threads));

            auto id = omp_get_thread_num();
            prefixsum[id].reserve(per_thread);

            EdgeID my_start = id * per_thread;
            EdgeID my_end = std::min((id + 1) * per_thread,
                                     G->number_of_edges());

            for (size_t m = my_start; m < my_end; ++m) {
                wgt += G->getEdgeWeight(m);
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

            double factor = configuration::getConfig()->contraction_factor;
            n_reduce = std::min(static_cast<NodeID>(
                                    static_cast<double>(G->number_of_nodes())
                                    * factor),
                                G->number_of_nodes() - 2);
            size_t contracted = 0;
            std::mt19937_64 m_mt(configuration::getConfig()->seed);
            n_reduce /= num_threads;

#pragma omp barrier
            while (contracted < n_reduce) {
                auto random = m_mt();
                EdgeWeight e_rand = (random % processor_prefixes[id])
                                    + my_prefix;

                auto edge = std::lower_bound(
                    prefixsum[id].begin(), prefixsum[id].end(), e_rand,
                    [](const auto& in1, const EdgeWeight& e) {
                        return in1 < e;
                    });

                EdgeID e = edge - prefixsum[id].begin() + (per_thread * id);

                NodeID src = G->getEdgeSource(e);
                NodeID tgt = G->getEdgeTarget(e);

                if (!uf->SameSet(src, tgt)) {
                    if (uf->Union(src, tgt)) {
                        ++contracted;
                    }
                }
            }
        }

        return G->number_of_nodes() - n_reduce * num_threads;
    }

    NodeID sample_contractible_separate(std::shared_ptr<graph_access> G,
                                        union_find* uf) {
        NodeID to_reduce = elementsToReduce(G);
        std::atomic<NodeID> reduced = 0;
        timer t;
        NodeID to_try = to_reduce / configuration::getConfig()->threads;
        const EdgeID my_range = G->number_of_edges()
                                / configuration::getConfig()->threads;

#pragma omp parallel
        {
            std::mt19937_64 m_mt(configuration::getConfig()->seed
                                 + omp_get_thread_num());
            EdgeID my_start = my_range * omp_get_thread_num();
            size_t tries = 0;
            NodeID my_reduced = 0;
            while (tries < to_try) {
                EdgeWeight e_rand = (m_mt() % my_range) + my_start;
                tries++;

                NodeID src = G->getEdgeSource(e_rand);
                NodeID tgt = G->getEdgeTarget(e_rand);

                if (uf->Union(src, tgt)) {
                    ++my_reduced;
                }
            }

            LOG1 << tries << " tries in " << t.elapsed()
                 << "s (success in " << my_reduced << ")";
        }

        return G->number_of_nodes() - reduced;
    }

    NodeID sample_geometric(std::shared_ptr<graph_access> G,
                            union_find* uf) {
        timer t;
    #pragma omp parallel
        {
            double n = G->number_of_nodes();
            double m = G->number_of_edges();
            double samples = n * configuration::getConfig()->contraction_factor;

            // as we do not want to sample an edge multiple times,
            // we sample without replacement
            // this is O(realized) by removing 'samples'
            // from the sampled range and adding +1 to e after every sample
            double prob =
                samples / (m * configuration::getConfig()->threads - samples);

            if (prob > 1) {
                LOG1 << "WARNING: Edge Sampling Probability larger than 1.";
                LOG1 << "Setting it to a valid value of 0.5.";
                prob = 0.5;
            }

            std::default_random_engine generator(
                configuration::getConfig()->seed + omp_get_thread_num());
            std::geometric_distribution<EdgeID> distribution(prob);

            EdgeID prev = 0;
            EdgeID e = prev;
            NodeID src = 0;
            while (e < G->number_of_edges()) {
                while (e >= G->get_first_edge(src+1)) {
                    ++src;
                }

                const NodeID tgt = G->getEdgeTarget(e);
                uf->Union(src, tgt);

                e += (distribution(generator) + 1);
            }
        }
        return uf->n();
    }

    NodeID sample_contractible(std::shared_ptr<graph_access> G,
                               union_find* uf) {
        timer t;
        NodeID to_reduce = elementsToReduce(G);
        NodeID to_try = to_reduce / configuration::getConfig()->threads;

#pragma omp parallel
        {
            NodeID n_reduce = elementsToReduce(G);
            std::mt19937_64 m_mt(configuration::getConfig()->seed
                                 + omp_get_thread_num());
            EdgeID num_edges = G->number_of_edges();
            size_t my_n = 0;
            size_t tries = 0;
            n_reduce = n_reduce / configuration::getConfig()->threads;

            while (tries < to_try) {
                EdgeWeight e_rand = m_mt() % num_edges;
                tries++;

                NodeID src = G->getEdgeSource(e_rand);
                NodeID tgt = G->getEdgeTarget(e_rand);

                if (!uf->SameSet(src, tgt)) {
                    if (uf->Union(src, tgt)) {
                        ++my_n;
                    }
                }
            }

            LOG1 << tries << " tries in "
                 << t.elapsed() << "s (success in " << my_n << ")";
        }

        return uf->n();
    }

    std::shared_ptr<graph_access> one_ks(
        std::shared_ptr<graph_access> G_in) {
        union_find uf(G_in->number_of_nodes());
        timer t;

        if (configuration::getConfig()->sampling_type == "geometric")
            sample_geometric(G_in, &uf);

        if (configuration::getConfig()->sampling_type == "random")
            sample_contractible(G_in, &uf);

        if (configuration::getConfig()->sampling_type == "separate")
            sample_contractible_separate(G_in, &uf);

        if (configuration::getConfig()->sampling_type == "weighted")
            sample_contractible_weighted(G_in, &uf);

        LOG1 << "t " << t.elapsed() << " sample";

        auto [map, rev_map] = createMappings(G_in, &uf);

        std::shared_ptr<graph_access> G2 =
            contraction::contractGraph(G_in, map, rev_map);

        LOG1 << "t " << t.elapsed() << " for contraction from "
             << G_in->number_of_nodes() << " to " << G2->number_of_nodes();
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
                } else {
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
            std::make_pair(-1, -1));

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
                    EdgeID coarseEdge =
                        G_out->new_edge(coarseNode, new_coarse_edge_target);
                    G_out->setEdgeWeight(coarseEdge, G_in->getEdgeWeight(e));
                    edge_positions[new_coarse_edge_target] =
                        std::make_pair(coarseNode, coarseEdge);
                } else {
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
                        G_out->setEdgeWeight(coarseEdge,
                                             G_in->getEdgeWeight(e));
                        edge_positions[new_coarse_edge_target] =
                            std::make_pair(coarseNode, coarseEdge);
                    } else {
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
        std::vector<NodeID>* ps) {
        std::vector<NodeID>& prefixsum = *ps;

        std::shared_ptr<graph_access> G_out = std::make_shared<graph_access>();

        EdgeWeight bound_deg = G_in->getPercentile(percentile);

        // new vertex ids in sparser graph
        prefixsum.resize(G_in->number_of_nodes());

        size_t ctr = 0;
        // this is an upper bound, as edges to dense vectors are also counted.
        // we can resize the edge vector afterwards though, no need to check
        // each edge (todo: if slow, this needs to be changed for parallelism)
        size_t existing_edges = 0;

        for (NodeID n : G_in->nodes()) {
            if (G_in->getWeightedNodeDegree(n) < bound_deg) {
                prefixsum[n] = ctr++;
                existing_edges += G_in->getNodeDegree(n);
            } else {
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
