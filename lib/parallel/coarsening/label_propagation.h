/******************************************************************************
 * label_propagation.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include "omp.h"

#include "data_structure/graph_access.h"
#include "definitions.h"
#include "tlx/logger.hpp"
#include "tools/random_functions.h"
#include "tools/string.h"

#include <cstdint>
#include <cstdlib>

#include <algorithm>
#include <deque>
#include <map>
#include <random>
#include <unordered_map>

class label_propagation
{
    static constexpr bool debug = false;
    static constexpr bool timing = true;

public:
    label_propagation() { }
    virtual ~label_propagation() { }

    std::vector<NodeID> propagate_labels(std::shared_ptr<graph_access> G) {

        std::vector<NodeID> cluster_mapping;

        cluster_mapping.resize(G->number_of_nodes());
        std::vector<NodeID> permutation(G->number_of_nodes());

        random_functions::permutate_vector_local(permutation, true);

        for (size_t i = 0; i < cluster_mapping.size(); ++i) {
            cluster_mapping[i] = i;
        }

        int iterations = 2;

        LOG << "Number of iterations: " << iterations;

        NodeID last_node = G->number_of_nodes();

#pragma omp parallel
        {
            std::mt19937 m_mt;

            m_mt.seed(random_functions::getSeed());// + omp_get_thread_num());

            std::vector<std::pair<NodeID, NodeID> > wgt(G->number_of_nodes(), std::make_pair(0, 0));
            timer t;
            for (int j = 0; j < iterations; j++) {
#pragma omp for schedule(dynamic,1024)
                for (NodeID node = 0; node < last_node; ++node) {
                    NodeID n = permutation[node];
                    PartitionID max_block = cluster_mapping[n];
                    EdgeWeight max_value = 0;

                    for (EdgeID e : G->edges_of(n)) {
                        NodeID target = G->getEdgeTarget(e);
                        PartitionID block = cluster_mapping[target];

                        if (wgt[block].first != n) {
                            wgt[block].first = n;
                            wgt[block].second = 0;
                        }

                        wgt[block].second += G->getEdgeWeight(e);

                        if (wgt[block].second > max_value ||
                            (wgt[block].second == max_value &&
                             m_mt() % 2)) {
                            // random_functions::next() % 2)) {
                            max_value = wgt[block].second;
                            max_block = block;
                        }
                    }

                    cluster_mapping[n] = max_block;
                }
                LOGC(timing && !omp_get_thread_num())
                    << "LP: Iteration " << j << " -  Timer: " << t.elapsedToZero();
            }
        }

        return cluster_mapping;
    }

    std::vector<std::vector<NodeID> > remap_cluster(
        std::shared_ptr<graph_access> G, std::vector<NodeID>& cluster_mapping, bool save_cut) {
        PartitionID cur_no_clusters = 0;
        std::unordered_map<PartitionID, PartitionID> remap;

        std::vector<std::vector<NodeID> > reverse_mapping;

        std::vector<NodeID> part(G->number_of_nodes(), UNDEFINED_NODE);

        for (NodeID node : G->nodes()) {
            PartitionID cur_cluster = cluster_mapping[node];
            // check wether we already had that
            if (part[cur_cluster] == UNDEFINED_NODE) {
                part[cur_cluster] = cur_no_clusters++;
                reverse_mapping.emplace_back();
            }

            cluster_mapping[node] = part[cur_cluster];
            if (save_cut) {
                G->setPartitionIndex(node, part[cur_cluster]);
            }
            reverse_mapping[part[cur_cluster]].push_back(node);
        }

        return reverse_mapping;
    }

    void certify_clusters(
        std::shared_ptr<graph_access> G,
        std::vector<NodeID>& cluster_id,
        std::vector<std::vector<NodeID> >& reverse_mapping) {

        std::vector<bool> found(G->number_of_nodes(), false);
        NodeID clusters = reverse_mapping.size();
        std::vector<std::vector<NodeID> > marked;
        {
            std::vector<long int> in (G->number_of_nodes(), 0);
            std::vector<long int> out(G->number_of_nodes(), 0);
            std::vector<int> cluster(G->number_of_nodes(), true);

            for (size_t p = 0; p < reverse_mapping.size(); ++p) {
                marked.emplace_back();
            }

#pragma omp parallel for schedule(dynamic)
            for (size_t p = 0; p < clusters; ++p) {

                bool change = false;

                for (size_t i = 0; i < reverse_mapping[p].size(); ++i) {
                    NodeID node = reverse_mapping[p][i];
                    for (EdgeID e : G->edges_of(node)) {
                        if (cluster_id[G->getEdgeTarget(e)] == p) {
                            in[node] += G->getEdgeWeight(e);
                        }
                        else {
                            out[node] += G->getEdgeWeight(e);
                        }
                    }

                    if (out[node] > in [node]) {
                        cluster[node] = false;
                        change = true;
                        for (EdgeID e : G->edges_of(node)) {
                            NodeID tgt = G->getEdgeTarget(e);
                            if (cluster_id[tgt] == p) {
                                in[tgt] -= G->getEdgeWeight(e);
                                out[tgt] += G->getEdgeWeight(e);
                            }
                        }

                        marked[p].push_back(i);
                    }
                }

                while (change) {
                    change = false;
                    for (size_t i = 0; i < reverse_mapping[p].size(); ++i) {
                        NodeID node = reverse_mapping[p][i];
                        if (cluster[node] && out[node] > in [node]) {
                            cluster[node] = false;
                            for (EdgeID e : G->edges_of(node)) {
                                NodeID tgt = G->getEdgeTarget(e);
                                if (cluster_id[tgt] == p) {
                                    in[tgt] -= G->getEdgeWeight(e);
                                    out[tgt] += G->getEdgeWeight(e);
                                }
                            }
                            change = true;

                            marked[p].push_back(i);
                        }
                    }
                }
            }

            for (size_t p = 0; p < clusters; ++p) {
                if (marked[p].size()) {
                    std::sort(marked[p].begin(), marked[p].end(), std::greater<NodeID>());
                    for (auto m : marked[p]) {
                        if (reverse_mapping[p].size() > 1) {
                            NodeID node = reverse_mapping[p][m];
                            std::iter_swap(reverse_mapping[p].begin() + m,
                                           reverse_mapping[p].end() - 1);
                            reverse_mapping[p].pop_back();
                            reverse_mapping.emplace_back();
                            reverse_mapping.back().push_back(node);
                            cluster_id[node] = reverse_mapping.size() - 1;
                        }
                    }
                }
            }
        }

        reverse_mapping.reserve(2 * reverse_mapping.capacity());

#pragma omp parallel for schedule(dynamic)
        for (size_t p = 0; p < clusters; ++p) {
            if (marked[p].empty())
                continue;
            std::deque<NodeID> bfs;
            size_t totalnodes = reverse_mapping[p].size();
            size_t foundnodes = 0;
            NodeID map = p;
            while (foundnodes + 1 < totalnodes) {
                NodeID n;
                bfs.push_back(reverse_mapping[map][0]);
                found[reverse_mapping[map][0]] = true;

                while (!bfs.empty()) {
                    foundnodes++;
                    n = bfs.front();
                    bfs.pop_front();
                    for (EdgeID e : G->edges_of(n)) {
                        NodeID tgt = G->getEdgeTarget(e);
                        if (!found[tgt] && cluster_id[tgt] == map) {
                            found[tgt] = true;
                            bfs.push_back(tgt);
                        }
                    }
                }

                if (foundnodes < totalnodes) {
                    // move not found nodes to new cluster
                    size_t lower = 0;
                    size_t upper = reverse_mapping[map].size() - 1;
                    while (lower < upper) {
                        while (found[reverse_mapping[map][lower]]) {
                            ++lower;
                        }

                        while (!found[reverse_mapping[map][upper]]) {
                            --upper;
                        }

                        if (lower < upper) {
                            std::iter_swap(reverse_mapping[map].begin() + lower,
                                           reverse_mapping[map].begin() + upper);
                        }
                    }

                    NodeID new_map;
#pragma omp critical
                    {
                        reverse_mapping.emplace_back();
                        new_map = reverse_mapping.size() - 1;
                    }

                    for (size_t i = lower; i < reverse_mapping[map].size(); ++i) {

                        NodeID el = reverse_mapping[map][i];
                        cluster_id[el] = new_map;
                        reverse_mapping[new_map].push_back(el);
                    }
                    reverse_mapping[map].erase(reverse_mapping[map].begin() + lower,
                                               reverse_mapping[map].end());
                    map = new_map;
                }
            }
        }
    }
};
