/******************************************************************************
 * graph_algorithms.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <memory>
#include <queue>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#pragma once

class graph_algorithms {
 public:
    static std::vector<NodeID> top_k_degrees(
        std::shared_ptr<graph_access> G, size_t k) {
        std::vector<std::pair<NodeID, EdgeWeight> > all_degrees;
        for (NodeID n : G->nodes()) {
            all_degrees.emplace_back(n, G->getNodeDegree(n));
        }

        return find_top_k(all_degrees, k);
    }

    static std::vector<NodeID> top_k_degrees(
        std::shared_ptr<mutable_graph> G, size_t k) {
        std::vector<std::pair<NodeID, EdgeWeight> > all_degrees;
        for (NodeID n : G->nodes()) {
            all_degrees.emplace_back(n, G->getWeightedNodeDegree(n));
        }

        return find_top_k(all_degrees, k);
    }

    static std::vector<NodeID> weighted_top_k_degrees(
        std::shared_ptr<graph_access> G, size_t k) {
        std::vector<std::pair<NodeID, EdgeWeight> > all_degrees;
        for (NodeID n : G->nodes()) {
            all_degrees.emplace_back(n, G->getWeightedNodeDegree(n));
        }

        return find_top_k(all_degrees, k);
    }

    static void checkGraphValidity(std::shared_ptr<graph_access> G) {
        std::unordered_map<uint64_t, EdgeWeight> weights;
        size_t reverse_found = 0;
        for (NodeID n : G->nodes()) {
            for (EdgeID e : G->edges_of(n)) {
                EdgeWeight w = G->getEdgeWeight(e);
                NodeID t = G->getEdgeTarget(e);

                if (t >= G->number_of_nodes()) {
                    LOG1 << "edge to element outside of graph";
                    exit(1);
                }

                uint64_t pair = uint64_from_pair(n, t);
                uint64_t reverse_pair = uint64_from_pair(t, n);
                
                if (weights.find(pair) != weights.end()) {
                        LOG1 << pair;
                    for (auto w : weights) {
                        LOG1 << w;
                    }

                    LOG1 << "ERROR: Multi-edge detected from "
                         << n << " to " << t;
                    exit(1);
                }
                auto el = weights.find(reverse_pair);
                if (el != weights.end()) {
                    if (el->second != w) {
                        LOG1 << "weight unequal!";
                        exit(1);
                    }
                    reverse_found++;
                } 
                weights.emplace(pair, w);               
            }
        }

        if (reverse_found * 2 != G->number_of_edges()) {
            LOG1 << "Edges not adding up! " << reverse_found << " " << G->number_of_edges();
            exit(1);
        }
        LOG1 << "Graph is valid!";
    }

    static void checkGraphValidity(std::shared_ptr<mutable_graph> G) {
        size_t edges = 0;
        for (NodeID n : G->nodes()) {
            for (NodeID v : G->containedVertices(n)) {
                if (G->getCurrentPosition(v) != n) {
                    LOG1 << "ERROR: Vertex pos of " << v << " is not " << n;
                    exit(1);
                }
            }

            EdgeWeight weight = 0;
            std::unordered_set<NodeID> targets;
            edges += G->get_first_invalid_edge(n);
            for (EdgeID e : G->edges_of(n)) {
                NodeID tgt = G->getEdgeTarget(n, e);
                NodeID rev = G->getReverseEdge(n, e);

                weight += G->getEdgeWeight(n, e);

                if (targets.count(tgt) > 0) {
                    LOG1 << "ERROR: Double edge from " << n << " to " << tgt;
                    exit(1);
                } else {
                    targets.insert(tgt);
                }

                if (tgt == n) {
                    LOG1 << "ERROR: Self edge from " << n << " to " << tgt;
                    exit(1);
                }

                if (tgt >= G->n()) {
                    LOG1 << "ERROR: Edge " << e << " of vertex "
                         << n << " points to " << tgt
                         << ", graph only has " << G->n() << "nodes!";
                    exit(1);
                }

                if (G->getEdgeTarget(tgt, rev) != n) {
                    LOG1 << "ERROR: Graph invalid: " << n
                         << "-" << e << " is not the correct edge target to "
                         << tgt << "-" << rev;
                    exit(1);
                }

                if (G->getEdgeWeight(tgt, rev) != G->getEdgeWeight(n, e)) {
                    LOG1 << "ERROR: Graph invalid: " << n << "-" << e
                         << " weight (" << G->getEdgeWeight(n, e)
                         << ") not equal to " << tgt
                         << "-" << rev << " ("
                         << G->getEdgeWeight(tgt, rev) << ")";
                    exit(1);
                }

                if (G->getReverseEdge(tgt, rev) != e) {
                    LOG1 << "ERROR: Graph invalid: " << n << "-"
                         << e << " is not the reverse edge to "
                         << tgt << "-" << rev;
                    exit(1);
                }
            }

            if (weight != G->getWeightedNodeDegree(n)) {
                LOG1 << "ERROR: Node weight in node " << n
                     << " invalid :" << weight
                     << " != " << G->getWeightedNodeDegree(n);
                exit(1);
            }
        }

        if (edges != G->m()) {
            LOG1 << "Error: Invalid number of edges!";
            exit(1);
        }

        LOG1 << "Graph is valid!";
    }

    static std::tuple<std::vector<NodeID>, std::vector<uint32_t>, NodeID>
    bfsDistances(std::shared_ptr<mutable_graph> G, NodeID start) {
        std::queue<NodeID> q;
        std::vector<NodeID> parent(G->n(), UNDEFINED_NODE);
        std::vector<uint32_t> distance(G->n(), 0);
        std::vector<bool> discovered(G->n(), false);
        NodeID vtcs = 1;
        discovered[start] = true;
        q.push(start);

        NodeID top = UNDEFINED_NODE;
        while (!q.empty()) {
            top = q.front();
            q.pop();
            for (EdgeID e : G->edges_of(top)) {
                NodeID t = G->getEdgeTarget(top, e);
                if (!discovered[t]) {
                    discovered[t] = true;
                    parent[t] = top;
                    distance[t] = distance[top] + 1;
                    vtcs++;
                    if (vtcs == G->n()) {
                        return std::make_tuple(parent, distance, t);
                    }
                    q.push(t);
                }
            }
        }

        return std::make_tuple(parent, distance, top);
    }

 private:
    static std::vector<NodeID> find_top_k(
        std::vector<std::pair<NodeID, EdgeWeight> > in, size_t k) {
        std::nth_element(in.begin(), in.end() - k, in.end(),
                         [](auto a1, auto a2) {
                             return a1.second < a2.second;
                         });

        std::vector<NodeID> out;
        for (size_t i = in.size() - k; i < in.size(); ++i) {
            out.emplace_back(in[i].first);
        }
        return out;
    }

    static inline uint64_t uint64_from_pair(NodeID cluster_a,
                                                NodeID cluster_b) {
        //if (cluster_a > cluster_b) {
        //    std::swap(cluster_a, cluster_b);
        //}
        return ((uint64_t)cluster_a << 32) | cluster_b;
    }

    static inline std::pair<NodeID, NodeID> pair_from_uint64(
        uint64_t data) {
        NodeID first = data >> 32;
        NodeID second = data;
        return std::make_pair(first, second);
    }
};
