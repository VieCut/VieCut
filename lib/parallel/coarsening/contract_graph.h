/******************************************************************************
 * contract_graph.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2017-2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <functional>
#include <memory>
#include <unordered_set>
#include <utility>
#include <vector>

#include "common/configuration.h"
#include "common/definitions.h"
#include "data-structures/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "parallel/data_structure/union_find.h"
#include "tlx/logger.hpp"
#include "tools/hash.h"
#include "tools/timer.h"

class contraction {
 public:
    static constexpr bool debug = false;

    static std::shared_ptr<graph_access> deleteEdge(
        std::shared_ptr<graph_access>, EdgeID) {
        LOG1 << "DELETE EDGE NOT IMPLEMENTED YET";
        exit(2);
    }

    static std::pair<std::shared_ptr<graph_access>, std::vector<NodeID> >
    contractEdge(std::shared_ptr<graph_access>,
                 std::vector<NodeID>,
                 EdgeID) {
        LOG1 << "CONTRACT EDGE NOT IMPLEMENTED YET";
        exit(2);
    }

    static inline uint64_t get_uint64_from_pair(NodeID cluster_a,
                                                NodeID cluster_b) {
        if (cluster_a > cluster_b) {
            std::swap(cluster_a, cluster_b);
        }
        return ((uint64_t)cluster_a << 32) | cluster_b;
    }

    static inline std::pair<NodeID, NodeID> get_pair_from_uint64(
        uint64_t data) {
        NodeID first = data >> 32;
        NodeID second = data;
        return std::make_pair(first, second);
    }

    static void findTrivialCuts(std::shared_ptr<graph_access> G,
                                std::vector<NodeID>* m,
                                std::vector<std::vector<NodeID> >* rm,
                                NodeWeight target_mindeg) {
        // create non-const references for better syntax
        std::vector<NodeID>& mapping = *m;
        std::vector<std::vector<NodeID> >& rev_mapping = *rm;

        LOG << "target min degree: " << target_mindeg;
#pragma omp parallel for schedule(dynamic, 1024)
        for (NodeID p = 0; p < rev_mapping.size(); ++p) {
            NodeID bestNode;
            NodeWeight improve = 0;
            NodeWeight node_degree = 0;
            NodeWeight block_degree = 0;
            if (rev_mapping[p].size() < std::log2(G->number_of_nodes())) {
                NodeID improve_idx;
                for (NodeID node = 0; node < rev_mapping[p].size(); ++node) {
                    for (EdgeID e : G->edges_of(rev_mapping[p][node])) {
                        NodeID contracted_target = mapping[G->getEdgeTarget(e)];

                        if (contracted_target == p) {
                            node_degree += G->getEdgeWeight(e);
                            continue;
                        }

                        node_degree -= G->getEdgeWeight(e);
                        block_degree += G->getEdgeWeight(e);
                    }

                    if (improve > node_degree) {
                        improve = node_degree;
                        bestNode = rev_mapping[p][node];
                        improve_idx = node;
                    }
                    node_degree = 0;
                }
                if (improve > 0 &&
                    block_degree + improve < target_mindeg &&
                    rev_mapping[p].size() > 1) {
                    target_mindeg = block_degree + improve;
                    rev_mapping[p].erase(rev_mapping[p].begin() +
                                         improve_idx);
                    VIECUT_ASSERT_LT(bestNode, G->number_of_nodes());
                    rev_mapping.push_back({ bestNode });
                    mapping[bestNode] = rev_mapping.size() - 1;
                }
            }
        }

        LOG << "target min degree now: " << target_mindeg;
    }

    // contraction global_mincut for small number of nodes in constructed graph,
    // we assume a full mesh and remove nonexistent edges afterwards.
    static std::shared_ptr<graph_access> contractGraphFullMesh(
        std::shared_ptr<graph_access> G,
        const std::vector<NodeID>& mapping,
        size_t num_nodes) {
        auto contracted = std::make_shared<graph_access>();

        std::vector<EdgeWeight> intermediate(num_nodes * (num_nodes - 1), 0);

#pragma omp parallel
        {
            std::vector<EdgeWeight> p_intermediate(
                num_nodes * (num_nodes - 1), 0);

#pragma omp for schedule(dynamic, 1024)
            for (NodeID n = 0; n < G->number_of_nodes(); ++n) {
                NodeID src = mapping[n];
                for (EdgeID e : G->edges_of(n)) {
                    NodeID tgt = mapping[G->getEdgeTarget(e)];

                    if (tgt != src) {
                        EdgeID edge_id =
                            src * (num_nodes - 1) + tgt - (tgt > src);
                        p_intermediate[edge_id] += G->getEdgeWeight(e);
                    }
                }
            }

#pragma omp critical
            {
                for (size_t i = 0; i < intermediate.size(); ++i) {
                    intermediate[i] += p_intermediate[i];
                }
            }
        }
        EdgeID existing_edges = intermediate.size();
        for (auto e : intermediate) {
            if (e == 0)
                --existing_edges;
        }

        contracted->start_construction(num_nodes, existing_edges);

        for (size_t i = 0; i < num_nodes; ++i) {
            contracted->new_node();
            for (size_t j = 0; j < num_nodes; ++j) {
                if (i == j)
                    continue;

                EdgeID edge_id = i * (num_nodes - 1) + j - (j > i);

                if (intermediate[edge_id] > 0) {
                    EdgeID edge = contracted->new_edge(i, j);
                    contracted->setEdgeWeight(edge, intermediate[edge_id]);
                }
            }
        }

        contracted->finish_construction();

        return contracted;
    }

    static std::shared_ptr<mutable_graph> fromUnionFind(
        std::shared_ptr<mutable_graph> G,
        union_find* uf) {
        std::vector<std::vector<NodeID> > rev_mapping(uf->n());

        std::vector<NodeID> part(G->number_of_nodes(), UNDEFINED_NODE);
        NodeID current_pid = 0;

        for (NodeID n : G->nodes()) {
            NodeID part_id = uf->Find(n);
            if (part[part_id] == UNDEFINED_NODE) {
                part[part_id] = current_pid++;
            }
            rev_mapping[part[part_id]].push_back(G->containedVertices(n)[0]);
        }

        for (size_t i = 0; i < rev_mapping.size(); ++i) {
            if (rev_mapping[i].size() > 1) {
                std::unordered_set<NodeID> vtx_to_ctr;
                for (auto v : rev_mapping[i]) {
                    vtx_to_ctr.emplace(G->getCurrentPosition(v));
                }
                G->contractVertexSet(vtx_to_ctr);
            }
        }

        return G;
    }

    static std::shared_ptr<graph_access> fromUnionFind(
        std::shared_ptr<graph_access> G,
        union_find* uf) {
        std::vector<std::vector<NodeID> > rev_mapping;

        std::vector<NodeID> mapping(G->number_of_nodes());
        std::vector<NodeID> part(G->number_of_nodes(), UNDEFINED_NODE);
        NodeID current_pid = 0;
        for (NodeID n : G->nodes()) {
            NodeID part_id = uf->Find(n);

            if (part[part_id] == UNDEFINED_NODE) {
                part[part_id] = current_pid++;
                rev_mapping.emplace_back();
            }

            mapping[n] = part[part_id];
            if (configuration::getConfig()->save_cut) {
                G->setPartitionIndex(n, part[part_id]);
            }
            rev_mapping[part[part_id]].push_back(n);
        }
        return contractGraph(G, mapping, rev_mapping.size());
    }

    static std::shared_ptr<graph_access>
    contractGraph(std::shared_ptr<graph_access> G,
                  const std::vector<NodeID>& mapping,
                  size_t num_nodes,
                  const std::vector<std::vector<NodeID> >& = { }) {
        if (num_nodes > std::sqrt(G->number_of_nodes())) {
            LOG << "SPARSE CONTRACT!";
            return contractGraphSparse(G, mapping, num_nodes);
        } else {
            LOG << "FULL MESH CONTRACT";
            return contractGraphFullMesh(G, mapping, num_nodes);
        }
    }

    // altered version of KaHiPs matching contraction
    static std::shared_ptr<graph_access>
    contractGraphSparse(std::shared_ptr<graph_access> G,
                        const std::vector<NodeID>& mapping,
                        size_t num_nodes) {
        // contested edge (both incident vertices have at least V/5 vertices)
        // compute value for this edge on every processor to allow parallelism
        timer t;
        EdgeID contested_edge = 0;
        NodeID block0 = 0;
        NodeID block1 = 0;

        if (G->number_of_edges() * 0.02
            < G->number_of_nodes() * G->number_of_nodes() &&
            G->number_of_nodes() > 100) {
            std::vector<uint32_t> el(num_nodes);
            for (size_t i = 0; i < mapping.size(); ++i) {
                ++el[mapping[i]];
            }

            std::vector<uint32_t> orig_el = el;
            std::nth_element(el.begin(), el.begin() + 1, el.end(),
                             std::greater<uint32_t>());

            if (el[1] > G->number_of_nodes() / 5) {
                block0 = std::distance(orig_el.begin(),
                                       std::find(orig_el.begin(),
                                                 orig_el.end(), el[0]));
                block1 = std::distance(orig_el.begin(),
                                       std::find(orig_el.begin(),
                                                 orig_el.end(), el[1]));

                contested_edge = get_uint64_from_pair(block1, block0);
            }
        }

        EdgeWeight sumweight_contested = 0;

        auto coarser = std::make_shared<graph_access>();
        std::vector<std::vector<std::pair<PartitionID, EdgeWeight> > >
        building_tool(num_nodes);
        std::vector<size_t> degrees(num_nodes);
        growt::uaGrow<xxhash<uint64_t> > new_edges(1024 * 1024);
        t.restart();
        std::vector<size_t> cur_degrees(num_nodes);
#pragma omp parallel
        {
            EdgeWeight contested_weight = 0;

            std::vector<uint64_t> my_keys;
            auto handle = new_edges.get_handle();
#pragma omp for schedule(guided)
            for (NodeID n = 0; n < G->number_of_nodes(); ++n) {
                NodeID p = mapping[n];
                for (EdgeID e : G->edges_of(n)) {
                    NodeID contracted_target = mapping[G->getEdgeTarget(e)];
                    if (contracted_target <= p) {
                        // self-loops are not in graph
                        // smaller do not need to be stored
                        // as their other side will be
                        continue;
                    }
                    EdgeWeight edge_weight = G->getEdgeWeight(e);
                    uint64_t key = get_uint64_from_pair(p, contracted_target);

                    if (key != contested_edge) {
                        if (handle.insert_or_update(key, edge_weight,
                                                    [](size_t& lhs,
                                                       const size_t& rhs) {
                                                        lhs += rhs;
                                                    }, edge_weight).second) {
#pragma omp atomic
                            ++degrees[p];
#pragma omp atomic
                            ++degrees[contracted_target];
                            my_keys.push_back(key);
                        }
                    } else {
                        contested_weight += edge_weight;
                    }
                }
            }

            if (contested_edge > 0) {
#pragma omp critical
                {
                    sumweight_contested += contested_weight;
                }
#pragma omp barrier
#pragma omp single
                {
                    if (sumweight_contested > 0) {
                        handle.insert_or_update(contested_edge,
                                                sumweight_contested,
                                                [](size_t& lhs,
                                                   const size_t& rhs) {
                                                    lhs += rhs;
                                                }, sumweight_contested);
                        my_keys.push_back(contested_edge);
                        ++degrees[block0];
                        ++degrees[block1];
                    }
                }
            }

#pragma omp single
            {
                size_t num_edges = 0;
                coarser->start_construction(num_nodes, 0);
                for (size_t i = 0; i < degrees.size(); ++i) {
                    cur_degrees[i] = num_edges;
                    num_edges += degrees[i];
                    coarser->new_node_hacky(num_edges);
                }
                coarser->resize_m(num_edges);
            }

            for (auto edge_uint : my_keys) {
                auto edge = get_pair_from_uint64(edge_uint);
                auto edge_weight = (*handle.find(edge_uint)).second;
                size_t firstdeg, seconddeg;

                while (true) {
                    firstdeg = cur_degrees[edge.first];
                    size_t plusone = cur_degrees[edge.first] + 1;
                    if (__sync_bool_compare_and_swap(&cur_degrees[edge.first],
                                                     firstdeg, plusone))
                        break;
                }

                while (true) {
                    seconddeg = cur_degrees[edge.second];
                    size_t plusone = cur_degrees[edge.second] + 1;
                    if (__sync_bool_compare_and_swap(&cur_degrees[edge.second],
                                                     seconddeg,
                                                     plusone))
                        break;
                }

                coarser->new_edge_and_reverse(edge.first, edge.second,
                                              firstdeg, seconddeg, edge_weight);
            }
        }
        coarser->finish_construction();
        return coarser;
    }

    static std::shared_ptr<graph_access>
    contractGraphSparseNoHash(std::shared_ptr<graph_access> G,
                              const std::vector<NodeID>& mapping,
                              const std::vector<std::vector<NodeID> >&
                              rev_mapping, size_t num_nodes) {
        std::vector<std::vector<NodeID> > rev_map;
        if (rev_mapping.size() == 0) {
            // create reverse mapping if it wasnt before
            rev_map.resize(num_nodes);
            for (size_t i = 0; i < mapping.size(); ++i) {
                rev_map[mapping[i]].push_back(i);
            }
        } else {
            rev_map = rev_mapping;
        }

        auto contracted = std::make_shared<graph_access>();

        std::vector<std::vector<std::pair<NodeID, EdgeWeight> > > edges;
        edges.resize(rev_map.size());
#pragma omp parallel
        {
#pragma omp single nowait
            {
                double average_degree =
                    static_cast<double>(G->number_of_edges()) /
                    static_cast<double>(G->number_of_nodes());

                EdgeID expected_edges = num_nodes * average_degree;
                // one worker can do this vector allocation while the others
                // build the contracted graph
                contracted->start_construction(num_nodes,
                                               std::min(G->number_of_edges(),
                                                        2 * expected_edges));
            }

            // first: coarse vertex which set this (to avoid total invalidation)
            // second: edge id in contracted graph
            std::vector<std::pair<NodeID, EdgeWeight> > edge_positions(
                num_nodes,
                std::make_pair(UNDEFINED_NODE, UNDEFINED_EDGE));

            std::vector<NodeID> non_null;

#pragma omp for schedule(dynamic)
            for (NodeID p = 0; p < num_nodes; ++p) {
                for (NodeID node = 0; node < rev_map[p].size(); ++node) {
                    for (EdgeID e : G->edges_of(rev_map[p][node])) {
                        NodeID contracted_target = mapping[G->getEdgeTarget(e)];

                        if (contracted_target == p)
                            continue;

                        NodeID last_use =
                            edge_positions[contracted_target].first;

                        if (last_use == p) {
                            edge_positions[contracted_target].second +=
                                G->getEdgeWeight(e);
                        } else {
                            edge_positions[contracted_target].first = p;
                            edge_positions[contracted_target].second =
                                G->getEdgeWeight(e);

                            non_null.push_back(contracted_target);
                        }
                    }
                }

                for (const auto& tgt : non_null) {
                    edges[p].emplace_back(tgt, edge_positions[tgt].second);
                }

                non_null.clear();
            }
        }

        for (const auto& vec : edges) {
            NodeID n = contracted->new_node();
            for (const auto& e : vec) {
                EdgeID e_new = contracted->new_edge(n, e.first);
                contracted->setEdgeWeight(e_new, e.second);
            }
        }

        contracted->finish_construction();

        return contracted;
    }
};
