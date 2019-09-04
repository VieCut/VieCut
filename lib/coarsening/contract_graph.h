/******************************************************************************
 * contract_graph.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/misc/graph_algorithms.h"
#include "common/configuration.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "data_structure/union_find.h"
#include "tlx/logger.hpp"
#include "tools/string.h"
#include "tools/vector.h"

class contraction {
 public:
    static constexpr bool debug = false;

    static std::shared_ptr<graph_access> deleteEdge(
        std::shared_ptr<graph_access> G, EdgeID e_del) {
        std::shared_ptr<graph_access> G_out = std::make_shared<graph_access>();

        G_out->start_construction(G->number_of_nodes(),
                                  G->number_of_edges() - 2);
        EdgeID e_del_rev = G->find_reverse_edge(e_del);

        for (NodeID n : G->nodes()) {
            for (EdgeID e : G->edges_of(n)) {
                if (e != e_del && e != e_del_rev && G->getEdgeWeight(e)) {
                    EdgeID new_e = G_out->new_edge(n, G->getEdgeTarget(e));
                    G_out->setEdgeWeight(new_e, G->getEdgeWeight(e));
                }
            }
            G_out->new_node();
        }
        G_out->finish_construction();
        return G_out;
    }

    static std::pair<std::shared_ptr<graph_access>, std::vector<NodeID> >
    contractEdge(std::shared_ptr<graph_access> G,
                 std::vector<NodeID> terminals,
                 EdgeID e_ctr) {
        NodeID src = G->getEdgeSource(e_ctr);
        NodeID tgt = G->getEdgeTarget(e_ctr);
        union_find uf(G->number_of_nodes());
        uf.Union(src, tgt);

        auto G_out = contraction::fromUnionFind(G, &uf);
        std::vector<NodeID> terminals_out;
        for (NodeID t : terminals) {
            terminals_out.emplace_back(G->getPartitionIndex(t));
        }

        return std::make_pair(G_out, terminals_out);
    }

    static void findTrivialCuts(std::shared_ptr<graph_access> G,
                                std::vector<NodeID>* m,
                                std::vector<std::vector<NodeID> >* rm,
                                int target_mindeg) {
        // non-const references for better syntax
        std::vector<NodeID>& mapping = *m;
        std::vector<std::vector<NodeID> > reverse_mapping = *rm;

        for (NodeID p = 0; p < reverse_mapping.size(); ++p) {
            NodeID bestNode;
            int64_t improve = 0;
            int64_t node_degree = 0;
            int64_t block_degree = 0;
            if (reverse_mapping[p].size() < std::log2(G->number_of_nodes())) {
                NodeID improve_idx;
                for (NodeID node = 0;
                     node < reverse_mapping[p].size(); ++node) {
                    for (EdgeID e : G->edges_of(reverse_mapping[p][node])) {
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
                        bestNode = reverse_mapping[p][node];
                        improve_idx = node;
                    }
                    node_degree = 0;
                }
                if (improve < 0 &&
                    block_degree + improve < target_mindeg &&
                    reverse_mapping[p].size() > 1) {
                    target_mindeg = block_degree + improve;
                    LOG << bestNode << " with " << block_degree +
                        improve << " possible degree in block " << p;
                    reverse_mapping[p].erase(reverse_mapping[p].begin() +
                                             improve_idx);
                    assert(bestNode < G->number_of_nodes());
                    reverse_mapping.push_back({ bestNode });
                    LOG << reverse_mapping.size() << " size";
                    mapping[bestNode] = reverse_mapping.size() - 1;
                    p--;
                }
            }
        }

        LOG << "target min degree now: " << target_mindeg;
    }

    // contraction global_mincut for small number of nodes in constructed graph
    // we assume a full mesh and remove nonexistent edges afterwards.
    static std::shared_ptr<graph_access> contractGraphFullMesh(
        std::shared_ptr<graph_access> G,
        const std::vector<NodeID>& mapping,
        const std::vector<std::vector<NodeID> >&
        reverse_mapping) {
        NodeID num_nodes = reverse_mapping.size();
        auto contracted = std::make_shared<graph_access>();
        std::vector<EdgeWeight> intermediate(num_nodes * (num_nodes - 1), 0);

        for (NodeID n = 0; n < G->number_of_nodes(); ++n) {
            NodeID src = mapping[n];
            for (EdgeID e : G->edges_of(n)) {
                NodeID tgt = mapping[G->getEdgeTarget(e)];
                if (tgt != src) {
                    EdgeID edge_id = src * (num_nodes - 1) + tgt - (tgt > src);
                    intermediate[edge_id] += G->getEdgeWeight(e);
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
        std::vector<std::vector<NodeID> > reverse_mapping(uf->n());

        std::vector<NodeID> part(G->number_of_nodes(), UNDEFINED_NODE);
        NodeID current_pid = 0;

        for (NodeID n : G->nodes()) {
            NodeID part_id = uf->Find(n);

            if (part[part_id] == UNDEFINED_NODE) {
                part[part_id] = current_pid++;
            }

            reverse_mapping[part[part_id]].push_back(
                G->containedVertices(n)[0]);
        }

        for (size_t i = 0; i < reverse_mapping.size(); ++i) {
            if (reverse_mapping[i].size() > 1) {
                std::unordered_set<NodeID> vtx_to_ctr;
                for (auto v : reverse_mapping[i]) {
                    vtx_to_ctr.emplace(G->getCurrentPosition(v));
                }

                G->contractVertexSet(vtx_to_ctr);
            }
        }

        if (debug) {
            graph_algorithms::checkGraphValidity(G);
        }

        return G;
    }

    static std::shared_ptr<graph_access> fromUnionFind(
        std::shared_ptr<graph_access> G,
        union_find* uf) {
        std::vector<std::vector<NodeID> > reverse_mapping(uf->n());

        std::vector<NodeID> mapping(G->number_of_nodes());
        std::vector<NodeID> part(G->number_of_nodes(), UNDEFINED_NODE);
        NodeID current_pid = 0;

        for (NodeID n : G->nodes()) {
            NodeID part_id = uf->Find(n);

            if (part[part_id] == UNDEFINED_NODE) {
                part[part_id] = current_pid++;
            }

            mapping[n] = part[part_id];
            if (configuration::getConfig()->save_cut) {
                G->setPartitionIndex(n, part[part_id]);
            }
            reverse_mapping[part[part_id]].push_back(n);
        }
        return contractGraph(G, mapping,
                             reverse_mapping.size(), reverse_mapping);
    }

    static std::shared_ptr<graph_access> contractGraph(
        std::shared_ptr<graph_access> G,
        const std::vector<NodeID>& mapping,
        size_t /*num_nodes*/,
        const std::vector<std::vector<NodeID> >& reverse_mapping) {
        if (reverse_mapping.size() > std::sqrt(G->number_of_nodes())) {
            return contractGraphSparse(G, mapping, reverse_mapping);
        } else {
            return contractGraphFullMesh(G, mapping, reverse_mapping);
        }
    }

    // altered version of KaHiPs matching contraction
    static std::shared_ptr<graph_access> contractGraphSparse(
        std::shared_ptr<graph_access> G,
        const std::vector<NodeID>& mapping,
        const std::vector<std::vector<NodeID> >&
        reverse_mapping) {
        // first: coarse vertex which set this (to avoid total invalidation)
        // second: edge id in contracted graph
        auto contracted = std::make_shared<graph_access>();
        contracted->start_construction(reverse_mapping.size(),
                                       G->number_of_edges());
        std::vector<std::pair<NodeID, EdgeID> > edge_positions(
            reverse_mapping.size(),
            std::make_pair(UNDEFINED_NODE, UNDEFINED_EDGE));

        for (NodeID p = 0; p < reverse_mapping.size(); ++p) {
            contracted->new_node();
            for (NodeID node = 0; node < reverse_mapping[p].size(); ++node) {
                for (EdgeID e : G->edges_of(reverse_mapping[p][node])) {
                    if (G->getEdgeWeight(e) == 0)
                        continue;

                    NodeID tgt = G->getEdgeTarget(e);

                    NodeID contracted_target = mapping[tgt];

                    if (contracted_target == p)
                        continue;

                    NodeID last_use = edge_positions[contracted_target].first;

                    if (last_use == p) {
                        EdgeWeight weight = G->getEdgeWeight(e);
                        EdgeID e_new = edge_positions[contracted_target].second;
                        contracted->setEdgeWeight(
                            e_new, contracted->getEdgeWeight(e_new) + weight);
                    } else {
                        EdgeID e_new =
                            contracted->new_edge(p, contracted_target);
                        edge_positions[contracted_target].first = p;
                        edge_positions[contracted_target].second = e_new;

                        contracted->setEdgeWeight(e_new, G->getEdgeWeight(e));
                    }
                }
            }
        }

        LOG << "Contracted from " << mapping.size() << " to " <<
            reverse_mapping.size() << "nodes!";

        contracted->finish_construction();

        return contracted;
    }
};
