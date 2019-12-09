/******************************************************************************
 * equal_neighborhood.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <memory>

#include "data_structure/union_find.h"
#include "data_structure/mutable_graph.h"

// from boost::hash
template <class T>
inline void hash_combine(std::size_t& seed, const T& v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

auto pairs = [](const std::pair<NodeID, EdgeWeight> p1,
                const std::pair<NodeID, EdgeWeight> p2) {
                    return p1.first < p2.first;
                };

class equal_neighborhood {
 public:
    //explicit equal_neighborhood() {};

    union_find findEqualNeighborhoods(
        std::shared_ptr<multicut_problem> problem) {
        std::shared_ptr<mutable_graph> G = problem->graph;
        union_find uf(G->n());

        std::unordered_set<NodeID> terminals;
        for (const auto& t: problem->terminals) {
            terminals.emplace(t.position);
        }

        for (NodeID n : G->nodes()) {
            if (terminals.count(n) > 0) {
                continue;
            }

            if (G->getNodeDegree(n) <= 10) {
                std::vector<NodeID> v;
                for (EdgeID e : G->edges_of(n)) {
                    NodeID tgt = G->getEdgeTarget(n, e);
                    v.emplace_back(tgt);
                }
                std::sort(v.begin(), v.end());
                size_t seed = 0;
                for (NodeID vtx : v) {
                    hash_combine(seed, vtx);
                }

                if (results.count(seed) > 0) {
                    NodeID o = results[seed];
                    std::vector<std::pair<NodeID, EdgeWeight>> neighbors_n;
                    std::vector<std::pair<NodeID, EdgeWeight>> neighbors_o;

                    for (EdgeID e : G->edges_of(n)) {
                        auto [tgt, wgt] = G->getEdge(n, e);
                        neighbors_n.emplace_back(tgt, wgt);
                    }

                    for (EdgeID e : G->edges_of(o)) {
                        auto [tgt, wgt] = G->getEdge(o, e);
                        neighbors_o.emplace_back(tgt, wgt);
                    }

                    std::sort(neighbors_n.begin(), neighbors_n.end(), pairs);
                    std::sort(neighbors_o.begin(), neighbors_o.end(), pairs);
                    bool difference_found = false;
                    for (size_t i = 0; i < neighbors_n.size(); ++i) {
                        if (neighbors_n[i] != neighbors_o[i]) {
                            difference_found = true;
                            break;
                        }
                    }

                    if (!difference_found) {
                        uf.Union(n, o);
                    }
                } else {
                    results.emplace(seed, n);
                }
            }
        }
        
        return uf;
    }

 private:
    std::unordered_map<size_t, NodeID> results;
};