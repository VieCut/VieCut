/******************************************************************************
 * find_bridges.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <data_structure/mutable_graph.h>

class find_bridges {
 public:
    find_bridges(std::shared_ptr<mutable_graph> G)
        : G(G),
          visited(G->n(), false),
          discovered(G->n(), UNDEFINED_NODE),
          parent(G->n(), UNDEFINED_NODE),
          lowest(G->n(), UNDEFINED_NODE),
          step(0) { };

    void visit(NodeID n) {
        // DFS, check whether alternative path exists when backtracking
        visited[n] = true;
        step++;
        discovered[n] = step;
        lowest[n] = step;
        for (EdgeID e : G->edges_of(n)) {
            NodeID t = G->getEdgeTarget(n, e);
            if (!visited[t]) {
                parent[t] = n;
                visit(t);
                lowest[n] = std::min(lowest[n], lowest[t]);
                if (lowest[t] > discovered[n]) {
                    LOG1 << "THERE IS A BRIDGE BETWEEN " << n << " and " << t
                    << "edge id " << e;
                }
            } else if (t != parent[n]) {
                lowest[n] = std::min(lowest[n], discovered[t]);
            }
        }
    }

    void find_all_bridges() {
        LOG1 << "trying to find bridges!";
        for (NodeID n : G->nodes()) {
            if (!visited[n]) {
                visit(n);
            }
        }
    }

 private:
    std::shared_ptr<mutable_graph> G;
    std::vector<bool> visited;
    std::vector<NodeID> discovered;
    std::vector<NodeID> parent;
    std::vector<NodeID> lowest;
    size_t step;
};
