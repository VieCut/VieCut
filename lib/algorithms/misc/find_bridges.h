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

#include <algorithms/multicut/multicut_problem.h>
#include <data_structure/mutable_graph.h>

class find_bridges {
 public:
    find_bridges(std::shared_ptr<mutable_graph> G)
        : G(G),
          visited(G->n(), false),
          discovered(G->n(), UNDEFINED_NODE),
          parent(G->n(), UNDEFINED_NODE),
          lowest(G->n(), UNDEFINED_NODE),
          step(0) { }

    void findAllBridgesInCC(NodeID vtx) {
        std::vector<NodeID> path;
        path.emplace_back(vtx);
        parent[vtx] = vtx;
        stack.push(vtx);
        while (!stack.empty()) {
            NodeID n = stack.top();
            stack.pop();
            if (visited[n])
                continue;

            while (path.size() > 1 && path.back() != parent[n]) {
                // explicit backtracking to avoid recursion
                NodeID b = path.back();
                path.pop_back();
                for (EdgeID e : G->edges_of(b)) {
                    NodeID t = G->getEdgeTarget(b, e);
                    if (parent[t] == b) {
                        lowest[b] = std::min(lowest[b], lowest[t]);
                        if (lowest[t] > discovered[b]) {
                            bridges.emplace_back(b, e);
                        }
                    } else {
                        if (t != parent[b]) {
                            lowest[b] = std::min(lowest[b], discovered[t]);
                        }
                    }
                }
            }

            path.emplace_back(n);
            // DFS, check whether alternative path exists when backtracking
            visited[n] = true;
            step++;
            discovered[n] = step;
            lowest[n] = step;
            for (EdgeID e : G->edges_of(n)) {
                NodeID t = G->getEdgeTarget(n, e);
                if (!visited[t]) {
                    parent[t] = n;
                    stack.push(t);
                }
            }
        }
    }

    void findAllBridges() {
        LOG1 << "trying to find bridges!";
        for (NodeID vtx : G->nodes()) {
            if (!visited[vtx]) {
                findAllBridgesInCC(vtx);
            }
        }
    }

    void checkTerminalSides(std::vector<terminal> terminals) {
        for (const auto& [n, e] : bridges) {
            LOG1 << "bridge " << n << " - " << e;
            for (const auto& t : terminals) {
                NodeID tpos = t.position;
                bool on_right = (lowest[tpos] <= discovered[n]);
                LOG1 << "terminal " << tpos << " is " << on_right;
            }
        }
    }

 private:
    std::shared_ptr<mutable_graph> G;
    std::vector<bool> visited;
    std::vector<NodeID> discovered;
    std::vector<NodeID> parent;
    std::vector<NodeID> lowest;

    std::stack<NodeID> stack;
    std::vector<std::pair<NodeID, EdgeID> > bridges;
    size_t step;
};
