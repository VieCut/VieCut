/******************************************************************************
 * most_balanced_minimum_cut.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <memory>
#include <queue>
#include <utility>
#include <vector>

#include "algorithms/global_mincut/cactus/balanced_cut_dfs.h"
#include "common/configuration.h"
#include "data_structure/mutable_graph.h"
#include "io/graph_io.h"

class most_balanced_minimum_cut {
 public:
    std::unordered_set<EdgeID> findCutFromCactus(
        std::shared_ptr<mutable_graph> G, EdgeWeight mincut,
        std::shared_ptr<graph_access> original_graph) {
        if (!configuration::getConfig()->save_cut) {
            LOG1 << "Error: can't find most balanced minimum cut "
                 << "when save_cut is not set";
            exit(1);
        }

        std::unordered_set<EdgeID> originalBestcutEdges;

        if (mincut == 0) {
            LOG1 << "G has multiple connected components and mincut is 0.";
            return originalBestcutEdges;
        }

        balanced_cut_dfs dfs(original_graph, G, mincut);
        auto [n1, e1, n2, e2, bestcutInCycle] = dfs.runDFS();

        NodeID rev_n1 = G->getEdgeTarget(n1, e1);
        NodeID rev_n2 = G->getEdgeTarget(n2, e2);

        for (NodeID n : original_graph->nodes()) {
            original_graph->setNodeInCut(n, false);
        }

        std::queue<NodeID> q;
        std::vector<bool> checked(G->n(), false);
        q.push(n1);
        checked[n1] = true;
        if (n1 != n2) {
            q.push(n2);
            checked[n2] = true;
        }
        // we want to set one side of most balanced cut to inCut
        // setting the targets of the edges to true makes sure the bfs doesn't
        // go into other side of most balanced cut
        checked[rev_n1] = true;
        checked[rev_n2] = true;

        while (!q.empty()) {
            NodeID top = q.front();
            q.pop();
            for (NodeID v : G->containedVertices(top)) {
                original_graph->setNodeInCut(v, true);
            }
            for (EdgeID e : G->edges_of(top)) {
                NodeID t = G->getEdgeTarget(top, e);
                if (!checked[t]) {
                    q.push(t);
                    checked[t] = true;
                }
            }
        }

        for (NodeID on : original_graph->nodes()) {
            for (EdgeID oe : original_graph->edges_of(on)) {
                NodeID ot = original_graph->getEdgeTarget(oe);
                if (original_graph->getNodeInCut(on)
                    != original_graph->getNodeInCut(ot)) {
                    originalBestcutEdges.insert(oe);
                }
            }
        }

        if (bestcutInCycle)
            LOG0 << "bestcut is in cycle";

        if (configuration::getConfig()->output_path != "") {
            LOG1 << "Printing output to file "
                 << configuration::getConfig()->output_path;

            graph_io::writeCut(original_graph,
                               configuration::getConfig()->output_path);
        }
        return originalBestcutEdges;
    }
};
