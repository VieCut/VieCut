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
#include <vector>

#include "algorithms/global_mincut/cactus/balanced_cut_dfs.h"
#include "common/configuration.h"
#include "data_structure/mutable_graph.h"
#include "io/graph_io.h"

class most_balanced_minimum_cut {
 public:
    void findCutFromCactus(std::shared_ptr<mutable_graph> G,
                           EdgeWeight mincut,
                           std::shared_ptr<graph_access> original_graph) {
        if (!configuration::getConfig()->save_cut) {
            LOG1 << "Error: can't find most balanced minimum cut "
                 << "when save_cut is not set";
            exit(1);
        }

        if (mincut == 0) {
            LOG1 << "G has multiple connected components and mincut is 0.";
            return;
        }

        balanced_cut_dfs dfs(G, mincut);
        auto [n1, e1, n2, e2] = dfs.runDFS();

        NodeID rev_n1 = G->getEdgeTarget(n1, e1);
        EdgeID rev_e1 = G->getReverseEdge(n1, e1);
        NodeID rev_n2 = G->getEdgeTarget(n2, e2);
        EdgeID rev_e2 = G->getReverseEdge(n2, e2);

        if (configuration::getConfig()->output_path != "") {
            LOG1 << "Printing output to file "
                 << configuration::getConfig()->output_path;

            std::queue<NodeID> q;
            q.push(0);

            for (NodeID vtx = 0; vtx < G->getOriginalNodes(); ++vtx) {
                original_graph->setNodeInCut(vtx, false);
            }

            std::vector<bool> checked(G->n(), false);

            checked[0] = true;

            while (!q.empty()) {
                NodeID v = q.front();
                q.pop();

                for (NodeID o : G->containedVertices(v)) {
                    original_graph->setNodeInCut(o, true);
                }

                for (EdgeID edge : G->edges_of(v)) {
                    if ((v == n1 && edge == e1) || (v == n2 && edge == e2)
                        || (v == rev_n1 && edge == rev_e1)
                        || (v == rev_n2 && edge == rev_e2)) {
                        // don't pass cut edges, we only want to find one side
                        continue;
                    }

                    NodeID tgt = G->getEdgeTarget(v, edge);
                    if (!checked[tgt]) {
                        q.push(tgt);
                        checked[tgt] = true;
                    }
                }
            }

            graph_io::writeCut(original_graph,
                               configuration::getConfig()->output_path);
        }
    }
};
