/******************************************************************************
 * noi_minimum_cut.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2017-2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/fifo_node_bucket_pq.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "data_structure/priority_queues/node_bucket_pq.h"
#include "data_structure/priority_queues/vecMaxNodeHeap.h"
#include "minimum_cut.h"
#include "tools/timer.h"

#ifdef PARALLEL
#include "parallel/coarsening/contract_graph.h"
#include "parallel/coarsening/contraction_tests.h"
#include "parallel/data_structure/union_find.h"
#else
#include "coarsening/contract_graph.h"
#include "coarsening/contraction_tests.h"
#include "data_structure/union_find.h"
#endif

#include "definitions.h"
#include "minimum_cut_helpers.h"
#include "tools/random_functions.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <data_structure/priority_queues/bucket_pq.h>
#include <functional>
#include <memory>
#include <unordered_map>

class noi_minimum_cut : public minimum_cut
{
public:
    noi_minimum_cut() { }
    ~noi_minimum_cut() { }
    static constexpr bool debug = false;
    static constexpr bool timing = true;

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G_ptr, bool save_cut) {
        return perform_minimum_cut(G_ptr, save_cut, false, "default");
    }

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G_ptr, bool save_cut, bool indirect) {
        return perform_minimum_cut(G_ptr, save_cut, indirect, "default");
    }

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G_ptr, bool save_cut, std::string pq) {
        return perform_minimum_cut(G_ptr, save_cut, false, pq);
    }

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G,
                                   bool save_cut,
                                   bool __attribute__ ((unused)) indirect,
                                   std::string pq) {

        if (!minimum_cut_helpers::graphValid(G))
            return -1;

        std::vector<std::shared_ptr<graph_access> > graphs;

        timer t;
        EdgeWeight global_mincut = G->getMinDegree();
        graphs.push_back(G);

        minimum_cut_helpers::setInitialCutValues(graphs, save_cut);

        while (graphs.back()->number_of_nodes() > 2 && global_mincut > 0) {

            std::vector<std::pair<NodeID, NodeID> > contractable;

            union_find uf(graphs.back()->number_of_nodes());
            timer time;
            global_mincut = std::min(global_mincut, modified_capforest(graphs.back(), global_mincut, uf, graphs, save_cut, pq));
            std::vector<std::vector<NodeID> > reverse_mapping(uf.n());
            std::vector<NodeID> mapping(graphs.back()->number_of_nodes());
            std::vector<NodeID> part(graphs.back()->number_of_nodes(), UNDEFINED_NODE);
            NodeID current_pid = 0;

            for (NodeID n : graphs.back()->nodes()) {
                NodeID part_id = uf.Find(n);
                if (part[part_id] == UNDEFINED_NODE) {
                    part[part_id] = current_pid++;
                }
                mapping[n] = part[part_id];
                graphs.back()->setPartitionIndex(n, part[part_id]);
                reverse_mapping[part[part_id]].push_back(n);
            }

            graphs.push_back(contraction::contractGraph(graphs.back(), mapping, reverse_mapping.size(), reverse_mapping));
            minimum_cut_helpers::updateCutValueAfterContraction(graphs, global_mincut, save_cut);
        }

        if (!indirect && save_cut)
            minimum_cut_helpers::retrieveMinimumCut(graphs);

        return global_mincut;
    }

    static priority_queue_interface * selectPq(std::string pq_str, EdgeWeight mincut, std::shared_ptr<graph_access> G, bool limit) {
        priority_queue_interface* pq;

        size_t max_deg = G->getMaxDegree();

        if (!limit) {
            mincut = max_deg;
        }

        if (pq_str == "default") {
            if (mincut > 10000 && mincut > G->number_of_nodes()) {
                pq = new vecMaxNodeHeap(G->number_of_nodes());
            }
            else {
                pq = new fifo_node_bucket_pq(G->number_of_nodes(), mincut);
            }
        }
        else {
            if (pq_str == "bqueue")
                pq = new fifo_node_bucket_pq(G->number_of_nodes(), mincut);
            else {
                if (pq_str == "heap")
                    pq = new vecMaxNodeHeap(G->number_of_nodes());
                else {
                    if (pq_str == "bstack") {
                        pq = new node_bucket_pq(G->number_of_nodes(), mincut);
                    }
                    else {
                        std::cerr << "unknown pq type " << pq_str << std::endl;
                        exit(1);
                        pq = new maxNodeHeap();
                    }
                }
            }
        }
        return pq;
    }

    EdgeWeight modified_capforest(
        std::shared_ptr<graph_access> G, EdgeWeight mincut,
        union_find& uf,
        std::vector<std::shared_ptr<graph_access> >& all_graphs,
        bool save_cut,
        std::string pq_str = "default",
        bool limit = true) {

        std::string aft;

        size_t alphamin = mincut;
        all_graphs.size();

        priority_queue_interface* pq = selectPq(pq_str, mincut, G, limit);

        std::vector<bool> visited(G->number_of_nodes(), false);
        std::vector<bool> seen(G->number_of_nodes(), false);
        std::vector<EdgeWeight> r_v(G->number_of_nodes(), 0);

        NodeID starting_node = random_functions::next() % G->number_of_nodes();

        NodeID current_node = starting_node;

        pq->insert(current_node, 0);
        EdgeWeight alpha = 0;

        timer t;
        size_t elements = 0;
        while (!pq->empty()) {
            elements++;
            current_node = pq->deleteMax();
            visited[current_node] = true;
            alpha += G->getWeightedNodeDegree(current_node) - (2 * r_v[current_node]);
            if (alpha < alphamin && (alpha > 0 || elements < G->number_of_nodes())) {

                alphamin = alpha;
                if (save_cut) {
                    for (NodeID idx : all_graphs[0]->nodes()) {
                        NodeID coarseID = idx;
                        for (size_t lv = 0; lv < all_graphs.size() - 1; ++lv) {
                            coarseID = all_graphs[lv]->getPartitionIndex(coarseID);
                        }

                        if (!visited[coarseID]) {
                            all_graphs[0]->setNodeInCut(idx, false);
                        }
                        else {
                            all_graphs[0]->setNodeInCut(idx, true);
                        }
                    }
                }
            }

            if (limit) {
                processVertexLimited(G, current_node, visited, r_v, seen, mincut, uf, pq);
            }
            else {
                processVertexUnlimited(G, current_node, visited, r_v, seen, mincut, uf, pq);
            }
        }
        delete pq;
        return alphamin;
    }

    void processVertexLimited(std::shared_ptr<graph_access> G,
                              NodeID current_node,
                              std::vector<bool>& visited,
                              std::vector<EdgeWeight>& r_v,
                              std::vector<bool>& seen,
                              EdgeWeight mincut,
                              union_find& uf,
                              priority_queue_interface* pq) {

        for (EdgeID e : G->edges_of(current_node)) {
            NodeID tgt = G->getEdgeTarget(e);

            if (!visited[tgt]) {

                bool increase = false;

                if (r_v[tgt] < mincut) {
                    increase = true;
                    if ((r_v[tgt] + G->getEdgeWeight(e)) >= mincut) {
                        uf.Union(current_node, tgt);
                    }
                }

                r_v[tgt] += G->getEdgeWeight(e);

                size_t new_rv = std::min(r_v[tgt], mincut);

                if (seen[tgt]) {
                    if (increase && !visited[tgt]) {
                        pq->increaseKey(tgt, new_rv);
                    }
                }
                else {
                    seen[tgt] = true;
                    pq->insert(tgt, new_rv);
                }
            }
        }
    }

    void processVertexUnlimited(std::shared_ptr<graph_access> G,
                                NodeID current_node,
                                std::vector<bool>& visited,
                                std::vector<EdgeWeight>& r_v,
                                std::vector<bool>& seen,
                                EdgeWeight mincut,
                                union_find& uf,
                                priority_queue_interface* pq) {

        for (EdgeID e : G->edges_of(current_node)) {
            NodeID tgt = G->getEdgeTarget(e);

            if (!visited[tgt]) {
                if (r_v[tgt] < mincut) {
                    if ((r_v[tgt] + G->getEdgeWeight(e)) >= mincut) {
                        uf.Union(current_node, tgt);
                    }
                }

                r_v[tgt] += G->getEdgeWeight(e);

                if (seen[tgt]) {
                    if (!visited[tgt]) {
                        pq->increaseKey(tgt, r_v[tgt]);
                    }
                }
                else {
                    seen[tgt] = true;
                    pq->insert(tgt, r_v[tgt]);
                }
            }
        }
    }
};
