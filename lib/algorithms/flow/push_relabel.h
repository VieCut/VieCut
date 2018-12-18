/******************************************************************************
 * push_relabel.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once

#include "data_structure/flow_graph.h"
#include "definitions.h"
#include "tools/timer.h"
#include <iostream>
#include <memory>

const int WORK_OP_RELABEL = 9;
const double GLOBAL_UPDATE_FRQ = 0.51;
const int WORK_NODE_TO_EDGES = 4;

class push_relabel
{
public:
    push_relabel() { }
    virtual ~push_relabel() { }

    void init(std::shared_ptr<flow_graph> G, std::vector<NodeID> sources, NodeID source) {

        m_Q = std::queue<NodeID>();
        m_excess.clear();
        m_distance.clear();
        m_active.clear();
        m_count.clear();
        m_bfstouched.clear();

        m_excess.resize(G->number_of_nodes(), 0);
        m_distance.resize(G->number_of_nodes(), 0);
        m_active.resize(G->number_of_nodes(), false);
        m_count.resize(2 * G->number_of_nodes(), 0);
        m_bfstouched.resize(G->number_of_nodes(), false);
        m_already_contracted.resize(G->number_of_nodes(), false);

        m_count[0] = G->number_of_nodes() - 1;
        m_count[G->number_of_nodes()] = 1;

        NodeID flow_source = sources[source];

        m_distance[flow_source] = G->number_of_nodes();

        for (NodeID n : sources) {
            m_active[n] = true;
        }

        for (EdgeID e : G->edges_of(flow_source)) {
            m_excess[flow_source] += G->getEdgeCapacity(flow_source, e);
            push(flow_source, e);
        }
    }

    // perform a backward bfs in the residual starting at the sink
    // to update distance labels
    void global_relabeling(std::vector<NodeID> sources, NodeID source) {

        std::queue<NodeID> Q;
        NodeID flow_source = sources[source];

        for (NodeID n : m_G->nodes()) {
            m_bfstouched[n] = false;
            m_distance[n] = std::max(m_distance[n], m_G->number_of_nodes());
        }

        for (NodeID sink : sources) {
            if (sink == flow_source)
                continue;

            Q.push(sink);
            m_bfstouched[sink] = true;
            m_distance[sink] = 0;
        }

        m_bfstouched[flow_source] = true;

        NodeID node = 0;
        while (!Q.empty()) {
            node = Q.front();
            Q.pop();

            for (EdgeID e : m_G->edges_of(node)) {
                NodeID target = m_G->getEdgeTarget(node, e);
                if (m_bfstouched[target]) continue;

                EdgeID rev_e = m_G->getReverseEdge(node, e);
                if (m_G->getEdgeCapacity(target, rev_e) - m_G->getEdgeFlow(target, rev_e) > 0) {
                    m_count[m_distance[target]]--;
                    m_distance[target] = m_distance[node] + 1;
                    m_count[m_distance[target]]++;
                    Q.push(target);
                    m_bfstouched[target] = true;
                }
            }
        }
    }

    // push flow from source to target if possible
    void push(NodeID source, EdgeID e) {
        m_pushes++;
        FlowType capacity = m_G->getEdgeCapacity(source, e);
        FlowType flow = m_G->getEdgeFlow(source, e);
        FlowType amount = std::min((long long)(capacity - flow), m_excess[source]);
        NodeID target = m_G->getEdgeTarget(source, e);

        if (m_distance[source] <= m_distance[target] || amount == 0) return;

        m_G->setEdgeFlow(source, e, flow + amount);

        EdgeID rev_e = m_G->getReverseEdge(source, e);
        FlowType rev_flow = m_G->getEdgeFlow(target, rev_e);
        m_G->setEdgeFlow(target, rev_e, rev_flow - amount);

        m_excess[source] -= amount;
        m_excess[target] += amount;

        enqueue(target);
    }

    // put a vertex in the FIFO queue of the global_mincut
    void enqueue(NodeID target) {
        if (m_active[target]) return;
        if (m_excess[target] > 0) {
            m_active[target] = true;
            // m_Q.push(target, m_distance[target]);
            m_Q.push(target);
        }
    }

    // try to push as much excess as possible out of the node node
    void discharge(NodeID node) {
        for (EdgeID e : m_G->edges_of(node)) {
            if (m_excess[node] == 0)
                break;

            push(node, e);
        }

        if (m_excess[node] > 0) {
            if (m_count[m_distance[node]] == 1 && m_distance[node] < m_G->number_of_nodes()) {
                // hence this layer will be empty after the relabel step
                gap_heuristic(m_distance[node]);
            }
            else {
                relabel(node);
            }
        }
    }

    // gap heuristic
    void gap_heuristic(NodeID level) {
        m_gaps++;
        for (NodeID node : m_G->nodes()) {
            if (m_distance[node] < level) continue;
            m_count[m_distance[node]]--;
            m_distance[node] = std::max(m_distance[node], m_G->number_of_nodes());
            m_count[m_distance[node]]++;
            enqueue(node);
        }
    }

    // relabel a node with respect to its
    // neighboring nodes
    void relabel(NodeID node) {
        m_work += WORK_OP_RELABEL;
        m_num_relabels++;

        m_count[m_distance[node]]--;
        m_distance[node] = 2 * m_G->number_of_nodes();

        for (EdgeID e : m_G->edges_of(node)) {
            if (m_G->getEdgeCapacity(node, e) - m_G->getEdgeFlow(node, e) > 0) {
                NodeID target = m_G->getEdgeTarget(node, e);
                m_distance[node] = std::min(m_distance[node], m_distance[target] + 1);
            }
            m_work++;
        }

        m_count[m_distance[node]]++;
        enqueue(node);
    }

    FlowType solve_max_flow_min_cut(std::shared_ptr<flow_graph> G,
                                    std::vector<NodeID> sources,
                                    NodeID curr_source,
                                    bool compute_source_set,
                                    std::vector<NodeID>& source_set) {

        for (NodeID s : sources) {
            if (s >= G->number_of_nodes()) {
                LOG1 << "source " << s << " is too large (only " << G->number_of_nodes() << " nodes)";
                return -1;
            }
        }

        m_G = G;
        m_work = 0;
        m_num_relabels = 0;
        m_gaps = 0;
        m_pushes = 0;
        m_global_updates = 1;
        NodeID src = sources[curr_source];

        init(G, sources, curr_source);
        global_relabeling(sources, curr_source);

        int work_todo = WORK_NODE_TO_EDGES * G->number_of_nodes() + G->number_of_edges();
        // main loop
        while (!m_Q.empty()) {
            NodeID v = m_Q.front();
            m_Q.pop();
            m_active[v] = false;
            discharge(v);

            if (m_work > GLOBAL_UPDATE_FRQ * work_todo) {
                // global_relabeling(sources, curr_source);
                m_work = 0;
                m_global_updates++;
            }
        }

        if (compute_source_set) {
            // perform bfs starting from source set
            source_set.clear();

            for (NodeID node : G->nodes()) {
                m_bfstouched[node] = false;
            }

            std::queue<NodeID> Q;
            for (NodeID tgt : sources) {
                if (tgt == src)
                    continue;

                Q.push(tgt);
                m_bfstouched[tgt] = true;
            }

            while (!Q.empty()) {
                NodeID node = Q.front();
                Q.pop();

                for (EdgeID e : G->edges_of(node)) {

                    EdgeID rev_e = G->getReverseEdge(node, e);

                    NodeID edge_source = G->getEdgeTarget(node, e);
                    FlowType resCap = G->getEdgeCapacity(edge_source, rev_e) - G->getEdgeFlow(edge_source, rev_e);
                    if (resCap > 0 && !m_bfstouched[edge_source]) {
                        Q.push(edge_source);
                        m_bfstouched[edge_source] = true;
                    }
                }
            }
        }

        FlowType total_flow = 0;
        // return value of flow
        for (NodeID n : sources) {
            if (n != src) {
                total_flow += m_excess[n];
            }
        }

        std::queue<NodeID> Q;
        Q.push(src);
        source_set.emplace_back(src);
        m_bfstouched[src] = true;
        m_already_contracted[src] = true;

        while (!Q.empty()) {
            NodeID node = Q.front();
            Q.pop();

            for (EdgeID e : G->edges_of(node)) {
                NodeID n = G->getEdgeTarget(node, e);
                if (!m_bfstouched[n] && !m_already_contracted[n]) {
                    source_set.emplace_back(n);
                    m_bfstouched[n] = true;
                    m_already_contracted[n] = true;
                    Q.push(n);
                }
            }
        }

        LOGC(extended_logs) << "updates " << m_global_updates << " relabel " << m_num_relabels << " pushes " << m_pushes;
        return total_flow;
    }

private:
    std::vector<long long> m_excess;
    std::vector<NodeID> m_distance;
    std::vector<bool> m_active;          // store which nodes are in the queue already
    std::vector<int> m_count;
    std::queue<NodeID> m_Q;
    std::vector<bool> m_bfstouched;
    std::vector<bool> m_already_contracted;
    // highest_label_queue    m_Q;
    int m_num_relabels;
    int m_gaps;
    int m_global_updates;
    int m_pushes;
    int m_work;
    std::shared_ptr<flow_graph> m_G;
    static const bool extended_logs = false;
};
