/******************************************************************************
 * push_relabel.h
 *
 * Source of VieCut.
 *
 * Adapted from KaHIP.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@univie.ac.at>
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <queue>
#include <utility>
#include <vector>

#include "common/definitions.h"
#include "data_structure/mutable_graph.h"
#include "tools/timer.h"

const int WORK_OP_RELABEL = 9;
const double GLOBAL_UPDATE_FRQ = 0.51;
const int WORK_NODE_TO_EDGES = 4;

class push_relabel {
 public:
    push_relabel() { }
    virtual ~push_relabel() { }

 private:
    void init(std::shared_ptr<mutable_graph> G,
              std::vector<NodeID> sources,
              NodeID source) {
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
            m_excess[flow_source] += G->getEdgeWeight(flow_source, e);
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
                if (m_G->getEdgeWeight(target, rev_e) -
                    getEdgeFlow(target, rev_e) > 0) {
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
        FlowType capacity = m_G->getEdgeWeight(source, e);
        FlowType flow = getEdgeFlow(source, e);

        FlowType amount = std::min(capacity - flow, m_excess[source]);
        NodeID target = m_G->getEdgeTarget(source, e);

        if (m_distance[source] <= m_distance[target] || amount == 0) return;

        setEdgeFlow(source, e, flow + amount);

        EdgeID rev_e = m_G->getReverseEdge(source, e);
        FlowType rev_flow = getEdgeFlow(target, rev_e);
        setEdgeFlow(target, rev_e, rev_flow - amount);

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
            if (m_count[m_distance[node]] == 1
                && m_distance[node] < m_G->number_of_nodes()) {
                // hence this layer will be empty after the relabel step
                gap_heuristic(m_distance[node]);
            } else {
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
            m_distance[node] = std::max(m_distance[node], m_G->n());
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
            if (m_G->getEdgeWeight(node, e) - getEdgeFlow(node, e) > 0) {
                NodeID target = m_G->getEdgeTarget(node, e);
                m_distance[node] =
                    std::min(m_distance[node], m_distance[target] + 1);
            }
            m_work++;
        }

        m_count[m_distance[node]]++;
        enqueue(node);
    }

    std::vector<NodeID> computeSourceSet(const std::vector<NodeID>& sources,
                                         NodeID curr_source) {
        std::vector<NodeID> source_set;
        // perform bfs starting from source set
        source_set.clear();
        NodeID src = sources[curr_source];

        for (NodeID node : m_G->nodes()) {
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

            for (EdgeID e : m_G->edges_of(node)) {
                EdgeID rev_e = m_G->getReverseEdge(node, e);

                NodeID edge_source = m_G->getEdgeTarget(node, e);
                FlowType resCap = m_G->getEdgeWeight(edge_source, rev_e)
                                  - getEdgeFlow(edge_source, rev_e);
                if (resCap > 0 && !m_bfstouched[edge_source]) {
                    Q.push(edge_source);
                    m_bfstouched[edge_source] = true;
                }
            }
        }

        std::queue<NodeID> Qsrc;
        Qsrc.push(src);
        source_set.emplace_back(src);
        m_bfstouched[src] = true;
        while (!Qsrc.empty()) {
            NodeID node = Qsrc.front();
            Qsrc.pop();

            for (EdgeID e : m_G->edges_of(node)) {
                NodeID n = m_G->getEdgeTarget(node, e);
                if (!m_bfstouched[n]) {
                    source_set.emplace_back(n);
                    m_bfstouched[n] = true;
                    Qsrc.push(n);
                }
            }
        }

        return source_set;
    }

    void setEdgeFlow(NodeID n, EdgeID e, FlowType f) {
        if (m_parallel_flows) {
            edge_flow[n][e] = f;
        } else {
            m_G->setEdgeFlow(n, e, f);
        }
    }

    FlowType getEdgeFlow(NodeID n, EdgeID e) {
        if (m_parallel_flows) {
            return edge_flow[n][e];
        } else {
            return m_G->getEdgeFlow(n, e);
        }
    }

 public:
    std::vector<NodeID> callable_max_flow(std::shared_ptr<mutable_graph> G,
                                          std::vector<NodeID> sources,
                                          NodeID curr_source,
                                          bool compute_source_set) {
        // this exists to be called by std::async.
        // thus, parallel_flows is set to true   
        auto source_set = solve_max_flow_min_cut(
            G, sources, curr_source, compute_source_set, true).second;

        return source_set;
    }

    std::pair<FlowType, std::vector<NodeID> > solve_max_flow_min_cut(
        std::shared_ptr<mutable_graph> G,
        std::vector<NodeID> sources,
        NodeID curr_source,
        bool compute_source_set,
        bool parallel_flows = false) {
        for (NodeID s : sources) {
            if (s >= G->number_of_nodes()) {
                LOG1 << "source " << s << " is too large (only "
                     << G->number_of_nodes() << " nodes)";
                return std::make_pair(-1, std::vector<NodeID>());
            }
        }

        if (parallel_flows) {
            // if we have parallel flows on the same graphs,
            // we can't use edge flows on the graph.
            // in sequential cases it's faster and flow values
            // might still be useful outside of this algorithm
            edge_flow.clear();
            for (NodeID n : G->nodes()) {
                edge_flow.emplace_back(G->get_first_invalid_edge(n), 0);
            }
        } else {
            for (NodeID n : G->nodes()) {
                for (EdgeID e : G->edges_of(n)) {
                    G->setEdgeFlow(n, e, 0);
                }
            }
        }

        m_parallel_flows = parallel_flows;
        m_G = G;
        m_work = 0;
        m_num_relabels = 0;
        m_gaps = 0;
        m_pushes = 0;
        m_global_updates = 1;
        NodeID src = sources[curr_source];

        init(G, sources, curr_source);
        global_relabeling(sources, curr_source);

        int work_todo = WORK_NODE_TO_EDGES * G->number_of_nodes()
                        + G->number_of_edges();
        // main loop
        while (!m_Q.empty()) {
            NodeID v = m_Q.front();
            m_Q.pop();
            m_active[v] = false;
            discharge(v);

            if (m_work > GLOBAL_UPDATE_FRQ * work_todo) {
                global_relabeling(sources, curr_source);
                m_work = 0;
                m_global_updates++;
            }
        }

        FlowType total_flow = 0;
        // return value of flow
        for (NodeID n : sources) {
            if (n != src) {
                total_flow += m_excess[n];
            }
        }

        std::vector<NodeID> source_set;

        if (compute_source_set) {
            source_set = computeSourceSet(sources, curr_source);
        }

        LOGC(extended_logs) << "updates " << m_global_updates
                            << " relabel " << m_num_relabels
                            << " pushes " << m_pushes;
        return std::make_pair(total_flow, source_set);
    }

 private:
    std::vector<FlowType> m_excess;
    std::vector<NodeID> m_distance;
    std::vector<bool> m_active;   // store which nodes are in the queue already
    std::vector<int> m_count;
    std::queue<NodeID> m_Q;
    std::vector<bool> m_bfstouched;
    std::vector<bool> m_already_contracted;
    std::vector<std::vector<FlowType> > edge_flow;
    int m_num_relabels;
    int m_gaps;
    int m_global_updates;
    int m_pushes;
    int m_work;
    std::shared_ptr<mutable_graph> m_G;
    static const bool extended_logs = false;
    bool m_parallel_flows;
};
