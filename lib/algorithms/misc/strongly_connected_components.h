/******************************************************************************
 * strongly_connected_components.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@univie.ac.at>
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#ifndef STRONGLY_CONNECTED_COMPONENTS_7ZJ8233R
#define STRONGLY_CONNECTED_COMPONENTS_7ZJ8233R

#include <memory>
#include <stack>
#include <vector>

#include "data_structure/graph_access.h"
#include "definitions.h"
#include "tools/graph_extractor.h"

class strongly_connected_components
{
public:
    const static bool debug = false;

    strongly_connected_components() { }
    virtual ~strongly_connected_components() { }

    size_t strong_components(graph_access& G, std::vector<int>& comp_num) {

        m_dfsnum.resize(G.number_of_nodes());
        m_comp_num.resize(G.number_of_nodes());
        m_dfscount = 0;
        m_comp_count = 0;

        for (NodeID node : G.nodes()) {
            // comp_num[node] = -1;
            m_comp_num[node] = -1;
            m_dfsnum[node] = -1;
        }

        for (NodeID node : G.nodes()) {
            if (m_dfsnum[node] == -1) {
                explicit_scc_dfs(node, G);
            }
        }

        for (NodeID node : G.nodes()) {
            comp_num[node] = m_comp_num[node];
        }

        return m_comp_count;
    }

    void explicit_scc_dfs(NodeID node, graph_access& G) {
        iteration_stack.push(std::pair<NodeID, EdgeID>(node, G.get_first_edge(node)));

        // make node a tentative scc of its own
        m_dfsnum[node] = m_dfscount++;
        m_unfinished.push(node);
        m_roots.push(node);

        while (!iteration_stack.empty()) {
            NodeID current_node = iteration_stack.top().first;
            EdgeID current_edge = iteration_stack.top().second;
            iteration_stack.pop();

            for (EdgeID e : G.edges_of_starting_at(current_node, current_edge)) {
                NodeID target = G.getEdgeTarget(e);
                // explore edge (node, target)
                if (m_dfsnum[target] == -1) {

                    iteration_stack.push(std::pair<NodeID, EdgeID>(current_node, e));
                    iteration_stack.push(std::pair<NodeID, EdgeID>(target, G.get_first_edge(target)));

                    m_dfsnum[target] = m_dfscount++;
                    m_unfinished.push(target);
                    m_roots.push(target);
                    break;
                }
                else if (m_comp_num[target] == -1) {
                    // merge scc's
                    while (m_dfsnum[m_roots.top()] > m_dfsnum[target]) m_roots.pop();
                }
            }

            // return from call of node node
            if (current_node == m_roots.top()) {
                NodeID w = 0;
                do {
                    w = m_unfinished.top();
                    m_unfinished.pop();
                    m_comp_num[w] = m_comp_count;
                } while (w != current_node);
                m_comp_count++;
                m_roots.pop();
            }
        }
    }

    std::shared_ptr<graph_access> largest_scc(std::shared_ptr<graph_access>& G) {
        std::vector<int> components(G->number_of_nodes());
        auto ct = strong_components(*G, components);
        LOG << "count of connected components: " << ct;

        std::vector<uint64_t> compsizes(static_cast<unsigned long>(ct));
        for (int component : components) {
            ++compsizes[component];
        }

        auto max_size = std::max_element(compsizes.begin(), compsizes.end());
        int max_comp = static_cast<int>(max_size - compsizes.begin());

        for (NodeID n : G->nodes()) {
            if (components[n] == max_comp) {
                G->setPartitionIndex(n, 0);
            }
            else {
                G->setPartitionIndex(n, 1);
            }
        }

        graph_extractor ge;
        return ge.extract_block(G, 0).first;
    }

private:
    int m_dfscount;
    size_t m_comp_count;

    std::vector<int> m_dfsnum;
    std::vector<int> m_comp_num;
    std::stack<NodeID> m_unfinished;
    std::stack<NodeID> m_roots;
    std::stack<std::pair<NodeID, EdgeID> > iteration_stack;
};

#endif /* end of include guard: STRONGLY_CONNECTED_COMPONENTS_7ZJ8233R */
