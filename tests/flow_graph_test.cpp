#include "io/graph_io.h"
#include <gtest/gtest.h>
#include <type_traits>

#include "data_structure/flow_graph.h"
#include "data_structure/graph_access.h"
#include "io/graph_io.h"

TEST(FlowGraphTest, EmptyGraph) {
    std::shared_ptr<graph_access> G = std::make_shared<graph_access>();
    std::shared_ptr<flow_graph> fG = graph_io::createFlowGraph(G);

    ASSERT_EQ(fG->number_of_nodes(), 0);
    ASSERT_EQ(fG->number_of_edges(), 0);
}

TEST(FlowGraphTest, NoEdges) {

    std::shared_ptr<graph_access> G = std::make_shared<graph_access>();
    G->start_construction(10, 0);
    for (size_t i = 0; i < 10; ++i) {
        G->new_node();
    }
    G->finish_construction();
    std::shared_ptr<flow_graph> fG = graph_io::createFlowGraph(G);

    ASSERT_EQ(fG->number_of_nodes(), 10);
    ASSERT_EQ(fG->number_of_edges(), 0);
}

TEST(FlowGraphTest, Clique) {

    std::shared_ptr<graph_access> G = std::make_shared<graph_access>();
    G->start_construction(10, 0);
    for (size_t i = 0; i < 10; ++i) {
        G->new_node();
        for (size_t j = 0; j < 10; ++j) {
            if (i != j)
                G->new_edge(i, j);
        }
    }
    G->finish_construction();
    std::shared_ptr<flow_graph> fG = graph_io::createFlowGraph(G);

    ASSERT_EQ(fG->number_of_nodes(), 10);
    ASSERT_EQ(fG->number_of_edges(), 180);

    for (size_t i = 0; i < 10; ++i) {
        ASSERT_EQ(fG->getCapacity(i), 9);
        for (size_t j = 0; j < 18; ++j) {

            NodeID tgt = 0;
            bool before = (j < i);
            bool behind = (j >= i + 9);
            bool forward = !before && !behind;

            if (before)
                tgt = j;
            if (forward)
                tgt = j - i + (j >= (i + i));
            if (behind)
                tgt = j - 8;

            ASSERT_EQ(fG->getEdgeTarget(i, j), tgt);
            ASSERT_EQ(fG->getEdgeCapacity(i, j), forward);
        }
    }
}

TEST(FlowGraphTest, SparseGraph) {
    std::shared_ptr<graph_access> G = std::make_shared<graph_access>();
    G->start_construction(100, 0);
    for (size_t i = 0; i < 100; ++i) {
        G->new_node();
        G->new_edge(i, (i + 50) % 100);
    }
    G->finish_construction();
    std::shared_ptr<flow_graph> fG = graph_io::createFlowGraph(G);
    ASSERT_EQ(fG->number_of_nodes(), 100);
    ASSERT_EQ(fG->number_of_edges(), 200);
    for (size_t i = 0; i < 100; ++i) {
        ASSERT_EQ(fG->getCapacity(i), 1);

        if (i < 50) {
            ASSERT_EQ(fG->getEdgeTarget(i, 0), i + 50);
            ASSERT_EQ(fG->getEdgeCapacity(i, 0), 1);
            ASSERT_EQ(fG->getEdgeTarget(i, 1), i + 50);
            ASSERT_EQ(fG->getEdgeCapacity(i, 1), 0);
        }
        else {
            ASSERT_EQ(fG->getEdgeTarget(i, 0), i - 50);
            ASSERT_EQ(fG->getEdgeCapacity(i, 0), 0);
            ASSERT_EQ(fG->getEdgeTarget(i, 1), i - 50);
            ASSERT_EQ(fG->getEdgeCapacity(i, 1), 1);
        }
    }
}
