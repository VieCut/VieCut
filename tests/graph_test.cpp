/******************************************************************************
 * graph_test.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <stdio.h>

#include <memory>
#include <string>
#include <vector>

#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "gtest/gtest_pred_impl.h"
#include "io/graph_io.h"
#include "tlx/logger.hpp"

graph_access make_circle() {
    graph_access G;
    G.start_construction(3, 3);
    G.new_node();
    G.new_edge(0, 1);
    G.new_edge(0, 2);
    G.new_node();
    G.new_edge(1, 0);
    G.new_edge(1, 2);
    G.new_node();
    G.new_edge(2, 0);
    G.new_edge(2, 1);
    G.finish_construction();
    return G;
}

TEST(Graph_Test, Create_Empty) {
    graph_access G;
    G.start_construction(0, 0);
    G.finish_construction();
    ASSERT_EQ(G.number_of_nodes(), 0);
    ASSERT_EQ(G.number_of_edges(), 0);
}

TEST(Graph_Test, AddVertices) {
    graph_access G;
    G.start_construction(3, 0);
    G.new_node();
    G.new_node();
    G.new_node();
    G.finish_construction();

    ASSERT_EQ(G.number_of_nodes(), 3);
    ASSERT_EQ(G.number_of_edges(), 0);
}

TEST(Graph_Test, AddEdges) {
    graph_access G = make_circle();
    ASSERT_EQ(G.number_of_edges(), 6);
    ASSERT_EQ(G.number_of_nodes(), 3);
}

TEST(Graph_Test, Degrees) {
    graph_access G = make_circle();
    ASSERT_EQ(G.getMaxDegree(), 2);
    ASSERT_EQ(G.getMinDegree(), 2);
    ASSERT_EQ(G.getMaxUnweightedDegree(), 2);
}

TEST(Graph_Test, WeightedGraph) {
    graph_access G = make_circle();

    // edge from 0 to 1
    G.setEdgeWeight(0, 3);
    G.setEdgeWeight(2, 3);

    // edge from 0 to 2
    G.setEdgeWeight(1, 10);
    G.setEdgeWeight(4, 10);

    // edge from 1 to 2
    G.setEdgeWeight(3, 20);
    G.setEdgeWeight(5, 20);

    ASSERT_EQ(G.getMaxDegree(), 30);
    ASSERT_EQ(G.getMinDegree(), 13);
    ASSERT_EQ(G.getMaxUnweightedDegree(), 2);
}

TEST(Graph_Test, Loops) {
    graph_access G = make_circle();
    size_t ed = 0, no = 0, ed2 = 0;

    for (EdgeID e : G.edges()) {
        LOG0 << e;
        ed++;
    }

    for (NodeID n : G.nodes()) {
        no++;
        for (EdgeID e : G.edges_of(n)) {
            LOG0 << n << e;
            ed2++;
        }
    }

    ASSERT_EQ(no, G.number_of_nodes());
    ASSERT_EQ(ed, G.number_of_edges());
    ASSERT_EQ(ed2, G.number_of_edges());
}

TEST(Graph_Test, GraphFromFile) {
    graphAccessPtr G =
        graph_io::readGraphWeighted(std::string(VIECUT_PATH)
                                    + "/graphs/small.metis");

    ASSERT_EQ(G->number_of_nodes(), 8);
    ASSERT_EQ(G->number_of_edges(), 28);
    ASSERT_EQ(G->getMinDegree(), 3);
    ASSERT_EQ(G->getMaxDegree(), 4);
}

TEST(Graph_Test, HackyGraphCreation) {
    graphAccessPtr G1 =
        graph_io::readGraphWeighted(std::string(VIECUT_PATH)
                                    + "/graphs/small.metis");

    graphAccessPtr G2 = std::make_shared<graph_access>();

    G2->start_construction(8, 28);
    G2->resize_m(28);
    G2->new_node_hacky(4);
    G2->new_node_hacky(7);
    G2->new_node_hacky(10);
    G2->new_node_hacky(14);
    G2->new_node_hacky(18);
    G2->new_node_hacky(21);
    G2->new_node_hacky(24);
    G2->new_node_hacky(28);
    G2->new_edge_and_reverse(0, 1, 0, 4, 1);
    G2->new_edge_and_reverse(0, 2, 1, 7, 1);
    G2->new_edge_and_reverse(0, 3, 2, 10, 1);
    G2->new_edge_and_reverse(0, 4, 3, 14, 1);
    G2->new_edge_and_reverse(1, 2, 5, 8);
    G2->new_edge_and_reverse(1, 3, 6, 11);
    G2->new_edge_and_reverse(2, 3, 9, 12);
    G2->new_edge_and_reverse(3, 7, 13, 24);
    G2->new_edge_and_reverse(4, 5, 15, 18);
    G2->new_edge_and_reverse(4, 6, 16, 21);
    G2->new_edge_and_reverse(4, 7, 17, 25);
    G2->new_edge_and_reverse(5, 6, 19, 22);
    G2->new_edge_and_reverse(5, 7, 20, 26);
    G2->new_edge_and_reverse(6, 7, 23, 27);
    G2->finish_construction();

    for (NodeID n : G1->nodes()) {
        for (EdgeID e : G1->edges_of(n)) {
            ASSERT_EQ(G1->getEdgeWeight(e), G2->getEdgeWeight(e));
            ASSERT_EQ(G1->getEdgeTarget(e), G2->getEdgeTarget(e));
        }
    }

    ASSERT_EQ(G1->getMaxDegree(), G2->getMaxDegree());
    ASSERT_EQ(G1->getMinDegree(), G2->getMinDegree());
}

TEST(Graph_Test, ReadWriteEqual) {
    std::vector<std::string> graphs = { "", "-wgt" };
    for (std::string graph : graphs) {
        std::string copystr = (std::string(VIECUT_PATH)
                               + "/graphs/copy" + graph + ".metis");

        graphAccessPtr G1 =
            graph_io::readGraphWeighted(std::string(VIECUT_PATH)
                                        + "/graphs/small" + graph + ".metis");
        graph_io::writeGraphWeighted(G1, copystr);
        graphAccessPtr G2 = graph_io::readGraphWeighted(copystr);
        remove(copystr.c_str());

        for (NodeID n : G1->nodes()) {
            for (EdgeID e : G1->edges_of(n)) {
                ASSERT_EQ(G1->getEdgeWeight(e), G2->getEdgeWeight(e));
                ASSERT_EQ(G1->getEdgeTarget(e), G2->getEdgeTarget(e));
            }
        }

        ASSERT_EQ(G1->getMaxDegree(), G2->getMaxDegree());
        ASSERT_EQ(G1->getMinDegree(), G2->getMinDegree());
    }
}

TEST(Graph_Test, ReadWeighted) {
    graphAccessPtr G =
        graph_io::readGraphWeighted(std::string(VIECUT_PATH)
                                    + "/graphs/small-wgt.metis");

    ASSERT_EQ(G->number_of_nodes(), 8);
    ASSERT_EQ(G->number_of_edges(), 28);
    ASSERT_EQ(G->getMinDegree(), 10);
    ASSERT_EQ(G->getMaxDegree(), 21);
    ASSERT_EQ(G->getMaxUnweightedDegree(), 4);
}
