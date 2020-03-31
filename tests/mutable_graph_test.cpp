/******************************************************************************
 * mutable_graph_test.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2018-2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <random>
#include <vector>

#include "data_structure/mutable_graph.h"
#include "gtest/gtest.h"
#include "io/graph_io.h"
#include "tlx/logger.hpp"

mutable_graph make_circle() {
    mutable_graph G;
    G.start_construction(3);
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

TEST(Mutable_Graph_Test, Create_Empty) {
    mutable_graph G;
    G.start_construction(0);
    G.finish_construction();
    ASSERT_EQ(G.number_of_nodes(), 0);
    ASSERT_EQ(G.number_of_edges(), 0);
}

TEST(Mutable_Graph_Test, AddVertices) {
    mutable_graph G;
    G.start_construction(3);
    G.new_node();
    G.new_node();
    G.new_node();
    G.finish_construction();

    ASSERT_EQ(G.number_of_nodes(), 3);
    ASSERT_EQ(G.number_of_edges(), 0);
}

TEST(Mutable_Graph_Test, AddEdges) {
    mutable_graph G = make_circle();
    ASSERT_EQ(G.number_of_edges(), 6);
    ASSERT_EQ(G.number_of_nodes(), 3);
}

TEST(Mutable_Graph_Test, Degrees) {
    mutable_graph G = make_circle();
    ASSERT_EQ(G.getMaxDegree(), 2);
    ASSERT_EQ(G.getMinDegree(), 2);
    ASSERT_EQ(G.getMaxUnweightedDegree(), 2);
}

TEST(Mutable_Graph_Test, WeightedGraph) {
    mutable_graph G = make_circle();

    // edge from 0 to 1
    G.setEdgeWeight(0, 0, 3);

    // edge from 0 to 2
    G.setEdgeWeight(0, 1, 10);

    // edge from 1 to 2
    G.setEdgeWeight(1, 1, 20);

    ASSERT_EQ(G.getMaxDegree(), 30);
    ASSERT_EQ(G.getMinDegree(), 13);
    ASSERT_EQ(G.getMaxUnweightedDegree(), 2);
}

TEST(Mutable_Graph_Test, Loops) {
    mutable_graph G = make_circle();
    size_t ed = 0, no = 0, ed2 = 0;

    for (NodeID n : G.nodes()) {
        for (EdgeID e : G.edges_of(n)) {
            LOG0 << e;
            ed++;
        }
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

TEST(Mutable_Graph_Test, FromFileThroughGraphAccess) {
    std::shared_ptr<graph_access> GA = graph_io::readGraphWeighted(
        std::string(VIECUT_PATH) + "/graphs/small.metis");

    mutableGraphPtr G =
        mutable_graph::from_graph_access(GA);

    ASSERT_EQ(G->number_of_nodes(), 8);
    ASSERT_EQ(G->number_of_edges(), 28);
    ASSERT_EQ(G->getMinDegree(), 3);
    ASSERT_EQ(G->getMaxDegree(), 4);
}

TEST(Mutable_Graph_Test, ToAndFromGraphAccess) {
    mutable_graph G_real = make_circle();
    auto G = std::make_shared<mutable_graph>(make_circle());

    auto GA = G->to_graph_access();

    auto G2 = mutable_graph::from_graph_access(GA);

    for (NodeID n : G2->nodes()) {
        for (EdgeID e : G2->edges_of(n)) {
            ASSERT_EQ(G->getEdgeWeight(n, e), G2->getEdgeWeight(n, e));
            ASSERT_EQ(G->getEdgeTarget(n, e), G2->getEdgeTarget(n, e));
            ASSERT_EQ(G->getReverseEdge(n, e), G2->getReverseEdge(n, e));
        }
    }
}

TEST(Mutable_Graph_Test, DeleteEdges) {
    mutable_graph G = make_circle();

    for (NodeID n : G.nodes()) {
        ASSERT_EQ(G.getUnweightedNodeDegree(n), 2);
    }

    G.deleteEdge(0, 0);

    ASSERT_EQ(G.getUnweightedNodeDegree(0), 1);
    ASSERT_EQ(G.getUnweightedNodeDegree(1), 1);
    ASSERT_EQ(G.getUnweightedNodeDegree(2), 2);

    G.deleteEdge(0, 0);

    ASSERT_EQ(G.getUnweightedNodeDegree(0), 0);
    ASSERT_EQ(G.getUnweightedNodeDegree(1), 1);
    ASSERT_EQ(G.getUnweightedNodeDegree(2), 1);

    G.deleteEdge(1, 0);

    ASSERT_EQ(G.getUnweightedNodeDegree(0), 0);
    ASSERT_EQ(G.getUnweightedNodeDegree(1), 0);
    ASSERT_EQ(G.getUnweightedNodeDegree(2), 0);
}

TEST(Mutable_Graph_Test, DeleteWeightedEdges) {
    mutable_graph G = make_circle();

    // edge from 0 to 1
    G.setEdgeWeight(0, 0, 3);

    // edge from 0 to 2
    G.setEdgeWeight(0, 1, 10);

    // edge from 1 to 2
    G.setEdgeWeight(1, 1, 20);

    G.deleteEdge(0, 0);

    ASSERT_EQ(G.getWeightedNodeDegree(0), 10);
    ASSERT_EQ(G.getWeightedNodeDegree(1), 20);
    ASSERT_EQ(G.getWeightedNodeDegree(2), 30);

    G.deleteEdge(0, 0);

    ASSERT_EQ(G.getWeightedNodeDegree(0), 0);
    ASSERT_EQ(G.getWeightedNodeDegree(1), 20);
    ASSERT_EQ(G.getWeightedNodeDegree(2), 20);

    G.deleteEdge(1, 0);

    ASSERT_EQ(G.getWeightedNodeDegree(0), 0);
    ASSERT_EQ(G.getWeightedNodeDegree(1), 0);
    ASSERT_EQ(G.getWeightedNodeDegree(2), 0);
}

TEST(Mutable_Graph_Test, ContractEdges) {
    mutable_graph G = make_circle();

    G.contractEdge(0, 0);

    ASSERT_EQ(G.getUnweightedNodeDegree(0), 1);
    ASSERT_EQ(G.getUnweightedNodeDegree(1), 1);
    ASSERT_EQ(G.getWeightedNodeDegree(0), 2);
    ASSERT_EQ(G.getWeightedNodeDegree(1), 2);
    ASSERT_EQ(G.number_of_nodes(), 2);
    ASSERT_EQ(G.number_of_edges(), 2);

    G.contractEdge(0, 0);

    ASSERT_EQ(G.getUnweightedNodeDegree(0), 0);
    ASSERT_EQ(G.getUnweightedNodeDegree(1), 0);
    ASSERT_EQ(G.number_of_nodes(), 1);
    ASSERT_EQ(G.number_of_edges(), 0);
}

TEST(Mutable_Graph_Test, ContractEdgesWeighted) {
    mutable_graph G = make_circle();

    // edge from 0 to 1
    G.setEdgeWeight(0, 0, 3);

    // edge from 0 to 2
    G.setEdgeWeight(0, 1, 10);

    // edge from 1 to 2
    G.setEdgeWeight(1, 1, 20);

    G.contractEdge(0, 0);

    ASSERT_EQ(G.getUnweightedNodeDegree(0), 1);
    ASSERT_EQ(G.getUnweightedNodeDegree(1), 1);
    ASSERT_EQ(G.getWeightedNodeDegree(0), 30);
    ASSERT_EQ(G.getWeightedNodeDegree(1), 30);
    ASSERT_EQ(G.number_of_nodes(), 2);
    ASSERT_EQ(G.number_of_edges(), 2);

    G.contractEdge(0, 0);

    ASSERT_EQ(G.getUnweightedNodeDegree(0), 0);
    ASSERT_EQ(G.getUnweightedNodeDegree(1), 0);
    ASSERT_EQ(G.number_of_nodes(), 1);
    ASSERT_EQ(G.number_of_edges(), 0);
}

TEST(Mutable_Graph_Test, LargerGraph) {
    for (size_t size : { 5, 10, 50, 100 }) {
        mutableGraphPtr G = std::make_shared<mutable_graph>();
        G->start_construction(size);
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
                if (i < j) {
                    G->new_edge(i, j, i + j);
                }
            }
        }

        EdgeWeight zero_weight = (size * (size - 1)) / 2;

        ASSERT_EQ(G->getUnweightedNodeDegree(0), size - 1);
        ASSERT_EQ(G->getWeightedNodeDegree(0), zero_weight);
        ASSERT_EQ(G->number_of_nodes(), size);
        ASSERT_EQ(G->number_of_edges(), size * (size - 1));

        for (size_t i = 0; i < size - 1; ++i) {
            NodeID tgt = G->getCurrentPosition(i + 1);
            for (EdgeID e : G->edges_of(0)) {
                if (tgt == G->getEdgeTarget(0, e)) {
                    G->contractEdge(0, e);
                    break;
                }
            }

            EdgeWeight rem = size - 2 - i;
            EdgeWeight added_weight = (rem * (i + 1));
            added_weight += ((size) * (size - 1)) / 2;
            added_weight -= ((i + 1) * (i + 2)) / 2;
            added_weight -= ((i + 1) * (i + 1) + ((i + 1) * i) / 2);
            zero_weight += added_weight;

            ASSERT_EQ(G->number_of_nodes(), rem + 1);
            ASSERT_EQ(G->number_of_edges(), (rem + 1) * (rem));
            ASSERT_EQ(G->getUnweightedNodeDegree(0), rem);
            ASSERT_EQ(G->getWeightedNodeDegree(0), zero_weight);
        }
    }
}

TEST(Mutable_Graph_Test, DeleteVertices) {
    for (size_t size : { 5, 10, 50, 100 }) {
        mutableGraphPtr G = std::make_shared<mutable_graph>();
        G->start_construction(size);
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
                if (i < j) {
                    G->new_edge(i, j, i + j);
                }
            }
        }

        for (size_t i = 0; i < size; ++i) {
            EdgeWeight deg = (size - 2) * i + ((size * (size - 1)) / 2);
            ASSERT_EQ(G->getWeightedNodeDegree(i), deg);
        }

        size_t graphsize = size;
        size_t deleted = 0;
        size_t deletedWeight = 0;
        std::random_device rd;
        std::mt19937 eng(rd());
        while (graphsize > 1) {
            std::uniform_int_distribution<> distribution(0, graphsize - 1);
            NodeID del = distribution(eng);
            deletedWeight += G->containedVertices(del)[0];
            G->deleteVertex(del);
            deleted++;
            graphsize--;

            for (size_t i = 0; i < graphsize; ++i) {
                size_t orig = G->containedVertices(i)[0];
                EdgeWeight deg = (size - 2) * orig + ((size * (size - 1)) / 2);
                ASSERT_EQ(G->getWeightedNodeDegree(i),
                          deg - deletedWeight - (deleted * orig));
            }

            for (NodeID n : G->nodes()) {
                for (EdgeID e : G->edges_of(n)) {
                    NodeID tgt = G->getEdgeTarget(n, e);
                    EdgeID rev = G->getReverseEdge(n, e);
                    EdgeWeight wgt = G->getEdgeWeight(n, e);

                    ASSERT_LE(tgt, G->number_of_nodes() - 1);

                    ASSERT_EQ(G->getEdgeTarget(tgt, rev), n);
                    ASSERT_EQ(G->getReverseEdge(tgt, rev), e);
                    ASSERT_EQ(G->getEdgeWeight(tgt, rev), wgt);
                }
            }
        }
    }
}

TEST(Mutable_Graph_Test, SerializeEmpty) {
    mutable_graph mG;
    mG.start_construction(0);
    mG.finish_construction();

    auto s = mG.serialize();
    std::vector<uint64_t> eq = { 0, 0, 0, 0, UNDEFINED_EDGE };

    ASSERT_EQ(s, eq);
}

TEST(Mutable_Graph_Test, SerializeCircle) {
    mutable_graph G = make_circle();

    auto s = G.serialize();
    std::vector<uint64_t> eq = { 3, 6, 0, 3,     // values
                                 6, 1, 1, 2, 1,  // n0
                                 7, 2, 1,        // n1
                                 8,              // n2
                                 UNDEFINED_EDGE,
                                 0, 0, 0,        // partition
                                 0, 1, 2,        // position
    };

    ASSERT_EQ(s, eq);
}

TEST(Mutable_Graph_Test, DeserializedCircleEqual) {
    mutable_graph G = make_circle();
    auto s = G.serialize();
    auto G2 = mutable_graph::deserialize(s);

    ASSERT_EQ(G.getOriginalNodes(), G2->getOriginalNodes());
    ASSERT_EQ(G.n(), G2->n());
    ASSERT_EQ(G.m(), G2->m());
    for (NodeID n : G.nodes()) {
        ASSERT_EQ(G.getWeightedNodeDegree(n), G2->getWeightedNodeDegree(n));
        ASSERT_EQ(G.getPartitionIndex(n), G2->getPartitionIndex(n));
        for (EdgeID e : G.edges_of(n)) {
            ASSERT_EQ(G.getEdgeTarget(n, e), G2->getEdgeTarget(n, e));
            ASSERT_EQ(G.getEdgeWeight(n, e), G2->getEdgeWeight(n, e));
        }
    }

    for (NodeID n = 0; n < G.getOriginalNodes(); ++n) {
        ASSERT_EQ(G.getCurrentPosition(n), G2->getCurrentPosition(n));
    }
}

TEST(Mutable_Graph_Test, DeserializedContractedEqual) {
    mutable_graph G = make_circle();
    G.contractEdge(0, 0);
    auto s = G.serialize();
    auto G2 = mutable_graph::deserialize(s);

    ASSERT_EQ(G.getOriginalNodes(), G2->getOriginalNodes());
    ASSERT_EQ(G.n(), G2->n());
    ASSERT_EQ(G.m(), G2->m());
    for (NodeID n : G.nodes()) {
        ASSERT_EQ(G.getWeightedNodeDegree(n), G2->getWeightedNodeDegree(n));
        ASSERT_EQ(G.getPartitionIndex(n), G2->getPartitionIndex(n));
        for (EdgeID e : G.edges_of(n)) {
            ASSERT_EQ(G.getEdgeTarget(n, e), G2->getEdgeTarget(n, e));
            ASSERT_EQ(G.getEdgeWeight(n, e), G2->getEdgeWeight(n, e));
        }
    }

    for (NodeID n = 0; n < G.getOriginalNodes(); ++n) {
        ASSERT_EQ(G.getCurrentPosition(n), G2->getCurrentPosition(n));
    }
}
