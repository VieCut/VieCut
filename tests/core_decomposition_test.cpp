#include "data_structure/graph_access.h"
#include "tlx/logger.hpp"
#include <algorithms/misc/core_decomposition.h>
#include <cstdio>
#include <gtest/gtest.h>
#include <io/graph_io.h>

std::shared_ptr<graph_access> make_G2() {
    auto G = std::make_shared<graph_access>();

    G->start_construction(5, 4);

    G->new_node();
    G->new_edge(0, 1);
    G->new_edge(0, 2);
    G->new_node();
    G->new_edge(1, 0);
    G->new_edge(1, 2);
    G->new_node();
    G->new_edge(2, 0);
    G->new_edge(2, 1);
    G->new_edge(2, 3);
    G->new_node();
    G->new_edge(3, 2);
    G->new_node();
    G->finish_construction();
    return G;
}

TEST(CoreDecompositionTest, CreateStruct) {
    k_cores kc = k_cores(100);

    ASSERT_EQ(kc.buckets.size(), 1);
    ASSERT_EQ(kc.degrees.size(), 100);
    ASSERT_EQ(kc.position.size(), 100);
    ASSERT_EQ(kc.vertices.size(), 100);
}

TEST(CoreDecompositionTest, AlgBZGraph1) {

    auto G = graph_io::readGraphWeighted(std::string(VIECUT_PATH) + "/graphs/small.metis");
    k_cores kCores = core_decomposition::batagelj_zaversnik(G);

    std::vector<NodeID> target_buckets = { 0, 0, 0, 0, 8 };
    std::vector<NodeID> target_degrees = { 3, 3, 3, 3, 3, 3, 3, 3 };
    std::vector<NodeID> target_position = { 0, 1, 2, 3, 4, 5, 6, 7 };
    std::vector<NodeID> target_vertices = { 0, 1, 2, 3, 4, 5, 6, 7 };

    ASSERT_EQ(kCores.buckets, target_buckets);
    ASSERT_EQ(kCores.degrees, target_degrees);
    // as all vertices have equal degree in core decomposition, these are non-deterministic.  therefore we only compare sizes
    ASSERT_EQ(kCores.position.size(), target_position.size());
    ASSERT_EQ(kCores.vertices.size(), target_vertices.size());
}

TEST(CoreDecomposisitonTest, AlgBZGraph2) {

    auto G = make_G2();

    k_cores kCores = core_decomposition::batagelj_zaversnik(G);

    std::vector<NodeID> target_buckets = { 0, 1, 2, 5 };
    std::vector<NodeID> target_degrees = { 2, 2, 2, 1, 0 };

    ASSERT_EQ(kCores.buckets, target_buckets);
    ASSERT_EQ(kCores.degrees, target_degrees);
    // as all vertices have equal degree in core decomposition, these are non-deterministic.
    // therefore we only compare sizes and the ones that are definitive
    ASSERT_EQ(kCores.position.size(), 5);
    ASSERT_EQ(kCores.position[3], 1);
    ASSERT_EQ(kCores.position[4], 0);
    ASSERT_EQ(kCores.vertices.size(), 5);
    ASSERT_EQ(kCores.vertices[1], 3);
    ASSERT_EQ(kCores.vertices[0], 4);
}

TEST(CoreDecompositionTest, CreateCoreGraph) {

    auto G = make_G2();

    k_cores kCores = core_decomposition::batagelj_zaversnik(G);

    for (size_t k : { 0, 1, 2 }) {
        auto new_graph = core_decomposition::createCoreGraph(kCores, k, G);
        ASSERT_EQ(new_graph->number_of_nodes(), std::min<EdgeID>(4, 5 - k));
        ASSERT_EQ(new_graph->number_of_edges(), std::min<EdgeID>(8, 10 - (2 * k)));
    }
}
