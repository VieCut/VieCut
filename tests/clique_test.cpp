/******************************************************************************
 * clique_test.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include "algorithms/misc/maximal_clique.h"
#include "algorithms/multicut/multicut_problem.h"
#include "common/configuration.h"
#include "data_structure/mutable_graph.h"
#include "gtest/gtest.h"
#include "io/graph_io.h"
#include "tools/random_functions.h"

TEST(CliqueCutTest, SingleClique) {
    std::shared_ptr<mutable_graph> G = std::make_shared<mutable_graph>();

    G->start_construction(4);

    for (NodeID i = 0; i < 4; ++i) {
        for (NodeID j = 0; j < 4; ++j) {
            G->new_edge(i, j);
        }
    }

    G->finish_construction();

    maximal_clique mc;

    auto r = mc.findCliques(G);

    ASSERT_EQ(r.size(), 1);
    ASSERT_EQ(r[0].size(), 4);
}

TEST(CliqueTest, MultipleCliques) {
    std::shared_ptr<mutable_graph> G = std::make_shared<mutable_graph>();

    G->start_construction(16);

    for (NodeID k = 0; k < 4; ++k) {
        for (NodeID i = 0; i < 4; ++i) {
            for (NodeID j = 0; j < 4; ++j) {
                G->new_edge((4 * k) + i, (4 * k) + j);
            }
        }
    }

    G->finish_construction();

    maximal_clique mc;

    auto r = mc.findCliques(G);

    ASSERT_EQ(r.size(), 4);
    for (size_t i = 0; i < 4; ++i) {
        ASSERT_EQ(r[i].size(), 4);
    }
}

TEST(CliqueTest, CliqueWithTriangleEars) {
    std::shared_ptr<mutable_graph> G = std::make_shared<mutable_graph>();

    G->start_construction(6);

    for (NodeID i = 0; i < 4; ++i) {
        for (NodeID j = 0; j < 4; ++j) {
            G->new_edge(i, j);
        }
    }

    G->new_edge(1, 4);
    G->new_edge(2, 4);
    G->new_edge(2, 5);
    G->new_edge(3, 5);

    G->finish_construction();

    LOG1 << G;

    maximal_clique mc;

    auto r = mc.findCliques(G);

    ASSERT_EQ(r.size(), 3);
    ASSERT_EQ(r[0].size(), 4);
    ASSERT_EQ(r[1].size(), 3);
    ASSERT_EQ(r[2].size(), 3);
}
