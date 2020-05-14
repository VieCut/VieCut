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

#include <stddef.h>

#include <gtest/gtest-message.h>
#include <gtest/gtest-test-part.h>
#include <memory>
#include <vector>

#include "algorithms/misc/maximal_clique.h"
#include "common/definitions.h"
#include "data_structure/mutable_graph.h"
#include "gtest/gtest_pred_impl.h"
#include "tlx/logger.hpp"

TEST(CliqueCutTest, SingleClique) {
    mutableGraphPtr G = std::make_shared<mutable_graph>();

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
    mutableGraphPtr G = std::make_shared<mutable_graph>();

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
    mutableGraphPtr G = std::make_shared<mutable_graph>();

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

    ASSERT_LE(r.size(), 3);

    for (const auto& c : r) {
        ASSERT_LE(c.size(), 4);
        ASSERT_GE(c.size(), 3);
    }
}
