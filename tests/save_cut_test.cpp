/******************************************************************************
 * save_cut_test.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include "data_structure/graph_access.h"
#include "tlx/logger.hpp"
#include <cstdio>
#include <gtest/gtest.h>
#include <io/graph_io.h>
#ifdef PARALLEL
#include "algorithms/global_mincut/viecut.h"
#include "parallel/algorithm/exact_parallel_minimum_cut.h"
#else
#include "algorithms/global_mincut/ks_minimum_cut.h"
#include "algorithms/global_mincut/matula_approx.h"
#include "algorithms/global_mincut/noi_minimum_cut.h"
#include "algorithms/global_mincut/padberg_rinaldi.h"
#include "algorithms/global_mincut/stoer_wagner_minimum_cut.h"
#include "algorithms/global_mincut/viecut.h"
#endif

template <typename T>
class SaveCutTest : public testing::Test { };

#ifdef PARALLEL
typedef testing::Types<viecut, exact_parallel_minimum_cut> MCAlgTypes;
#else
typedef testing::Types<viecut, noi_minimum_cut, padberg_rinaldi, matula_approx, ks_minimum_cut> MCAlgTypes;// maxNodeHeap, vecMaxNodeHeap> PQTypes;
#endif

TYPED_TEST_CASE(SaveCutTest, MCAlgTypes);

TYPED_TEST(SaveCutTest, UnweightedGraph) {
    std::shared_ptr<graph_access> G = graph_io::readGraphWeighted(std::string(VIECUT_PATH) + "/graphs/small.metis");
    TypeParam mc;

    EdgeWeight cut = mc.perform_minimum_cut(G, true);

#ifdef PARALLEL
    if (std::is_same<TypeParam, exact_parallel_minimum_cut>::value) {
#else
    if (std::is_same<TypeParam, noi_minimum_cut>::value) {
#endif
        ASSERT_EQ(cut, (EdgeWeight)2);
        std::vector<EdgeID> e;
        if (G->getNodeInCut(0))
            e = { 1, 1, 1, 1, 0, 0, 0, 0 };
        else
            e = { 0, 0, 0, 0, 1, 1, 1, 1 };

        for (size_t i = 0; i < e.size(); ++i) {
            ASSERT_EQ(G->getNodeInCut(i), e[i]);
        }
    }
    else {
        // we cannot guarantee that cut is correct in exact algorithms, however, not all vertices can be on same side of cut

        ASSERT_LE(cut, (EdgeWeight)3);
        ASSERT_GE(cut, (EdgeWeight)2);
        size_t cut_sum = 0;
        for (NodeID n : G->nodes()) {
            cut_sum += G->getNodeInCut(n);
        }

        ASSERT_LE(cut_sum, (EdgeWeight)7);
        ASSERT_GE(cut_sum, (EdgeWeight)1);
    }
}

TYPED_TEST(SaveCutTest, WeightedGraph) {
    std::shared_ptr<graph_access> G = graph_io::readGraphWeighted(std::string(VIECUT_PATH) + "/graphs/small-wgt.metis");
    TypeParam mc;

    EdgeWeight cut = mc.perform_minimum_cut(G, true);

#ifdef PARALLEL
    if (std::is_same<TypeParam, exact_parallel_minimum_cut>::value) {
#else
    if (std::is_same<TypeParam, noi_minimum_cut>::value) {
#endif
        ASSERT_EQ(cut, (EdgeWeight)3);
        std::vector<EdgeID> e;
        if (G->getNodeInCut(0))
            e = { 1, 1, 1, 1, 0, 0, 0, 0 };
        else
            e = { 0, 0, 0, 0, 1, 1, 1, 1 };

        for (size_t i = 0; i < e.size(); ++i) {
            ASSERT_EQ(G->getNodeInCut(i), e[i]);
        }
    }
    else {
        // we cannot guarantee that cut is correct in exact algorithms, however, not all vertices can be on same side of cut

        ASSERT_LE(cut, (EdgeWeight)10);
        ASSERT_GE(cut, (EdgeWeight)3);
        size_t cut_sum = 0;
        for (NodeID n : G->nodes()) {
            cut_sum += G->getNodeInCut(n);
        }

        ASSERT_LE(cut_sum, (EdgeWeight)7);
        ASSERT_GE(cut_sum, (EdgeWeight)1);
    }
}
