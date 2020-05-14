/******************************************************************************
 * mincut_algo_test.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <string>
#include <type_traits>

#ifdef PARALLEL
#include "algorithms/global_mincut/viecut.h"
#include "parallel/algorithm/exact_parallel_minimum_cut.h"
#include "parallel/algorithm/parallel_cactus.h"
#else
#include "algorithms/global_mincut/cactus/cactus_mincut.h"
#include "algorithms/global_mincut/ks_minimum_cut.h"
#include "algorithms/global_mincut/matula_approx.h"
#include "algorithms/global_mincut/noi_minimum_cut.h"
#include "algorithms/global_mincut/padberg_rinaldi.h"
#include "algorithms/global_mincut/stoer_wagner_minimum_cut.h"
#include "algorithms/global_mincut/viecut.h"
#endif
#include "common/definitions.h"
#include "gtest/gtest_pred_impl.h"
#include "io/graph_io.h"

template <class GraphPtr>
class exact_parallel_minimum_cut;
template <class GraphPtr>
class parallel_cactus;
template <class GraphPtr>
class viecut;
class ks_minimum_cut;
template <class GraphPtr>
class cactus_mincut;
template <class GraphPtr>
class matula_approx;
template <class GraphPtr>
class noi_minimum_cut;
template <class GraphPtr>
class padberg_rinaldi;

template <typename T>
class MincutAlgoTest : public testing::Test { };

#ifdef PARALLEL
typedef testing::Types<viecut<graphAccessPtr>,
                       exact_parallel_minimum_cut<graphAccessPtr>,
                       parallel_cactus<graphAccessPtr>,
                       viecut<mutableGraphPtr>,
                       exact_parallel_minimum_cut<mutableGraphPtr>,
                       parallel_cactus<mutableGraphPtr> >
    MCAlgTypes;
#else
typedef testing::Types<viecut<graphAccessPtr>,
                       noi_minimum_cut<graphAccessPtr>,
                       padberg_rinaldi<graphAccessPtr>,
                       matula_approx<graphAccessPtr>,
                       ks_minimum_cut,
                       cactus_mincut<graphAccessPtr>,
                       viecut<mutableGraphPtr>,
                       noi_minimum_cut<mutableGraphPtr>,
                       padberg_rinaldi<mutableGraphPtr>,
                       matula_approx<mutableGraphPtr>,
                       cactus_mincut<mutableGraphPtr> >
    MCAlgTypes;
#endif

TYPED_TEST_CASE(MincutAlgoTest, MCAlgTypes);

TYPED_TEST(MincutAlgoTest, NoGraph) {
    typename TypeParam::GraphPtrType G;
    TypeParam mc;

    EdgeWeight cut = mc.perform_minimum_cut(G);
    ASSERT_EQ(cut, (EdgeWeight) - 1);
}

TYPED_TEST(MincutAlgoTest, UnweightedGraphFromFile) {
    typename TypeParam::GraphPtrType G =
        graph_io::readGraphWeighted<
            typename TypeParam::GraphPtrType::element_type>(
            std::string(VIECUT_PATH) + "/graphs/small.metis");
    TypeParam mc;
    EdgeWeight cut = mc.perform_minimum_cut(G);

#ifdef PARALLEL
    if (std::is_same<TypeParam,
                     exact_parallel_minimum_cut<graphAccessPtr> >::value ||
        std::is_same<TypeParam,
                     exact_parallel_minimum_cut<mutableGraphPtr> >::value) {
#else
    if (std::is_same<TypeParam, noi_minimum_cut<graphAccessPtr> >::value ||
        std::is_same<TypeParam, noi_minimum_cut<mutableGraphPtr> >::value) {
#endif
        ASSERT_EQ(cut, 2);
    } else {
        // inexact, we can only guarantee that minimum cut is
        // between minimum degree and minimum cut
        ASSERT_LE(cut, 3);
        ASSERT_GE(cut, 2);
    }
}

TYPED_TEST(MincutAlgoTest, WeightedGraphFromFile) {
    typename TypeParam::GraphPtrType G =
        graph_io::readGraphWeighted<
            typename TypeParam::GraphPtrType::element_type>(
            std::string(VIECUT_PATH) + "/graphs/small-wgt.metis");
    TypeParam mc;
    EdgeWeight cut = mc.perform_minimum_cut(G);

#ifdef PARALLEL
    if (std::is_same<TypeParam,
                     exact_parallel_minimum_cut<graphAccessPtr> >::value ||
        std::is_same<TypeParam,
                     exact_parallel_minimum_cut<mutableGraphPtr> >::value) {
#else
    if (std::is_same<TypeParam, noi_minimum_cut<graphAccessPtr> >::value ||
        std::is_same<TypeParam, noi_minimum_cut<mutableGraphPtr> >::value) {
#endif
        ASSERT_EQ(cut, 3);
    } else {
        // inexact, we can only guarantee that minimum cut is
        // between minimum degree and minimum cut
        ASSERT_LE(cut, 10);
        ASSERT_GE(cut, 3);
    }
}
