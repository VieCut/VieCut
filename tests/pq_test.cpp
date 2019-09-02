/******************************************************************************
 * pq_test.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <random>
#include <type_traits>

#include "data_structure/priority_queues/bucket_pq.h"
#include "data_structure/priority_queues/fifo_node_bucket_pq.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "data_structure/priority_queues/node_bucket_pq.h"
#include "data_structure/priority_queues/priority_queue_interface.h"
#include "data_structure/priority_queues/vecMaxNodeHeap.h"
#include "gtest/gtest.h"
#include "tlx/logger.hpp"

template <typename T>
class PQTest : public testing::Test { };

typedef testing::Types<vecMaxNodeHeap, maxNodeHeap, node_bucket_pq,
                       fifo_node_bucket_pq, bucket_pq> PQTypes;
TYPED_TEST_CASE(PQTest, PQTypes);

TYPED_TEST(PQTest, EmptyAtStart) {
    auto priority_queue = new TypeParam(10, 10);

    ASSERT_EQ(priority_queue->size(), 0);
    ASSERT_TRUE(priority_queue->empty());
}

TYPED_TEST(PQTest, AddElements) {
    auto priority_queue = new TypeParam(10, 10);

    priority_queue->insert(5, 3);

    ASSERT_EQ(priority_queue->size(), 1);
    ASSERT_EQ(priority_queue->maxValue(), 3);
    ASSERT_EQ(priority_queue->maxElement(), 5);

    priority_queue->insert(9, 1);

    ASSERT_EQ(priority_queue->size(), 2);
    ASSERT_EQ(priority_queue->maxValue(), 3);
    ASSERT_EQ(priority_queue->maxElement(), 5);

    priority_queue->insert(1, 4);

    ASSERT_EQ(priority_queue->size(), 3);
    ASSERT_EQ(priority_queue->maxValue(), 4);
    ASSERT_EQ(priority_queue->maxElement(), 1);
}

TYPED_TEST(PQTest, DeleteMax) {
    auto priority_queue = new TypeParam(10, 10);

    priority_queue->insert(5, 3);
    priority_queue->insert(1, 4);
    priority_queue->insert(9, 1);

    priority_queue->deleteMax();
    ASSERT_EQ(priority_queue->size(), 2);
    ASSERT_EQ(priority_queue->maxValue(), 3);
    ASSERT_EQ(priority_queue->maxElement(), 5);

    priority_queue->deleteMax();
    ASSERT_EQ(priority_queue->size(), 1);
    ASSERT_EQ(priority_queue->maxValue(), 1);
    ASSERT_EQ(priority_queue->maxElement(), 9);

    priority_queue->deleteMax();
    ASSERT_TRUE(priority_queue->empty());
}

TYPED_TEST(PQTest, ChangeKeys) {
    auto priority_queue = new TypeParam(10, 10);

    priority_queue->insert(5, 3);
    priority_queue->insert(1, 4);
    priority_queue->insert(9, 1);

    // increase key to new maximum
    priority_queue->changeKey(9, 5);
    ASSERT_EQ(priority_queue->size(), 3);
    ASSERT_EQ(priority_queue->maxValue(), 5);
    ASSERT_EQ(priority_queue->maxElement(), 9);

    // lower keys so (3,5) is maximum
    priority_queue->changeKey(9, 2);
    priority_queue->changeKey(1, 2);
    ASSERT_EQ(priority_queue->size(), 3);
    ASSERT_EQ(priority_queue->maxValue(), 3);
    ASSERT_EQ(priority_queue->maxElement(), 5);
}

TYPED_TEST(PQTest, HeapSortRandom) {
    size_t num_el = 1000;
    auto priority_queue = new TypeParam(num_el, num_el);

    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_int_distribution<> distribution(0, 999);

    for (size_t i = 0; i < num_el; ++i) {
        priority_queue->insert(i, distribution(eng));
    }

    std::vector<size_t> vec;
    for (size_t i = 0; i < num_el; ++i) {
        EdgeWeight max_prio = priority_queue->maxValue();
        priority_queue->deleteMax();
        vec.push_back(max_prio);
    }

    for (size_t i = 0; i + 1 < num_el; ++i) {
        ASSERT_GE(vec[i], vec[i + 1]);
    }
}
