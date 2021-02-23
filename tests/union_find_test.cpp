/******************************************************************************
 * union_find_test.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <omp.h>
#include <stddef.h>

#include <algorithm>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include "gtest/gtest.h"

#ifdef PARALLEL
#include "parallel/data_structure/union_find.h"
#else
#include "data_structure/union_find.h"
#endif

TEST(UnionFindTest, CreateEmpty) {
    union_find uf(0);
    ASSERT_EQ(uf.n(), 0);
}

TEST(UnionFindTest, CreateWithSize) {
    union_find uf(200);
    ASSERT_EQ(uf.n(), 200);
}

TEST(UnionFindTest, FindInitial) {
    union_find uf(200);

    for (size_t i = 0; i < 200; ++i) {
        ASSERT_EQ(uf.Find(i), i);
    }
}

TEST(UnionFindTest, UnionTwoToOne) {
    union_find uf(1000);

#ifdef PARALLEL
    omp_set_num_threads(4);
#pragma omp parallel for
#endif
    for (size_t i = 0; i < 500; ++i) {
        uf.Union(2 * i, (2 * i) + 1);
    }

    ASSERT_EQ(uf.n(), 500);

    for (size_t i = 0; i < 500; ++i) {
        ASSERT_EQ(uf.Find(2 * i), uf.Find((2 * i) + 1));
        ASSERT_TRUE((uf.Find(2 * i) >= 2 * i)
                    && (uf.Find((2 * i) + 1) <= (2 * i) + 1));
    }
}

TEST(UnionFindTest, UnionBlocks) {
    // In this test we want to Union() blocks of [i,i+50)
    // we randomly swap the operations to ensure
    // that this works correctly for all orderings.
    size_t num_blocks = 50;
    size_t blocksize_outer = 5;
    size_t blocksize_inner = 10;
    size_t blocksize = blocksize_inner * blocksize_outer;
    union_find uf(num_blocks* blocksize);
    std::random_device rd;
    std::mt19937 eng(rd());
    std::uniform_int_distribution<> distribution(0, blocksize_inner - 1);

    std::vector<std::pair<size_t, size_t> > union_ops;

    for (size_t i = 0; i < num_blocks; ++i) {
        size_t block_start = i * blocksize;
        for (size_t j = 0; j < blocksize_outer; ++j) {
            size_t inner_start = block_start + j * blocksize_inner;
            for (size_t k = 0; k < blocksize_inner - 1; ++k) {
                // union all in inner block
                union_ops.emplace_back(inner_start + k, inner_start + k + 1);
            }
            // union random element of inner block to
            // random element of next block
            if (j + 1 < blocksize_outer)
                union_ops.emplace_back(inner_start + distribution(eng),
                                       inner_start + blocksize_inner
                                       + distribution(eng));
        }
    }

    std::random_shuffle(union_ops.begin(), union_ops.end());
#ifdef PARALLEL
    omp_set_num_threads(4);
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (size_t i = 0; i < union_ops.size(); ++i) {
        uf.Union(union_ops[i].first, union_ops[i].second);
    }

    ASSERT_EQ(uf.n(), num_blocks);

    for (size_t i = 0; i < num_blocks; ++i) {
        size_t block_start = i * blocksize;
        size_t blockval = uf.Find(block_start);
        // value is one of elements in block...
        ASSERT_TRUE((blockval >= block_start)
                    && (blockval < block_start + blocksize));
        // ..and all elements in block have same value
        for (size_t j = 0; j < blocksize; ++j) {
            ASSERT_EQ(blockval, uf.Find(block_start + j));
        }
    }
}
