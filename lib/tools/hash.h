//
// Created by noea45cs on 28.08.18.
//

#pragma once
#define XXH_PRIVATE_API
#include "misc/submodules/xxhash/xxhash.h"

template <typename T>
struct xxhash {
    using h = size_t;
    static const size_t significant_digits = 64;

    inline h operator () (const T& t) const {
        return XXH64(&t, sizeof(t), 0);
    }
};
