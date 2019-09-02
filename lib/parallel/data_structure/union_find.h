/******************************************************************************
 * union_find.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "tlx/logger.hpp"

// A simple Union-Find datastructure implementation.
// This is sometimes also caled "disjoint sets datastructure.
class union_find {
 public:
    explicit union_find(unsigned n) : m_parent(n), m_rank(n) {
        for (unsigned i = 0; i < m_parent.size(); i++) {
            m_parent[i] = i;
            m_rank[i] = 0;
        }
    }

    union_find(const union_find& uf) = default;

    inline bool Union(unsigned lhs, unsigned rhs) {
        unsigned r_lhs, r_rhs;

        while (true) {
            lhs = Find(lhs);
            rhs = Find(rhs);

            if (lhs == rhs) {
                return false;
            }

            r_lhs = m_rank[lhs];
            r_rhs = m_rank[rhs];

            if (r_lhs > r_rhs || (r_lhs == r_rhs && lhs > rhs)) {
                std::swap(lhs, rhs);
                std::swap(r_lhs, r_rhs);
            }

            if (UpdateRoot(lhs, r_lhs, rhs, r_rhs))
                break;
        }

        if (r_lhs == r_rhs) {
            UpdateRoot(rhs, r_rhs, rhs, r_rhs + 1);
        }

        return true;
    }

    inline bool SameSet(unsigned x, unsigned y) {
        while (true) {
            x = Find(x);
            y = Find(y);
            if (x == y)
                return true;
            if (m_parent[x] == x) {
                return false;
            }
        }
    }

    inline unsigned Find(unsigned element) const {
        while (m_parent[element] != element) {
            unsigned next = m_parent[element];
            element = m_parent[next];
        }
        return element;
    }

    inline unsigned Find(unsigned element) {
        while (m_parent[element] != element) {
            unsigned next = m_parent[element];
            // CAS path halving
            __sync_bool_compare_and_swap(&m_parent[element], next,
                                         m_parent[next]);
            element = m_parent[next];
        }
        return element;
    }

    inline unsigned n() const {
        std::vector<bool> found(m_parent.size(), false);
        unsigned n = 0;
        for (auto& e : m_parent) {
            if (!found[Find(e)]) {
                ++n;
                found[Find(e)] = true;
            }
        }
        return n;
    }

 private:
    inline bool UpdateRoot(unsigned x, unsigned xr, unsigned y, unsigned yr) {
        unsigned old = x;
        if (m_parent[old] != x || m_rank[old] != xr) {
            return false;
        }

        if (__sync_bool_compare_and_swap(&m_parent[x], old, y)) {
            if (__sync_bool_compare_and_swap(&m_rank[x], m_rank[old], yr)) {
                return true;
            }

            __sync_bool_compare_and_swap(&m_parent[x], y, old);
        }
        return false;
    }

    std::vector<unsigned> m_parent;
    std::vector<unsigned> m_rank;
};
