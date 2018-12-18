/******************************************************************************
 * test_wrapper.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#pragma once

#include "data_structure/graph_access.h"
#include "definitions.h"
#include "tlx/logger.hpp"

#ifdef PARALLEL
#include "parallel/coarsening/contraction_tests.h"
#include "parallel/data_structure/union_find.h"
#else
#include "coarsening/contraction_tests.h"
#include "data_structure/union_find.h"
#endif

#include <cstdint>
#include <cstdlib>
#include <memory>

class test_wrapper
{

public:
    const static bool timing = false;
    const static bool debug = false;

    static std::pair<std::shared_ptr<graph_access>, EdgeWeight> run_pr_12(std::vector<std::shared_ptr<graph_access> > graphs, EdgeWeight mincut) {
        timer t;
        auto uf = tests::prTests12(graphs.back(), mincut);
        std::shared_ptr<graph_access> G_out = contraction::contractFromUnionFind(graphs.back(), uf);
        mincut = minimum_cut_helpers::updateCutValueAfterContraction(graphs, mincut);
        LOGC(timing) << "Padberg-Rinaldi Tests 1-2 (to "
                     << G_out->number_of_nodes() << "nodes): " << t.elapsedToZero();

        return std::make_pair(G_out, mincut);
    }

    static std::shared_ptr<graph_access> run_pr_34(std::vector<std::shared_ptr<graph_access> > graphs, EdgeWeight mincut) {
        timer t;
        auto uf = tests::prTests34(graphs.back(), mincut);
        std::shared_ptr<graph_access> G_out = contraction::contractFromUnionFind(graphs.back(), uf);
        mincut = minimum_cut_helpers::updateCutValueAfterContraction(graphs, mincut);
        LOGC(timing) << "Padberg-Rinaldi Tests 3-4 (to "
                     << G_out->number_of_nodes() << "nodes): " << t.elapsedToZero();

        return G_out;
    }
};
