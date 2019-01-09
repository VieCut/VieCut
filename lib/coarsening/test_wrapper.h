/******************************************************************************
 * test_wrapper.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
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

    static std::pair<std::shared_ptr<graph_access>, EdgeWeight> run_pr_12(std::vector<std::shared_ptr<graph_access> > graphs, EdgeWeight mincut, bool save_cut) {
        timer t;
        auto uf = tests::prTests12(graphs.back(), mincut);
        std::shared_ptr<graph_access> G_out = contraction::contractFromUnionFind(graphs.back(), uf, save_cut);
        mincut = minimum_cut_helpers::updateCutValueAfterContraction(graphs, mincut, save_cut);
        LOGC(timing) << "Padberg-Rinaldi Tests 1-2 (to "
                     << G_out->number_of_nodes() << "nodes): " << t.elapsedToZero();

        return std::make_pair(G_out, mincut);
    }

    static std::shared_ptr<graph_access> run_pr_34(std::vector<std::shared_ptr<graph_access> > graphs, EdgeWeight mincut, bool save_cut) {
        timer t;
        auto uf = tests::prTests34(graphs.back(), mincut);
        std::shared_ptr<graph_access> G_out = contraction::contractFromUnionFind(graphs.back(), uf, save_cut);
        mincut = minimum_cut_helpers::updateCutValueAfterContraction(graphs, mincut, save_cut);
        LOGC(timing) << "Padberg-Rinaldi Tests 3-4 (to "
                     << G_out->number_of_nodes() << "nodes): " << t.elapsedToZero();

        return G_out;
    }
};
