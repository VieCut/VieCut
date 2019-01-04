/******************************************************************************
 * padberg_rinaldi.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include "minimum_cut.h"
#include "minimum_cut_helpers.h"

#ifdef PARALLEL
#include "parallel/coarsening/contract_graph.h"
#include "parallel/coarsening/contraction_tests.h"
#else // PARALLEL
#include "coarsening/contract_graph.h"
#include "coarsening/contraction_tests.h"
#endif // PARALLEL

#include "definitions.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <unordered_map>

#include <cstdint>
#include <cstdlib>
class padberg_rinaldi : public minimum_cut
{
public:
    padberg_rinaldi() { }
    virtual ~padberg_rinaldi() { }
    static constexpr bool debug = false;

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G) {

        if (!minimum_cut_helpers::graphValid(G))
            return -1;
        EdgeWeight cut = G->getMinDegree();
        std::vector<std::shared_ptr<graph_access> > graphs;
        graphs.push_back(G);
        NodeID last_nodes = G->number_of_nodes() + 1;
        timer t;
        minimum_cut_helpers::setInitialCutValues(graphs);

        while (graphs.back()->number_of_nodes() > 2
               && graphs.back()->number_of_nodes() < last_nodes) {
            last_nodes = graphs.back()->number_of_nodes();
            union_find uf_34 = tests::prTests34(graphs.back(), cut);
            graphs.push_back(contraction::contractFromUnionFind(graphs.back(), uf_34));
            cut = minimum_cut_helpers::updateCutValueAfterContraction(graphs, cut);
            union_find uf_12 = tests::prTests12(graphs.back(), cut);
            graphs.push_back(contraction::contractFromUnionFind(graphs.back(), uf_12));
            cut = minimum_cut_helpers::updateCutValueAfterContraction(graphs, cut);
        }

        minimum_cut_helpers::retrieveMinimumCut(graphs);

        LOG << "nodesleft=" << graphs.back()->number_of_nodes();

        return cut;
    }
};
