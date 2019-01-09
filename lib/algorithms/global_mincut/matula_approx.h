/******************************************************************************
 * matula_approx.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "minimum_cut.h"
#include "noi_minimum_cut.h"

#include "tools/timer.h"

#ifdef PARALLEL
#include "parallel/coarsening/contract_graph.h"
#else
#include "coarsening/contract_graph.h"
#endif

#include "definitions.h"
#include "minimum_cut.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <memory>
#include <unordered_map>

#include <cstdint>
#include <cstdlib>
class matula_approx : public minimum_cut
{
public:
    matula_approx() { }
    virtual ~matula_approx() { }
    static constexpr bool debug = false;

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G, bool save_cut) {

        if (!minimum_cut_helpers::graphValid(G))
            return -1;

        std::vector<std::shared_ptr<graph_access> > graphs;

        timer t;
        EdgeWeight global_mincut = G->getMinDegree();
        graphs.push_back(G);

        minimum_cut_helpers::setInitialCutValues(graphs, save_cut);

        while (graphs.back()->number_of_nodes() > 2 && global_mincut > 0) {

            std::vector<std::pair<NodeID, NodeID> > contractable;

            union_find uf(graphs.back()->number_of_nodes());
            timer time;
            noi_minimum_cut noi;
            noi.modified_capforest(graphs.back(), global_mincut / 2, uf, graphs, save_cut);
            graphs.emplace_back(contraction::contractFromUnionFind(graphs.back(), uf, save_cut));
            global_mincut = minimum_cut_helpers::updateCutValueAfterContraction(graphs, global_mincut, save_cut);
        }

        if (save_cut)
            minimum_cut_helpers::retrieveMinimumCut(graphs);

        return global_mincut;
    }
};
