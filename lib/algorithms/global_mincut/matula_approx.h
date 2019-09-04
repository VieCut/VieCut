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

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include "algorithms/global_mincut/minimum_cut.h"
#include "algorithms/global_mincut/noi_minimum_cut.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/priority_queues/maxNodeHeap.h"
#include "tools/timer.h"

#ifdef PARALLEL
#include "parallel/coarsening/contract_graph.h"
#else
#include "coarsening/contract_graph.h"
#endif

class matula_approx : public minimum_cut {
 public:
    matula_approx() { }
    virtual ~matula_approx() { }
    static constexpr bool debug = false;

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G) {
        if (!minimum_cut_helpers::graphValid(G))
            return -1;

        std::vector<std::shared_ptr<graph_access> > graphs;
        EdgeWeight mincut = G->getMinDegree();
        graphs.push_back(G);
        minimum_cut_helpers::setInitialCutValues(graphs);

        while (graphs.back()->number_of_nodes() > 2 && mincut > 0) {
            std::vector<std::pair<NodeID, NodeID> > contractable;
            timer time;
            noi_minimum_cut noi;
            auto uf = noi.modified_capforest(graphs.back(),
                                             std::max(mincut / 2, 1UL));
            graphs.emplace_back(contraction::fromUnionFind(graphs.back(), &uf));
            mincut = minimum_cut_helpers::updateCut(graphs, mincut);
        }

        if (configuration::getConfig()->save_cut)
            minimum_cut_helpers::retrieveMinimumCut(graphs);

        return mincut;
    }
};
