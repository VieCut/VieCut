/******************************************************************************
 * padberg_rinaldi.h
 *
 * Source of VieCut
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
