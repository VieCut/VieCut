/******************************************************************************
 * viecut.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2017-2018 Alexander Noe <alexander.noe@univie.ac.at>
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

#include "algorithms/flow/excess_scaling.h"
#include "algorithms/misc/strongly_connected_components.h"
#include "data_structure/flow_graph.h"
#include "data_structure/graph_access.h"
#include "definitions.h"
#include "minimum_cut.h"
#include "minimum_cut_helpers.h"
#include "noi_minimum_cut.h"
#include "tlx/logger.hpp"
#include "tools/graph_extractor.h"
#include "tools/timer.h"

#ifdef PARALLEL
#include "parallel/coarsening/contract_graph.h"
#include "parallel/coarsening/contraction_tests.h"
#include "parallel/coarsening/label_propagation.h"
#else // PARALLEL
#include "coarsening/contract_graph.h"
#include "coarsening/contraction_tests.h"
#include "coarsening/label_propagation.h"
#endif // PARALLEL

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <unordered_map>

#include <cstdint>
#include <cstdlib>

class viecut : public minimum_cut
{
public:
    static constexpr bool debug = false;
    static constexpr bool timing = true;
    viecut() { }

    virtual ~viecut() { }

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G) {
        return perform_minimum_cut(G, false);
    }

    EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G,
                                   bool indirect) {

        if (!minimum_cut_helpers::graphValid(G))
            return -1;
        EdgeWeight cut = G->getMinDegree();

        std::vector<std::shared_ptr<graph_access> > graphs;
        graphs.push_back(G);

        minimum_cut_helpers::setInitialCutValues(graphs);

        while (graphs.back()->number_of_nodes() > 10000 &&
               (graphs.size() == 1 ||
                (graphs.back()->number_of_nodes() <
                 graphs[graphs.size() - 2]->number_of_nodes()))) {

            timer t;
            G = graphs.back();

            label_propagation lp;
            std::vector<NodeID> cluster_mapping = lp.propagate_labels(G);
            std::vector<std::vector<NodeID> > reverse_mapping = minimum_cut_helpers::remap_cluster(G, cluster_mapping);
            LOGC(timing) << "LP (total): " << t.elapsedToZero();

            contraction::findTrivialCuts(G, cluster_mapping, reverse_mapping, cut);
            LOGC(timing) << "Trivial Cut Local Search: " << t.elapsedToZero();

            graphs.push_back(contraction::contractGraph(G, cluster_mapping, reverse_mapping.size(), reverse_mapping));
            cut = minimum_cut_helpers::updateCutValueAfterContraction(graphs, cut);
            LOGC(timing) << "Graph Contraction (to " << graphs.back()->number_of_nodes() << " nodes): " << t.elapsedToZero();

            union_find uf = tests::prTests12(graphs.back(), cut);
            graphs.push_back(contraction::contractFromUnionFind(graphs.back(), uf));
            cut = minimum_cut_helpers::updateCutValueAfterContraction(graphs, cut);
            union_find uf2 = tests::prTests34(graphs.back(), cut);
            graphs.push_back(contraction::contractFromUnionFind(graphs.back(), uf2));
            cut = minimum_cut_helpers::updateCutValueAfterContraction(graphs, cut);
            LOGC(timing) << "Padberg-Rinaldi Tests (to " << graphs.back()->number_of_nodes() << " nodes): " << t.elapsedToZero();
        }

        if (graphs.back()->number_of_nodes() > 1) {
            timer t;
            noi_minimum_cut noi;
            cut = std::min(cut, noi.perform_minimum_cut(graphs.back(), true));

            LOGC(timing) << "Exact Algorithm:" << t.elapsedToZero() << " deg: " << graphs.back()->getMinDegree();
        }

        if (!indirect)
            minimum_cut_helpers::retrieveMinimumCut(graphs);

        return cut;
    }
};
