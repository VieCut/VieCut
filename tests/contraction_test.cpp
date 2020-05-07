/******************************************************************************
 * contraction_test.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <omp.h>

#ifdef PARALLEL
#include "parallel/coarsening/contract_graph.h"
#else
#include "coarsening/contract_graph.h"
#endif
#include "data_structure/graph_access.h"
#include "gtest/gtest.h"
#include "io/graph_io.h"
#include "tlx/logger.hpp"

TEST(ContractionTest, NoContr) {
#ifdef PARALLEL
    omp_set_num_threads(4);
#endif
    graphAccessPtr G = graph_io::readGraphWeighted(
        std::string(VIECUT_PATH) + "/graphs/small.metis");

    std::vector<NodeID> mapping;
    std::vector<std::vector<NodeID> > reverse_mapping;

    for (NodeID n : G->nodes()) {
        mapping.push_back(n);
        reverse_mapping.emplace_back();
        reverse_mapping.back().emplace_back(n);
    }

    graphAccessPtr cntr = contraction::contractGraph(
        G, mapping, reverse_mapping);

    ASSERT_EQ(G->number_of_nodes(), cntr->number_of_nodes());
    ASSERT_EQ(G->number_of_edges(), cntr->number_of_edges());
}

TEST(ContractionTest, ContrBlock) {
#ifdef PARALLEL
    omp_set_num_threads(4);
#endif
    std::vector<std::string> graphs = { "", "-wgt" };
    for (std::string graph : graphs) {
        graphAccessPtr G = graph_io::readGraphWeighted(
            std::string(VIECUT_PATH) + "/graphs/small" + graph + ".metis");

        std::vector<NodeID> mapping;
        std::vector<std::vector<NodeID> > reverse_mapping;

        reverse_mapping.emplace_back();
        reverse_mapping.emplace_back();

        for (NodeID n : G->nodes()) {
            mapping.push_back(n / 4);
            reverse_mapping[n / 4].emplace_back(n);
        }

        graphAccessPtr cntr = contraction::contractGraph(
            G, mapping, reverse_mapping);

        ASSERT_EQ(cntr->number_of_edges(), 2);
        ASSERT_EQ(cntr->number_of_nodes(), 2);

        ASSERT_EQ(cntr->getEdgeTarget(0), 1);
        ASSERT_EQ(cntr->getEdgeTarget(1), 0);

        // in weighted graph edge weight is 3, in unweighted 2
        if (graph == "") {
            ASSERT_EQ(cntr->getEdgeWeight(0), 2);
            ASSERT_EQ(cntr->getEdgeWeight(1), 2);
        } else {
            ASSERT_EQ(cntr->getEdgeWeight(0), 3);
            ASSERT_EQ(cntr->getEdgeWeight(1), 3);
        }
    }
}
