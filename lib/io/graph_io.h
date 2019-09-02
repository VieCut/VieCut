/******************************************************************************
 * graph_io.h
 *
 * Source of VieCut.
 *
 * Adapted from KaHIP.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@univie.ac.at>
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "common/definitions.h"
#include "data_structure/flow_graph.h"
#include "data_structure/graph_access.h"
#include "tools/string.h"

class graph_io {
 public:
    graph_io();

    virtual ~graph_io();

    static
    std::shared_ptr<graph_access> readGraphWeighted(std::string filename);

    static
    int writeGraphWeighted(std::shared_ptr<graph_access> G,
                           std::string filename);

    static
    int writeGraph(std::shared_ptr<graph_access> G, std::string filename);

    static
    int writeGraphDimacsKS(std::shared_ptr<graph_access> G,
                           std::string filename,
                           std::string format = "FORMAT");

    static
    int writeGraphDimacsJP(std::shared_ptr<graph_access> G,
                           std::string filename);

    static void writeCut(std::shared_ptr<graph_access> G,
                         std::string filename) {
        std::ofstream f(filename.c_str());
        LOG1 << "writing partition to " << filename << " ... ";

        for (NodeID node : G->nodes()) {
            f << G->getNodeInCut(node) << std::endl;
        }

        f.close();
    }

    static std::shared_ptr<flow_graph> createFlowGraph(
        std::shared_ptr<graph_access> G) {
        std::shared_ptr<flow_graph> fg = std::make_shared<flow_graph>();
        fg->start_construction(G->number_of_nodes());

        for (NodeID n : G->nodes()) {
            for (EdgeID e : G->edges_of(n)) {
                NodeID tgt = G->getEdgeTarget(e);
                fg->new_edge(n, tgt, G->getEdgeWeight(e));
            }
        }

        fg->finish_construction();

        VIECUT_ASSERT_EQ(fg->number_of_nodes(), G->number_of_nodes());
        VIECUT_ASSERT_EQ(fg->number_of_edges(), 2 * G->number_of_edges());

        return fg;
    }

    template <typename vectortype>
    static std::vector<vectortype> readVector(std::string filename) {
        std::vector<vectortype> vec;
        std::string line;
        // open file for reading
        std::ifstream instream(filename.c_str());
        if (!instream) {
            std::cerr << "Error opening vectorfile" << filename << std::endl;
            exit(5);
        }

        std::getline(instream, line);
        while (!instream.eof()) {
            if (line[0] == '%') {         // Comment
                continue;
            }

            vectortype value = (vectortype)atof(line.c_str());
            vec.emplace_back(value);
            std::getline(instream, line);
        }

        instream.close();
        return vec;
    }

    template <typename vectortype>
    void writeVector(const std::vector<vectortype>& vec, std::string filename) {
        std::ofstream f(filename.c_str());
        for (unsigned i = 0; i < vec.size(); ++i) {
            f << vec[i] << std::endl;
        }
        f.close();
    }
};
