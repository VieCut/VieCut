/******************************************************************************
 * graph_io.h
 *
 * Source of KaHIP -- Karlsruhe High Quality Partitioning.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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

#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "data_structure/flow_graph.h"
#include "data_structure/graph_access.h"
#include "definitions.h"
#include "tools/string.h"

class graph_io
{
public:
    graph_io();

    virtual ~graph_io();

    static
    std::shared_ptr<graph_access> readGraphWeighted(std::string filename);

    static
    int writeGraphWeighted(std::shared_ptr<graph_access> G, std::string filename);

    static
    int writeGraph(std::shared_ptr<graph_access> G, std::string filename);

    static
    int writeGraphDimacsKS(std::shared_ptr<graph_access> G, std::string filename, std::string format = "FORMAT");

    static
    int writeGraphDimacsJP(std::shared_ptr<graph_access> G, std::string filename);

    static void writeCut(std::shared_ptr<graph_access> G, std::string filename) {
        std::ofstream f(filename.c_str());
        std::cout << "writing partition to " << filename << " ... " << std::endl;

        for (NodeID node : G->nodes()) {
            f << G->getNodeInCut(node) << std::endl;
        }

        f.close();
    }

    static std::shared_ptr<flow_graph> createFlowGraph(std::shared_ptr<graph_access> G) {
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
    std::vector<vectortype> readVector(std::string filename);

    template <typename vectortype>
    void writeVector(std::vector<vectortype>& vec, std::string filename);
};
