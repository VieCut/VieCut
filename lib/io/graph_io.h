/******************************************************************************
 * graph_io.h
 *
 * Source of VieCut.
 *
 * Adapted from KaHIP.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@univie.ac.at>
 * Copyright (C) 2017-2020 Alexander Noe <alexander.noe@univie.ac.at>
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
    graph_io() { }

    virtual ~graph_io() { }

    template <class Graph = graph_access>
    static std::shared_ptr<Graph> readGraphWeighted(std::string file) {
        std::shared_ptr<Graph> G = std::make_shared<Graph>();
        std::string line;
        // open file for reading
        std::ifstream instream(file.c_str());
        if (!instream) {
            std::cerr << "Error opening " << file << std::endl;
            exit(2);
        }

        uint64_t nmbNodes;
        uint64_t nmbEdges;

        std::getline(instream, line);
        // skip comments
        while (line[0] == '%') {
            std::getline(instream, line);
        }

        int ew = 0;
        std::stringstream ss(line);
        ss >> nmbNodes;
        ss >> nmbEdges;
        ss >> ew;

        bool read_ew = false;

        if (ew == 1) {
            read_ew = true;
        } else if (ew == 11) {
            read_ew = true;
        }
        nmbEdges *= 2;  // since we have forward and backward edges

        NodeID node_counter = 0;
        EdgeID edge_counter = 0;

        G->start_construction(nmbNodes, nmbEdges);

        while (std::getline(instream, line)) {
            if (line[0] == '%') {     // a comment in the file
                continue;
            }

            NodeID node = G->new_node();
            node_counter++;
            G->setPartitionIndex(node, 0);

            size_t line_ptr = 0;

            // remove leading whitespaces
            while (line[line_ptr] == ' ')
                ++line_ptr;

            while (line_ptr < line.size()) {
                NodeID target = naive(line, &line_ptr);
                if (!target) break;
                // check for self-loops
                if (target - 1 == node) {
                    LOG0 << "The graph file contains self-loops. "
                         << "This is not supported. "
                         << "Please remove them from the file.";
                    continue;
                }

                EdgeWeight edge_weight = 1;
                if (read_ew) {
                    edge_weight = naive(line, &line_ptr);
                }
                edge_counter++;
                G->new_edge(node, target - 1, edge_weight);
            }

            if (instream.eof()) {
                break;
            }
        }

        if (edge_counter != (EdgeID)nmbEdges) {
            std::cerr << "number of specified edges mismatch" << std::endl;
            std::cerr << edge_counter << " " << nmbEdges << std::endl;
            // As we discard self-loops, this might happen. Thus, no exiting!
            // Node counter not matching is however probably a problem with G.
        }

        if (node_counter != (NodeID)nmbNodes) {
            std::cerr << "number of specified nodes mismatch" << std::endl;
            std::cerr << node_counter << " " << nmbNodes << std::endl;
            exit(4);
        }

        G->finish_construction();
        G->computeDegrees();
        return G;
    }

    static int writeGraphWeighted(graphAccessPtr G,
                                  std::string filename) {
        std::ofstream f(filename.c_str());
        f << G->number_of_nodes() << " "
          << G->number_of_edges() / 2 << " 1" << std::endl;

        for (NodeID node : G->nodes()) {
            for (EdgeID e : G->edges_of(node)) {
                f << (G->getEdgeTarget(e) + 1) << " "
                  << G->getEdgeWeight(e) << " ";
            }
            f << std::endl;
        }

        f.close();
        return 0;
    }

    static
    int writeGraph(graphAccessPtr G, std::string filename) {
        std::ofstream f(filename.c_str());
        f << G->number_of_nodes() << " "
          << G->number_of_edges() / 2 << " 0" << std::endl;

        for (NodeID node : G->nodes()) {
            for (EdgeID e : G->edges_of(node)) {
                f << (G->getEdgeTarget(e) + 1) << " ";
            }
            f << std::endl;
        }

        f.close();
        return 0;
    }

    static
    int writeGraphDimacsKS(graphAccessPtr G,
                           std::string filename,
                           std::string format = "FORMAT") {
        std::ofstream f(filename.c_str());
        f << "p " << format << " " << G->number_of_nodes() << " "
          << G->number_of_edges() / 2 << std::endl;

        for (NodeID node : G->nodes()) {
            for (EdgeID e : G->edges_of(node)) {
                if (G->getEdgeTarget(e) > node) {
                    f << "a " << node + 1 << " " << G->getEdgeTarget(e) + 1
                      << " " << G->getEdgeWeight(e) << std::endl;
                }
            }
        }

        f.close();
        return 0;
    }

    static void writeCut(graphAccessPtr G,
                         std::string filename) {
        std::ofstream f(filename.c_str());
        LOG1 << "writing partition to " << filename << " ... ";

        for (NodeID node : G->nodes()) {
            f << G->getNodeInCut(node) << std::endl;
        }

        f.close();
    }

    static std::shared_ptr<flow_graph> createFlowGraph(
        graphAccessPtr G) {
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

 private:
// orig. from http://tinodidriksen.com/uploads/code/cpp/speed-string-to-int.cpp
    static int naive(const std::string& str, size_t* line_ptr) {
        int x = 0;

        while (str[*line_ptr] >= '0' && str[*line_ptr] <= '9') {
            x = (x * 10) + (str[*line_ptr] - '0');
            ++(*line_ptr);
        }
        ++(*line_ptr);
        return x;
    }
};
