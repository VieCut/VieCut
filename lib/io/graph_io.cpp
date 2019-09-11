/******************************************************************************
 * graph_io.cpp
 *
 * Source of VieCut.
 *
 * Adapted from KaHIP.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <sstream>
#include <string>

#include "io/graph_io.h"
#include "tlx/logger.hpp"

graph_io::graph_io() { }

graph_io::~graph_io() { }

int graph_io::writeGraphWeighted(std::shared_ptr<graph_access> G,
                                 std::string filename) {
    std::ofstream f(filename.c_str());
    f << G->number_of_nodes() << " "
      << G->number_of_edges() / 2 << " 1" << std::endl;

    for (NodeID node : G->nodes()) {
        for (EdgeID e : G->edges_of(node)) {
            f << (G->getEdgeTarget(e) + 1) << " " << G->getEdgeWeight(e) << " ";
        }
        f << std::endl;
    }

    f.close();
    return 0;
}

int graph_io::writeGraph(std::shared_ptr<graph_access> G,
                         std::string filename) {
    std::ofstream f(filename.c_str());
    f << G->number_of_nodes() << " " << G->number_of_edges() / 2 << std::endl;

    for (NodeID node : G->nodes()) {
        for (EdgeID e : G->edges_of(node)) {
            f << (G->getEdgeTarget(e) + 1) << " ";
        }
        f << std::endl;
    }

    f.close();
    return 0;
}

int graph_io::writeGraphDimacsJP(std::shared_ptr<graph_access> G,
                                 std::string filename) {
    std::ofstream f(filename.c_str());
    f << G->number_of_nodes() << " " << G->number_of_edges() / 2 << std::endl;

    for (NodeID node : G->nodes()) {
        for (EdgeID e : G->edges_of(node)) {
            if (G->getEdgeTarget(e) > node) {
                f << node << " " << G->getEdgeTarget(e) << " "
                  << G->getEdgeWeight(e) << std::endl;
            }
        }
    }
    f << G->number_of_nodes() << std::endl;

    for (NodeID node : G->nodes()) {
        f << node << std::endl;
    }

    f.close();
    return 0;
}

int graph_io::writeGraphDimacsKS(std::shared_ptr<graph_access> G,
                                 std::string filename,
                                 std::string format) {
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

// orig. from http://tinodidriksen.com/uploads/code/cpp/speed-string-to-int.cpp
int naive(const std::string& str, size_t* line_ptr) {
    int x = 0;

    while (str[*line_ptr] >= '0' && str[*line_ptr] <= '9') {
        x = (x * 10) + (str[*line_ptr] - '0');
        ++(*line_ptr);
    }
    ++(*line_ptr);
    return x;
}

std::shared_ptr<graph_access> graph_io::readGraphWeighted(std::string file) {
    std::shared_ptr<graph_access> G = std::make_shared<graph_access>();
    std::string line;

    // open file for reading
    std::ifstream in(file.c_str());
    if (!in) {
        std::cerr << "Error opening " << file << std::endl;
        exit(2);
    }

    uint64_t nmbNodes;
    uint64_t nmbEdges;

    std::getline(in, line);
    // skip comments
    while (line[0] == '%') {
        std::getline(in, line);
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
    nmbEdges *= 2;     // since we have forward and backward edges

    NodeID node_counter = 0;
    EdgeID edge_counter = 0;

    G->start_construction(nmbNodes, nmbEdges);

    while (std::getline(in, line)) {
        if (line[0] == '%') {         // a comment in the file
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
                LOG1 << "The graph file contains self-loops. "
                     << "This is not supported. "
                     << "Please remove them from the file.";
            }

            EdgeWeight edge_weight = 1;
            if (read_ew) {
                edge_weight = naive(line, &line_ptr);
            }
            edge_counter++;
            EdgeID e = G->new_edge(node, target - 1);

            G->setEdgeWeight(e, edge_weight);
        }

        if (in.eof()) {
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
