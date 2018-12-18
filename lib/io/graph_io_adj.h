/******************************************************************************
 * graph_io_adj.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
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

#ifndef GRAPHIOADJ_H_
#define GRAPHIOADJ_H_

#include <fstream>
#include <iostream>
#include <limits>
#include <ostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "data_structure/adjlist_graph.h"
#include "definitions.h"

class graph_io_adj
{
public:
    graph_io_adj();
    virtual ~graph_io_adj();

    static int readGraphWeightedDimacs(adjlist_graph& G, std::string filename) {
        std::string line;

        // open file for reading
        std::ifstream in (filename.c_str());
        if (!in) {
            std::cerr << "Error opening " << filename << std::endl;
            return 1;
        }

        long nmbNodes;
        long nmbEdges;

        std::getline(in, line);
        // skip comments
        while (line[0] == '%') {
            std::getline(in, line);
        }

        std::string devnull;
        std::stringstream ss(line);
        ss >> devnull >> devnull;
        ss >> nmbNodes;
        ss >> nmbEdges;

        nmbEdges *= 2;
        EdgeID edge_counter = 0;

        std::cout << "NODES : " << nmbNodes << " AND EDGYNESS: " << nmbEdges << std::endl;

        G.start_construction(nmbNodes, nmbEdges);

        for (NodeID i = 0; i < nmbNodes; ++i) {
            G.new_node();
        }

        while (std::getline(in, line)) {

            if (line[0] == '%') { // a comment in the file
                continue;
            }

            std::stringstream ss(line);

            char a;
            NodeID source;
            NodeID target;
            EdgeWeight weight;

            while (ss >> a >> source >> target >> weight) {

                // check for self-loops
                if (target - 1 == source - 1) {
                    std::cerr << "The graph file contains self-loops. This is not supported. Please remove them from the file." << std::endl;
                }

                if (a != 'a') {
                    std::cerr << "file is not in format" << std::endl;
                }

                edge_counter += 2;
                G.new_edge(source - 1, target - 1, weight);
                G.new_edge(target - 1, source - 1, weight);
            }

            if (in.eof()) {
                break;
            }
        }

        if (edge_counter != (EdgeID)nmbEdges) {
            std::cerr << "number of specified edges mismatch" << std::endl;
            std::cerr << edge_counter << " " << nmbEdges << std::endl;
            exit(0);
        }

        G.finish_construction();
        return 0;
    }

    static
    int readGraphWeighted(adjlist_graph& G, std::string filename) {
        std::string line;

        // open file for reading
        std::ifstream in (filename.c_str());
        if (!in) {
            std::cerr << "Error opening " << filename << std::endl;
            return 1;
        }

        long nmbNodes;
        long nmbEdges;

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

        /* if( 2*nmbEdges > std::numeric_limits<int>::max() || nmbNodes > std::numeric_limits<int>::max()) {
           std::cerr <<  "The graph is too large. Currently only 32bit supported!"  << std::endl;
           exit(0);
           }*/

        bool read_ew = false;
        bool read_nw = false;

        if (ew == 1) {
            read_ew = true;
        }
        else if (ew == 11) {
            read_ew = true;
            read_nw = true;
        }
        else if (ew == 10) {
            read_nw = true;
        }
        nmbEdges *= 2; // since we have forward and backward edges

        NodeID node_counter = 0;
        EdgeID edge_counter = 0;
        long long total_nodeweight = 0;

        G.start_construction(nmbNodes, nmbEdges);

        while (std::getline(in, line)) {

            if (line[0] == '%') { // a comment in the file
                continue;
            }

            node_counter++;

            std::stringstream ss(line);

            NodeWeight weight = 1;
            if (read_nw) {
                ss >> weight;
                total_nodeweight += weight;
                if (total_nodeweight > (long long)std::numeric_limits<NodeWeight>::max()) {
                    std::cerr << "The sum of the node weights is too large (it exceeds the node weight type)." << std::endl;
                    std::cerr << "Currently not supported. Please scale your node weights." << std::endl;
                    exit(0);
                }
            }
            NodeID node = G.new_node(weight);
            NodeID target;
            while (ss >> target) {
                // check for self-loops
                if (target - 1 == node) {
                    std::cerr << "The graph file contains self-loops. This is not supported. Please remove them from the file." << std::endl;
                }

                EdgeWeight edge_weight = 1;
                if (read_ew) {
                    ss >> edge_weight;
                }
                edge_counter++;
                G.new_edge(node, target - 1, edge_weight);
            }

            if (in.eof()) {
                break;
            }
        }

        if (edge_counter != (EdgeID)nmbEdges) {
            std::cerr << "number of specified edges mismatch" << std::endl;
            std::cerr << edge_counter << " " << nmbEdges << std::endl;
            exit(0);
        }

        if (node_counter != (NodeID)nmbNodes) {
            std::cerr << "number of specified nodes mismatch" << std::endl;
            std::cerr << node_counter << " " << nmbNodes << std::endl;
            exit(0);
        }

        G.finish_construction();
        return 0;
    }

    static
    int writeGraphWeighted(adjlist_graph& G, std::string filename) {
        return 0;
    }
};

#endif
