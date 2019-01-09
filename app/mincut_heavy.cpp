/******************************************************************************
 * mincut_heavy.cpp
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2017-2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <math.h>
#include <memory>
#include <omp.h>
#include <regex.h>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>

#include "algorithms/global_mincut/algorithms.h"
#include "data_structure/graph_access.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tlx/cmdline_parser.hpp"
#include "tlx/logger.hpp"
#include "tools/quality_metrics.h"
#include "tools/string.h"
#include "tools/timer.h"
#ifdef PARALLEL
#include "parallel/coarsening/sparsify.h"
#else
#include "coarsening/sparsify.h"
#endif

int main(int argn, char** argv) {

    static constexpr bool debug = false;

    tlx::CmdlineParser cmdl;

    std::string graph_filename;
    std::string algo;
    std::string pq = "default";
    double contraction_factor = 0.9;
    int iter = 1;
    bool save_cut = false;

    cmdl.add_param_string("graph", graph_filename, "path to graph file");

#ifdef PARALLEL
    std::vector<std::string> procs;
    cmdl.add_stringlist('p', "proc", procs, "number of processes");
    cmdl.add_param_string("algo", algo, "algorithm name ('vc', 'exact')");
#else
    cmdl.add_param_string("algo", algo, "algorithm name ('vc', 'noi', 'matula', 'pr', 'ks')");
#endif

    cmdl.add_string('q', "pq", pq, "name of priority queue implementation");
    cmdl.add_int('i', "iter", iter, "number of iterations");
    cmdl.add_double('c', "contraction factor", contraction_factor, "contract until only n/(1-contraction_factor) vertices are left");
    cmdl.add_bool('s', "save_cut", save_cut, "compute and save minimum cut");

    if (!cmdl.process(argn, argv))
        return -1;

    std::vector<int> numthreads;

    std::shared_ptr<graph_access> G = graph_io::readGraphWeighted(graph_filename);
    timer t;

    LOG << "io time: " << t.elapsed();
    t.restart();

    sparsify sf;

    std::vector<NodeID> indices;
    auto G2 = sf.remove_heavy_vertices(G, contraction_factor, indices);
    LOG << "sparsification time: " << t.elapsed();

#ifdef PARALLEL
    LOG1 << "PARALLEL DEFINED";
    size_t i;
    try {
        for (i = 0; i < procs.size(); ++i) {
            numthreads.emplace_back(std::stoi(procs[i]));
        }
    } catch (...) {
        LOG1 << procs[i] << " is not a valid number of workers! Continuing without.";
    }
#else
    LOG1 << "PARALLEL NOT DEFINED";
#endif

    if (numthreads.empty())
        numthreads.emplace_back(1);

    for (int i = 0; i < iter; ++i) {
        for (size_t numthread : numthreads) {
            random_functions::setSeed(i);
            t.restart();
            minimum_cut* mc = selectMincutAlgorithm(algo);
            mc->perform_minimum_cut(G2, save_cut);
            omp_set_num_threads(numthread);

            for (NodeID n : G->nodes()) {
                if (indices[n] < G->number_of_nodes()) {
                    G->setNodeInCut(n, G2->getNodeInCut(indices[n]));
                }
                else {
                    size_t in = 0;
                    size_t out = 0;
                    for (EdgeID e : G->edges_of(n)) {
                        NodeID target = G->getEdgeTarget(e);
                        PartitionID pid = indices[target];
                        if (pid < G->number_of_nodes()) {
                            if (G2->getNodeInCut(pid)) {
                                in++;
                            }
                            else {
                                out++;
                            }
                        }
                    }
                    if (in > out) {
                        G->setNodeInCut(n, 1);
                    }
                    else if (out > in) {
                        G->setNodeInCut(n, 0);
                    }
                    else {
                        LOG1 << "TIE - THIS SHOULD NOT HAPPEN";
                        G->setNodeInCut(n, 1);
                    }
                }
            }

            size_t inside = 0, outside = 0;
            for (NodeID n : G->nodes()) {
                if (G->getNodeInCut(n)) {
                    inside++;
                    G->setPartitionIndex(n, 0);
                }
                else {
                    outside++;
                    G->setPartitionIndex(n, 1);
                }
            }

            LOG1 << "smaller side of cut has " << std::min(inside, outside) << " nodes.";

            quality_metrics qm;

            if (save_cut) {
                graph_io::writeCut(G, graph_filename + ".cut");
            }

            std::string graphname = string::basename(graph_filename);
            std::cout << "RESULT source=taa algo=viecut graph="
                      << graphname << " time=" << t.elapsed()
                      << " cut=" << qm.edge_cut(*G)
                      << " n=" << G->number_of_nodes()
                      << " m=" << G->number_of_edges() << " processes="
                      << omp_get_max_threads() << std::endl;
        }
    }
}
