/******************************************************************************
 * mincut_recursive.cpp
 *
 * Source of VieCut.
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

#include <algorithm>
#include <argtable2.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <math.h>
#include <memory>
#include <regex.h>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>

#ifdef PARALLEL
#include "parallel/algorithm/exact_parallel_minimum_cut.h"
#else
#include "algorithms/global_mincut/noi_minimum_cut.h"
#endif
#include "algorithms/global_mincut/viecut.h"
#include "algorithms/misc/strongly_connected_components.h"
#include "data_structure/graph_access.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tlx/logger.hpp"
#include "tools/graph_extractor.h"
#include "tools/macros_assertions.h"
#include "tools/timer.h"

int main(int argn, char** argv) {

    static constexpr bool debug = true;
    timer t;
    tlx::CmdlineParser cmdl;
    std::string graph_filename;
    bool output = false;
    cmdl.add_param_string("graph", graph_filename, "path to graph file");
    cmdl.add_bool('o', "output", output, "Write intermediate graphs to disk");

    if (!cmdl.process(argn, argv))
        return -1;

    std::vector<std::shared_ptr<graph_access> > graphs;

    std::shared_ptr<graph_access> G = graph_io::readGraphWeighted(graph_filename);
    graphs.push_back(G);

#ifdef PARALLEL
    exact_parallel_minimum_cut mc;
#else
    noi_minimum_cut mc;
#endif
    LOG1 << "io time: " << t.elapsed();

    EdgeWeight minimum_degree = G->getMinDegree();
    EdgeWeight cut = 0;
    EdgeWeight last_cut = 0;
    timer t_this;

    while (cut < minimum_degree) {
        t_this.restart();
        LOG << G->number_of_nodes() << " " << G->number_of_edges();
        cut = mc.perform_minimum_cut(G, false);
        size_t inside = 0;

        for (NodeID n : G->nodes()) {
            NodeWeight inCut = G->getNodeInCut(n) ? 1 : 0;
            G->setPartitionIndex(n, inCut);
            inside += inCut;
        }

        size_t larger = (G->number_of_nodes() < inside * 2) ? 1 : 0;

        graph_extractor ge;
        strongly_connected_components scc;

        auto new_g = graphs.emplace_back(ge.extract_block(G, larger).first);

        graphs.emplace_back(scc.largest_scc(new_g));
        if (cut > last_cut && last_cut != 0 && cut < G->getMinDegree() && output) {
            graph_io::writeGraph(G, graph_filename + "_" + std::to_string(cut));
        }
        G = graphs.back();

        LOG1 << "New graph: (" << G->number_of_nodes() << ";" << G->number_of_edges() << ") " << cut << " < " << G->getMinDegree();

        last_cut = cut;
    }
}
