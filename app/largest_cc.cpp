/******************************************************************************
 * largest_cc.cpp
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2017-2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <sstream>

#include "algorithms/global_mincut/minimum_cut_helpers.h"
#include "algorithms/misc/strongly_connected_components.h"
#include "data_structure/graph_access.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tlx/string.hpp"
#include "tools/random_functions.h"
#include "tools/string.h"
#include "tools/timer.h"

int main(int argn, char** argv) {
    static constexpr bool debug = false;

    LOG << "Find and save largest CC!";

    tlx::CmdlineParser cmdl;
    std::string path;
    bool weighted = false;
    cmdl.add_param_string("graph", path, "path to graph file");
    cmdl.add_flag('w', "weighted", weighted, "weighted graph");
    configuration::getConfig()->graph_filename = path;

    if (!cmdl.process(argn, argv))
        return -1;

    std::shared_ptr<graph_access> G =
        graph_io::readGraphWeighted(path);

    strongly_connected_components scc;

    auto out = scc.largest_scc(G);
    if (weighted) {
        graph_io::writeGraphWeighted(out, tlx::split(".", path)[0] + ".cc");
    } else {
        graph_io::writeGraph(out, tlx::split(".", path)[0] + ".cc");
    }
}
