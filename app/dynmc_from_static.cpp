/******************************************************************************
 * dynmc_from_static.cpp
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2020 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <stddef.h>

#include <ext/alloc_traits.h>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#ifdef PARALLEL
#include "parallel/algorithm/exact_parallel_minimum_cut.h"
#include "parallel/algorithm/parallel_cactus.h"
#else
#include "algorithms/global_mincut/cactus/cactus_mincut.h"
#include "algorithms/global_mincut/noi_minimum_cut.h"
#endif

#include "algorithms/global_mincut/dynamic/dynamic_mincut.h"
#include "common/configuration.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tlx/logger.hpp"
#include "tlx/string.hpp"
#include "tools/random_functions.h"

int main(int argn, char** argv) {
    tlx::CmdlineParser cmdl;
    auto cfg = configuration::getConfig();
    std::string init_graph = "";
    bool run_static = false;
    double insert_factor = 0.0;
    double remove_factor = 0.0;
    size_t batch_size = 1;
    cmdl.add_string('g', "initial_graph", init_graph, "path to graph file");
    cmdl.add_double('i', "insert_factor", insert_factor,
                    "factor of edges that are added dynamically");
    cmdl.add_double('r', "remove_factor", remove_factor,
                    "factor of edges that are removed dynamically");
    cmdl.add_size_t('b', "batch_size", batch_size, "batch size");
    cmdl.add_size_t('s', "seed", configuration::getConfig()->seed, "rnd seed");

    if (!cmdl.process(argn, argv))
        return -1;

    random_functions::setSeed(configuration::getConfig()->seed);
    mutableGraphPtr G = graph_io::readGraphWeighted<mutable_graph>(init_graph);

    size_t insert_edges = std::floor(insert_factor * G->m());
    size_t remove_edges = std::floor(remove_factor * G->m());

    LOG1 << "insert " << insert_edges << " remove " << remove_edges; 

}