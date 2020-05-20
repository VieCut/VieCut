/******************************************************************************
 * create_augment_edges.cpp
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/


#include "algorithms/global_mincut/dynamic/dynamic_mincut.h"
#include "common/configuration.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tlx/logger.hpp"
#include "tools/random_functions.h"

int main(int argn, char** argv) {
    tlx::CmdlineParser cmdl;
    auto cfg = configuration::getConfig();
    cmdl.add_param_string("graph", cfg->graph_filename, "path to graph file");
    size_t insert_edges = 0;
    size_t delete_edges = 0;
#ifdef PARALLEL
    size_t procs = 1;
    cmdl.add_size_t('p', "proc", procs, "number of processes");
#endif
    cmdl.add_size_t('d', "delete", delete_edges, "number of edges to delete");
    cmdl.add_size_t('i', "insert", insert_edges, "number of edges to insert");
    cmdl.add_flag('v', "verbose", cfg->verbose, "more verbose logs");
    cmdl.add_size_t('r', "seed", cfg->seed, "random seed");

    if (!cmdl.process(argn, argv))
        return -1;

    if (delete_edges > insert_edges) {
        LOG1 << "Error: Trying to delete more edges that were inserted!";
        exit(1);
    }

    
    // TODO(anoe): run cactus algorithm, find set of edges that augment
    // connectivity, repeat until number of edges is reached


}
