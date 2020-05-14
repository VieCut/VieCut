/******************************************************************************
 * dynamic_mincut.cpp
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2020 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <omp.h>

#include "algorithms/global_mincut/algorithms.h"
#include "algorithms/global_mincut/minimum_cut.h"
#include "common/configuration.h"
#include "common/definitions.h"
#include "io/graph_io.h"
#include "tools/timer.h"

int main(int argn, char** argv) {
    static constexpr bool debug = false;
    tlx::CmdlineParser cmdl;
    int num_iterations = 1;
    auto cfg = configuration::getConfig();
    cmdl.add_param_string("graph", cfg->graph_filename, "path to graph file");
#ifdef PARALLEL
    std::vector<std::string> procs;
    cmdl.add_stringlist('p', "proc", procs, "number of processes");
#endif
    cmdl.add_param_string("algo", cfg->algorithm, "algorithm name");
    cmdl.add_size_t('i', "iter", num_iterations, "number of iterations");
    cmdl.add_flag('v', "verbose", cfg->verbose, "more verbose logs");
    cmdl.add_size_t('r', "seed", cfg->seed, "random seed");

    if (!cmdl.process(argn, argv))
        return -1;

    cfg->save_cut = true;

    // TODO(anoe): make this actually do something :)
}
