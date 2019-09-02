/******************************************************************************
 * most_balanced_minimum_cut.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <memory>

#include "algorithms/global_mincut/cactus/balanced_cut_dfs.h"
#include "common/configuration.h"
#include "data_structure/mutable_graph.h"

class most_balanced_minimum_cut {
 public:
    void findCutFromCactus(std::shared_ptr<mutable_graph> G,
                           EdgeWeight mincut) {
        if (!configuration::getConfig()->save_cut) {
            LOG1 << "Error: can't find most balanced minimum cut "
                 << "when save_cut is not set";
            exit(1);
        }
        balanced_cut_dfs dfs(G, mincut);
        dfs.runDFS();
    }
};
