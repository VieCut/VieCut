/******************************************************************************
 * minimum_cut.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <memory>

#include "common/definitions.h"
#include "data_structure/graph_access.h"

class minimum_cut {
 public:
    virtual ~minimum_cut() { }

    virtual EdgeWeight perform_minimum_cut(graphAccessPtr) {
        return perform_minimum_cut();
    }

    virtual EdgeWeight perform_minimum_cut(mutableGraphPtr) {
        return perform_minimum_cut();
    }

    virtual EdgeWeight perform_minimum_cut() {
#ifdef PARALLEL
        LOG1 << "Please select a parallel minimum cut"
             << " algorithm [inexact, exact, cactus]!";
        LOG1 << "inexact - Run heuristic VieCut algorithm";
        LOG1 << "exact - Run shared-memory exact algorithm";
        LOG1 << "cactus - Find all minimum cuts and build cactus graph!";
#else
        LOG1 << "Please select a minimum cut global_mincut"
             << " [vc, noi, pr, matula, ks, cactus]!";
        LOG1 << "vc - Run heuristic VieCut algorithm";
        LOG1 << "noi - Run algorithm of Nagamochi, Ono and Ibaraki";
        LOG1 << "pr - Repeated run of routines of Padberg and Rinaldi";
        LOG1 << "matula - Run algorithm of Matula";
        LOG1 << "ks - Run algorithm of Karger and Stein";
        LOG1 << "cactus - Find all minimum cuts and build cactus graph!";
#endif
        exit(1);
        return 42;
    }
};
