/******************************************************************************
 * algorithms.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include "ks_minimum_cut.h"
#include "matula_approx.h"
#include "minimum_cut.h"
#include "noi_minimum_cut.h"
#include "padberg_rinaldi.h"
#include "parallel/algorithm/exact_parallel_minimum_cut.h"
#include "stoer_wagner_minimum_cut.h"
#include "viecut.h"

static minimum_cut * selectMincutAlgorithm(std::string argv_str) {
#ifndef PARALLEL
    if (argv_str == "ks")
        return new ks_minimum_cut();
    if (argv_str == "noi")
        return new noi_minimum_cut();
    if (argv_str == "matula")
        return new matula_approx();
    if (argv_str == "vc")
        return new viecut();
    if (argv_str == "pr")
        return new padberg_rinaldi();
#endif
#ifdef PARALLEL
    if (argv_str == "vc")
        return new viecut();
    if (argv_str == "exact")
        return new exact_parallel_minimum_cut();
#endif
    return new minimum_cut();
}
