/******************************************************************************
 * algorithms.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
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
