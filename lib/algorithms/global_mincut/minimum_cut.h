/******************************************************************************
 * minimum_cut.h
 *
 * Source of VieCut
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

#include "data_structure/graph_access.h"
#include "definitions.h"
#include <memory>

class minimum_cut
{
public:
    virtual ~minimum_cut() { }

    virtual EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G, std::string, bool) {
        return perform_minimum_cut(G);
    }

    virtual EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access>) {
#ifdef PARALLEL
        std::cout << "Please select a parallel minimum cut algorithm [vc, exact]!" << std::endl;
#else
        std::cout << "Please select a minimum cut global_mincut [vc, noi, pr, matula, ks]!" << std::endl;
#endif
        exit(1);
        return 42;
    }
};
