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
