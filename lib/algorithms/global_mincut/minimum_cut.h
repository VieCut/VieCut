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

    virtual EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access> G, bool save_cut, std::string, bool) {
        return perform_minimum_cut(G, save_cut);
    }

    virtual EdgeWeight perform_minimum_cut(std::shared_ptr<graph_access>, bool) {
#ifdef PARALLEL
        std::cout << "Please select a parallel minimum cut algorithm [vc, exact]!" << std::endl;
        std::cout << "vc - Run heuristic VieCut algorithm" << std::endl;
        std::cout << "extact - Run shared-memory exact algorithm" << std::endl;
#else
        std::cout << "Please select a minimum cut global_mincut [vc, noi, pr, matula, ks]!" << std::endl;
        std::cout << "vc - Run heuristic VieCut algorithm" << std::endl;
        std::cout << "noi - Run algorithm of Nagamochi, Ono and Ibaraki" << std::endl;
        std::cout << "pr - Repeated run of local contraction routines of Padberg and Rinaldi" << std::endl;
        std::cout << "matula - Run algorithm of Matula" << std::endl;
        std::cout << "ks - Run algorithm of Karger and Stein" << std::endl;
#endif
        exit(1);
        return 42;
    }
};
