/******************************************************************************
 * gnp_torus.cpp
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

#ifdef PARALLEL
#include "parallel/algorithm/parallel_cactus.h"
#else
#include "algorithms/global_mincut/cactus/cactus_mincut.h"
#endif
#include "algorithms/global_mincut/viecut.h"
#include "algorithms/misc/strongly_connected_components.h"
#include "data_structure/graph_access.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tlx/logger.hpp"
#include "tools/graph_extractor.h"
#include "tools/macros_assertions.h"
#include "tools/random_functions.h"
#include "tools/timer.h"

int main(int argn, char** argv) {
    timer t;
    tlx::CmdlineParser cmdl;
    auto cfg = configuration::getConfig();

    size_t vertices_per_block = 0;
    size_t cut = 0;
    size_t blocks_per_dimension = 0;
    size_t seed = 0;
    std::string output_folder = "";

    cmdl.add_size_t('n', "vertices", vertices_per_block,
                    "number of vertices in block");
    cmdl.add_size_t('c', "cut", cut, "cut edges between gnp blocks");
    cmdl.add_size_t('b', "blocks", blocks_per_dimension,
                    "number of blocks per dimension");
    cmdl.add_size_t('s', "seed", seed, "random seed");
    cmdl.add_string('o', "output", output_folder, "folder to write output to");

    if (!cmdl.process(argn, argv))
        return -1;

    random_functions::setSeed(seed);
    std::shared_ptr<graph_access> G = std::make_shared<graph_access>();

    // blocks * blocks gnp graphs with vertices vertices each
    size_t blocks = blocks_per_dimension * blocks_per_dimension;
    size_t num_vertices = blocks * vertices_per_block;

    LOG1 << "Generating inter-block edges...";

    std::vector<std::vector<NodeID> > interblockedges(num_vertices);

    for (size_t i = 0; i < blocks; ++i) {
        NodeID my_start = vertices_per_block * i;
        NodeID x = i % blocks_per_dimension;
        NodeID y = i / blocks_per_dimension;
        // right
        for (size_t k = 0; k < cut; ++k) {
            size_t r1 = random_functions::nextInt(0, vertices_per_block - 1);
            size_t r2 = random_functions::nextInt(0, vertices_per_block - 1);
            NodeID v1 = my_start + r1;
            NodeID v2 = ((y * blocks_per_dimension)
                         + ((x + 1) % blocks_per_dimension))
                        * vertices_per_block + r2;

            interblockedges[v1].emplace_back(v2);
            interblockedges[v2].emplace_back(v1);
        }

        // down
        for (size_t k = 0; k < cut; ++k) {
            size_t r1 = random_functions::nextInt(0, vertices_per_block - 1);
            size_t r2 = random_functions::nextInt(0, vertices_per_block - 1);
            NodeID v1 = my_start + r1;
            NodeID v2 = ((((y + 1) % blocks_per_dimension)
                          * blocks_per_dimension) + x)
                        * vertices_per_block + r2;

            interblockedges[v1].emplace_back(v2);
            interblockedges[v2].emplace_back(v1);
        }
    }

    LOG1 << "...done!";

    LOG1 << "Memory allocation...";

    G->start_construction(num_vertices, num_vertices * 6 * cut);
    std::vector<std::vector<EdgeWeight> > inblock(vertices_per_block);

    for (auto& i : inblock) {
        i.resize(vertices_per_block, 0);
    }

    LOG1 << "...done";
    LOG1 << "Setting in-block edges...";

    for (size_t b = 0; b < blocks; ++b) {
        LOG1 << "Generating edges for block " << b << "...";
        NodeID my_start = b * vertices_per_block;
        for (auto& i : inblock) {
            for (auto& j : i) {
                j = 0;
            }
        }

        for (size_t i = 0; i < inblock.size(); ++i) {
            for (size_t m = 0; m < 4 * cut + 1; ++m) {
                size_t r = random_functions::nextInt(0, vertices_per_block - 1);
                inblock[i][r] += 1;
                inblock[r][i] += 1;
            }
        }

        LOG1 << "Writing block " << b << "...";
        for (size_t i = 0; i < inblock.size(); ++i) {
            G->new_node();
            for (size_t j = 0; j < inblock[i].size(); ++j) {
                if (inblock[i][j] > 0)
                    G->new_edge(i + my_start, j + my_start, inblock[i][j]);
            }

            if (interblockedges[i + my_start].size() > 0) {
                std::map<size_t, size_t> intermap;
                std::for_each(interblockedges[i + my_start].begin(),
                              interblockedges[i + my_start].end(),
                              [&intermap](NodeID val) { intermap[val]++; });

                for (const auto& p : intermap) {
                    G->new_edge(i + my_start, p.first, p.second);
                }
            }
        }
    }

    G->finish_construction();

    std::string name = "gnp_" + std::to_string(vertices_per_block) + "_"
                       + std::to_string(blocks_per_dimension) + "_"
                       + std::to_string(cut) + "_" + std::to_string(seed);

    if (output_folder != "") {
        name = output_folder + "/" + name;
    }

    LOG1 << "Writing to file " << name << "...";
    graph_io::writeGraphWeighted(G, name);
    LOG1 << "Done!";
}
