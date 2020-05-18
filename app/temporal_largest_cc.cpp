/******************************************************************************
 * temporal_largest_cc.cpp
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2020 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include <ext/alloc_traits.h>

#include <memory>
#include <string>
#include <vector>

#include "algorithms/misc/strongly_connected_components.h"
#include "common/configuration.h"
#include "common/definitions.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tlx/logger.hpp"
#include "tlx/string/split.hpp"

int main(int argn, char** argv) {
    static constexpr bool debug = false;

    LOG << "Find and save largest CC!";

    tlx::CmdlineParser cmdl;
    std::string path;
    cmdl.add_param_string("graph", path, "path to graph file");
    configuration::getConfig()->graph_filename = path;

    if (!cmdl.process(argn, argv))
        return -1;

    auto [numV, tempEdges] = graph_io::readTemporalGraph(path);
    LOG << "Creating graph with " << numV << " vertices!";
    mutableGraphPtr G = std::make_shared<mutable_graph>();
    G->start_construction(numV);
    for (auto [a, b, c, d] : tempEdges) {
        G->new_edge_order(a, b, c);
    }
    G->finish_construction();
    graphAccessPtr GA = G->to_graph_access();

    strongly_connected_components scc;
    std::vector<int32_t> components(GA->number_of_nodes());
    auto ct = scc.strong_components(GA, &components);
    LOG << "count of connected components: " << ct;

    std::vector<uint64_t> compsizes(static_cast<uint64_t>(ct));
    for (int32_t component : components) {
        ++compsizes[component];
    }

    auto max_size = std::max_element(compsizes.begin(), compsizes.end());
    int max_comp = static_cast<int>(max_size - compsizes.begin());

    std::vector<NodeID> maxCompID(GA->number_of_nodes());

    NodeID mCCtr = 1;

    for (size_t i = 0; i < components.size(); ++i) {
        if (components[i] == max_comp) {
            maxCompID[i] = mCCtr;
            mCCtr++;
        }
    }

    std::string outpath = path + ".cc";
    std::ofstream f(outpath.c_str());

    for (auto [a, b, c, d] : tempEdges) {
        if (components[a] == max_comp && components[b] == max_comp) {
            NodeID newa = maxCompID[a];
            NodeID newb = maxCompID[b];
            f << newa << " " << newb << " " << c << " " << d << "\n";
        }
    }

    f.close();
}
