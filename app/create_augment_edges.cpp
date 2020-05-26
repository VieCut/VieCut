/******************************************************************************
 * create_augment_edges.cpp
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#include "algorithms/global_mincut/dynamic/dynamic_mincut.h"
#include "common/configuration.h"
#include "common/definitions.h"
#include "data_structure/graph_access.h"
#include "data_structure/mutable_graph.h"
#include "io/graph_io.h"
#include "tlx/cmdline_parser.hpp"
#include "tlx/logger.hpp"
#include "tools/random_functions.h"

int main(int argn, char** argv) {
    tlx::CmdlineParser cmdl;
    auto cfg = configuration::getConfig();
    cmdl.add_param_string("graph", cfg->graph_filename, "path to graph file");

    size_t insert_edges = 0;
    size_t delete_edges = 0;
#ifdef PARALLEL
    size_t procs = 1;
    cmdl.add_size_t('p', "proc", procs, "number of processes");
#endif
    cmdl.add_size_t('d', "delete", delete_edges, "number of edges to delete");
    cmdl.add_size_t('i', "insert", insert_edges, "number of edges to insert");
    cmdl.add_flag('v', "verbose", cfg->verbose, "more verbose logs");
    cmdl.add_size_t('r', "seed", cfg->seed, "random seed");

    if (!cmdl.process(argn, argv))
        return -1;

    cfg->save_cut = true;

    auto output = cfg->graph_filename + "." + std::to_string(insert_edges)
                  + "." + std::to_string(delete_edges);

    std::ofstream f(output, std::ofstream::trunc);

    if (delete_edges > insert_edges) {
        LOG1 << "Error: Trying to delete more edges that were inserted!";
        exit(1);
    }

    auto G = graph_io::readGraphWeighted<mutable_graph>(cfg->graph_filename);

    dynamic_mincut dynmc;
    dynmc.initialize(G);
    std::vector<std::pair<NodeID, NodeID> > insEdges;

    for (size_t ins = 0; ins < insert_edges; ++ins) {
        auto curr = dynmc.getCurrentCactus();
        auto original_graph = dynmc.getOriginalGraph();
        NodeID s = UNDEFINED_NODE;
        NodeID t = UNDEFINED_NODE;

        while (s == UNDEFINED_NODE) {
            NodeID r = random_functions::nextInt(0, curr->n() - 1);
            NodeID size = curr->numContainedVertices(r);
            if (size > 0) {
                NodeID r2 = random_functions::nextInt(0, size - 1);
                s = curr->containedVertices(r)[r2];
            }
        }

        while (t == UNDEFINED_NODE) {
            NodeID r = random_functions::nextInt(0, curr->n() - 1);
            NodeID size = curr->containedVertices(r).size();
            if (size > 0) {
                NodeID r2 = random_functions::nextInt(0, size - 1);
                NodeID preliminary_t = curr->containedVertices(r)[r2];
                if (s == preliminary_t)
                    continue;
                bool neighbors = false;
                for (EdgeID e : original_graph->edges_of(s)) {
                    // only create edge if it doesn't exist yet
                    NodeID tgt = original_graph->getEdgeTarget(s, e);
                    if (tgt == preliminary_t) {
                        neighbors = true;
                        break;
                    }
                }

                if (!neighbors) {
                    t = preliminary_t;
                }
            }
        }

        f << s << " " << t << " +1\n";
        dynmc.addEdge(s, t, 1);
        insEdges.emplace_back(s, t);
    }

    std::vector<bool> alreadyDeleted(insert_edges, false);
    size_t del = 0;
    while (del < delete_edges) {
        size_t r = random_functions::nextInt(0, insert_edges - 1);
        if (!alreadyDeleted[r]) {
            del++;
            alreadyDeleted[r] = true;
            auto [s, t] = insEdges[r];
            f << s << " " << t << " -1\n";
        }
    }
    f.close();
}
