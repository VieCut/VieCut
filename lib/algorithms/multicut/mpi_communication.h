/******************************************************************************
 * mpi_communication.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2020 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/
#include <mpi.h>

#include <vector>

#include "multicut_problem.h"

class mpi_communication {

 public:
    mpi_communication() {}

    static void sendProblem(
        std::shared_ptr<multicut_problem> problem, size_t tgt) {
        MPI_Send(&problem->lower_bound, 1, MPI_LONG, tgt, 523, MPI_COMM_WORLD);
        MPI_Send(&problem->upper_bound, 1, MPI_LONG, tgt, 0, MPI_COMM_WORLD);
        MPI_Send(&problem->deleted_weight, 1, MPI_LONG, tgt, 0, MPI_COMM_WORLD);

        // terminals
        size_t termsize = problem->terminals.size();
        std::vector<NodeID> termpos;
        std::vector<NodeID> origid;
        for (const auto& [t, o, b] : problem->terminals) {
            termpos.emplace_back(t);
            origid.emplace_back(o);
        }

        MPI_Send(&termsize, 1, MPI_LONG,
                 tgt, 0, MPI_COMM_WORLD);
        MPI_Send(&termpos.front(), termpos.size(),
                 MPI_INT, tgt, 0, MPI_COMM_WORLD);
        MPI_Send(&origid.front(), origid.size(),
                 MPI_INT, tgt, 0, MPI_COMM_WORLD);

        // mappings
        size_t mapsize = problem->mappings.size();
        MPI_Send(&mapsize, 1, MPI_LONG, tgt, 0, MPI_COMM_WORLD);
        for (size_t i = 0; i < problem->mappings.size(); ++i) {
            size_t mapisize = problem->mappings[i]->size();
            MPI_Send(&mapisize, 1, MPI_LONG, tgt, 0, MPI_COMM_WORLD);
            MPI_Send(&problem->mappings[i]->front(),
                     problem->mappings[i]->size(),
                     MPI_INT, tgt, 0, MPI_COMM_WORLD);
        }

        // graph
        auto serial = problem->graph->serialize();
        size_t serialsize = serial.size();
        MPI_Send(&serialsize, 1, MPI_LONG, tgt, 0, MPI_COMM_WORLD);
        MPI_Send(&serial.front(), serial.size(), MPI_LONG,
                 tgt, 0, MPI_COMM_WORLD);
        return;
    }

    static std::shared_ptr<multicut_problem> recvProblem() {
        FlowType lower_bound;
        FlowType upper_bound;
        EdgeWeight deleted_wgt;
        MPI_Status status;
        MPI_Recv(&lower_bound, 1, MPI_LONG, MPI_ANY_SOURCE,
                 523, MPI_COMM_WORLD, &status);
        MPI_Recv(&upper_bound, 1, MPI_LONG, status.MPI_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
        MPI_Recv(&deleted_wgt, 1, MPI_LONG, status.MPI_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, NULL);

        // terminals
        size_t termsize = 0;
        MPI_Recv(&termsize, 1, MPI_LONG, status.MPI_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
        std::vector<NodeID> termids;
        std::vector<NodeID> origids;
        termids.resize(termsize);
        origids.resize(termsize);
        MPI_Recv(&termids.front(), termsize, MPI_INT, status.MPI_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
        MPI_Recv(&origids.front(), termsize, MPI_INT, status.MPI_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
        std::vector<terminal> terminals;
        for (size_t i = 0; i < termids.size(); ++i) {
            terminals.emplace_back(termids[i], origids[i], true);
        }

        // mappings
        size_t mappingsize = 0;
        MPI_Recv(&mappingsize, 1, MPI_LONG, status.MPI_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
        std::vector<std::shared_ptr<std::vector<NodeID> > > mappings;
        for (size_t i = 0; i < mappingsize; ++i) {
            auto map = std::make_shared<std::vector<NodeID> >();
            size_t currmapsize = 0;
            MPI_Recv(&currmapsize, 1, MPI_LONG, status.MPI_SOURCE,
                     MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
            map->resize(currmapsize);
            MPI_Recv(&map->front(), currmapsize, MPI_INT, status.MPI_SOURCE,
                     MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
            mappings.emplace_back(map);
        }

        size_t serialsize = 0;
        MPI_Recv(&serialsize, 1, MPI_LONG, status.MPI_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
        std::vector<uint64_t> serial(serialsize);
        MPI_Recv(&serial.front(), serialsize, MPI_LONG, status.MPI_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
        auto G = mutable_graph::deserialize(serial);

        auto problem = std::make_shared<multicut_problem>();
        problem->graph = G;
        problem->mappings = mappings;
        problem->terminals = terminals;
        problem->lower_bound = lower_bound;
        problem->upper_bound = upper_bound;
        problem->deleted_weight = deleted_wgt;
        return problem;
    }
};