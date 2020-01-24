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
#pragma once

#include <mpi.h>

#include <memory>
#include <vector>

#include "algorithms/multicut/multicut_problem.h"

class mpi_communication {
 public:
    static const bool debug = true;

    mpi_communication() : bestSolutionLocal(false),
                          best_solution(UNDEFINED_NODE) {
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        req.resize(mpi_size);
    }

    FlowType getGlobalBestSolution() {
        int result = 1;
        while (result > 0) {
            MPI_Status status;
            MPI_Iprobe(MPI_ANY_SOURCE, 15123, MPI_COMM_WORLD, &result, &status);
            if (result > 0) {
                FlowType recv;
                MPI_Recv(&recv, 1, MPI_LONG, status.MPI_SOURCE,
                         15123, MPI_COMM_WORLD, NULL);
                if (best_solution > recv) {
                    best_solution = recv;
                    bestSolutionLocal = false;
                    LOG1 << "Received better solution " << best_solution;
                }
            }
        }
        return best_solution;
    }

    void broadcastImprovedSolution(FlowType solution) {
        if (solution < best_solution) {
            best_solution = solution;
            bestSolutionLocal = true;
        }

        for (int i = 0; i < mpi_size; ++i) {
            if (i != mpi_rank) {
                LOG << mpi_rank << "sending " << best_solution << " to " << i;
                MPI_Isend(&best_solution, 1, MPI_LONG,
                          i, 15123, MPI_COMM_WORLD, &req[i]);
            }
        }
    }



    void sendProblem(std::shared_ptr<multicut_problem> problem, size_t tgt) {
        LOG << mpi_rank << " send problem to " << tgt;
        MPI_Send(&problem->lower_bound, 1, MPI_LONG, tgt, 1010, MPI_COMM_WORLD);
        MPI_Send(&problem->upper_bound, 1, MPI_LONG, tgt, 1020, MPI_COMM_WORLD);
        MPI_Send(&problem->deleted_weight, 1, MPI_LONG,
                 tgt, 1030, MPI_COMM_WORLD);

        // terminals
        size_t termsize = problem->terminals.size();
        std::vector<NodeID> termpos;
        std::vector<NodeID> origid;
        for (const auto& [t, o, b] : problem->terminals) {
            termpos.emplace_back(t);
            origid.emplace_back(o);
        }

        MPI_Send(&termsize, 1, MPI_LONG,
                 tgt, 1040, MPI_COMM_WORLD);
        MPI_Send(&termpos.front(), termpos.size(),
                 MPI_INT, tgt, 1050, MPI_COMM_WORLD);
        MPI_Send(&origid.front(), origid.size(),
                 MPI_INT, tgt, 1060, MPI_COMM_WORLD);

        // mappings
        size_t mapsize = problem->mappings.size();
        MPI_Send(&mapsize, 1, MPI_LONG, tgt, 1070, MPI_COMM_WORLD);
        for (size_t i = 0; i < problem->mappings.size(); ++i) {
            size_t mapisize = problem->mappings[i]->size();
            MPI_Send(&mapisize, 1, MPI_LONG, tgt, 1080, MPI_COMM_WORLD);
            MPI_Send(&problem->mappings[i]->front(),
                     problem->mappings[i]->size(),
                     MPI_INT, tgt, 1090, MPI_COMM_WORLD);
        }

        // graph
        auto serial = problem->graph->serialize();
        size_t serialsize = serial.size();
        MPI_Send(&serialsize, 1, MPI_LONG, tgt, 1100, MPI_COMM_WORLD);
        MPI_Send(&serial.front(), serial.size(), MPI_LONG,
                 tgt, 1110, MPI_COMM_WORLD);
    }

    bool checkForReceiver(std::shared_ptr<multicut_problem> problem) {
        int dest = (mpi_rank == mpi_size - 1) ? 0 : mpi_rank + 1;
        int result;
        MPI_Iprobe(dest, 2000, MPI_COMM_WORLD, &result, MPI_STATUS_IGNORE);
        if (result > 0) {
            LOG << mpi_rank << " found probe from " << dest;
            MessageStatus message = haveProblem;
            MPI_Send(&message, 1, MPI_INT, dest, 3000, MPI_COMM_WORLD);
            sendProblem(problem, dest);
            return true;
        }
        return false;
    }

    std::shared_ptr<multicut_problem> waitForProblem() {
        MessageStatus message = needProblem;
        int src = (mpi_rank == 0) ? mpi_size - 1 : mpi_rank - 1;
        LOG << mpi_rank << " waiting for " << src;
        MPI_Request request;
        MPI_Isend(&message, 1, MPI_INT, src, 2000, MPI_COMM_WORLD, &request);
        MPI_Probe(src, 3000, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        LOG << mpi_rank << " heard answer from " << src;
        return recvProblem();
    }

    std::shared_ptr<multicut_problem> recvProblem() {
        LOG1 << mpi_rank << " receives problem";
        FlowType lower_bound;
        FlowType upper_bound;
        EdgeWeight deleted_wgt;
        MPI_Status status;
        MPI_Recv(&lower_bound, 1, MPI_LONG, MPI_ANY_SOURCE,
                 1010, MPI_COMM_WORLD, &status);
        MPI_Recv(&upper_bound, 1, MPI_LONG, status.MPI_SOURCE,
                 1020, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&deleted_wgt, 1, MPI_LONG, status.MPI_SOURCE,
                 1030, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        LOG1 << mpi_rank << " currently receiving from " << status.MPI_SOURCE;

        // terminals
        size_t termsize = 0;
        MPI_Recv(&termsize, 1, MPI_LONG, status.MPI_SOURCE,
                 1040, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::vector<NodeID> termids;
        std::vector<NodeID> origids;
        termids.resize(termsize);
        origids.resize(termsize);
        MPI_Recv(&termids.front(), termsize, MPI_INT, status.MPI_SOURCE,
                 1050, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&origids.front(), termsize, MPI_INT, status.MPI_SOURCE,
                 1060, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::vector<terminal> terminals;
        for (size_t i = 0; i < termids.size(); ++i) {
            terminals.emplace_back(termids[i], origids[i], true);
        }

        // mappings
        size_t mappingsize = 0;
        MPI_Recv(&mappingsize, 1, MPI_LONG, status.MPI_SOURCE,
                 1070, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::vector<std::shared_ptr<std::vector<NodeID> > > mappings;
        for (size_t i = 0; i < mappingsize; ++i) {
            auto map = std::make_shared<std::vector<NodeID> >();
            size_t currmapsize = 0;
            MPI_Recv(&currmapsize, 1, MPI_LONG, status.MPI_SOURCE,
                     1080, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            map->resize(currmapsize);
            MPI_Recv(&map->front(), currmapsize, MPI_INT, status.MPI_SOURCE,
                     1090, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            mappings.emplace_back(map);
        }

        size_t serialsize = 0;
        MPI_Recv(&serialsize, 1, MPI_LONG, status.MPI_SOURCE,
                 1100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::vector<uint64_t> serial(serialsize);
        MPI_Recv(&serial.front(), serialsize, MPI_LONG, status.MPI_SOURCE,
                 1110, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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

    bool bestSolutionIsLocal() {
        return bestSolutionLocal;
    }

 private:
    int mpi_size;
    int mpi_rank;
    std::vector<MPI_Request> req;
    bool bestSolutionLocal;
    FlowType best_solution;

    enum MessageStatus { needProblem, haveProblem };
};
