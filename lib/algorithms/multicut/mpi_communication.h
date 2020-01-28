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
#include "common/definitions.h"

class mpi_communication {
 public:
    static const bool debug = false;

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
                    LOG1 << mpi_rank << " received better solution "
                         << best_solution << " from " << status.MPI_SOURCE;
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
                LOG0 << mpi_rank << " sending " << best_solution << " to " << i;
                MPI_Isend(&best_solution, 1, MPI_LONG,
                          i, 15123, MPI_COMM_WORLD, &req[i]);
            }
        }
    }

    void sendProblem(std::shared_ptr<multicut_problem> problem, size_t tgt) {
        LOG1 << mpi_rank << " sending problem to " << tgt;
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
        int result;
        LOG << mpi_rank << "probing";
        MPI_Iprobe(MPI_ANY_SOURCE, 2000, MPI_COMM_WORLD,
                   &result, MPI_STATUS_IGNORE);
        if (result > 0) {
            MPI_Status st;
            MPI_Recv(&result, 1, MPI_LONG, MPI_ANY_SOURCE, 2000,
                     MPI_COMM_WORLD, &st);
            LOG << mpi_rank << " found probe from " << st.MPI_SOURCE;
            MessageStatus message = haveProblem;
            MPI_Send(&message, 1, MPI_INT, st.MPI_SOURCE, 3000, MPI_COMM_WORLD);
            sendProblem(problem, st.MPI_SOURCE);
            return true;
        }
        LOG << mpi_rank << " no receiver found";
        return false;
    }

    std::variant<int, bool> waitForProblem() {
        if (mpi_size == 1) {
            return true;
        }

        MessageStatus message = needProblem;
        int src = 0;

        do {
            src = random_functions::nextInt(0, mpi_size - 1);
        } while (src == mpi_rank);

        LOG << mpi_rank << " waiting for " << src;
        MPI_Request request;
        MPI_Isend(&message, 1, MPI_INT, src, 2000, MPI_COMM_WORLD, &request);

        int incoming = 0;

        while (incoming == 0) {
            int result = 1;
            while (result > 0) {
                int ms;
                MPI_Status stat;
                MPI_Iprobe(MPI_ANY_SOURCE, 2000,
                           MPI_COMM_WORLD, &result, &stat);
                if (result > 0) {
                    MPI_Recv(&ms, 1, MPI_LONG, stat.MPI_SOURCE, 2000,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    LOG << mpi_rank << " found probe from "
                        << stat.MPI_SOURCE << " but is empty";
                    MessageStatus answer = emptyAsWell;
                    MPI_Request rq;
                    MPI_Isend(&answer, 1, MPI_INT, stat.MPI_SOURCE,
                              3000, MPI_COMM_WORLD, &rq);
                }
            }

            // LOG << mpi_rank << " PROBING ON 3000";
            MPI_Iprobe(src, 3000, MPI_COMM_WORLD, &incoming, MPI_STATUS_IGNORE);
            /*if (!incoming) {
                LOG << mpi_rank << " sleeping for 1 second!";
                std::this_thread::sleep_for(std::chrono::seconds(1));
            }*/
        }
        MessageStatus re;
        MPI_Recv(&re, 1, MPI_INT, src, 3000, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        LOG << mpi_rank << " received " << re;

        if (re == MessageStatus::allEmpty) {
            return true;
        }

        if (re == MessageStatus::emptyAsWell) {
            return false;
        }

        if (re == MessageStatus::needProblem) {
            LOG1 << "Error: Invalid MPI message!";
            exit(1);
        }

        LOG << mpi_rank << " heard answer from " << src;
        return src;
    }

    std::shared_ptr<multicut_problem> recvProblem(int src) {
        LOG << mpi_rank << " receives problem";
        FlowType lower_bound;
        FlowType upper_bound;
        EdgeWeight deleted_wgt;
        MPI_Status status;
        MPI_Recv(&lower_bound, 1, MPI_LONG, src, 1010, MPI_COMM_WORLD, &status);
        MPI_Recv(&upper_bound, 1, MPI_LONG, src,
                 1020, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&deleted_wgt, 1, MPI_LONG, src,
                 1030, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        LOG << mpi_rank << " currently receiving from " << status.MPI_SOURCE;

        // terminals
        size_t termsize = 0;
        MPI_Recv(&termsize, 1, MPI_LONG, src,
                 1040, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::vector<NodeID> termids;
        std::vector<NodeID> origids;
        termids.resize(termsize);
        origids.resize(termsize);
        MPI_Recv(&termids.front(), termsize, MPI_INT, src,
                 1050, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&origids.front(), termsize, MPI_INT, src,
                 1060, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::vector<terminal> terminals;
        for (size_t i = 0; i < termids.size(); ++i) {
            terminals.emplace_back(termids[i], origids[i], true);
        }

        // mappings
        size_t mappingsize = 0;
        MPI_Recv(&mappingsize, 1, MPI_LONG, src,
                 1070, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::vector<std::shared_ptr<std::vector<NodeID> > > mappings;
        for (size_t i = 0; i < mappingsize; ++i) {
            auto map = std::make_shared<std::vector<NodeID> >();
            size_t currmapsize = 0;
            MPI_Recv(&currmapsize, 1, MPI_LONG, src,
                     1080, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            map->resize(currmapsize);
            MPI_Recv(&map->front(), currmapsize, MPI_INT, src,
                     1090, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            mappings.emplace_back(map);
        }

        size_t serialsize = 0;
        MPI_Recv(&serialsize, 1, MPI_LONG, src,
                 1100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::vector<uint64_t> serial;
        serial.resize(serialsize);
        MPI_Recv(&serial.front(), serialsize, MPI_LONG, src,
                 1110, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        auto G = mutable_graph::deserialize(serial);

        auto problem = std::make_shared<multicut_problem>();
        problem->graph = G;
        problem->mappings = mappings;
        problem->terminals = terminals;
        problem->lower_bound = lower_bound;
        problem->upper_bound = upper_bound;
        problem->deleted_weight = deleted_wgt;
        LOG << mpi_rank << " returns new problem";
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
};
