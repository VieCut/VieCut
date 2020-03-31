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

#include <chrono>
#include <memory>
#include <thread>
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

    void sendProblem(problemPointer problem, size_t tgt) {
        LOG1 << mpi_rank << " sends problem to " << tgt;
        std::vector<uint64_t> data;
        data.emplace_back(problem->lower_bound);
        data.emplace_back(problem->upper_bound);
        data.emplace_back(problem->deleted_weight);
        data.emplace_back(problem->terminals.size());

        for (const auto& [t, o, b] : problem->terminals) {
            (void)b;
            data.emplace_back(t);
            data.emplace_back(o);
        }

        data.emplace_back(problem->mappings.size());
        for (size_t i = 0; i < problem->mappings.size(); ++i) {
            data.emplace_back(problem->mappings[i]->size());
            data.insert(data.end(), problem->mappings[i]->begin(),
                        problem->mappings[i]->end());
        }

        auto serial = problem->graph->serialize();
        data.emplace_back(serial.size());
        data.insert(data.end(), serial.begin(), serial.end());

        size_t datasize = data.size();
        MPI_Send(&datasize, 1, MPI_LONG, tgt, 1010, MPI_COMM_WORLD);
        MPI_Send(&data.front(), datasize, MPI_LONG, tgt, 1020, MPI_COMM_WORLD);
    }

    problemPointer recvProblem(size_t src) {
        size_t datasize = 0;
        MPI_Recv(&datasize, 1, MPI_LONG, src, 1010,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::vector<uint64_t> data(datasize);
        MPI_Recv(&data.front(), datasize, MPI_LONG, src, 1020,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        auto problem = std::make_shared<multicut_problem>();
        problem->lower_bound = data[0];
        problem->upper_bound = data[1];
        problem->deleted_weight = data[2];
        size_t num_terminals = data[3];
        size_t next_index = 4;
        for (size_t i = 0; i < num_terminals; ++i) {
            NodeID t = data[next_index++];
            NodeID o = data[next_index++];
            problem->terminals.emplace_back(t, o, true);
        }

        size_t num_mappings = data[next_index++];
        for (size_t i = 0; i < num_mappings; ++i) {
            size_t map_size_i = data[next_index++];
            problem->mappings.emplace_back(
                std::make_shared<std::vector<NodeID> >());
            problem->mappings[i]->insert(problem->mappings[i]->end(),
                                         &data[next_index],
                                         &data[next_index + map_size_i]);
            next_index += map_size_i;
        }

        size_t serialsize = data[next_index++];
        std::vector<uint64_t> serial;
        serial.insert(serial.end(), &data[next_index],
                      &data[next_index + serialsize]);

        problem->graph = mutable_graph::deserialize(serial);
        LOG << mpi_rank << " returns new problem";
        return problem;
    }

    std::optional<int> checkForReceiver() {
        int result;
        LOG0 << mpi_rank << "probing";
        MPI_Iprobe(MPI_ANY_SOURCE, 3000, MPI_COMM_WORLD,
                   &result, MPI_STATUS_IGNORE);
        if (result > 0) {
            MPI_Status st;
            MPI_Recv(&result, 1, MPI_LONG, MPI_ANY_SOURCE, 3000,
                     MPI_COMM_WORLD, &st);
            LOG << mpi_rank << " found probe from " << st.MPI_SOURCE;
            MessageStatus message = haveProblem;
            MPI_Send(&message, 1, MPI_INT, st.MPI_SOURCE, 3000, MPI_COMM_WORLD);
            return st.MPI_SOURCE;
        }
        LOG0 << mpi_rank << " no receiver found";
        return std::nullopt;
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
        MPI_Isend(&message, 1, MPI_INT, src, 3000, MPI_COMM_WORLD, &request);

        while (true) {
            int ms;
            MPI_Status stat;

            MPI_Recv(&ms, 1, MPI_LONG, MPI_ANY_SOURCE, 3000,
                     MPI_COMM_WORLD, &stat);

            if (ms == MessageStatus::allEmpty) {
                return true;
            }

            if (ms == MessageStatus::emptyAsWell) {
                if (stat.MPI_SOURCE != src) {
                    LOG1 << mpi_rank << "WRONG SENDER "
                         << stat.MPI_SOURCE << "EMPTY AS WELL !";
                    exit(1);
                }
                return false;
            }

            if (ms == MessageStatus::needProblem) {
                MessageStatus answer = emptyAsWell;
                MPI_Request rq;
                MPI_Isend(&answer, 1, MPI_INT, stat.MPI_SOURCE,
                          3000, MPI_COMM_WORLD, &rq);
            }

            if (ms == MessageStatus::haveProblem && stat.MPI_SOURCE == src) {
                if (stat.MPI_SOURCE != src) {
                    LOG1 << mpi_rank << "WRONG SENDER "
                         << stat.MPI_SOURCE << "HAS PROBLEM !";
                    exit(1);
                }
                LOG << mpi_rank << " heard answer from " << src;
                return src;
            }
        }
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
