/******************************************************************************
 * branch_multicut.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2018-2020 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <mpi.h>

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdio>
#include <iterator>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/global_mincut/noi_minimum_cut.h"
#include "algorithms/misc/graph_algorithms.h"
#include "algorithms/multicut/edge_selection.h"
#include "algorithms/multicut/graph_contraction.h"
#include "algorithms/multicut/kernelization_criteria.h"
#include "algorithms/multicut/local_search.h"
#include "algorithms/multicut/maximum_flow.h"
#include "algorithms/multicut/measurements.h"
#include "algorithms/multicut/mpi_communication.h"
#include "algorithms/multicut/multicut_problem.h"
#include "algorithms/multicut/problem_management.h"
#include "coarsening/contract_graph.h"
#include "common/configuration.h"
#include "data_structure/union_find.h"
#include "io/graph_io.h"
#include "tools/timer.h"
#include "tools/vector.h"

#ifdef USE_GUROBI
#include "algorithms/multicut/ilp_model.h"
#endif

#ifdef USE_TCMALLOC
#include "gperftools/malloc_extension.h"
#endif

class branch_multicut {
 public:
#ifndef NDEBUG
    static const bool debug = true;
#else
    static const bool debug = false;
#endif
    static const bool testing = true;

    branch_multicut(const mutable_graph& original_graph,
                    std::vector<NodeID> original_terminals,
                    std::vector<bool> fixed_vertex)
        : original_graph(original_graph),
          original_terminals(original_terminals),
          fixed_vertex(fixed_vertex),
          total_time(),
          q_mutex(configuration::getConfig()->threads),
          num_threads(configuration::getConfig()->threads),
          branch_invalid(configuration::getConfig()->threads, 0),
          kc(original_terminals),
          mf(original_terminals),
          pm(this->original_graph, this->original_terminals,
             this->fixed_vertex),
          msm(this->original_graph, this->original_terminals),
          last_sent_flow(UNDEFINED_FLOW),
          log_timer(0),
          finished(false) {
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        mpi_num_done = 0;
        mpi_done.resize(mpi_size, 0);
        sent_done = false;
    }

    ~branch_multicut() { }

    std::pair<std::vector<NodeID>, size_t> find_multiterminal_cut(
        problemPointer problem) {
        std::vector<NodeID> sol;
        size_t numTerminals = problem->terminals.size();
        if (mpi_rank == 0) {
            mf.maximumIsolatingFlow(problem, 0, /* parallel */ true);
            sol = msm.getSolution(problem);
            pm.addProblem(problem, 0, false);
            pm.updateBound(problem->upper_bound);
        }

        std::vector<std::thread> threads;
        for (size_t i = 0; i < num_threads; ++i) {
            threads.emplace_back(
                std::thread(&branch_multicut::pollWork, this, i));
            if (!configuration::getConfig()->disable_cpu_affinity) {
                cpu_set_t cpuset;
                CPU_ZERO(&cpuset);
                CPU_SET(i, &cpuset);
                pthread_setaffinity_np(
                    threads[i].native_handle(), sizeof(cpu_set_t), &cpuset);
            }
        }

        if (mpi_rank == 0) {
            updateBestSolution(&sol, numTerminals);
        }

        for (auto& t : threads) {
            pm.notifyAllThreads();
            t.join();
        }

        MPI_Request all_done;
        MPI_Ibarrier(MPI_COMM_WORLD, &all_done);

        int finished = 0;
        while (finished == 0) {
            MPI_Test(&all_done, &finished, MPI_STATUS_IGNORE);
            int ms;
            MPI_Status status;

            MPI_Iprobe(MPI_ANY_SOURCE, 3000, MPI_COMM_WORLD, &ms, &status);

            if (ms > 0) {
                int re;
                MPI_Recv(&re, 1, MPI_LONG, status.MPI_SOURCE, 3000,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                MessageStatus answer = allEmpty;
                MPI_Request rq;
                MPI_Isend(&answer, 1, MPI_INT, status.MPI_SOURCE,
                          3000, MPI_COMM_WORLD, &rq);
            }
        }

        FlowType total_weight = pm.bestCut();

        FlowType local_weight = total_weight;
        MPI_Allreduce(&local_weight, &total_weight, 1, MPI_LONG,
                      MPI_MIN, MPI_COMM_WORLD);

        size_t local_bcast_id = mpi_size;
        if (mpic.bestSolutionIsLocal() && local_weight == total_weight) {
            // find processor which has best solution to broadcast
            local_bcast_id = mpi_rank;
        }

        size_t global_bcast_id = local_bcast_id;
        MPI_Allreduce(&local_bcast_id, &global_bcast_id, 1, MPI_LONG,
                      MPI_MIN, MPI_COMM_WORLD);

        std::vector<NodeID> best_solution;
        if (mpi_rank == static_cast<int>(global_bcast_id)) {
            best_solution = pm.getBestSolution();
        } else {
            best_solution = std::vector<NodeID>(original_graph.n());
        }

        size_t bsize = best_solution.size();

        MPI_Bcast(&best_solution.front(), bsize, MPI_INT,
                  global_bcast_id, MPI_COMM_WORLD);

        total_weight = msm.flowValue(false, best_solution);
        VIECUT_ASSERT_LEQ(total_weight, pm.bestCut());

        return std::make_pair(best_solution, total_weight);
    }

 private:
    void pollWork(size_t thread_id) {
        bool im_idle = false;
        while (!pm.finished()) {
            if (last_sent_flow > pm.bestCut()) {
                last_sent_flow = pm.bestCut();
                mpic.broadcastImprovedSolution(last_sent_flow);
            }
            pm.prepareQueue(thread_id);
            if (!pm.queueEmpty(thread_id) || pm.haveASendProblem()) {
                if (thread_id == 0) {
                    pm.updateBound(mpic.getGlobalBestSolution());
                }

                std::optional<int> sending = std::nullopt;
                if (mpi_size > 1 && thread_id == 0 && pm.numProblems() > 1) {
                    sending = mpic.checkForReceiver();
                }

                auto problem = pm.pullProblem(thread_id, sending.has_value());
                if (!problem.has_value()) {
                    continue;
                }

                if (sending.has_value()) {
                    // forget this problem if it was sent to another worker
                    mpic.sendProblem(problem.value(), sending.value());
                } else {
                    solveProblem(problem.value(), thread_id);
                }
            } else {
                if (!im_idle) {
                    pm.incrementIdleThreads();
                    im_idle = true;
                }

                if (pm.allThreadsIdle() && pm.allEmpty()) {
                    if (mpi_mutex.try_lock()) {
                        auto src = mpic.waitForProblem();
                        if (std::holds_alternative<int>(src)) {
                            if (sent_done) {
                                int done = 0;
                                sent_done = false;
                                MPI_Request rq;
                                MPI_Isend(&done, 1, MPI_INT, mpi_size - 1,
                                          4000, MPI_COMM_WORLD, &rq);
                            }
                            auto p = mpic.recvProblem(std::get<int>(src));
                            pm.addProblem(p, thread_id, true);
                        } else {
                            if (std::get<bool>(src)) {
                                pm.setFinish();
                                pm.notifyAllThreads();
                                return;
                            }

                            if (mpi_rank == mpi_size - 1) {
                                int result = 1;
                                while (result > 0) {
                                    MPI_Status stat;
                                    MPI_Iprobe(MPI_ANY_SOURCE, 4000,
                                               MPI_COMM_WORLD, &result, &stat);
                                    if (result > 0) {
                                        int done;
                                        int source = stat.MPI_SOURCE;
                                        MPI_Recv(&done, 1, MPI_INT, source,
                                                 4000, MPI_COMM_WORLD,
                                                 MPI_STATUS_IGNORE);
                                        int done_before = mpi_done[source];
                                        mpi_done[source] = done;
                                        mpi_num_done += (done - done_before);
                                    }
                                }

                                if (mpi_num_done == mpi_size - 1) {
                                    for (int i = 0; i < mpi_size - 1; ++i) {
                                        MessageStatus fnl =
                                            MessageStatus::allEmpty;
                                        MPI_Request fr;
                                        MPI_Isend(&fnl, 1, MPI_INT, i, 3000,
                                                  MPI_COMM_WORLD, &fr);
                                        pm.setFinish();
                                        pm.notifyAllThreads();
                                        return;
                                    }
                                }
                            } else {
                                if (!sent_done) {
                                    int done = 1;
                                    sent_done = true;
                                    MPI_Request rq;
                                    MPI_Isend(&done, 1, MPI_INT, mpi_size - 1,
                                              4000, MPI_COMM_WORLD, &rq);
                                }
                            }
                        }
                        mpi_mutex.unlock();
                    }
                }

                pm.waitForProblem(thread_id);

                if (im_idle) {
                    pm.decrementIdleThreads();
                    im_idle = false;
                }
            }
        }
    }

    bool outOfMemory() {
#ifdef USE_TCMALLOC
        uint64_t heapsize = 0;
        MallocExtension::instance()->GetNumericProperty(
            "generic.heap_size", &heapsize);

        uint64_t max_size = 32UL * 1024UL * 1024UL * 1024UL;
        if (heapsize > max_size) {
            LOG1 << "Memoryout!";
            finished = true;
            return true;
        }
#endif
        return false;
    }

    void solveProblem(problemPointer problem,
                      size_t thread_id) {
        if (finished)
            return;

        auto c = configuration::getConfig();
        if (problem == NULL) {
            LOG1 << "ERROR: Problem is NULL. This should not happen!";
            exit(1);
        }

        if (!pm.checkProblem(problem)) {
            return;
        }

        if (outOfMemory())
            return;

        if (total_time.elapsed() > configuration::getConfig()->timeoutSeconds) {
            LOG1 << "Timeout!";
            finished = true;
            return;
        }

        if (total_time.elapsed() > log_timer) {
            double logs_per_second = 2;
            double time_added = 1.0 / logs_per_second;

            log_timer = time_added + log_timer;

            LOGC(testing) << "Rank " << mpi_rank << " - "
                          << "After " << total_time.elapsed()
                          << " - terminals:" << problem->terminals.size()
                          << " vertices:" << problem->graph->n()
                          << " deleted:" << problem->deleted_weight
                          << " lower:" << problem->lower_bound
                          << " upper:" << problem->upper_bound
                          << " global_upper:" << pm.bestCut()
                          << " queue.size:" << pm.numProblems();
        }

        if (problem->graph->n() * 2 < problem->graph->getOriginalNodes()) {
            auto map = std::make_shared<std::vector<NodeID> >();
            for (size_t i = 0; i < problem->graph->getOriginalNodes(); ++i) {
                map->emplace_back(problem->graph->getCurrentPosition(i));
            }
            problem->mappings.emplace_back(map);
            problem->graph = problem->graph->simplify();
        }
        graph_contraction::setTerminals(problem, original_terminals);
        nonBranchingContraction(problem);
        if (outOfMemory())
            return;

        if (problem->deleted_weight > static_cast<EdgeWeight>(pm.bestCut())) {
            return;
        }

        bool branchHere = true;
#ifdef USE_GUROBI
        size_t ilpLimit = 50000;
        branchHere = problem->graph->m() > ilpLimit || (!c->use_ilp);
#endif

        auto path = c->first_branch_path;
        if (path != "") {
            multicut_problem::writeGraph(problem, path);
            exit(1);
        }

        if (outOfMemory())
            return;

        if (branchHere) {
            branchOnEdge(problem, thread_id);
        } else {
            solve_with_ilp(problem, thread_id);
        }
    }

    void contractHeaviestTerminal(problemPointer problem) {
        auto c = configuration::getConfig();
        NodeID heaviest_t = 0;
        EdgeWeight heaviest_weight = 0;
        for (auto t : problem->terminals) {
            NodeID pos = t.position;
            auto deg = problem->graph->getWeightedNodeDegree(pos);
            if (deg > heaviest_weight) {
                heaviest_weight = deg;
                heaviest_t = t.position;
            }
        }
        std::vector<bool> found(problem->graph->n(), false);
        std::unordered_set<NodeID> contractSet;
        std::unordered_set<NodeID> isOtherTerminal;
        std::queue<std::pair<NodeID, size_t> > nodeAndDistance;
        nodeAndDistance.emplace(heaviest_t, 0);
        for (auto t : problem->terminals) {
            if (t.position != heaviest_t) {
                isOtherTerminal.insert(t.position);
            }
        }

        while (!nodeAndDistance.empty()) {
            auto [n, d] = nodeAndDistance.front();
            nodeAndDistance.pop();
            bool neighboringOtherTerminal = false;
            for (EdgeID e : problem->graph->edges_of(n)) {
                NodeID nbr = problem->graph->getEdgeTarget(n, e);
                if (isOtherTerminal.count(nbr) > 0) {
                    neighboringOtherTerminal = true;
                    break;
                }

                if (d < c->contractionDepthAroundTerminal && !found[nbr]) {
                    found[nbr] = true;
                    nodeAndDistance.emplace(nbr, d + 1);
                }
            }
            if (!neighboringOtherTerminal) {
                contractSet.insert(n);
            }
        }
        if (contractSet.size() > 0) {
            problem->graph->contractVertexSet(contractSet);
        }
    }

    void deleteLowestDegTerminals(problemPointer problem) {
        size_t lowestTerminals =
            std::ceil(
                static_cast<double>(problem->terminals.size()) *
                configuration::getConfig()->removeTerminalsBeforeBranch);

        for (size_t i = 0; i < lowestTerminals; ++i) {
            NodeID lightest_t = 0;
            NodeID lightest_oid = 0;
            EdgeWeight lightest_weight = UNDEFINED_FLOW;

            for (auto t : problem->terminals) {
                NodeID pos = t.position;
                auto deg = problem->graph->getWeightedNodeDegree(pos);
                if (deg < lightest_weight && deg > 0) {
                    lightest_weight = deg;
                    lightest_t = t.position;
                    lightest_oid = t.original_id;
                }
            }

            for (size_t i = 0; i < original_terminals.size(); ++i) {
                if (i != lightest_oid) {
                    problem->addFinishedPair(
                        i, lightest_oid, original_terminals.size());
                }
            }

            std::unordered_set<NodeID> contractIntoTerminal;
            contractIntoTerminal.insert(lightest_t);
            std::vector<NodeID> contractVertices(problem->graph->n(),
                                                 UNDEFINED_NODE);

            if (pm.isBestSolutionInitialized()) {
                auto best_solution = pm.getBestSolution();
                for (size_t i = 0; i < best_solution.size(); ++i) {
                    NodeID map = problem->mapped(i);
                    NodeID cp = problem->graph->getCurrentPosition(map);
                    if (best_solution[i] != lightest_oid) {
                        contractVertices[cp] = best_solution[i];
                    } else {
                        if (contractVertices[cp] != UNDEFINED_NODE) {
                            problem->removeFinishedPair(
                                lightest_oid, contractVertices[cp],
                                original_terminals.size());
                        }
                    }
                }

                for (size_t i = 0; i < contractVertices.size(); ++i) {
                    if (contractVertices[i] == UNDEFINED_NODE) {
                        contractIntoTerminal.insert(i);
                    }
                }
                NodeID invtx = problem->graph->containedVertices(lightest_t)[0];
                problem->graph->contractVertexSet(contractIntoTerminal);
                lightest_t = problem->graph->getCurrentPosition(invtx);
                problem->graph->setPartitionIndex(lightest_t, lightest_oid);
            }

            graph_contraction::deleteTermEdges(problem, original_terminals);
            EdgeID e1 = problem->graph->get_first_invalid_edge(lightest_t);

            for (EdgeID e = e1; e-- != 0; ) {
                auto wgt = problem->graph->getEdgeWeight(lightest_t, e);
                problem->graph->deleteEdge(lightest_t, e);
                problem->deleted_weight += wgt;
            }
            graph_contraction::setTerminals(problem, original_terminals);
        }
    }

    void updateBestSolution(std::vector<NodeID>* sol, size_t numTerminals) {
        auto s = pm.findBestSolution(sol, numTerminals);
        if (s.has_value()) {
            mpic.broadcastImprovedSolution(s.value());
        }
    }

    void printBoundaries() {
        std::vector<NodeID> boundary;
        for (NodeID n : original_graph.nodes()) {
            for (EdgeID e : original_graph.edges_of(n)) {
                NodeID tgt = original_graph.getEdgeTarget(n, e);
                if (original_graph.getPartitionIndex(n)
                    != original_graph.getPartitionIndex(tgt)) {
                    boundary.emplace_back(n);
                    break;
                }
            }
        }

        graph_io gio;
        std::string sol_str = "sol_" + std::to_string(print_index);
        LOG1 << "PRINTING TO " << sol_str;
        gio.writeVector(boundary, sol_str);
        print_index++;
    }

    void branchOnEdge(problemPointer problem,
                      size_t thread_id) {
        if (configuration::getConfig()->inexact) {
            deleteLowestDegTerminals(problem);
            contractHeaviestTerminal(problem);
        }

        graph_contraction::deleteTermEdges(problem, original_terminals);
        if (!problem->graph->number_of_edges()) {
            problem->upper_bound = problem->deleted_weight;
            auto sol = msm.getSolution(problem);
            updateBestSolution(&sol, problem->terminals.size());
            return;
        }

        if (outOfMemory())
            return;
        pm.branch(problem, thread_id);

        if (outOfMemory())
            return;
    }

    void nonBranchingContraction(problemPointer problem) {
        auto pe = kc.kernelization(problem, pm.bestCut(),
                                   pm.numProblems() == 0
                                   && configuration::getConfig()->threads > 1);
        if (pe.has_value()) {
            problem->priority_edge = *pe;
        }
    }

#ifdef USE_GUROBI
    void solve_with_ilp(problemPointer problem,
                        size_t thread_id) {
        graph_contraction::deleteTermEdges(problem, original_terminals);
        std::vector<NodeID> presets(problem->graph->n(),
                                    original_terminals.size());

        for (size_t i = 0; i < original_terminals.size(); ++i) {
            NodeID map = problem->mapped(original_terminals[i]);
            NodeID pos = problem->graph->getCurrentPosition(map);
            presets[pos] = i;
        }

        LOG1 << "start ILP with deleted " << problem->deleted_weight;

        auto [result, wgt, reIntroduce] =
            ilp.computeIlp(problem, presets, original_terminals.size(),
                           pm.numProblems() == 0, thread_id);
        problem->upper_bound = problem->deleted_weight + wgt;

        if (pm.runLocalSearch(problem)) {
            for (const auto& n : problem->graph->nodes()) {
                problem->graph->setPartitionIndex(n, result[n]);
            }
            auto sol = msm.getSolution(problem);
            updateBestSolution(&sol, problem->terminals.size());
        }

        if (reIntroduce) {
            branchOnEdge(problem, thread_id);
        }
    }
#else
    void solve_with_ilp(problemPointer, size_t) {
        LOG1 << "Error: Code not compiled with option -DUSE_GUROBI,"
             << " but called ILP solver. Exiting!";
        exit(1);
    }
#endif

    mutable_graph original_graph;
    std::vector<NodeID> original_terminals;
    std::vector<bool> fixed_vertex;
    timer total_time;
    std::pair<NodeID, EdgeID> priority_edge;

    // parallel
    std::vector<std::mutex> q_mutex;
    size_t num_threads;
    std::vector<size_t> branch_invalid;

    kernelization_criteria kc;
    maximum_flow mf;
    problem_management pm;
    measurements msm;
    FlowType last_sent_flow;
    std::atomic<double> log_timer;
    bool finished;

#ifdef USE_GUROBI
    ilp_model ilp;
#endif

    // MPI
    int mpi_size;
    int mpi_rank;
    mpi_communication mpic;
    int mpi_num_done;
    // only used in rank mpi_size - 1
    std::vector<int> mpi_done;
    bool sent_done;
    std::mutex mpi_mutex;

    size_t print_index = 0;
};
