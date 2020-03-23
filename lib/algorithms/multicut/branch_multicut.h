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
#include <chrono>
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
#include "algorithms/multicut/mpi_communication.h"
#include "algorithms/multicut/multicut_problem.h"
#include "algorithms/multicut/problem_queues/per_thread_problem_queue.h"
#include "algorithms/multicut/problem_queues/single_problem_queue.h"
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

using namespace std::chrono_literals;

class branch_multicut {
 public:
#ifndef NDEBUG
    static const bool debug = true;
#else
    static const bool debug = false;
#endif
    static const bool testing = true;

    branch_multicut(mutable_graph original_graph,
                    std::vector<NodeID> original_terminals,
                    std::vector<bool> fixed_vertex)
        : original_graph(original_graph),
          original_terminals(original_terminals),
          fixed_vertex(fixed_vertex),
          global_upper_bound(std::numeric_limits<FlowType>::max()),
          non_ls_global_upper_bound(std::numeric_limits<FlowType>::max()),
          problems(configuration::getConfig()->threads,
                   configuration::getConfig()->queue_type),
          total_time(),
          q_cv(configuration::getConfig()->threads),
          q_mutex(configuration::getConfig()->threads),
          num_threads(configuration::getConfig()->threads),
          branch_invalid(configuration::getConfig()->threads, 0),
          kc(original_terminals),
          mf(original_terminals),
          log_timer(0) {
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        mpi_num_done = 0;
        mpi_done.resize(mpi_size, 0);
        sent_done = false;
    }

    ~branch_multicut() { }

    std::pair<std::vector<NodeID>, size_t> find_multiterminal_cut(
        std::shared_ptr<multicut_problem> problem) {
        best_solution.resize(original_graph.number_of_nodes());
        if (mpi_rank == 0) {
            mf.maximumIsolatingFlow(problem, 0, /* parallel */ true);
            problems.addProblem(problem, 0);
        }
        updateBestSolution(problem);

        std::vector<std::thread> threads;
        is_finished = false;
        idle_threads = 0;
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

        for (auto& t : threads) {
            for (size_t j = 0; j < num_threads; ++j) {
                q_cv[j].notify_all();
            }
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

        FlowType total_weight = global_upper_bound;

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
        size_t bsize = best_solution.size();
        MPI_Bcast(&best_solution.front(), bsize, MPI_INT,
                  global_bcast_id, MPI_COMM_WORLD);

        total_weight = flowValue(false, best_solution);
        VIECUT_ASSERT_LEQ(total_weight, global_upper_bound);

        return std::make_pair(best_solution, total_weight);
    }

 private:
    bool queueNotEmpty(size_t thread_id) {
        return !problems.empty(thread_id) || is_finished
               || (idle_threads == num_threads);
    }

    void pollWork(size_t thread_id) {
        bool im_idle = false;
        while (!is_finished) {
            problems.prepareQueue(thread_id, global_upper_bound);
            if (!problems.empty(thread_id) || problems.haveASendProblem()) {
                if (thread_id == 0) {
                    FlowType update = mpic.getGlobalBestSolution();
                    global_upper_bound = std::min(global_upper_bound, update);
                }

                std::optional<int> sending = std::nullopt;
                if (mpi_size > 1 && thread_id == 0 && problems.size() > 1) {
                    sending = mpic.checkForReceiver();
                }

                auto problem = problems.pullProblem(thread_id,
                                                    sending.has_value());
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
                    ++idle_threads;
                    im_idle = true;
                }

                if (idle_threads == num_threads && problems.all_empty()) {
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
                            problems.addProblem(p, thread_id);
                        } else {
                            if (std::get<bool>(src)) {
                                is_finished = true;
                                for (size_t j = 0; j < num_threads; ++j) {
                                    q_cv[j].notify_all();
                                }
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
                                        is_finished = true;
                                        for (size_t j = 0; j < num_threads;
                                             ++j) {
                                            q_cv[j].notify_all();
                                        }
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

                std::unique_lock<std::mutex> lck(q_mutex[thread_id]);
                q_cv[thread_id].wait_for(
                    lck, 1000ms,
                    [this, thread_id] {
                        return queueNotEmpty(thread_id);
                    });

                if (im_idle) {
                    idle_threads--;
                    im_idle = false;
                }
            }
        }
    }

    void solveProblem(std::shared_ptr<multicut_problem> problem,
                      size_t thread_id) {
        if (problem == NULL) {
            LOG1 << "ERROR: Problem is NULL. This should not happen!";
            exit(1);
        }

        if (problem->lower_bound >= global_upper_bound)
            return;

#ifdef USE_TCMALLOC
        uint64_t heapsize = 0;
        MallocExtension::instance()->GetNumericProperty(
            "generic.heap_size", &heapsize);

        uint64_t max_size = 4UL * 1024UL * 1024UL * 1024UL;
        if (heapsize > max_size) {
            LOG1 << "RESULT Memoryout";
            exit(1);
        }
#endif

        if (total_time.elapsed() > 3600) {
            LOG1 << "RESULT Timeout!";
            exit(1);
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
                          << " global_upper:" << global_upper_bound
                          << " non_ls_upper:" << non_ls_global_upper_bound
                          << " queue.size:" << problems.size();
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
        if (problem->deleted_weight >
            static_cast<EdgeWeight>(global_upper_bound)) {
            return;
        }

        bool branchOnCurrentInstance = true;
        auto c = configuration::getConfig();
#ifdef USE_GUROBI
        branchOnCurrentInstance =
            problem->graph->n() > 1000 || (!c->use_ilp);
        if (!c->differences_set) {
            c->bound_difference = problem->upper_bound
                                  - problem->lower_bound;
            c->n = problem->graph->n();
            c->m = problem->graph->m();
            c->differences_set = true;
        }
#endif

        auto path = c->first_branch_path;
        if (path != "") {
            multicut_problem::writeGraph(problem, path);
            exit(1);
        }

        if (c->inexact) {
            size_t lowestTerminals =
                std::ceil(static_cast<double>(problem->terminals.size()) *
                          c->removeTerminalsBeforeBranch);

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

                std::unordered_set<NodeID> terminalpositions;
                for (auto t : problem->terminals) {
                    terminalpositions.insert(t.position);
                }

                std::unordered_set<NodeID> contractIntoTerminal;
                contractIntoTerminal.insert(lightest_t);
                for (size_t i = 0; i < best_solution.size(); ++i) {
                    NodeID cp = problem->graph->getCurrentPosition(i);
                    if (terminalpositions.count(cp) > 0
                        || cp >= problem->graph->n())
                        continue;
                    if (best_solution[i] == lightest_oid && cp != lightest_t) {
                        contractIntoTerminal.insert(cp);
                        for (auto v : problem->graph->containedVertices(cp)) {
                            best_solution[v] = lightest_oid;
                        }
                    }
                }

                NodeID invtx = problem->graph->containedVertices(lightest_t)[0];
                problem->graph->contractVertexSet(contractIntoTerminal);
                lightest_t = problem->graph->getCurrentPosition(invtx);
                size_t zero = 0;
                for (auto v : best_solution) {
                    if (v == 0) {
                        zero++;
                    }
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

        if (branchOnCurrentInstance) {
            branchOnEdge(problem, thread_id);
        } else {
            solve_with_ilp(problem, thread_id);
        }
    }

    void updateBestSolution(std::shared_ptr<multicut_problem> problem) {
        if (problem->upper_bound < non_ls_global_upper_bound) {
            auto cfg = configuration::getConfig();
            std::vector<NodeID> current_solution(original_graph.n(), 0);
            std::vector<size_t> blocksize(original_terminals.size(), 0);
            for (NodeID n = 0; n < current_solution.size(); ++n) {
                NodeID n_coarse = problem->mapped(n);
                auto t = problem->graph->getCurrentPosition(n_coarse);
                current_solution[n] = problem->graph->getPartitionIndex(t);
                blocksize[current_solution[n]]++;
            }

            if (debug) {
                std::vector<size_t> term_in_block(original_terminals.size(), 0);
                for (const auto& t : original_terminals) {
                    ++term_in_block[current_solution[t]];
                }

                for (size_t i = 0; i < term_in_block.size(); ++i) {
                    if (term_in_block[i] != 1) {
                        LOG1 << term_in_block[i] << " terminals on block " << i;
                        exit(1);
                    }
                }
            }
            std::unordered_set<NodeID> origtermset;
            for (NodeID t : original_terminals) {
                origtermset.insert(t);
            }

            FlowType prev_gub = flowValue(false, current_solution);
            local_search ls(problem, original_graph, original_terminals,
                            fixed_vertex, &current_solution);
            FlowType total_improvement = ls.improveSolution();
            FlowType ls_bound = prev_gub - total_improvement;

            bestsol_mutex.lock();

            if (ls_bound < global_upper_bound) {
                LOGC(testing) << "Improvement after time="
                              << total_time.elapsed() << " upper_bound="
                              << ls_bound;
            } else {
                LOGC(testing) << "No improvement, as " << ls_bound << " >= "
                              << global_upper_bound << " even though it was"
                              << " better before local search "
                              << "(" << prev_gub << " < "
                              << non_ls_global_upper_bound << ")";
            }

            if (ls_bound <= global_upper_bound) {
                global_upper_bound = ls_bound;
                non_ls_global_upper_bound = prev_gub;
                for (size_t i = 0; i < current_solution.size(); ++i) {
                    best_solution[i] = current_solution[i];
                }
            }
            bestsol_mutex.unlock();

            mpic.broadcastImprovedSolution(global_upper_bound);
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

    FlowType flowValue(bool verbose, const std::vector<NodeID>& sol) {
        std::vector<size_t> block_sizes(original_terminals.size(), 0);

        EdgeWeight total_weight = 0;
        for (NodeID n : original_graph.nodes()) {
            block_sizes[sol[n]]++;
            for (EdgeID e : original_graph.edges_of(n)) {
                NodeID tgt = original_graph.getEdgeTarget(n, e);
                EdgeWeight wgt = original_graph.getEdgeWeight(n, e);
                if (sol[n] != sol[tgt]) {
                    total_weight += wgt;
                }
            }
        }
        total_weight /= 2;

        for (size_t i = 0; i < block_sizes.size(); ++i) {
            size_t terminals = 0;
            LOGC(verbose) << "Block " << i << " has size " << block_sizes[i];
            if (block_sizes[i] == 0) {
                LOG1 << "ERROR: Empty block " << i;
                exit(1);
            }

            for (NodeID ot : original_terminals) {
                if (sol[ot] == i) {
                    terminals++;
                }
            }

            if (terminals != 1) {
                LOG1 << terminals << " terminals in block " << i;
                exit(1);
            }
        }
        return total_weight;
    }

    void branchOnEdge(std::shared_ptr<multicut_problem> problem,
                      size_t thread_id) {
        graph_contraction::deleteTermEdges(problem, original_terminals);
        if (!problem->graph->number_of_edges()) {
            problem->upper_bound = problem->deleted_weight;
            if (problem->upper_bound < non_ls_global_upper_bound) {
                updateBestSolution(problem);
            }
            return;
        }

        if (configuration::getConfig()->multibranch) {
            multiBranch(problem, thread_id);
        } else {
            singleBranch(problem, thread_id);
        }
    }

    void multiBranch(std::shared_ptr<multicut_problem> problem,
                     size_t thread_id) {
        auto [vertex, terminal_ids] = findEdgeMultiBranch(problem);
        std::unordered_set<NodeID> terminals;
        for (const auto& t : problem->terminals) {
            terminals.emplace(t.position);
        }

        for (size_t i = 0; i < terminal_ids.size(); ++i) {
            NodeID ctr_terminal = terminal_ids[i];
            std::shared_ptr<multicut_problem> new_p;
            if (i < terminal_ids.size() - 1) {
                new_p = std::make_shared<multicut_problem>();
                new_p->graph = std::make_shared<mutable_graph>(*problem->graph);
                new_p->terminals = problem->terminals;
                new_p->mappings = problem->mappings;
                new_p->priority_edge = { UNDEFINED_NODE, UNDEFINED_EDGE };
                new_p->lower_bound = problem->lower_bound;
                new_p->upper_bound = problem->upper_bound;
                new_p->deleted_weight = problem->deleted_weight;
                new_p->finished_blockpairs = problem->finished_blockpairs;
                for (size_t j = 0; j < new_p->terminals.size(); ++j) {
                    new_p->terminals[j].invalid_flow = true;
                }
            } else {
                new_p = problem;
            }

            NodeID coarse_vtx = new_p->graph->containedVertices(vertex)[0];
            bool finished = false;
            // first delete edges to terminals not picked
            while (!finished) {
                finished = true;
                vertex = new_p->graph->getCurrentPosition(coarse_vtx);
                for (size_t e = 0; e <
                     new_p->graph->get_first_invalid_edge(vertex); ++e) {
                    auto [tgt, wgt] = new_p->graph->getEdge(vertex, e);
                    if (terminals.count(tgt) > 0 && tgt != ctr_terminal) {
                        new_p->graph->deleteEdge(vertex, e);
                        new_p->deleted_weight += wgt;
                        auto p = new_p->graph->getCurrentPosition(coarse_vtx);
                        if (p != vertex) {
                            vertex = p;
                            finished = false;
                            break;
                        }
                        --e;
                    }
                }
            }

            for (EdgeID e : new_p->graph->edges_of(vertex)) {
                NodeID tgt = new_p->graph->getEdgeTarget(vertex, e);
                if (tgt == ctr_terminal) {
                    new_p->graph->contractEdge(vertex, e);
                    break;
                }
            }

            graph_contraction::deleteTermEdges(new_p, original_terminals);
            processNewProblem(new_p, thread_id);
        }
    }

    void processNewProblem(std::shared_ptr<multicut_problem> new_p,
                           size_t thread_id) {
        if (new_p->terminals.size() < 2
            && new_p->upper_bound < non_ls_global_upper_bound) {
            updateBestSolution(new_p);
        }

        if (new_p->terminals.size() == 2) {
            mf.maximumSTFlow(new_p);
            if (new_p->upper_bound < non_ls_global_upper_bound) {
                updateBestSolution(new_p);
            }
        } else {
            mf.maximumIsolatingFlow(new_p, thread_id, problems.size() == 0);
            graph_contraction::deleteTermEdges(new_p, original_terminals);
            if (new_p->graph->m() == 0) {
                new_p->upper_bound = new_p->deleted_weight;
                if (new_p->upper_bound < non_ls_global_upper_bound) {
                    updateBestSolution(new_p);
                }
                return;
            }

            if (new_p->lower_bound < global_upper_bound) {
                if (new_p->upper_bound < non_ls_global_upper_bound) {
                    updateBestSolution(new_p);
                }

                size_t thr = problems.addProblem(new_p, thread_id);
                q_cv[thr].notify_all();
            }
        }
    }

    void singleBranch(std::shared_ptr<multicut_problem> problem,
                      size_t thread_id) {
        auto [b_vtx, b_edge] = findEdgeSingleBranch(problem);
        NodeID target_terminal_id = UNDEFINED_NODE;
        auto b = problem->graph->getEdgeTarget(b_vtx, b_edge);
        for (size_t i = 0; i < problem->terminals.size(); ++i) {
            if (problem->terminals[i].position == b ||
                problem->terminals[i].position == b_vtx) {
                target_terminal_id = i;
            }
        }
        size_t max_wgt = problem->graph->getEdgeWeight(b_vtx, b_edge);
        // |-> edge in multicut
        if (max_wgt + problem->deleted_weight
            < static_cast<EdgeWeight>(global_upper_bound)) {
            // ^- if this is not true, there can not be a better cut
            // where the deleted edge is in the optimal multicut
            auto del_p = std::make_shared<multicut_problem>();
            del_p->graph = std::make_shared<mutable_graph>(*problem->graph);
            del_p->graph->deleteEdge(b_vtx, b_edge);
            del_p->terminals = problem->terminals;

            for (size_t i = 0; i < del_p->terminals.size(); ++i) {
                if (i != target_terminal_id) {
                    del_p->terminals[i].invalid_flow = true;
                }
            }

            del_p->mappings = problem->mappings;
            del_p->deleted_weight = problem->deleted_weight + max_wgt;
            del_p->priority_edge = { UNDEFINED_NODE, UNDEFINED_EDGE };
            del_p->lower_bound = problem->lower_bound;
            del_p->upper_bound = problem->upper_bound + max_wgt;
            processNewProblem(del_p, thread_id);
        }
        // |-> edge not in multicut
        if (!degreeThreeContraction(problem, b_vtx, b_edge)) {
            problem->graph->contractEdge(b_vtx, b_edge);
            problem->priority_edge = { UNDEFINED_NODE, UNDEFINED_EDGE };
            for (size_t i = 0; i < problem->terminals.size(); ++i) {
                if (i == target_terminal_id) {
                    problem->terminals[i].invalid_flow = true;
                }
            }
        }

        graph_contraction::deleteTermEdges(problem, original_terminals);
        processNewProblem(problem, thread_id);
    }

    bool degreeThreeContraction(std::shared_ptr<multicut_problem> problem,
                                NodeID b_vtx, EdgeID b_edge) {
        auto& g = problem->graph;
        auto c = g->getEdgeWeight(b_vtx, b_edge);
        NodeID term = g->getEdgeTarget(b_vtx, b_edge);
        EdgeWeight max_incident_wgt = 0;
        if (!(g->getUnweightedNodeDegree(b_vtx) == 3))
            return false;

        for (EdgeID e : g->edges_of(b_vtx)) {
            auto wgt = g->getEdgeWeight(b_vtx, e);
            if (wgt > max_incident_wgt) {
                max_incident_wgt = wgt;
            }
        }

        if (max_incident_wgt != c)
            return false;

        for (EdgeID e = 0; e < g->getUnweightedNodeDegree(b_vtx); ++e) {
            if (g->getEdgeTarget(b_vtx, e) != term) {
                problem->deleted_weight += g->getEdgeWeight(b_vtx, e);
                g->deleteEdge(b_vtx, e);
                --e;
            }
        }
        g->contractEdge(b_vtx, 0);

        for (auto& t : problem->terminals) {
            t.invalid_flow = true;
        }
        return true;
    }

    void nonBranchingContraction(std::shared_ptr<multicut_problem> problem) {
        auto pe = kc.kernelization(problem, global_upper_bound,
                                   problems.size() == 0);
        if (pe.has_value()) {
            problem->priority_edge = *pe;
        }
    }

#ifdef USE_GUROBI
    void solve_with_ilp(std::shared_ptr<multicut_problem> problem,
                        size_t thread_id) {
        graph_contraction::deleteTermEdges(problem, original_terminals);
        std::vector<NodeID> presets(problem->graph->n(),
                                    original_terminals.size());

        for (size_t i = 0; i < original_terminals.size(); ++i) {
            NodeID map = problem->mapped(original_terminals[i]);
            NodeID pos = problem->graph->getCurrentPosition(map);
            presets[pos] = i;
        }

        LOG1 << "start with deleted " << problem->deleted_weight;

        auto [result, wgt, reIntroduce] =
            ilp.computeIlp(problem, presets, original_terminals.size(),
                           problems.size() == 0, thread_id);
        problem->upper_bound = problem->deleted_weight + wgt;

        if (problem->upper_bound < non_ls_global_upper_bound) {
            for (const auto& n : problem->graph->nodes()) {
                problem->graph->setPartitionIndex(n, result[n]);
            }
            updateBestSolution(problem);
        }
        if (reIntroduce) {
            problems.addProblem(problem, thread_id);
        }
    }
#else
    void solve_with_ilp(std::shared_ptr<multicut_problem>, size_t) {
        LOG1 << "Error: Code not compiled with option -DUSE_GUROBI,"
             << " but called ILP solver. Exiting!";
        exit(1);
    }
#endif

    void findSubproblems(std::shared_ptr<multicut_problem> problem) {
        std::queue<NodeID> block;
        std::vector<bool> checked(problem->graph->number_of_nodes(), false);
        std::vector<bool> terminal(problem->graph->number_of_nodes(), false);
        NodeID n = 0;

        for (const auto& m : problem->terminals) {
            checked[m.position] = true;
            terminal[m.position] = true;
        }

        while (n < problem->graph->number_of_nodes()) {
            while (checked[n]) {
                n++;
            }

            if (n >= problem->graph->number_of_nodes()) {
                break;
            }

            std::vector<bool> t_this(problem->graph->number_of_nodes(), false);
            size_t terms = 0;
            size_t vtcs = 1;

            block.emplace(n);
            checked[n++] = true;

            while (!block.empty()) {
                NodeID node = block.front();
                block.pop();

                for (EdgeID e : problem->graph->edges_of(node)) {
                    NodeID tgt = problem->graph->getEdgeTarget(node, e);
                    if (!checked[tgt]) {
                        checked[tgt] = true;
                        block.emplace(tgt);
                        ++vtcs;
                    }

                    if (terminal[tgt]) {
                        if (!t_this[tgt]) {
                            t_this[tgt] = true;
                            ++terms;
                        }
                    }
                }
            }

            if (terms)
                LOG << "Number of terminals: " << terms
                    << " with vertices " << vtcs;
        }
    }

    mutable_graph original_graph;
    std::vector<NodeID> original_terminals;
    std::vector<bool> fixed_vertex;
    FlowType global_upper_bound;
    FlowType non_ls_global_upper_bound;
    per_thread_problem_queue problems;
    std::vector<NodeID> best_solution;
    timer total_time;
    std::pair<NodeID, EdgeID> priority_edge;

    // parallel
    std::vector<std::condition_variable> q_cv;
    std::vector<std::mutex> q_mutex;
    std::atomic<uint> idle_threads;
    size_t num_threads;
    std::vector<size_t> branch_invalid;
    bool is_finished;

    kernelization_criteria kc;
    maximum_flow mf;
    std::atomic<double> log_timer;
    std::mutex bestsol_mutex;

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
