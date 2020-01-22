/******************************************************************************
 * branch_multicut.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2018 Alexander Noe <alexander.noe@univie.ac.at>
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
#include <unordered_set>
#include <utility>
#include <vector>

#include "algorithms/global_mincut/noi_minimum_cut.h"
#include "algorithms/misc/graph_algorithms.h"
#include "algorithms/multicut/edge_selection.h"
#include "algorithms/multicut/graph_contraction.h"
#include "algorithms/multicut/kernelization_criteria.h"
#include "algorithms/multicut/maximum_flow.h"
#include "algorithms/multicut/multicut_problem.h"
#include "algorithms/multicut/problem_queues/per_thread_problem_queue.h"
#include "algorithms/multicut/problem_queues/single_problem_queue.h"
#include "coarsening/contract_graph.h"
#include "common/configuration.h"
#include "data_structure/union_find.h"
#include "gperftools/malloc_extension.h"
#include "io/graph_io.h"
#include "tools/timer.h"
#include "tools/vector.h"

#ifdef USE_GUROBI
#include "algorithms/multicut/ilp_model.h"
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

    branch_multicut(std::shared_ptr<mutable_graph> original_graph,
                    std::vector<NodeID> original_terminals)
        : original_graph(*original_graph),
          original_terminals(original_terminals),
          global_upper_bound(std::numeric_limits<FlowType>::max()),
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
        mpi_finished.resize(mpi_size, true);
        mpi_finished[0] = false;
    }

    ~branch_multicut() { }

    size_t find_multiterminal_cut(std::shared_ptr<multicut_problem> problem) {
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
            cpu_set_t cpuset;
            CPU_ZERO(&cpuset);
            CPU_SET(i, &cpuset);
            pthread_setaffinity_np(threads[i].native_handle(),
                                   sizeof(cpu_set_t), &cpuset);
        }

        for (auto& t : threads) {
            for (size_t j = 0; j < num_threads; ++j) {
                q_cv[j].notify_all();
            }
            t.join();
        }
        FlowType total_weight = flowValue(false, best_solution);

        VIECUT_ASSERT_EQ(total_weight, global_upper_bound);

        return (EdgeWeight)total_weight;
    }

 private:
    bool queueNotEmpty(size_t thread_id) {
        return !problems.empty(thread_id) || is_finished;
    }

    void pollWork(size_t thread_id) {
        bool im_idle = false;
        while (!is_finished) {
            if (!problems.empty(thread_id)) {
                std::shared_ptr<multicut_problem> problem =
                    problems.pullProblem(thread_id);

                // Print small graph to file, left in as used in testing

                /*if (problem->graph->n() < 2000) {
                    bestsol_mutex.lock();
                    LOG1 << "Writing...";
                    std::string gid = "small_graph";
                    graph_io::writeGraphWeighted(
                        problem->graph->to_graph_access(), gid);
                    for (auto t : problem->terminals) {
                        LOG1 << t.position;
                    }
                    LOG1 << "...done!";
                    exit(1);
                }*/
                for (size_t i = 0; i < mpi_finished.size(); ++i) {
                    if (mpi_finished[i]) {
                        // sendProblem(problem, i);
                        break;
                    }
                }
                solveProblem(problem, thread_id);
            } else {
                if (!im_idle) {
                    idle_threads++;
                    im_idle = true;
                }
                if (idle_threads == num_threads && problems.all_empty()) {
                    is_finished = true;
                    // recvProblem();
                    for (size_t j = 0; j < num_threads; ++j) {
                        q_cv[j].notify_all();
                    }

                    LOG1 << "IM FINISHED";
                    // MPI_Barrier(MPI_COMM_WORLD);

                    return;
                }

                std::unique_lock<std::mutex> lck(q_mutex[thread_id]);
                q_cv[thread_id].wait_for(
                    lck, 1000ms,
                    [this, thread_id] { return queueNotEmpty(thread_id); });

                if (im_idle) {
                    idle_threads--;
                    im_idle = false;
                }
            }
        }
    }

    void sendProblem(std::shared_ptr<multicut_problem> problem, size_t tgt) {
        MPI_Send(&problem->lower_bound, 1, MPI_LONG, tgt, 523, MPI_COMM_WORLD);
        MPI_Send(&problem->upper_bound, 1, MPI_LONG, tgt, 0, MPI_COMM_WORLD);
        MPI_Send(&problem->deleted_weight, 1, MPI_LONG, tgt, 0, MPI_COMM_WORLD);

        // terminals
        size_t termsize = problem->terminals.size();
        MPI_Send(&termsize, 1, MPI_LONG,
                 tgt, 0, MPI_COMM_WORLD);
        MPI_Send(&problem->terminals.front(), problem->terminals.size(),
                 MPI_LONG, tgt, 0, MPI_COMM_WORLD);

        // mappings
        size_t mapsize = problem->mappings.size();
        MPI_Send(&mapsize, 1, MPI_LONG, tgt, 0, MPI_COMM_WORLD);
        for (size_t i = 0; i < problem->mappings.size(); ++i) {
            size_t mapisize = problem->mappings[i]->size();
            MPI_Send(&mapisize, 1, MPI_LONG, tgt, 0, MPI_COMM_WORLD);
            MPI_Send(&problem->mappings[i]->front(),
                     problem->mappings[i]->size(),
                     MPI_LONG, tgt, 0, MPI_COMM_WORLD);
        }
    }

    void recvProblem() {
        std::vector<NodeID> terminals;
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
        terminals.resize(termsize);
        MPI_Recv(&terminals.front(), termsize, MPI_LONG, status.MPI_SOURCE,
                 MPI_ANY_TAG, MPI_COMM_WORLD, NULL);

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
            MPI_Recv(&map->front(), currmapsize, MPI_LONG, status.MPI_SOURCE,
                     MPI_ANY_TAG, MPI_COMM_WORLD, NULL);
            mappings.emplace_back(map);
        }
    }

    void solveProblem(std::shared_ptr<multicut_problem> problem,
                      size_t thread_id) {
        if (problem->lower_bound >= global_upper_bound)
            return;

#ifdef USE_TCMALLOC
        uint64_t heapsize = 0;
        MallocExtension::instance()->GetNumericProperty(
            "generic.heap_size", &heapsize);

        uint64_t max_size = 50UL * 1024UL * 1024UL * 1024UL;
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
            double logs_per_second = 2.0;
            double time_added = 1.0 / logs_per_second;

            log_timer = time_added + log_timer;

            LOGC(testing) << "After " << total_time.elapsed()
                          << " - terminals:" << problem->terminals.size()
                          << " vertices:" << problem->graph->n()
                          << " deleted:" << problem->deleted_weight
                          << " lower:" << problem->lower_bound
                          << " upper:" << problem->upper_bound
                          << " global_upper:" << global_upper_bound
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
        bool branchOnCurrentInstance = true;
#ifdef USE_GUROBI
        auto c = configuration::getConfig();
        branchOnCurrentInstance = problem->graph->m() > 100000;
        if (!c->differences_set) {
            c->bound_difference = problem->upper_bound
                                  - problem->lower_bound;
            c->n = problem->graph->n();
            c->m = problem->graph->m();
            c->differences_set = true;
        }
#endif
        if (branchOnCurrentInstance) {
            branchOnEdge(problem, thread_id);
        } else {
            solve_with_ilp(problem, thread_id);
        }
    }

    void updateBestSolution(std::shared_ptr<multicut_problem> problem) {
        bestsol_mutex.lock();
        if (problem->upper_bound < global_upper_bound) {
            global_upper_bound = problem->upper_bound;
            for (NodeID n = 0; n < best_solution.size(); ++n) {
                NodeID n_coarse = problem->mapped(n);
                auto t = problem->graph->getCurrentPosition(n_coarse);
                best_solution[n] = problem->graph->getPartitionIndex(t);
            }

            if (debug) {
                std::vector<size_t> term_in_block(original_terminals.size(), 0);
                for (const auto& t : original_terminals) {
                    ++term_in_block[best_solution[t]];
                }

                for (size_t i = 0; i < term_in_block.size(); ++i) {
                    if (term_in_block[i] != 1) {
                        LOG1 << term_in_block[i] << " terminals on block " << i;
                        exit(1);
                    }
                }
            }

            global_upper_bound = flowValue(false, best_solution);

            LOGC(testing) << "Improvement after time="
                          << total_time.elapsed() << " upper_bound="
                          << global_upper_bound;
        }
        bestsol_mutex.unlock();
    }

    FlowType flowValue(bool verbose, const std::vector<NodeID>& sol) {
        std::vector<size_t> block_sizes(original_terminals.size(), 0);
        for (NodeID n : original_graph.nodes()) {
            original_graph.setPartitionIndex(n, sol[n]);
            block_sizes[sol[n]]++;
        }

        EdgeWeight total_weight = 0;
        for (NodeID n : original_graph.nodes()) {
            for (EdgeID e : original_graph.edges_of(n)) {
                NodeID tgt = original_graph.getEdgeTarget(n, e);
                EdgeWeight wgt = original_graph.getEdgeWeight(n, e);
                if (original_graph.getPartitionIndex(n)
                    != original_graph.getPartitionIndex(tgt)) {
                    LOG << "PARTITION INDEX OF " << n << " IS "
                        << original_graph.getPartitionIndex(n)
                        << " WHILE " << tgt << " IS "
                        << original_graph.getPartitionIndex(tgt);

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
            if (problem->upper_bound < global_upper_bound) {
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
            if (new_p->terminals.size() == 2) {
                mf.maximumSTFlow(new_p);
                if (new_p->upper_bound < global_upper_bound) {
                    updateBestSolution(new_p);
                }
            } else if (new_p->terminals.size() > 2) {
                mf.maximumIsolatingFlow(new_p, thread_id, problems.size() == 0);
                graph_contraction::deleteTermEdges(new_p, original_terminals);
                if (new_p->graph->m() == 0) {
                    new_p->upper_bound = new_p->deleted_weight;
                    if (new_p->upper_bound < global_upper_bound) {
                        updateBestSolution(new_p);
                    }
                    continue;
                }

                if (new_p->lower_bound < global_upper_bound) {
                    if (new_p->upper_bound < global_upper_bound) {
                        updateBestSolution(new_p);
                    }

                    size_t thr = problems.addProblem(new_p, thread_id);
                    q_cv[thr].notify_all();
                }
            } else {
                if (new_p->upper_bound < global_upper_bound) {
                    updateBestSolution(new_p);
                }
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
            if (del_p->lower_bound < global_upper_bound) {
                size_t thr = problems.addProblem(del_p, thread_id);
                q_cv[thr].notify_all();
            }
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
        if (problem->lower_bound < global_upper_bound) {
            size_t thr = problems.addProblem(problem, thread_id);
            q_cv[thr].notify_all();
        }
    }

    bool degreeThreeContraction(std::shared_ptr<multicut_problem> problem,
                                NodeID b_vtx, EdgeID b_edge) {
        auto& g = problem->graph;
        auto c = g->getEdgeWeight(b_vtx, b_edge);
        NodeID term = g->getEdgeTarget(b_vtx, b_edge);
        EdgeWeight max_incident_wgt = 0;
        if (!(g->getNodeDegree(b_vtx) == 3))
            return false;

        for (EdgeID e : g->edges_of(b_vtx)) {
            auto wgt = g->getEdgeWeight(b_vtx, e);
            if (wgt > max_incident_wgt) {
                max_incident_wgt = wgt;
            }
        }

        if (max_incident_wgt != c)
            return false;

        for (EdgeID e = 0; e < g->getNodeDegree(b_vtx); ++e) {
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
        auto pe = kc.kernelization(problem, global_upper_bound);
        if (pe.has_value()) {
            problem->priority_edge = *pe;
        }
    }

#ifdef USE_GUROBI
    void solve_with_ilp(std::shared_ptr<multicut_problem> problem,
                        size_t thread_id) {
        std::vector<NodeID> presets(problem->graph->n(),
                                    original_terminals.size());

        for (size_t i = 0; i < original_terminals.size(); ++i) {
            NodeID map = problem->mapped(original_terminals[i]);
            NodeID pos = problem->graph->getCurrentPosition(map);
            presets[pos] = i;
        }

        LOG1 << "start with deleted " << problem->deleted_weight;

        auto [result, wgt] = ilp_model::computeIlp(problem, presets,
                                                   original_terminals.size(),
                                                   problems.size() == 0,
                                                   thread_id);
        problem->upper_bound = problem->deleted_weight + wgt;

        if (problem->upper_bound < global_upper_bound) {
            for (const auto& n : problem->graph->nodes()) {
                problem->graph->setPartitionIndex(n, result[n]);
            }
            updateBestSolution(problem);
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
    FlowType global_upper_bound;
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

    // MPI
    int mpi_size;
    int mpi_rank;
    std::vector<bool> mpi_finished;
};
