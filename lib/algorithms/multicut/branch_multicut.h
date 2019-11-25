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

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <future>
#include <iterator>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>

#include "algorithms/flow/push_relabel.h"
#include "algorithms/global_mincut/noi_minimum_cut.h"
#include "algorithms/misc/graph_algorithms.h"
#include "algorithms/multicut/edge_selection.h"
#include "algorithms/multicut/graph_contraction.h"
#include "algorithms/multicut/kernelization_criteria.h"
#include "algorithms/multicut/multicut_problem.h"
#include "algorithms/multicut/problem_queues/per_thread_problem_queue.h"
#include "algorithms/multicut/problem_queues/single_problem_queue.h"
#include "coarsening/contract_graph.h"
#include "common/configuration.h"
#include "data_structure/union_find.h"
#include "gperftools/malloc_extension.h"
#include "io/graph_io.h"
#include "tlx/math/div_ceil.hpp"
#include "tools/timer.h"
#include "tools/vector.h"

#ifdef USE_GUROBI
#include "algorithms/multicut/ilp_model.h"
#endif

using namespace std::chrono_literals;

class branch_multicut {
 public:
    static const bool debug = false;
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
          kc(configuration::getConfig()->contraction_type, original_terminals),
          log_timer(0) { }

    ~branch_multicut() { }

    size_t find_multiterminal_cut(std::shared_ptr<multicut_problem> mcp) {
        best_solution.resize(original_graph.number_of_nodes());
        auto cfg = configuration::getConfig();
        problems.addProblem(mcp, 0);
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
                std::shared_ptr<multicut_problem> mcp =
                    problems.pullProblem(thread_id);
                solveProblem(mcp, thread_id);
            } else {
                if (!im_idle) {
                    idle_threads++;
                    im_idle = true;
                }
                if (idle_threads == num_threads && problems.all_empty()) {
                    is_finished = true;
                    for (size_t j = 0; j < num_threads; ++j) {
                        q_cv[j].notify_all();
                    }
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

    void solveProblem(std::shared_ptr<multicut_problem> current_problem,
                      size_t thread_id) {
        if (current_problem->lower_bound >= global_upper_bound)
            return;

        uint64_t heapsize;
        MallocExtension::instance()->GetNumericProperty(
            "generic.heap_size", &heapsize);

        uint64_t max_size = 250UL * 1024UL * 1024UL * 1024UL;
        if (heapsize > max_size) {
            LOG1 << "RESULT Memoryout";
            exit(1);
        }

        if (total_time.elapsed() > 3600) {
            LOG1 << "RESULT Timeout!";
            exit(1);
        }

        if (total_time.elapsed() > log_timer) {
            double logs_per_second = 2.0;
            double time_added = 1.0 / logs_per_second;

            log_timer = time_added + log_timer;

            LOGC(testing) << "After " << total_time.elapsed()
                          << " - terminals:"
                          << current_problem->terminals.size()
                          << " vertices:"
                          << current_problem->graph->number_of_nodes()
                          << " deleted:" << current_problem->deleted_weight
                          << " lower:" << current_problem->lower_bound
                          << " upper:" << current_problem->upper_bound
                          << " global_upper:" << global_upper_bound
                          << " queue.size:" << problems.size();

            // LOG1 << current_problem.path;
        }

        if (current_problem->graph->n() * 2
            < current_problem->graph->getOriginalNodes()) {
            auto map = std::make_shared<std::vector<NodeID> >();

            for (size_t i = 0;
                 i < current_problem->graph->getOriginalNodes(); ++i) {
                map->emplace_back(
                    current_problem->graph->getCurrentPosition(i));
            }

            current_problem->mappings.emplace_back(map);
            current_problem->graph = current_problem->graph->simplify();
        }

        graph_contraction::setTerminals(current_problem, original_terminals);
        NodeID edges_before = current_problem->graph->m();

        if (current_problem->terminals.size() == 2) {
            FlowType max =
                maximumFlow(current_problem) + current_problem->deleted_weight;

            if (max < global_upper_bound) {
                bestsol_mutex.lock();
                if (max < global_upper_bound) {
                    LOGC(testing) << "Improvement after time="
                                  << total_time.elapsed() << " upper_bound="
                                  << current_problem->upper_bound;
                    LOGC(testing) << current_problem->path;

                    global_upper_bound = max;
                    for (NodeID n = 0; n < best_solution.size(); ++n) {
                        NodeID n_coarse = current_problem->mapped(n);
                        auto t = current_problem->graph->getCurrentPosition(
                            n_coarse);
                        best_solution[n] =
                            current_problem->graph->getPartitionIndex(t);
                    }
                }
                bestsol_mutex.unlock();
            }
        } else {
            nonBranchingContraction(current_problem, thread_id);
            if (current_problem->graph->m() >= edges_before
                && current_problem->terminals.size() > 1) {
#ifdef USE_GUROBI
                auto c = configuration::getConfig();
                bool branchOnCurrentInstance = !c->use_ilp;
                if (!c->differences_set) {
                    c->bound_difference = current_problem->upper_bound
                                          - current_problem->lower_bound;
                    c->n = current_problem->graph->n();
                    c->m = current_problem->graph->m();
                    c->differences_set = true;
                }

                if (branchOnCurrentInstance) {
#endif
                branchOnEdge(current_problem, thread_id);
#ifdef USE_GUROBI
            } else {
                solve_with_ilp(current_problem);
            }
#endif
            } else {
                if (current_problem->lower_bound < global_upper_bound) {
                    size_t thr =
                        problems.addProblem(current_problem, thread_id);
                    q_cv[thr].notify_all();
                } else {
                    // notify all sleeping threads, as we might be finished
                    for (size_t j = 0; j < num_threads; ++j) {
                        q_cv[j].notify_all();
                    }
                }
            }
        }
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

    void branchOnEdge(std::shared_ptr<multicut_problem> current_problem,
                      size_t thread_id) {
        graph_contraction::deleteEdgesBetweenTerminals(
            current_problem, original_terminals);

        if (!current_problem->graph->number_of_edges()) {
            if (static_cast<FlowType>(current_problem->deleted_weight)
                < global_upper_bound) {
                bestsol_mutex.lock();
                if (static_cast<FlowType>(current_problem->deleted_weight)
                    < global_upper_bound) {
                    global_upper_bound =
                        std::min(static_cast<FlowType>(
                                     current_problem->deleted_weight),
                                 global_upper_bound);
                    LOGC(testing) << "noedges - Improvement after time="
                                  << total_time.elapsed() << " upper_bound="
                                  << current_problem->upper_bound;
                    LOGC(testing) << current_problem->path;

                    for (NodeID n = 0; n < best_solution.size(); ++n) {
                        NodeID n_coarse = current_problem->mapped(n);
                        auto t = current_problem->graph->getCurrentPosition(
                            n_coarse);
                        best_solution[n] =
                            current_problem->graph->getPartitionIndex(t);
                    }
                }
                bestsol_mutex.unlock();
            }
            return;
        }

        auto [branch_vtx, branch_edge] =
            findEdge(current_problem, edge_selection);

        size_t max_wgt =
            current_problem->graph->getEdgeWeight(branch_vtx, branch_edge);
        // |-> edge in multicut
        if (max_wgt + current_problem->deleted_weight
            < (EdgeWeight)global_upper_bound) {
            // ^- if this is not true, there can not be a better cut
            // where the deleted edge is in the optimal multicut
            auto delete_problem = std::make_shared<multicut_problem>();

            delete_problem->graph = std::make_shared<mutable_graph>(
                *current_problem->graph);
            delete_problem->graph->deleteEdge(branch_vtx, branch_edge);
            delete_problem->terminals = current_problem->terminals;

            for (auto& t : delete_problem->terminals) {
                if (t.position != branch_vtx) {
                    t.invalid_flow = true;
                }
            }

            delete_problem->mappings = current_problem->mappings;
            delete_problem->lower_bound = current_problem->lower_bound;
            delete_problem->deleted_weight =
                current_problem->deleted_weight + max_wgt;

            delete_problem->lower_bound = current_problem->lower_bound;
            delete_problem->upper_bound =
                current_problem->upper_bound + max_wgt;

            if (delete_problem->lower_bound < global_upper_bound) {
                size_t thr = problems.addProblem(delete_problem, thread_id);
                q_cv[thr].notify_all();
            }
        }

        // |-> edge not in multicut
        current_problem->graph->contractEdge(branch_vtx, branch_edge);

        for (auto& t : current_problem->terminals) {
            if (t.position == branch_vtx) {
                t.invalid_flow = true;
            }
        }

        graph_contraction::deleteEdgesBetweenTerminals(current_problem,
                                                       original_terminals);

        if (current_problem->lower_bound < global_upper_bound) {
            size_t thr = problems.addProblem(current_problem, thread_id);
            q_cv[thr].notify_all();
        }
    }

    void nonBranchingContraction(std::shared_ptr<multicut_problem> mcp,
                                 size_t thread_id) {
        FlowType contracting_flow = maximumIsolatingFlow(mcp, thread_id);
        if (mcp->upper_bound < global_upper_bound) {
            bestsol_mutex.lock();
            // check again inside mutex as global_upper_bound might have changed
            if (mcp->upper_bound < global_upper_bound) {
                global_upper_bound = mcp->upper_bound;

                for (NodeID n = 0; n < best_solution.size(); ++n) {
                    NodeID n_coarse = mcp->mapped(n);
                    auto t = mcp->graph->getCurrentPosition(n_coarse);
                    best_solution[n] = mcp->graph->getPartitionIndex(t);
                }

                global_upper_bound = flowValue(false, best_solution);
                LOGC(testing) << "Improvement after time="
                              << total_time.elapsed() << " upper_bound="
                              << global_upper_bound;
            }
            bestsol_mutex.unlock();
        }
        kc.perform_kernelization(mcp, global_upper_bound, contracting_flow);
    }

    FlowType maximumFlow(std::shared_ptr<multicut_problem> problem) {
        push_relabel pr;
        auto G = problem->graph;

        std::vector<NodeID> current_terminals;
        for (const auto& t : problem->terminals) {
            current_terminals.emplace_back(t.position);
        }

        auto [flow, isolating_block] =
            pr.solve_max_flow_min_cut(G, current_terminals, 0, true);

        NodeID term0 = problem->terminals[0].original_id;
        NodeID term1 = problem->terminals[1].original_id;

        for (NodeID n : problem->graph->nodes()) {
            G->setPartitionIndex(n, term1);
        }

        for (size_t i = 0; i < original_terminals.size(); ++i) {
            NodeID c =
                G->getCurrentPosition(problem->mapped(original_terminals[i]));
            G->setPartitionIndex(c, i);
        }

        for (NodeID n : isolating_block) {
            G->setPartitionIndex(n, term0);
        }

        return flow;
    }

    FlowType maximumIsolatingFlow(std::shared_ptr<multicut_problem> problem,
                                  size_t thread_id) {
        graph_contraction::setTerminals(problem, original_terminals);
        std::vector<FlowType> isolating_flow;
        std::vector<std::vector<NodeID> > maxVolIsoBlock;
        std::vector<NodeID> curr_terminals;
        std::vector<NodeID> orig_index;
        for (const auto& t : problem->terminals) {
            curr_terminals.emplace_back(t.position);
            orig_index.emplace_back(t.original_id);
        }

        bool parallel_flows = false;
        std::vector<std::future<std::vector<NodeID> > > futures;
        // so futures don't lose their object :)
        std::vector<push_relabel> prs(problem->terminals.size());
        if (problems.size() == 0) {
            // in the beginning when we don't have many problems
            // already (but big graphs), we can start a thread per flow.
            // later on, we have a problem for each processor to work on,
            // so we run sequential flows
            parallel_flows = true;
        }

        for (NodeID i = 0; i < problem->terminals.size(); ++i) {
            if (problem->terminals[i].invalid_flow) {
                if (parallel_flows) {
                    maxVolIsoBlock.emplace_back();
                    cpu_set_t all_cores;
                    CPU_ZERO(&all_cores);
                    for (size_t i = 0; i < num_threads; ++i) {
                        CPU_SET(i, &all_cores);
                    }

                    sched_setaffinity(0, sizeof(cpu_set_t), &all_cores);
                    futures.emplace_back(
                        std::async(&push_relabel::callable_max_flow,
                                   &prs[i],
                                   problem->graph, curr_terminals, i, true));
                } else {
                    push_relabel pr;
                    maxVolIsoBlock.emplace_back(
                        pr.solve_max_flow_min_cut(problem->graph,
                                                  curr_terminals,
                                                  i,
                                                  true).second);
                }

                problem->terminals[i].invalid_flow = false;
            } else {
                maxVolIsoBlock.emplace_back();
                maxVolIsoBlock.back().emplace_back(
                    problem->terminals[i].position);
            }
        }

        if (parallel_flows) {
            for (size_t i = 0; i < futures.size(); ++i) {
                auto& t = futures[i];
                maxVolIsoBlock[i] = t.get();
            }

            cpu_set_t my_id;
            CPU_ZERO(&my_id);
            CPU_SET(thread_id, &my_id);
            sched_setaffinity(0, sizeof(cpu_set_t), &my_id);
        }

        graph_contraction::contractIsolatingBlocks(problem, maxVolIsoBlock);

        EdgeWeight maximum = 0;
        FlowType sum = 0;
        size_t max_index = 0;

        for (size_t i = 0; i < problem->terminals.size(); ++i) {
            NodeID orig_id = problem->terminals[i].original_id;
            NodeID term_id = problem->mapped(original_terminals[orig_id]);
            NodeID pos = problem->graph->getCurrentPosition(term_id);
            EdgeWeight deg = problem->graph->getWeightedNodeDegree(pos);

            if (deg > maximum) {
                max_index = orig_id;
                maximum = deg;
            }
            sum += deg;

            isolating_flow.emplace_back(deg);
        }

        std::sort(isolating_flow.begin(), isolating_flow.end());

        FlowType contr_flow = 0;

        if (isolating_flow.size() >= 2)
            contr_flow = sum - isolating_flow[isolating_flow.size() - 1]
                         - isolating_flow[isolating_flow.size() - 2];

        for (NodeID n : problem->graph->nodes()) {
            problem->graph->setPartitionIndex(n, max_index);
        }

        for (size_t l = 0; l < original_terminals.size(); ++l) {
            NodeID l_coarse = problem->mapped(original_terminals[l]);
            NodeID v = problem->graph->getCurrentPosition(l_coarse);
            problem->graph->setPartitionIndex(v, l);
        }

        problem->upper_bound = problem->deleted_weight + sum - maximum;
        problem->lower_bound = problem->deleted_weight + tlx::div_ceil(sum, 2);

        // findSubproblems(problem);
        if (debug) {
            graph_algorithms::checkGraphValidity(problem->graph);
        }

        return contr_flow;
    }

#ifdef USE_GUROBI
    void solve_with_ilp(std::shared_ptr<multicut_problem> mcp) {
        std::vector<NodeID> presets(mcp->graph->n(), mcp->terminals.size());

        for (size_t i = 0; i < original_terminals.size(); ++i) {
            NodeID map = mcp->mapped(original_terminals[i]);
            presets[mcp->graph->getCurrentPosition(map)] = i;
        }

        auto [result, wgt] = ilp_model::computeIlp(mcp->graph, presets,
                                                   original_terminals.size());

        if (mcp->deleted_weight + wgt <
            static_cast<EdgeWeight>(global_upper_bound)) {
            for (NodeID n = 0; n < best_solution.size(); ++n) {
                NodeID n_coarse = mcp->mapped(n);
                auto t = mcp->graph->getCurrentPosition(n_coarse);
                best_solution[n] = result[t];
            }
            global_upper_bound = wgt + mcp->deleted_weight;
            LOGC(testing) << "Improvement (IN ILP) after time="
                          << total_time.elapsed()
                          << " upper_bound=" << global_upper_bound;
        }
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

    // parallel
    std::vector<std::condition_variable> q_cv;
    std::vector<std::mutex> q_mutex;
    std::atomic<uint> idle_threads;
    size_t num_threads;
    std::vector<size_t> branch_invalid;
    bool is_finished;

    std::string edge_selection;
    kernelization_criteria kc;
    std::atomic<double> log_timer;
    std::mutex bestsol_mutex;
};
