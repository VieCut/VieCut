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
#include <unordered_set>
#include <utility>
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
          kc(original_terminals),
          log_timer(0) { }

    ~branch_multicut() { }

    size_t find_multiterminal_cut(std::shared_ptr<multicut_problem> problem) {
        best_solution.resize(original_graph.number_of_nodes());
        auto cfg = configuration::getConfig();
        problems.addProblem(problem, 0);
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

                /*if (problem->graph->n() < 2000) {
                    graphs++;
                    if (graphs == 100) {
                    LOG1 << "Writing...";

                    std::string gid = "graph100";
                    graph_io::writeGraphWeighted(problem->graph->to_graph_access(), gid);

                    for (auto t : problem->terminals) {
                        LOG1 << t.position;
                    }
                    LOG1 << "...done!";
                    exit(1);
                    }
                }*/

                solveProblem(problem, thread_id);
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

    void solveProblem(std::shared_ptr<multicut_problem> problem,
                      size_t thread_id) {
        if (problem->lower_bound >= global_upper_bound)
            return;

        uint64_t heapsize;
        MallocExtension::instance()->GetNumericProperty(
            "generic.heap_size", &heapsize);

        uint64_t max_size = 200UL * 1024UL * 1024UL * 1024UL;
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
                          << " - terminals:" << problem->terminals.size()
                          << " vertices:" << problem->graph->n()
                          << " deleted:" << problem->deleted_weight
                          << " lower:" << problem->lower_bound
                          << " upper:" << problem->upper_bound
                          << " global_upper:" << global_upper_bound
                          << " queue.size:" << problems.size();
        }

        if (problem->graph->n() * 2
            < problem->graph->getOriginalNodes()) {
            auto map = std::make_shared<std::vector<NodeID> >();

            for (size_t i = 0;
                 i < problem->graph->getOriginalNodes(); ++i) {
                map->emplace_back(
                    problem->graph->getCurrentPosition(i));
            }

            problem->mappings.emplace_back(map);
            problem->graph = problem->graph->simplify();
        }

        graph_contraction::setTerminals(problem, original_terminals);
        NodeID edges_before = problem->graph->m();

        if (problem->terminals.size() == 2) {
            maximumFlow(problem);
            if (problem->upper_bound < global_upper_bound) {
                updateBestSolution(problem);
            }
        } else {
            nonBranchingContraction(problem, thread_id);
            if (problem->graph->m() >= edges_before
                && problem->terminals.size() > 1) {
#ifdef USE_GUROBI
                auto c = configuration::getConfig();
                NodeID r = random_functions::nextInt(0, 300000);
                bool branchOnCurrentInstance = problem->graph->m() > r;
                if (!c->differences_set) {
                    c->bound_difference = problem->upper_bound
                                          - problem->lower_bound;
                    c->n = problem->graph->n();
                    c->m = problem->graph->m();
                    c->differences_set = true;
                }

                if (branchOnCurrentInstance) {
#endif
                branchOnEdge(problem, thread_id);
#ifdef USE_GUROBI
            } else {
                solve_with_ilp(problem);
            }
#endif
            } else {
                if (problem->lower_bound < global_upper_bound) {
                    size_t thr =
                        problems.addProblem(problem, thread_id);
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

    void updateBestSolution(std::shared_ptr<multicut_problem> problem) {
        bestsol_mutex.lock();
        if (problem->upper_bound < global_upper_bound) {
            LOGC(testing) << "Improvement after time="
                            << total_time.elapsed() << " upper_bound="
                            << problem->upper_bound;

            global_upper_bound = problem->upper_bound;
            for (NodeID n = 0; n < best_solution.size(); ++n) {
                NodeID n_coarse = problem->mapped(n);
                auto t = problem->graph->getCurrentPosition(n_coarse);
                best_solution[n] = problem->graph->getPartitionIndex(t);
            }
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
                for (size_t e = 0; e <
                     new_p->graph->get_first_invalid_edge(vertex); ++e) {
                    auto [tgt, wgt] = new_p->graph->getEdge(vertex, e);

                    if (terminals.count(tgt) > 0 && tgt != ctr_terminal) {
                        new_p->graph->deleteEdge(vertex, e);
                        auto p = new_p->graph->getCurrentPosition(coarse_vtx);
                        new_p->deleted_weight += wgt;
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
            graph_contraction::setTerminals(new_p, original_terminals);

            if (new_p->lower_bound < global_upper_bound) {
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

    void nonBranchingContraction(std::shared_ptr<multicut_problem> problem,
                                 size_t thread_id) {
        FlowType flow = maximumIsolatingFlow(problem, thread_id);
        if (problem->upper_bound < global_upper_bound) {
            updateBestSolution(problem);
        }
        auto pe = kc.kernelization(problem, global_upper_bound, flow);
        if (pe.has_value()) {
            problem->priority_edge = *pe;
        }
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

        problem->upper_bound = flow + problem->deleted_weight;
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
    void solve_with_ilp(std::shared_ptr<multicut_problem> problem) {
        std::vector<NodeID> presets(problem->graph->n(),
                                    problem->terminals.size());

        for (size_t i = 0; i < original_terminals.size(); ++i) {
            NodeID map = problem->mapped(original_terminals[i]);
            presets[problem->graph->getCurrentPosition(map)] = i;
        }

        auto [result, wgt] = ilp_model::computeIlp(problem->graph, presets,
                                                   original_terminals.size(),
                                                   problem->terminals.size());

        problem->upper_bound = problem->deleted_weight + wgt;


        if (problem->upper_bound < global_upper_bound) {
            updateBestSolution(problem);
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
    std::pair<NodeID, EdgeID> priority_edge;

    // parallel
    std::vector<std::condition_variable> q_cv;
    std::vector<std::mutex> q_mutex;
    std::atomic<uint> idle_threads;
    size_t num_threads;
    std::vector<size_t> branch_invalid;
    bool is_finished;

    kernelization_criteria kc;
    std::atomic<double> log_timer;
    std::mutex bestsol_mutex;

    std::atomic<size_t> graphs = 0;
};
