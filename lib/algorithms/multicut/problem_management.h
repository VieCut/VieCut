/******************************************************************************
 * problem_management.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2020 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/
#pragma once

#include <algorithm>
#include <chrono>
#include <limits>
#include <memory>
#include <unordered_set>
#include <vector>

#include "algorithms/multicut/measurements.h"
#include "algorithms/multicut/multicut_problem.h"
#include "algorithms/multicut/problem_queues/per_thread_problem_queue.h"
#include "algorithms/multicut/problem_queues/single_problem_queue.h"
#include "common/configuration.h"
#include "data_structure/mutable_graph.h"

using namespace std::chrono_literals;

class problem_management {
 private:
    static const bool testing = false;
#ifndef NDEBUG
    static const bool debug = true;
#else
    static const bool debug = false;
#endif

    const std::vector<NodeID>& original_terminals;
    const mutable_graph& original_graph;
    const std::vector<bool>& fixed_vertex;
    per_thread_problem_queue problems;
    maximum_flow mf;
    std::mutex bestsol_mutex;
    size_t num_threads;
    std::vector<std::mutex> q_mutex;
    std::vector<std::condition_variable> q_cv;
    bool is_finished;
    std::atomic<uint> idle_threads;

    FlowType global_upper_bound;
    std::vector<FlowType> terminalGUB;
    std::vector<FlowType> beforeLSGUB;

    std::vector<NodeID> best_solution;
    bool bestSolutionInitialized;
    measurements msm;

 public:
    problem_management(const mutable_graph& original_graph,
                       const std::vector<NodeID>& original_terminals,
                       const std::vector<bool>& fixed_vertex)
        : original_terminals(original_terminals),
          original_graph(original_graph),
          fixed_vertex(fixed_vertex),
          problems(configuration::getConfig()->threads,
                   configuration::getConfig()->queue_type),
          mf(original_terminals),
          num_threads(configuration::getConfig()->threads),
          q_mutex(configuration::getConfig()->threads),
          q_cv(configuration::getConfig()->threads),
          is_finished(false),
          idle_threads(0),
          global_upper_bound(UNDEFINED_FLOW),
          terminalGUB(original_terminals.size(), UNDEFINED_FLOW),
          beforeLSGUB(original_terminals.size(), UNDEFINED_FLOW),
          bestSolutionInitialized(false),
          msm(this->original_graph, this->original_terminals) {
        best_solution.resize(original_graph.number_of_nodes());
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

    std::optional<std::shared_ptr<multicut_problem> > pullProblem(
        size_t thread_id, bool send) {
        return problems.pullProblem(thread_id, send);
    }

    void branch(std::shared_ptr<multicut_problem> problem, size_t thread_id) {
        if (configuration::getConfig()->multibranch) {
            multiBranch(problem, thread_id);
        } else {
            singleBranch(problem, thread_id);
        }
    }

    std::optional<std::vector<NodeID> > getBestSolution() {
        if (!bestSolutionInitialized)
            return std::nullopt;
        bestsol_mutex.lock();
        std::vector<NodeID> ret = best_solution;
        bestsol_mutex.unlock();
        return ret;
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

    std::optional<FlowType> processNewProblem(
        std::shared_ptr<multicut_problem> new_p, size_t thread_id) {
        size_t numTerminals = new_p->terminals.size();
        if (numTerminals < 2 && runLocalSearch(new_p)) {
            auto sol = msm.getSolution(new_p);
            return findBestSolution(&sol, numTerminals);
        }

        if (numTerminals == 2) {
            mf.maximumSTFlow(new_p);
            if (runLocalSearch(new_p)) {
                auto sol = msm.getSolution(new_p);
                return findBestSolution(&sol, numTerminals);
            }
            return std::nullopt;
        }

        mf.maximumIsolatingFlow(new_p, thread_id, problems.size() == 0);
        graph_contraction::deleteTermEdges(new_p, original_terminals);
        if (new_p->graph->m() == 0) {
            new_p->upper_bound = new_p->deleted_weight;
            if (runLocalSearch(new_p)) {
                auto sol = msm.getSolution(new_p);
                return findBestSolution(&sol, numTerminals);
            }
            return std::nullopt;
        }

        if (checkProblem(new_p)) {
            std::vector<NodeID> sol;
            bool runLS = false;
            if (runLocalSearch(new_p)) {
                sol = msm.getSolution(new_p);
                runLS = true;
            }
            size_t thr = problems.addProblem(new_p, thread_id,
                                             runLocalSearch(new_p));
            q_cv[thr].notify_all();
            if (runLS) {
                return findBestSolution(&sol, numTerminals);
            }
        }
        return std::nullopt;
    }

    std::optional<FlowType> findBestSolution(
        std::vector<NodeID>* current_solution, NodeID numTerminals) {
        local_search ls(original_graph, original_terminals,
                        fixed_vertex, current_solution);
        FlowType prev_gub = msm.flowValue(false, *current_solution);
        if (prev_gub > beforeLSGUB[numTerminals])
            return std::nullopt;
        FlowType total_improvement = ls.improveSolution();
        FlowType ls_bound = prev_gub - total_improvement;

        if (ls_bound <= terminalGUB[numTerminals]) {
            terminalGUB[numTerminals] = ls_bound;
            beforeLSGUB[numTerminals] = prev_gub;
        }

        if (ls_bound <= global_upper_bound) {
            global_upper_bound = ls_bound;

            bestsol_mutex.lock();
            for (size_t i = 0; i < current_solution->size(); ++i) {
                best_solution[i] = (*current_solution)[i];
            }
            initalizeBestSolution();
            bestsol_mutex.unlock();
            return ls_bound;
        }
        return std::nullopt;
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

    bool allEmpty() {
        return problems.all_empty();
    }

    bool queueEmpty(size_t thread_id) {
        return problems.empty(thread_id);
    }

    bool haveASendProblem() {
        return problems.haveASendProblem();
    }

    size_t numProblems() {
        return problems.size();
    }

    void prepareQueue(size_t thread_id) {
        problems.prepareQueue(thread_id, global_upper_bound);
    }

    void addProblem(std::shared_ptr<multicut_problem> p,
                    size_t thread_id, bool preferLocal) {
        problems.addProblem(p, thread_id, preferLocal);
    }

    bool checkProblem(std::shared_ptr<multicut_problem> problem) {
        return problem->lower_bound < global_upper_bound;
    }

    bool runLocalSearch(std::shared_ptr<multicut_problem> problem) {
        return problem->upper_bound < beforeLSGUB[problem->terminals.size()];
    }

    void updateBound(FlowType newSolution) {
        global_upper_bound = std::min(newSolution, global_upper_bound);
    }

    FlowType bestCut() {
        return global_upper_bound;
    }

    void notifyThread(size_t thread_id) {
        q_cv[thread_id].notify_all();
    }

    void initalizeBestSolution() {
        bestSolutionInitialized = true;
    }

    void notifyAllThreads() {
        for (size_t i = 0; i < num_threads; ++i) {
            notifyThread(i);
        }
    }

    bool leaveWaitState(size_t thread_id) {
        return !queueEmpty(thread_id) || is_finished || (allThreadsIdle());
    }

    void waitForProblem(size_t thread_id) {
        std::unique_lock<std::mutex> lck(q_mutex[thread_id]);
        q_cv[thread_id].wait_for(
            lck, 1000ms,
            [this, thread_id] {
                return leaveWaitState(thread_id);
            });
    }

    bool allThreadsIdle() {
        return idle_threads == num_threads;
    }

    void incrementIdleThreads() {
        ++idle_threads;
    }

    void decrementIdleThreads() {
        --idle_threads;
    }

    bool finished() {
        return is_finished;
    }

    void setFinish() {
        is_finished = true;
    }
};
