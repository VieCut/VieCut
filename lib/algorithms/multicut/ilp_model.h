/******************************************************************************
 * ilp_model.h
 *
 * Source of VieCut
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexandra Henzinger <ahenz@stanford.edu>
 * Copyright (C) 2017-2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/

#pragma once

#include <memory>
#include <utility>
#include <vector>

#include "gurobi_c++.h" // NOLINT

#include "algorithms/multicut/multicut_problem.h"

class ilp_model {
 public:
    inline static GRBEnv env;

    static std::pair<std::vector<NodeID>, EdgeWeight> computeIlp(
        std::shared_ptr<multicut_problem> problem,
        const std::vector<NodeID>& presets,
        size_t num_terminals,
        bool parallel,
        size_t thread_id) {
        try {
            std::shared_ptr<mutable_graph> graph = problem->graph;
            timer ilp_timer;
            LOG1 << "starting ilp on graph with " << graph->n() << " vertices "
                 << "and " << graph->m() << " edges!";
            GRBModel model = GRBModel(ilp_model::env);
            NodeID numNodes = graph->number_of_nodes();
            NodeID numEdges = graph->number_of_edges() / 2;

            NodeID max_weight = 0;
            NodeID max_id = 0;
            for (size_t i = 0; i < problem->terminals.size(); ++i) {
                NodeID v = problem->terminals[i].position;
                if (graph->getWeightedNodeDegree(v) > max_weight) {
                    max_weight = graph->getWeightedNodeDegree(v);
                    max_id = i;
                }
            }

            std::vector<GRBVar> nodes(numNodes * num_terminals);
            std::vector<GRBVar> edges(numEdges);

            model.set(GRB_IntParam_Threads, 1);
            model.set(GRB_StringAttr_ModelName, "Partition");
            model.set(GRB_DoubleParam_MIPGap, 0);
            if (!parallel) {
                model.set(GRB_IntParam_Threads, 1);
                model.set(GRB_IntParam_LogToConsole, 0);
            } else {
                size_t threads = configuration::getConfig()->threads;
                model.set(GRB_IntParam_Threads, threads);
                LOG1 << "Running parallel ILP with "
                     << model.get(GRB_IntParam_Threads) << " threads";
                cpu_set_t all_cores;
                CPU_ZERO(&all_cores);
                for (size_t i = 0; i < threads; ++i) {
                    CPU_SET(i, &all_cores);
                }
                model.set(GRB_DoubleParam_Heuristics, 0.2);

                sched_setaffinity(0, sizeof(cpu_set_t), &all_cores);
            }
            model.set(GRB_IntParam_PoolSearchMode, 0);
            model.set(GRB_DoubleParam_TimeLimit, 3600.0);
            // Set decision variables for nodes
            for (size_t q = 0; q < num_terminals; q++) {
                GRBLinExpr nodeTot = 0;
                for (NodeID i = 0; i < numNodes; i++) {
                    if (presets[i] < num_terminals) {
                        bool isCurrent = (presets[i] == q);
                        double f = isCurrent ? 1.0 : 0.0;
                        nodes[i + q * numNodes] =
                            model.addVar(f, f, 0, GRB_BINARY);
                        nodes[i + q * numNodes].set(GRB_DoubleAttr_Start, f);
                    } else {
                        nodes[i + q * numNodes] =
                            model.addVar(0.0, 1.0, 0, GRB_BINARY);
                        nodes[i + q * numNodes].set(GRB_DoubleAttr_Start,
                                                    (max_id == q));
                    }
                }
            }

            size_t j = 0;
            // Decision variables for edges
            GRBLinExpr edgeTot = 0;
            for (NodeID n : graph->nodes()) {
                for (EdgeID e : graph->edges_of(n)) {
                    NodeID t = graph->getEdgeTarget(n, e);
                    if (n > graph->getEdgeTarget(n, e)) {
                        edges[j] = model.addVar(
                            0.0, 1.0, graph->getEdgeWeight(n, e), GRB_BINARY);
                        // there is no edge between terminals, we mark edges
                        // that are incident to non-maximal weight terminal
                        double start = ((presets[n] < num_terminals
                                         || presets[t] < num_terminals)
                                        && presets[n] != max_id
                                        && presets[t] != max_id)
                            ? 1.0 : 0.0;
                        edges[j].set(GRB_DoubleAttr_Start, start);
                        for (size_t q = 0; q < num_terminals; q++) {
                            GRBLinExpr cons = nodes[n + q * numNodes]
                                              - nodes[t + q * numNodes];
                            // Add constraint: valid partiton
                            model.addConstr(edges[j] >= cons, "valid part.");
                            model.addConstr(edges[j] >= -cons, "valid part.");
                        }
                        j++;
                    }
                }
            }

            // Add constraint: sum of all decision variables for 1 node is 1
            for (size_t i = 0; i < numNodes; i++) {
                GRBLinExpr sumCons = 0;
                for (size_t q = 0; q < num_terminals; q++) {
                    sumCons += nodes[i + q * numNodes];
                }
                model.addConstr(sumCons == 1);
            }

            model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
            // Optimize model
            model.optimize();

            std::vector<NodeID> result(numNodes);
            // if solution is found
            if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL ||
                model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
                // set partition
                for (PartitionID q = 0; q < num_terminals; q++) {
                    for (size_t i = 0; i < numNodes; i++) {
                        auto v = nodes[i + q * numNodes].get(GRB_DoubleAttr_X);
                        if (v == 1) {
                            result[i] = q;
                        }
                    }
                }
            } else {
                LOG1 << "No solution";
            }
            EdgeWeight wgt =
                static_cast<EdgeWeight>(model.get(GRB_DoubleAttr_ObjVal));

            LOG1 << "SOLUTION time=" << ilp_timer.elapsed()
                 << " graph=" << configuration::getConfig()->graph_filename
                 << " n=" << graph->n()
                 << " m=" << graph->m()
                 << " wgt=" << wgt
                 << " total=" << wgt + problem->deleted_weight
                 << " original_terminals=" << num_terminals
                 << " current_terminals=" << problem->terminals.size();

            if (parallel) {
                cpu_set_t my_id;
                CPU_ZERO(&my_id);
                CPU_SET(thread_id, &my_id);
                sched_setaffinity(0, sizeof(cpu_set_t), &my_id);
            }

            return std::make_pair(result, wgt);
        } catch (GRBException e) {
            LOG1 << e.getErrorCode() << " Message: " << e.getMessage();
            exit(1);
        }
    }
};
