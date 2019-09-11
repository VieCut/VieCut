/******************************************************************************
 * configuration.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 *****************************************************************************/
#pragma once

#include <memory>
#include <string>
#include <vector>

class configuration {
 public:
    configuration(configuration const&) = delete;
    void operator = (configuration const&) = delete;

    static std::shared_ptr<configuration> getConfig() {
        static std::shared_ptr<configuration> instance{ new configuration };
        return instance;
    }

    ~configuration() { }

    // Settings - these are public for ease of use
    // don't change in program when not necessary
    std::string graph_filename;
    std::string partition_file = "";
    std::string output_path = "";
    size_t seed = 0;
    bool verbose = false;

    // multiterminal cut parameters
    std::string edge_selection = "heavy_vertex";
    std::string queue_type = "bound_sum";
    std::vector<std::string> term_strings;
    int top_k = 0;
    int random_k = 0;
    size_t bfs_size = 0;
    size_t threads = 1;
    size_t contraction_type = 0;
    size_t preset_percentage = 0;
    size_t total_terminals = 0;
    bool noBranching = false;
    size_t print_cc = 0;

    // minimum cut parameters
    bool save_cut = false;
    std::string algorithm;
    std::string sampling_type = "geometric";
    std::string pq = "default";
    size_t num_iterations = 1;
    bool disable_limiting = false;
    double contraction_factor = 0.0;
    bool find_most_balanced_cut = false;

    // karger-stein:
    size_t optimal = 0;

 private:
    configuration() { }

    // configuration(configuration const&);
    // void operator=(configuration const&);
};
