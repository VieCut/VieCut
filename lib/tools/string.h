/*******************************************************************************
 * string.h
 *
 * Source of VieCut.
 *
 * Copyright (C) 2015 Timo Bingmann <tb@panthema.net>
 * Copyright (C) 2017-2018 Alexander Noe <alexander.noe@univie.ac.at>
 ******************************************************************************/

#pragma once

#include <data_structure/graph_access.h>
#include <definitions.h>
#include <memory>
#include <vector>

class string
{
public:
//! Logging helper to print vectors as [a1,a2,a3,...]
    template <typename T>
    static std::string VecToStr(const std::vector<T>& data) {
        std::ostringstream oss;
        oss << '[';
        for (typename std::vector<T>::const_iterator it = data.begin();
             it != data.end(); ++it)
        {
            if (it != data.begin()) oss << ',';
            oss << *it;
        }
        oss << ']';
        return oss.str();
    }

    template <typename T>
    static typename std::enable_if<std::is_integral<T>::value, std::string>::type
    NonZeroVecToStr(const std::vector<T>& data) {
        std::ostringstream oss;
        oss << '[';

        bool first = true;
        for (typename std::vector<T>::const_iterator it = data.begin();
             it != data.end(); ++it)
        {
            size_t i = it - data.begin();
            if (*it != 0) {
                if (first) {
                    first = false;
                }
                else {
                    oss << ',';
                }
                oss << i << ":";
                oss << *it;
            }
        }
        oss << ']';
        return oss.str();
    }

    template <typename T>
    static typename std::enable_if<!std::is_integral<T>::value, std::string>::type
    NonZeroVecToStr(const std::vector<T>& data) {
        return "This method only works for integral types. Use VecToStr()";
    }

    static std::string graphToString [[maybe_unused]] (std::shared_ptr<graph_access> G) {
        std::ostringstream oss;

        for (NodeID n : G->nodes()) {
            oss << n << ": [";
            for (EdgeID e : G->edges_of(n)) {
                oss << G->getEdgeTarget(e);
                if (e < G->get_first_invalid_edge(n) - 1) {
                    oss << ",";
                }
            }
            oss << "]\n";
        }
        return oss.str();
    }

    static std::string weightedGraphToString [[maybe_unused]] (std::shared_ptr<graph_access> G) {
        std::ostringstream oss;

        for (NodeID n : G->nodes()) {
            oss << n << ": [";
            for (EdgeID e : G->edges_of(n)) {
                oss << G->getEdgeTarget(e) << "(" << G->getEdgeWeight(e) << ")";
                if (e < G->get_first_invalid_edge(n) - 1) {
                    oss << ",";
                }
            }
            oss << "]\n";
        }
        return oss.str();
    }

    static std::string basename(std::string filename) {
        if (filename.find_last_of("\\/") == std::string::npos) {
            return filename;
        }
        else {
            return filename.substr(filename.find_last_of("\\/") + 1);
        }
    }
};
