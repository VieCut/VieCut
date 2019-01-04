/******************************************************************************
 * string.h
 *
 * Source of VieCut.
 *
 ******************************************************************************
 * Copyright (C) 2017 Alexander Noe <alexander.noe@univie.ac.at>
 * 
 * Published under the MIT license in the LICENSE file.
 *****************************************************************************/
#pragma once

#include <string>

class string_tools
{
public:
    static std::string basename(std::string filename) {
        if (filename.find_last_of("\\/") == std::string::npos) {
            return filename;
        }
        else {
            return filename.substr(filename.find_last_of("\\/") + 1);
        }
    }
};
