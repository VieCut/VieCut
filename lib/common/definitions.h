/******************************************************************************
 * definitions.h
 *
 * Source of VieCut.
 *
 * Adapted from KaHIP.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@univie.ac.at>
 * Copyright (C) 2017-2019 Alexander Noe <alexander.noe@univie.ac.at>
 *
 *****************************************************************************/

#pragma once

#include <limits>
#include <queue>
#include <vector>

#include "tools/macros_assertions.h"

// allows us to disable most of the output during partitioning
#ifdef KAFFPAOUTPUT
        #define PRINT(x) x
#else
        #define PRINT(x) do { } while (false);
#endif

/**********************************************
 * Constants
 * ********************************************/
// Types needed for the graph ds
typedef uint32_t NodeID;
typedef double EdgeRatingType;
typedef uint64_t EdgeID;
typedef uint64_t PathID;
typedef uint32_t PartitionID;
typedef uint32_t NodeWeight;
typedef uint64_t EdgeWeight;
typedef NodeWeight Gain;
typedef int32_t Color;
typedef uint64_t Count;
typedef int64_t FlowType;

const EdgeID UNDEFINED_EDGE = std::numeric_limits<EdgeID>::max();
const EdgeID NOTMAPPED = std::numeric_limits<EdgeID>::max();
const NodeID UNDEFINED_NODE = std::numeric_limits<NodeID>::max();
const NodeID UNASSIGNED = std::numeric_limits<NodeID>::max();
const NodeID ASSIGNED = std::numeric_limits<NodeID>::max() - 1;
const Count UNDEFINED_COUNT = std::numeric_limits<Count>::max();
const int NOTINQUEUE = std::numeric_limits<int>::max();
const int ROOT = 0;

typedef enum {
    UNDISCOVERED,
    ACTIVE,
    CYCLE,
    FINISHED
} DFSVertexStatus;

enum MessageStatus { needProblem, haveProblem, emptyAsWell, allEmpty };
