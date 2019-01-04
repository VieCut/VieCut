/******************************************************************************
 * definitions.h
 *
 * Source of VieCut.
 * 
 * Adapted from KaHIP.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@univie.ac.at>
 * Copyright (C) 2017-2018 Alexander Noe <alexander.noe@univie.ac.at>
 *
 *****************************************************************************/

#pragma once

#include <limits>
#include <queue>
#include <vector>

#include "limits.h"
#include "stdio.h"
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
typedef unsigned int NodeID;
typedef double EdgeRatingType;
typedef unsigned long int EdgeID;
typedef unsigned long int PathID;
typedef unsigned int PartitionID;
typedef unsigned int NodeWeight;
typedef unsigned long int EdgeWeight;
typedef EdgeWeight Gain;
typedef int Color;
typedef unsigned long int Count;
typedef std::vector<NodeID> boundary_starting_nodes;
typedef long FlowType;

const EdgeID UNDEFINED_EDGE = std::numeric_limits<EdgeID>::max();
const EdgeID NOTMAPPED = std::numeric_limits<EdgeID>::max();
const NodeID UNDEFINED_NODE = std::numeric_limits<NodeID>::max();
const NodeID UNASSIGNED = std::numeric_limits<NodeID>::max();
const NodeID ASSIGNED = std::numeric_limits<NodeID>::max() - 1;
const PartitionID INVALID_PARTITION = std::numeric_limits<PartitionID>::max();
const PartitionID BOUNDARY_STRIPE_NODE = std::numeric_limits<PartitionID>::max();
const Count UNDEFINED_COUNT = std::numeric_limits<Count>::max();
const int NOTINQUEUE = std::numeric_limits<int>::max();
const int ROOT = 0;

// for the gpa global_mincut
struct edge_source_pair {
    EdgeID e;
    NodeID source;
};

struct source_target_pair {
    NodeID source;
    NodeID target;
};

// matching array has size (no_of_nodes), so for entry in this table we get the matched neighbor
typedef std::vector<NodeID> CoarseMapping;
typedef std::vector<NodeID> Matching;
typedef std::vector<NodeID> NodePermutationMap;

typedef double ImbalanceType;

typedef enum {
    PERMUTATION_QUALITY_NONE,
    PERMUTATION_QUALITY_FAST,
    PERMUTATION_QUALITY_GOOD
} PermutationQuality;

