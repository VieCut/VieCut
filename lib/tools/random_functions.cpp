/******************************************************************************
 * random_functions.cpp
 *
 * Source of VieCut.
 * 
 * Adapted from KaHIP.
 *
 ******************************************************************************
 * Copyright (C) 2013-2015 Christian Schulz <christian.schulz@kit.edu>
 * 
 *****************************************************************************/

#include "random_functions.h"

MersenneTwister random_functions::m_mt;
int random_functions::m_seed = 0;

random_functions::random_functions() { }

random_functions::~random_functions() { }
