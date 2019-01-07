cmake_minimum_required(VERSION 3.5)
project(vc)

macro(build_and_link TARGETNAME)
    add_executable(${TARGETNAME} app/${TARGETNAME}.cpp)
    target_link_libraries(${TARGETNAME} ${LIBS})
endmacro(build_and_link)

add_definitions("-DVIECUT_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}\"")

set(CMAKE_CXX_STANDARD 17)

OPTION(SAVE_CUTS "Save minimum cut to file" OFF)
OPTION(ENABLE_PARALLEL "Shared-memory parallelism enabled" OFF)
OPTION(RUN_TESTS "Compile and run tests" ON)

# prohibit in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(SEND_ERROR "In-source builds are not allowed.")
endif()

IF (ENABLE_PARALLEL)
    add_definitions(-DPARALLEL)
ENDIF(ENABLE_PARALLEL)


IF(SAVE_CUTS)
    add_definitions(-DSAVECUT)
ENDIF(SAVE_CUTS)

MESSAGE(STATUS "Option: SAVE_CUTS " ${SAVE_CUTS})
MESSAGE(STATUS "Option: ENABLE_PARALLEL " ${ENABLE_PARALLEL})
MESSAGE(STATUS "Option: RUN_TESTS " ${RUN_TESTS})

set(default_build_type "Release")

include_directories(
        app
        extlib/growt
        extlib/tlx
        lib
        /usr/include/mpi)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pipe -msse4.2 -mcx16 -Wall -Wextra")

set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

set (CMAKE_CXX_FLAGS_DEBUG "-g3")

add_subdirectory(extlib/tlx)

set(SOURCE_FILES
        lib/algorithms/global_mincut/algorithms.h
        lib/algorithms/global_mincut/ks_minimum_cut.h
        lib/algorithms/global_mincut/matula_approx.h
        lib/algorithms/global_mincut/viecut.h
        lib/algorithms/global_mincut/minimum_cut.h
        lib/algorithms/global_mincut/minimum_cut_helpers.h
        lib/algorithms/global_mincut/noi_minimum_cut.h
        lib/algorithms/global_mincut/padberg_rinaldi.h
        lib/algorithms/global_mincut/stoer_wagner_minimum_cut.h
        lib/algorithms/misc/strongly_connected_components.h
        lib/coarsening/contract_graph.h
        lib/coarsening/contraction_tests.h
        lib/algorithms/flow/excess_scaling.h
        lib/coarsening/label_propagation.h
        lib/coarsening/sparsify.h
        lib/coarsening/test_wrapper.h
        lib/algorithms/flow/unit_flow.h
        lib/tools/string.h
        lib/data_structure/adjlist_graph.h
        lib/data_structure/graph_access.h
        lib/io/graph_io_adj.h
        lib/parallel/algorithm/exact_parallel_minimum_cut.h
        lib/parallel/coarsening/contract_graph.h
        lib/parallel/coarsening/contraction_tests.h
        lib/parallel/coarsening/label_propagation.h
        lib/parallel/coarsening/sparsify.h
        lib/parallel/data_structure/union_find.h
        lib/algorithms/global_mincut/minimum_cut.h
        lib/tools/random_functions.cpp
        lib/tools/graph_extractor.h
        lib/tools/quality_metrics.h
        lib/io/graph_io.cpp
        lib/definitions.h
        lib/algorithms/misc/core_decomposition.h
        lib/algorithms/flow/push_relabel.h
        lib/data_structure/flow_graph.h
        lib/tools/vector.h 
        lib/algorithms/misc/graph_algorithms.h)

SET(ARGTABLE2_SEARCH_PATHS
        /usr/local
        /usr
        /opt
        )

FIND_PACKAGE(OpenMP REQUIRED)
#include_directories(${OPENMP_INCLUDE_DIRS})

add_library(viecut ${SOURCE_FILES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

set(EXTLIBS tlx)
set(LIBS viecut ${EXTLIBS})

build_and_link(kcore)
build_and_link(mincut_contract)
build_and_link(mincut_heavy)
build_and_link(mincut_recursive)
build_and_link(mincut)


if(RUN_TESTS)
    enable_testing()
    add_subdirectory(extlib/googletest)
    add_subdirectory(tests)
endif(RUN_TESTS)