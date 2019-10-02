cmake_minimum_required(VERSION 3.9)
project(vc)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

find_package(Tcmalloc REQUIRED)
find_package(Profiler)
find_package(OpenMP REQUIRED)

macro(bal_seq TARGETNAME)
    set (SEQ_NAME "${TARGETNAME}") 
    add_executable(${SEQ_NAME} app/${TARGETNAME}.cpp)
    target_link_libraries(${SEQ_NAME} PUBLIC ${LIBS} OpenMP::OpenMP_CXX)
endmacro(bal_seq)

macro(bal_par TARGETNAME)
    set (PAR_NAME "${TARGETNAME}_parallel")
    add_executable(${PAR_NAME} app/${TARGETNAME}.cpp)
    target_link_libraries(${PAR_NAME} PUBLIC ${LIBS} OpenMP::OpenMP_CXX)
    target_compile_definitions(${PAR_NAME} PUBLIC -DPARALLEL)
endmacro(bal_par)

macro(build_and_link TARGETNAME)
    bal_seq(${TARGETNAME})
    bal_par(${TARGETNAME})
endmacro(build_and_link)

if(NOT CMAKE_BUILD_TYPE)
  message("Defaulting CMAKE_BUILD_TYPE to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()


add_definitions("-DVIECUT_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}\"")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

OPTION(RUN_TESTS "Compile and run tests" ON)

# prohibit in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(SEND_ERROR "In-source builds are not allowed.")
endif()

MESSAGE(STATUS "Option: RUN_TESTS " ${RUN_TESTS})

add_subdirectory(extlib/tlx)

include_directories(
    app 
    lib             
    /usr/include/mpi
    extlib/growt
    extlib/tlx
)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fopenmp -g")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread")
set (CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -DNDEBUG")
set (CMAKE_CXX_FLAGS_DEBUG "-g3 -march=native")

#TCMalloc
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free")

set(SOURCE_FILES

    lib/algorithms/flow/excess_scaling.h
    lib/algorithms/flow/push_relabel.h
    lib/algorithms/flow/unit_flow.h

    lib/algorithms/global_mincut/algorithms.h
    lib/algorithms/global_mincut/ks_minimum_cut.h
    lib/algorithms/global_mincut/matula_approx.h
    lib/algorithms/global_mincut/minimum_cut.h
    lib/algorithms/global_mincut/minimum_cut_helpers.h
    lib/algorithms/global_mincut/noi_minimum_cut.h
    lib/algorithms/global_mincut/padberg_rinaldi.h
    lib/algorithms/global_mincut/stoer_wagner_minimum_cut.h
    lib/algorithms/global_mincut/viecut.h

    lib/algorithms/global_mincut/cactus/balanced_cut_dfs.h
    lib/algorithms/global_mincut/cactus/cactus_mincut.h
    lib/algorithms/global_mincut/cactus/most_balanced_minimum_cut.h
    lib/algorithms/global_mincut/cactus/recursive_cactus.h
    
    lib/algorithms/misc/core_decomposition.h
    lib/algorithms/misc/graph_algorithms.h
    lib/algorithms/misc/strongly_connected_components.h

    lib/algorithms/multicut/branch_multicut.h
    lib/algorithms/multicut/edge_selection.h
    lib/algorithms/multicut/graph_contraction.h
    lib/algorithms/multicut/kernelization_criteria.h
    lib/algorithms/multicut/multicut_problem.h
    lib/algorithms/multicut/multiterminal_cut.h

    lib/algorithms/multicut/problem_queues/per_thread_problem_queue.h
    lib/algorithms/multicut/problem_queues/single_problem_queue.h

    lib/coarsening/contract_graph.h
    lib/coarsening/contraction_tests.h
    lib/coarsening/label_propagation.h
    lib/coarsening/sparsify.h
    lib/coarsening/test_wrapper.h
    
    lib/common/configuration.h
    lib/common/definitions.h

    lib/data_structure/adjlist_graph.h
    lib/data_structure/flow_graph.h
    lib/data_structure/graph_access.h
    lib/data_structure/mutable_graph.h
    lib/data_structure/union_find.h

    lib/data_structure/priority_queues/bucket_pq.h
    lib/data_structure/priority_queues/fifo_node_bucket_pq.h
    lib/data_structure/priority_queues/maxNodeHeap.h
    lib/data_structure/priority_queues/node_bucket_pq.h
    lib/data_structure/priority_queues/priority_queue_interface.h
    lib/data_structure/priority_queues/vecMaxNodeHeap.h

    lib/io/graph_io.cpp
    lib/io/graph_io.h

    lib/parallel/algorithm/exact_parallel_minimum_cut.h
    lib/parallel/algorithm/parallel_cactus.h

    lib/parallel/coarsening/contract_graph.h
    lib/parallel/coarsening/contraction_tests.h
    lib/parallel/coarsening/label_propagation.h
    lib/parallel/coarsening/sparsify.h

    lib/parallel/data_structure/union_find.h

    lib/tools/graph_extractor.h
    lib/tools/hash.h
    lib/tools/macros_assertions.h
    lib/tools/macros_common.h
    lib/tools/quality_metrics.h
    lib/tools/random_functions.cpp
    lib/tools/random_functions.h
    lib/tools/string.h
    lib/tools/timer.h
    lib/tools/vector.h 

)

add_library(viecut ${SOURCE_FILES})

set(EXTLIBS tlx)
if(USE_PROFILER)
set(LIBS viecut ${EXTLIBS} ${Profiler_LIBRARIES})
else()
set(LIBS viecut ${EXTLIBS} ${Tcmalloc_LIBRARIES})
endif()

build_and_link(kcore)
build_and_link(mincut_contract)
build_and_link(mincut_heavy)
build_and_link(mincut_recursive)
build_and_link(mincut)
build_and_link(balanced_small_cut)
bal_seq(multiterminal_cut)
bal_seq(largest_cc)
bal_seq(make_graph_weighted)
bal_seq(gnp_torus)

message(STATUS "CXX Flags ${CMAKE_CXX_FLAGS}")
message(STATUS "Linker flags ${CMAKE_EXE_LINKER_FLAGS}")

if (CMAKE_BUILD_TYPE STREQUAL Release)
    message (STATUS "Release build!")
    message (STATUS "Release CXX Flags: ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

if (CMAKE_BUILD_TYPE STREQUAL Debug)
    message (STATUS "Debug build!")
    message (STATUS "Debug CXX Flags: ${CMAKE_CXX_FLAGS_DEBUG}")
endif()

if(RUN_TESTS)
    enable_testing()
    add_subdirectory(extlib/googletest)
    add_subdirectory(tests)
endif(RUN_TESTS)
