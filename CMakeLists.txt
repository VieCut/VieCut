cmake_minimum_required(VERSION 3.9)
project(vc)

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

INCLUDE(CheckCXXCompilerFlag)

check_cxx_compiler_flag("-std=c++17" HAS_MODERN_COMPILER)
MESSAGE(Status: " Compiler supports c++17: " ${HAS_MODERN_COMPILER})

IF (HAS_MODERN_COMPILER)
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -std=c++17)
ELSE (HAS_MODERN_COMPILER)
    set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -std=c++14)
ENDIF (HAS_MODERN_COMPILER)

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


set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -fopenmp")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread -ltcmalloc")
set (CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -g -fno-omit-frame-pointer -DNDEBUG -pthread")

set (CMAKE_CXX_FLAGS_DEBUG "-g3 -march=native")

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
        lib/common/definitions.h
        lib/algorithms/misc/core_decomposition.h
        lib/algorithms/flow/push_relabel.h
        lib/data_structure/flow_graph.h
        lib/tools/vector.h 
        lib/algorithms/misc/graph_algorithms.h)

add_library(viecut ${SOURCE_FILES})

set(EXTLIBS tlx)
set(LIBS viecut ${EXTLIBS})

build_and_link(kcore)
build_and_link(mincut_contract)
build_and_link(mincut_heavy)
build_and_link(mincut_recursive)
build_and_link(mincut)
bal_seq(multiterminal_cut)
bal_seq(largest_cc)

target_link_libraries(mincut PUBLIC /usr/lib/x86_64-linux-gnu/libtcmalloc_and_profiler.so)
target_link_libraries(mincut_parallel PUBLIC /usr/lib/x86_64-linux-gnu/libtcmalloc_and_profiler.so)
target_link_libraries(multiterminal_cut PUBLIC /usr/lib/x86_64-linux-gnu/libtcmalloc_and_profiler.so)

if(RUN_TESTS)
    enable_testing()
    add_subdirectory(extlib/googletest)
    add_subdirectory(tests)
endif(RUN_TESTS)
