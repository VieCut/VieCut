cmake_minimum_required(VERSION 3.9)
project(vc)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

OPTION(RUN_TESTS "Compile and run tests" ON)
OPTION(USE_TCMALLOC "Replace builtin malloc with TCMalloc" ON)
OPTION(USE_PROFILER "Use TCMalloc profiler (requires TCMalloc)" OFF)
OPTION(USE_GUROBI "Use Gurobi for ILP Solving in multiterminal cut" OFF)

if (USE_TCMALLOC)
    find_package(Tcmalloc REQUIRED)
endif()

if (USE_PROFILER)
    find_package(Tcmalloc REQUIRED)
    find_package(Profiler REQUIRED)
endif()

if (USE_GUROBI)
    find_package(Gurobi REQUIRED)
    MESSAGE("Using Gurobi for ILP Solving in multiterminal cut...")
    add_definitions(-DUSE_GUROBI)
endif()

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

macro(bal_mpi TARGETNAME)
    set (SEQ_NAME "${TARGETNAME}") 
    add_executable(${SEQ_NAME} app/${TARGETNAME}.cpp)
    target_link_libraries(${SEQ_NAME} PUBLIC ${LIBS} OpenMP::OpenMP_CXX ${MPI_LIBRARIES})
endmacro(bal_mpi)

macro(build_and_link TARGETNAME)
    bal_seq(${TARGETNAME})
    bal_par(${TARGETNAME})
endmacro(build_and_link)

if(NOT CMAKE_BUILD_TYPE)
  message("Defaulting CMAKE_BUILD_TYPE to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()


add_definitions("-DVIECUT_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}\"")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# prohibit in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
    message(SEND_ERROR "In-source builds are not allowed.")
endif()

MESSAGE(STATUS "Option: RUN_TESTS " ${RUN_TESTS})
MESSAGE(STATUS "Option: USE_TCMALLOC " ${USE_TCMALLOC})
MESSAGE(STATUS "Option: USE_PROFILER " ${USE_PROFILER})
MESSAGE(STATUS "Option: USE_GUROBI " ${USE_GUROBI})

MESSAGE(STATUS ${GUROBI_INCLUDE_DIR})

add_subdirectory(extlib/tlx)

include_directories(
    app 
    lib             
    /usr/include/mpi
    extlib/growt
    extlib/tlx
    ${GUROBI_INCLUDE_DIR}
)


set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fopenmp")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lpthread")
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set (CMAKE_CXX_FLAGS_DEBUG "-g3 -march=native")

#TCMalloc
if (USE_TCMALLOC)
    add_definitions(-DUSE_TCMALLOC)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free")
endif()


set(EXTLIBS tlx)

if(USE_PROFILER)
set(LIBS ${EXTLIBS} ${Profiler_LIBRARIES} ${GUROBI_LIBRARIES})
else()
set(LIBS ${EXTLIBS} ${Tcmalloc_LIBRARIES} ${GUROBI_LIBRARIES})
endif()

build_and_link(kcore)
build_and_link(mincut_contract)
build_and_link(mincut_heavy)
build_and_link(mincut_recursive)
build_and_link(mincut)
build_and_link(balanced_small_cut)
bal_mpi(multiterminal_cut)
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
