
macro(build_and_test TESTNAME PARALLEL)

set(COMPLETENAME ${TESTNAME})

if (${PARALLEL} STREQUAL "TRUE")
    set(COMPLETENAME ${COMPLETENAME}_par)
endif()

add_executable(${COMPLETENAME} ${TESTNAME}.cpp)
target_link_libraries(${COMPLETENAME} ${TESTLIBS} ${Tcmalloc_LIBRARIES})
if (${PARALLEL} STREQUAL "TRUE")
    target_compile_definitions(${COMPLETENAME} PUBLIC -DPARALLEL)
endif()
run_test(${COMPLETENAME})
endmacro(build_and_test)

macro(run_test TESTNAME)
add_test(
        NAME ${TESTNAME}
	COMMAND ${TESTNAME}
)
set_tests_properties(${TESTNAME} PROPERTIES ENVIRONMENT OMP_NUM_THREADS=4)
endmacro(run_test)

set(TESTLIBS ${LIBS} gtest gtest_main)

build_and_test(graph_test FALSE)
build_and_test(mutable_graph_test FALSE)
build_and_test(pq_test FALSE)
build_and_test(union_find_test FALSE)
build_and_test(union_find_test TRUE)
build_and_test(contraction_test FALSE)
build_and_test(contraction_test TRUE)
build_and_test(mincut_algo_test TRUE)
build_and_test(mincut_algo_test FALSE)
build_and_test(core_decomposition_test FALSE)
build_and_test(save_cut_test FALSE)
build_and_test(save_cut_test TRUE)
build_and_test(clique_test FALSE)
build_and_test(flow_graph_test FALSE)
build_and_test(push_relabel_test FALSE)
build_and_test(multiterminal_cut_test FALSE)
build_and_test(cactus_cut_test FALSE)
build_and_test(cactus_cut_test TRUE)

target_link_libraries(multiterminal_cut_test -lpthread ${MPI_LIBRARIES})

