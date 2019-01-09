
macro(build_and_test TESTNAME PARALLEL)

    set(COMPLETENAME ${TESTNAME})

    if (${PARALLEL} STREQUAL "TRUE")
        set(COMPLETENAME ${COMPLETENAME}_par)
    endif()

    add_executable(${COMPLETENAME} ${TESTNAME}.cpp)
    target_link_libraries(${COMPLETENAME} ${TESTLIBS})
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
endmacro(run_test)

set(TESTLIBS ${LIBS} gtest gtest_main)

#build viecut library with all parameters enabled and disabled for testing purposes
remove_definitions(-DPARALLEL)

build_and_test(graph_test FALSE)
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
build_and_test(flow_graph_test FALSE)