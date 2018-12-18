
macro(build_and_test TESTNAME PARALLEL SAVECUT)

    set(COMPLETENAME ${TESTNAME})

    if (${PARALLEL} STREQUAL "TRUE")
        set(COMPLETENAME ${COMPLETENAME}_par)
    endif()
    if (${SAVECUT} STREQUAL "TRUE")
        set(COMPLETENAME ${COMPLETENAME}_save)
    endif()

    add_executable(${COMPLETENAME} ${TESTNAME}.cpp)
    target_link_libraries(${COMPLETENAME} ${TESTLIBS})
    if (${PARALLEL} STREQUAL "TRUE")
        target_compile_definitions(${COMPLETENAME} PUBLIC -DPARALLEL)
    endif()
    if (${SAVECUT} STREQUAL "TRUE")
        target_compile_definitions(${COMPLETENAME} PUBLIC -DSAVECUT)
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
remove_definitions(-DSAVECUT)

build_and_test(graph_test FALSE FALSE)
build_and_test(pq_test FALSE FALSE)
build_and_test(union_find_test FALSE FALSE)
build_and_test(union_find_test TRUE FALSE)
build_and_test(contraction_test FALSE FALSE)
build_and_test(contraction_test TRUE FALSE)
build_and_test(mincut_algo_test TRUE FALSE)
build_and_test(mincut_algo_test FALSE FALSE)
build_and_test(core_decomposition_test FALSE FALSE)
build_and_test(save_cut_test FALSE TRUE)
build_and_test(save_cut_test TRUE TRUE)
build_and_test(flow_graph_test FALSE FALSE)