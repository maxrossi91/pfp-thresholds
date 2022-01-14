add_executable(pfp_thresholds pfp_thresholds.cpp)
target_link_libraries(pfp_thresholds common pfp gsacak sdsl malloc_count)

add_executable(pfp_thresholds64 pfp_thresholds.cpp)
target_link_libraries(pfp_thresholds64 common pfp gsacak64 sdsl malloc_count)
target_compile_options(pfp_thresholds64 PUBLIC -DM64)

add_executable(pfp_lcp pfp_lcp.cpp)
target_link_libraries(pfp_lcp common pfp gsacak sdsl malloc_count)

add_executable(pfp_lcp64 pfp_lcp.cpp)
target_link_libraries(pfp_lcp64 common pfp gsacak64 sdsl malloc_count)
target_compile_options(pfp_lcp64 PUBLIC -DM64)

add_executable(sdsl_thresholds sdsl_thresholds.cpp)
target_link_libraries(sdsl_thresholds common gsacak sdsl malloc_count divsufsort divsufsort64)

if(${COMPILE_BWT2LCP})
    add_executable(bwt_lcp_thresholds bwt_lcp_thresholds.cpp)
    target_link_libraries(bwt_lcp_thresholds common gsacak sdsl malloc_count divsufsort divsufsort64)
    target_include_directories(bwt_lcp_thresholds PUBLIC "${rlbwt2lcp_SOURCE_DIR}/internal")

    add_executable(bwt_lcp bwt_lcp.cpp)
    target_link_libraries(bwt_lcp common gsacak sdsl malloc_count divsufsort divsufsort64)
    target_include_directories(bwt_lcp PUBLIC "${bwt2lcp_SOURCE_DIR}/internal")
endif()

if(${COMPILE_GSACAK})
    add_executable(gsacak_thresholds gsacak_thresholds.cpp)
    target_link_libraries(gsacak_thresholds common gsacak sdsl malloc_count)

    add_executable(gsacak_thresholds64 gsacak_thresholds.cpp)
    target_link_libraries(gsacak_thresholds64 common gsacak64 sdsl malloc_count)
    target_compile_options(gsacak_thresholds64 PUBLIC -DM64)

    add_executable(gsacak_lcp gsacak_lcp.cpp)
    target_link_libraries(gsacak_lcp common gsacak sdsl malloc_count)

    add_executable(gsacak_lcp64 gsacak_lcp.cpp)
    target_link_libraries(gsacak_lcp64 common gsacak64 sdsl malloc_count)
    target_compile_options(gsacak_lcp64 PUBLIC -DM64)
endif()

FetchContent_GetProperties(r-index)
add_executable(matching_statistics_pointers matching_statistics_pointers.cpp)
target_link_libraries(matching_statistics_pointers common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(matching_statistics_pointers PUBLIC "../../include/ms")
# target_include_directories(matching_statistics_pointers PUBLIC "${r-index_SOURCE_DIR}/internal")

add_executable(matching_statistics matching_statistics.cpp)
target_link_libraries(matching_statistics common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(matching_statistics PUBLIC "../../include/ms")
target_include_directories(matching_statistics PUBLIC "../../include/pfp")

add_executable(pfp_ms_build_only pfp_ms_build_only.cpp)
target_link_libraries(pfp_ms_build_only common sdsl divsufsort divsufsort64 malloc_count ri)
target_include_directories(pfp_ms_build_only PUBLIC "../../include/ms")
target_include_directories(pfp_ms_build_only PUBLIC "../../include/pfp")

add_executable(sdsl_matching_statistics sdsl_matching_statistics.cpp)
target_link_libraries(sdsl_matching_statistics common sdsl divsufsort divsufsort64 malloc_count)

if(${COMPILE_TESTS})
    add_executable(pfp_ms_test pfp_ms_test.cpp)
    target_link_libraries(pfp_ms_test common gsacak sdsl divsufsort divsufsort64 malloc_count ri benchmark pthread gtest_main)
    target_include_directories(pfp_ms_test PUBLIC "../../include/ms")
    target_include_directories(pfp_ms_test PUBLIC "../../include/pfp")
    target_include_directories(pfp_ms_test PUBLIC "../../benchmarks/src")

    add_executable(ms_rle_string_test ms_rle_string_test.cpp)
    target_link_libraries(ms_rle_string_test common gsacak sdsl divsufsort divsufsort64 malloc_count ri benchmark pthread gtest_main)
    target_include_directories(ms_rle_string_test PUBLIC "../../include/ms")
    target_include_directories(ms_rle_string_test PUBLIC "../../include/pfp")
    target_include_directories(ms_rle_string_test PUBLIC "../../benchmarks/src")
endif()