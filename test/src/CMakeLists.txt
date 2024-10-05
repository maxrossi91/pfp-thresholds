add_executable(pfp_thresholds pfp_thresholds.cpp)
target_link_libraries(pfp_thresholds common pfp gsacak sdsl malloc_count z)

add_executable(pfp_thresholds64 pfp_thresholds.cpp)
target_link_libraries(pfp_thresholds64 common pfp gsacak64 sdsl malloc_count z)
target_compile_options(pfp_thresholds64 PUBLIC -DM64)

add_executable(pfp_lcp pfp_lcp.cpp)
target_link_libraries(pfp_lcp common pfp gsacak sdsl malloc_count z)

add_executable(pfp_lcp64 pfp_lcp.cpp)
target_link_libraries(pfp_lcp64 common pfp gsacak64 sdsl malloc_count z)
target_compile_options(pfp_lcp64 PUBLIC -DM64)

add_executable(sdsl_thresholds sdsl_thresholds.cpp)
target_link_libraries(sdsl_thresholds common gsacak sdsl malloc_count divsufsort divsufsort64 z)

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

