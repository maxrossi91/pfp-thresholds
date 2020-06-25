add_executable(pfp_thresholds pfp_thresholds.cpp)
target_link_libraries(pfp_thresholds common pfp gsacak sdsl malloc_count)

add_executable(sdsl_thresholds sdsl_thresholds.cpp)
target_link_libraries(sdsl_thresholds common gsacak sdsl malloc_count divsufsort divsufsort64)

add_executable(bwt_lcp_thresholds bwt_lcp_thresholds.cpp)
target_link_libraries(bwt_lcp_thresholds common gsacak sdsl malloc_count divsufsort divsufsort64)
target_include_directories(bwt_lcp_thresholds PUBLIC "${rlbwt2lcp_SOURCE_DIR}/internal")