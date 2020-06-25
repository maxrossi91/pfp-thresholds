add_executable(pfp_thresholds pfp_thresholds.cpp)
target_link_libraries(pfp_thresholds common strdup pfp gsacak sdsl malloc_count)

add_executable(sdsl_thresholds sdsl_thresholds.cpp)
target_link_libraries(sdsl_thresholds common strdup gsacak sdsl malloc_count divsufsort divsufsort64)