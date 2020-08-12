
add_executable(pfp_ms_benchmarks pfp_ms_benchmarks.cpp)
target_link_libraries(pfp_ms_benchmarks common pfp ri gsacak sdsl divsufsort divsufsort64 malloc_count benchmark pthread)
target_include_directories(pfp_ms_benchmarks PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_include_directories(pfp_ms_benchmarks PUBLIC "../../include/ms")
target_include_directories(pfp_ms_benchmarks PUBLIC "../../include/pfp")
target_compile_options(pfp_ms_benchmarks PUBLIC "-std=c++14")

# target_compile_options(ds_building_wt_sdsl_264 PUBLIC -DM64)