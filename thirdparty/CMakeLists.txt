# Thirdparty libraries
include(FetchContent)

# Log directory for external content builds
make_directory(${CMAKE_CURRENT_BINARY_DIR}/logs)

# Set environment

list(APPEND CMAKE_LIBRARY_PATH "${CMAKE_CURRENT_BINARY_DIR}/lib ")
list(APPEND CMAKE_INCLUDE_PATH "${CMAKE_CURRENT_BINARY_DIR}/include ")
list(APPEND CMAKE_PREFIX_PATH  "${CMAKE_CURRENT_BINARY_DIR}")

## Add malloc_count
FetchContent_Declare(
  malloc_count
  GIT_REPOSITORY https://github.com/bingmann/malloc_count
  )
  
FetchContent_GetProperties(malloc_count)
if(NOT malloc_count_POPULATED)
  FetchContent_Populate(malloc_count)

  add_library(malloc_count OBJECT ${malloc_count_SOURCE_DIR}/malloc_count.c ${malloc_count_SOURCE_DIR}/malloc_count.h)
  target_link_libraries(malloc_count dl)
  target_include_directories(malloc_count PUBLIC "${malloc_count_SOURCE_DIR}")

  add_library(memprofile OBJECT ${malloc_count_SOURCE_DIR}/memprofile.h)
  target_include_directories(memprofile PUBLIC "${malloc_count_SOURCE_DIR}")
endif()


## Add Big-BWT
FetchContent_Declare(
  bigbwt
  GIT_REPOSITORY https://github.com/alshai/Big-BWT.git
  )
  
FetchContent_GetProperties(bigbwt)
if(NOT bigbwt_POPULATED)
  FetchContent_Populate(bigbwt)
  add_subdirectory(${bigbwt_SOURCE_DIR} ${bigbwt_BINARY_DIR})

endif()
  
## Add gsacak
FetchContent_Declare(
  gsacak
  GIT_REPOSITORY https://github.com/felipelouza/gsa-is.git
)

FetchContent_GetProperties(gsacak)
if(NOT gsacak_POPULATED)
    FetchContent_Populate(gsacak)
    add_library(gsacak OBJECT ${gsacak_SOURCE_DIR}/gsacak.c ${gsacak_SOURCE_DIR}/gsacak.h)
    target_include_directories(gsacak PUBLIC "${gsacak_SOURCE_DIR}")
    
    add_library(gsacak64 OBJECT ${gsacak_SOURCE_DIR}/gsacak.c ${gsacak_SOURCE_DIR}/gsacak.h)
    target_include_directories(gsacak64 PUBLIC "${gsacak_SOURCE_DIR}")
    target_compile_options(gsacak64 PUBLIC -DM64)
endif()



# DIVSUFSORT
find_library(DIVSUFSORT_LIB divsufsort HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
if(NOT DIVSUFSORT_LIB)
    message(STATUS "divsufsort library not found. Building as an external content")
else()
    message(STATUS "divsufsort library found at ${DIVSUFSORT_LIB}.")
endif()

# DIVSUFSORT64
find_library(DIVSUFSORT64_LIB divsufsort64 HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
if(NOT DIVSUFSORT64_LIB)
    message(STATUS "divsufsort64 library not found. Building as an external content")
else()
    message(STATUS "divsufsort64 library found at ${DIVSUFSORT64_LIB}.")
endif()

# SDSL
find_library(SDSL_LIB sdsl HINTS ${CMAKE_CURRENT_BINARY_DIR}/lib PATHS ${CONDA_PREFIX}/lib)
find_path(SDSL_SRC sdsl HINTS ${CMAKE_CURRENT_BINARY_DIR}/include PATHS ${CONDA_PREFIX}/include)
if(NOT SDSL_LIB)
    message(STATUS "sdsl library not found. Building as an external content")
    FetchContent_Declare(sdsl
            GIT_REPOSITORY "https://github.com/simongog/sdsl-lite.git"
            GIT_TAG "v2.1.1"
    )

    if(NOT sdsl_POPULATED)
        FetchContent_Populate(sdsl)
        # Install in the current binary dir
        execute_process(COMMAND ./install.sh ${CMAKE_CURRENT_BINARY_DIR}
                RESULT_VARIABLE result
                WORKING_DIRECTORY ${sdsl_SOURCE_DIR} 
                OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log
                ERROR_FILE ${CMAKE_CURRENT_BINARY_DIR}/logs/sdsl.log)
        if(result)
            message(FATAL_ERROR "CMake step for sdsl-lite failed: ${result}")
        endif()
    endif()

    set(SDSL_SRC ${CMAKE_CURRENT_BINARY_DIR}/include)
    set(DIVSUFSORT_SRC ${CMAKE_CURRENT_BINARY_DIR}/include)
    set(SDSL_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libsdsl.a)
    set(DIVSUFSORT_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libdivsufsort.a)
    set(DIVSUFSORT64_LIB ${CMAKE_CURRENT_BINARY_DIR}/lib/libdivsufsort64.a)

else()
    message(STATUS "sdsl library found at ${SDSL_LIB}.")
    message(STATUS "sdsl sources found at ${SDSL_SRC}.")
endif()

add_library(sdsl STATIC IMPORTED GLOBAL)
set_property(TARGET sdsl PROPERTY IMPORTED_LOCATION ${SDSL_LIB})
target_include_directories(sdsl INTERFACE ${SDSL_SRC})

add_library(divsufsort STATIC IMPORTED GLOBAL)
set_property(TARGET divsufsort PROPERTY IMPORTED_LOCATION ${DIVSUFSORT_LIB})
    
add_library(divsufsort64 STATIC IMPORTED GLOBAL)
set_property(TARGET divsufsort64 PROPERTY IMPORTED_LOCATION ${DIVSUFSORT64_LIB})


if(${COMPILE_BENCHMARKS})

  ## Google benchmark
  FetchContent_Declare(
    benchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG        main
    )
    
  FetchContent_GetProperties(benchmark)
  if(NOT benchmark_POPULATED)
    FetchContent_Populate(benchmark)
    set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Enable testing of the benchmark library." FORCE)
    add_subdirectory(${benchmark_SOURCE_DIR} ${benchmark_BINARY_DIR} EXCLUDE_FROM_ALL)
  endif()

endif()



if(${COMPILE_BWT2LCP})
  ## Add rlbwt2lcp
  FetchContent_Declare(
    rlbwt2lcp
    GIT_REPOSITORY https://github.com/nicolaprezza/rlbwt2lcp
  )

  FetchContent_GetProperties(rlbwt2lcp)
  if(NOT rlbwt2lcp_POPULATED)
    FetchContent_Populate(rlbwt2lcp)

    add_subdirectory(${rlbwt2lcp_SOURCE_DIR} ${rlbwt2lcp_BINARY_DIR} )
  endif()

  ## Add bwt2lcp
  FetchContent_Declare(
    bwt2lcp
    GIT_REPOSITORY https://github.com/nicolaprezza/bwt2lcp
  )

  FetchContent_GetProperties(bwt2lcp)
  if(NOT bwt2lcp_POPULATED)
    FetchContent_Populate(bwt2lcp)

    add_subdirectory(${bwt2lcp_SOURCE_DIR} ${bwt2lcp_BINARY_DIR} )
  endif()
endif()
