add_executable(fasta2plain fasta2plain.cpp)
target_link_libraries(fasta2plain common sdsl malloc_count)

add_executable(sample_lcp sample_lcp.cpp)
target_link_libraries(sample_lcp common sdsl malloc_count)
