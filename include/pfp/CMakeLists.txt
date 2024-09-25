set(PFP_SOURCES dictionary.hpp
                parse.hpp
                pfp.hpp)

if (NOT TARGET pfp)
    add_library(pfp OBJECT ${PFP_SOURCES})
    target_link_libraries(pfp common sdsl divsufsort divsufsort64)
endif()
