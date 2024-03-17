set ( CACHING_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
message( "Software-assisted Memory Buffer - CACHING_DIR: " ${CACHING_DIR} )

set (CACHING_DIR_INC ${CACHING_DIR}/include)
include_directories( ${CACHING_DIR_INC})
message( "Software-assisted Memory Buffer - Included DIR: " ${CACHING_DIR_INC} )

include_directories(${PARALIA_INSTALL_INC})
message( "Software-assisted Memory Buffer - Included DIR: " ${PARALIA_INSTALL_INC} )

set (CACHING_DIR_SRC ${CACHING_DIR}/src)
message( "Software-assisted Memory Buffer - Source DIR: " ${CACHING_DIR_SRC} )

## Use continuous allocation for buffer memory blocks in all locations
add_definitions(-DENABLE_BUFFER_CONTINUOUS_ALLOC)

add_library(caching SHARED
${CACHING_DIR_SRC}/DataCaching.cpp
)

target_link_libraries(caching ${INTERNAL_LINK_LIB} chl_grid_amalgamation)
configure_file(${CACHING_DIR_INC}/DataCaching.hpp ${PARALIA_INSTALL_INC}/DataCaching.hpp COPYONLY)
message( "------------------------------------------------------------------------------------------------")
