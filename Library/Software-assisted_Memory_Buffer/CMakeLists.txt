cmake_minimum_required (VERSION 3.10)
### Initialize the project
project (Software-assisted_Memory_Buffer VERSION ${PROJECT_VERSION} DESCRIPTION "Software-assisted Memory Buffer implementation for automated multi-GPU data movement" LANGUAGES CXX CUDA)
message( "------------------------------------------------------------------------------------------------")
message( "Software-assisted Memory Buffer - PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )
#message( "Library Scheduler - CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )
#message( "Library Scheduler - CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )

set ( CACHING_DIR ${PROJECT_SOURCE_DIR} )
message( "Software-assisted Memory Buffer - CACHING_DIR: " ${CACHING_DIR} )

set (CACHING_DIR_INC ${CACHING_DIR}/include)
include_directories( ${CACHING_DIR_INC})
message( "Software-assisted Memory Buffer - Included DIR: " ${CACHING_DIR_INC} )

set (CACHING_DIR_TEST ${PROJECT_SOURCE_DIR}/Testing)
#include_directories( ${CACHING_DIR_TEST})
message( "Software-assisted Memory Buffer - Included DIR: " ${CACHING_DIR_TEST} )

include_directories(${CHANGELINK_INSTALL_INC})
message( "Software-assisted Memory Buffer - Included DIR: " ${CHANGELINK_INSTALL_INC} )

include_directories(${LIBSC_INSTALL_INC})
message( "Software-assisted Memory Buffer - Included DIR: " ${LIBSC_INSTALL_INC} )

set (CACHING_DIR_SRC ${CACHING_DIR}/src)
message( "Software-assisted Memory Buffer - Source DIR: " ${CACHING_DIR_SRC} )

set (CACHING_INSTALL ${CMAKE_INSTALL_PROJECT}/Library_scheduler_exp/Software-assisted_Memory_Buffer)
file(MAKE_DIRECTORY ${CACHING_INSTALL})

set (CACHING_INSTALL_INC ${CACHING_INSTALL}/include)
file(MAKE_DIRECTORY ${CACHING_INSTALL_INC})

set (CACHING_INSTALL_LIB ${CACHING_INSTALL}/lib)
file(MAKE_DIRECTORY ${CACHING_INSTALL_LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CACHING_INSTALL_LIB})

set (TEST_INSTALL_BIN ${CACHING_INSTALL}/test-bin)
file(MAKE_DIRECTORY ${TEST_INSTALL_BIN})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TEST_INSTALL_BIN})

## Use continuous allocation for buffer memory blocks in all locations
add_definitions(-DENABLE_BUFFER_CONTINUOUS_ALLOC)

add_library(caching SHARED
${CACHING_DIR_SRC}/DataCaching.cpp
)

set (BUFFER_VAL 'N')
add_definitions(-DBUFFER_SCHEDULING_POLICY=${BUFFER_VAL})
set (BUFFER_VAL ${BUFFER_VAL} PARENT_SCOPE)

target_link_libraries(caching ${INTERNAL_LINK_LIB} chl_grid_amalgamation)
configure_file(${CACHING_DIR_INC}/DataCaching.hpp ${CACHING_INSTALL_INC}/DataCaching.hpp COPYONLY)

add_executable (functionTest ${CACHING_DIR_TEST}/functionTest.cpp )
target_link_libraries(functionTest ${INTERNAL_LINK_LIB} caching)

set (CACHING_INSTALL_INC ${CACHING_INSTALL_INC} PARENT_SCOPE)
message( "->\nSoftware-assisted Memory Buffer - Shared Library: " ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} )
message( "Software-assisted Memory Buffer - Shared Library headers: " ${CACHING_INSTALL_INC} )
message( "------------------------------------------------------------------------------------------------")
