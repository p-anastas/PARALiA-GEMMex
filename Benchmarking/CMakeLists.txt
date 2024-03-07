cmake_minimum_required (VERSION 3.10)
### Initialize the project
project (PARALiA-Benchmarks VERSION ${PROJECT_VERSION} DESCRIPTION "Library wrappers for performance comparisson with PARALiA" LANGUAGES CXX CUDA)
message( "------------------------------------------------------------------------------------------------")
message( "Benchmarking - PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )
#message( "Library Scheduler - CMAKE_CURRENT_SOURCE_DIR: " ${CMAKE_CURRENT_SOURCE_DIR} )
#message( "Library Scheduler - CMAKE_CURRENT_BINARY_DIR: " ${CMAKE_CURRENT_BINARY_DIR} )

set ( BENCH_DIR ${PROJECT_SOURCE_DIR}/${BACKEND} )
message( "Benchmarking - LIBSC_DIR: " ${BENCH_DIR} )

set (BENCH_DIR_INC ${BENCH_DIR}/include)
include_directories( ${BENCH_DIR_INC})
message( "Benchmarking - Included DIR: " ${BENCH_DIR_INC} )

include_directories(${PARALIA_INSTALL_INC})
message( "Benchmarking - Included DIR: " ${PARALIA_INSTALL_INC} )

include_directories(${AUTOTUNE_INSTALL_INC})
message( "Benchmarking - Included DIR: " ${AUTOTUNE_INSTALL_INC} )

include_directories(${LIBSC_INSTALL_INC})
message( "Benchmarking - Included DIR: " ${LIBSC_INSTALL_INC} )

include_directories(${CACHING_INSTALL_INC})
message( "Benchmarking - Included DIR: " ${CACHING_INSTALL_INC} )

set (BENCH_DIR_SRC ${BENCH_DIR}/src)
message( "Benchmarking - Source DIR: " ${BENCH_DIR_SRC} )

set (BENCH_INSTALL ${CMAKE_INSTALL_PROJECT}/Benchmarking)
file(MAKE_DIRECTORY ${BENCH_INSTALL})

set (BENCH_INSTALL_INC ${BENCH_INSTALL}/include)
file(MAKE_DIRECTORY ${BENCH_INSTALL_INC})

#set (BENCH_INSTALL_LIB ${BENCH_INSTALL}/lib)
#file(MAKE_DIRECTORY ${BENCH_INSTALL_LIB})
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIBSC_INSTALL_LIB})

set (BENCH_INSTALL_BIN ${BENCH_INSTALL}/testing-bin)
file(MAKE_DIRECTORY ${BENCH_INSTALL_BIN})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BENCH_INSTALL_BIN})

file(MAKE_DIRECTORY ${BENCH_INSTALL}/testLogs)
add_definitions(-DTESTLIBDIR="${BENCH_INSTALL}/testLogs")

add_library(benchparalia SHARED
${BENCH_DIR_SRC}/cuBLASXt/cuBLASXtWrappers.cpp
${BENCH_DIR_SRC}/testing.cpp
)
target_link_libraries(benchparalia ${INTERNAL_LINK_LIB} ${OPENBLAS_LINK} paralia autotuner chl_grid_amalgamation )

#Enable results validation in "Runner" scripts (adds extra overhead, "Tester" scripts should be used for result validation)
add_definitions(-DRUNVALIDATION)
#Enable checking the existing logs for each benchmark before executing it, to avoid running the same exact test multiple times.
#add_definitions(-DCHECKLOG)

## Internal testers for PARALia components
#add_executable (testLinkOverlap ${BENCH_DIR_SRC}/Component_testers/testLinkOverlap.cpp)
#add_executable (testFasTransfer ${BENCH_DIR_SRC}/Component_testers/testFasTransfer.cpp)

#target_link_libraries(testLinkOverlap ${INTERNAL_LINK_LIB} benchparalia )
#target_link_libraries(testFasTransfer ${INTERNAL_LINK_LIB} benchparalia )

add_executable (cuBLASXt_dgemm_runner ${BENCH_DIR_SRC}/cuBLASXt/cuBLASXt_dgemm_runner.cpp )
## PARALia Runners and Testers
add_executable (dgemm_tester ${BENCH_DIR_SRC}/PARALia/dgemm_tester.cpp )
#add_executable (sgemm_tester ${BENCH_DIR_SRC}/PARALia/sgemm_tester.cpp )
add_executable (dgemm_runner ${BENCH_DIR_SRC}/PARALia/dgemm_runner.cpp )
#add_executable (dgemm_runner_hetero ${BENCH_DIR_SRC}/PARALia/dgemm_runner_hetero.cu )
add_executable (dgemm_predictor ${BENCH_DIR_SRC}/PARALia/dgemm_predictor.cpp )

target_link_libraries(cuBLASXt_dgemm_runner ${INTERNAL_LINK_LIB} benchparalia )
target_link_libraries(dgemm_tester ${INTERNAL_LINK_LIB} benchparalia )
target_link_libraries(dgemm_runner ${INTERNAL_LINK_LIB} benchparalia )
target_link_libraries(dgemm_predictor ${INTERNAL_LINK_LIB} benchparalia )

set (BENCH_INSTALL_INC ${BENCH_INSTALL_INC} PARENT_SCOPE)
#message( "->\nBenchmarking - Shared Library: " ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} )
#message( "Benchmarking - Shared Library headers: " ${BENCH_INSTALL_INC} )
message( "Benchmarking - Library testing bin DIR: " ${BENCH_INSTALL_BIN} )
message( "------------------------------------------------------------------------------------------------")
