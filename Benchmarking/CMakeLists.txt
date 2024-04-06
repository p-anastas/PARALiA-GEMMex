set ( BENCH_DIR_BASE ${CMAKE_CURRENT_SOURCE_DIR}/ )
message( "Benchmarking - BENCH_DIR_BASE: " ${BENCH_DIR_BASE} )
set (BENCH_DIR_BASE ${BENCH_DIR_BASE} PARENT_SCOPE)

set ( BENCH_DIR ${BENCH_DIR_BASE}/${BACKEND} )
message( "Benchmarking - BENCH_DIR: " ${BENCH_DIR} )

set (BENCH_DIR_INC ${BENCH_DIR}/include)
include_directories( ${BENCH_DIR_INC})
message( "Benchmarking - Included DIR: " ${BENCH_DIR_INC} )

include_directories(${PARALIA_INSTALL_INC})
message( "Benchmarking - Included DIR: " ${PARALIA_INSTALL_INC} )

set (BENCH_DIR_SRC ${BENCH_DIR}/src)
message( "Benchmarking - Source DIR: " ${BENCH_DIR_SRC} )

file(MAKE_DIRECTORY ${PARALIA_INSTALL}/testLogs)
add_definitions(-DTESTLIBDIR="${PARALIA_INSTALL}/testLogs")

set (BENCHMARK_INSTALL_BIN ${PARALIA_INSTALL}/testing-bin)
file(MAKE_DIRECTORY ${BENCHMARK_INSTALL_BIN})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BENCHMARK_INSTALL_BIN})

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
add_executable (cuBLASXt_sgemm_runner ${BENCH_DIR_SRC}/cuBLASXt/cuBLASXt_sgemm_runner.cpp )

## PARALia Runners and Testers
add_executable (dgemm_tester ${BENCH_DIR_SRC}/PARALia/dgemm_tester.cpp )
add_executable (dgemm_runner ${BENCH_DIR_SRC}/PARALia/dgemm_runner.cpp )

add_executable (sgemm_tester ${BENCH_DIR_SRC}/PARALia/sgemm_tester.cpp )
add_executable (sgemm_runner ${BENCH_DIR_SRC}/PARALia/sgemm_runner.cpp )
add_executable (hgemm_runner ${BENCH_DIR_SRC}/PARALia/hgemm_runner.cpp )

#add_executable (dgemm_runner_hetero ${BENCH_DIR_SRC}/PARALia/dgemm_runner_hetero.cu )
add_executable (dgemm_predictor ${BENCH_DIR_SRC}/PARALia/dgemm_predictor.cpp )

target_link_libraries(cuBLASXt_dgemm_runner ${INTERNAL_LINK_LIB} benchparalia )
target_link_libraries(cuBLASXt_sgemm_runner ${INTERNAL_LINK_LIB} benchparalia )

#target_link_libraries(dgemm_tester_B-${MAX_BACKEND_L}_${OUTPUT_ALGO_MODE}_S-${STREAMING_BUFFER_OVERLAP}_${DISTRIBUTION}-${ORDER_2DBC}_${TASK_ORDER}_F-${FETCH_ROUTING}_WB-${WB_ROUTING}_${PREDICT_OPTIMIZE_TARGET} ${INTERNAL_LINK_LIB} benchparalia )
#target_link_libraries(dgemm_runner_B-${MAX_BACKEND_L}_${OUTPUT_ALGO_MODE}_S-${STREAMING_BUFFER_OVERLAP}_${DISTRIBUTION}-${ORDER_2DBC}_${TASK_ORDER}_F-${FETCH_ROUTING}_WB-${WB_ROUTING}_${PREDICT_OPTIMIZE_TARGET}  ${INTERNAL_LINK_LIB} benchparalia )
target_link_libraries(dgemm_tester ${INTERNAL_LINK_LIB} benchparalia )
target_link_libraries(dgemm_runner ${INTERNAL_LINK_LIB} benchparalia )

target_link_libraries(sgemm_tester ${INTERNAL_LINK_LIB} benchparalia )
target_link_libraries(sgemm_runner ${INTERNAL_LINK_LIB} benchparalia )
target_link_libraries(hgemm_runner ${INTERNAL_LINK_LIB} benchparalia )

target_link_libraries(dgemm_predictor ${INTERNAL_LINK_LIB} benchparalia )

message( "------------------------------------------------------------------------------------------------")