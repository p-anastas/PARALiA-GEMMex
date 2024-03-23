set (PARALIA_BACKEND_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
message( "PARALIA - Backend - PARALIA_BACKEND_SOURCE_DIR = " ${PARALIA_BACKEND_SOURCE_DIR} )

message( "PARALIA - Backend - Included smart_wrappers DIR: " ${PARALIA_BACKEND_SOURCE_DIR}/smart_wrappers )
include_directories(${PARALIA_BACKEND_SOURCE_DIR}/smart_wrappers)

message( "PARALIA - Backend - Included grid_amalgamation DIR: " ${PARALIA_BACKEND_SOURCE_DIR}/grid_amalgamation )
include_directories(${PARALIA_BACKEND_SOURCE_DIR}/grid_amalgamation)

add_library(chl_smart_wrappers SHARED
  ${PARALIA_BACKEND_SOURCE_DIR}/smart_wrappers/general_utils.cpp
  ${PARALIA_BACKEND_SOURCE_DIR}/smart_wrappers/data_helpers.cpp
  #${PARALIA_BACKEND_SOURCE_DIR}/smart_wrappers/numa_awareness.cpp # This is deprecated
  ${PARALIA_BACKEND_SOURCE_DIR}/smart_wrappers/${BACKEND}/backend_util_wrappers.cu
  ${PARALIA_BACKEND_SOURCE_DIR}/smart_wrappers/${BACKEND}/backend_class_wrappers.cu
  ${PARALIA_BACKEND_SOURCE_DIR}/smart_wrappers/${BACKEND}/backend_data_wrappers.cu
  ${PARALIA_BACKEND_SOURCE_DIR}/smart_wrappers/${BACKEND}/backend_lib_wrappers.cu
  ${PARALIA_BACKEND_SOURCE_DIR}/smart_wrappers/${BACKEND}/backend_routine_wrappers.cpp
)

target_link_libraries(chl_smart_wrappers ${INTERNAL_LINK_LIB} ${OPENBLAS_LINK})

add_library(chl_grid_amalgamation SHARED
  ${PARALIA_BACKEND_SOURCE_DIR}/grid_amalgamation/chain_transfer.cpp
  ${PARALIA_BACKEND_SOURCE_DIR}/grid_amalgamation/grid_amalgamation.cpp
  ${PARALIA_BACKEND_SOURCE_DIR}/grid_amalgamation/problem_comm_translate.cpp
)

target_link_libraries(chl_grid_amalgamation ${INTERNAL_LINK_LIB} chl_smart_wrappers)

configure_file(${PARALIA_BACKEND_SOURCE_DIR}/smart_wrappers/smart_wrappers.hpp ${PARALIA_INSTALL_INC}/chl_smart_wrappers.hpp COPYONLY)
configure_file(${PARALIA_BACKEND_SOURCE_DIR}/grid_amalgamation/grid_amalgamation.hpp ${PARALIA_INSTALL_INC}/chl_grid_amalgamation.hpp COPYONLY)
configure_file(${PARALIA_BACKEND_SOURCE_DIR}/smart_wrappers/${BACKEND}/backend_wrappers.hpp ${PARALIA_INSTALL_INC}/backend_wrappers.hpp COPYONLY)

###########################################################################################
#------------------------------------microbenchmarks--------------------------------------#
###########################################################################################

## TODO: Choose between ->
## i) Automatic micro-benchmarking (Population checked within 95% CI) using boost (recommended if boost available).
message("PARALIA - Backend - Database Builder - Using Boost for microbenchmark statistics")
add_definitions(-DAUTO_BENCH_USE_BOOST)
set(BOOSTO_PREFIX $ENV{PARALIA_BOOST_PREFIX})
set (BOOSTO_INCLUDE_DIRS ${BOOSTO_PREFIX}/include)
include_directories(${BOOSTO_INCLUDE_DIRS})
message("PARALIA - Backend - Database Builder - Included (Boost) DIR: " ${BOOSTO_INCLUDE_DIRS})
## OR
## ii) Benchmark with a predifined itteration number for each problem size.
#add_definitions(-DITER=100)
#message("Database Builder - Performing predefined microbenchmark itterations")

include_directories(${PARALIA_INSTALL_INC})
message( "PARALIA - Backend - Database Builder - Included DIR: " ${PARALIA_INSTALL_INC} )

include_directories(${PARALIA_BACKEND_SOURCE_DIR}/microbenchmarks)
message( "PARALIA - Backend - Included microbenchmarks DIR: " ${PARALIA_BACKEND_SOURCE_DIR}/microbenchmarks)

set(NVEMAPI_PREFIX "${PARALIA_BACKEND_SOURCE_DIR}/nvidia-energy-measure")
set (NVEMAPI_INCLUDE_DIRS ${NVEMAPI_PREFIX}/include) # -L${CUDA_LD}
include_directories(${NVEMAPI_INCLUDE_DIRS} )
message( "PARALIA - Backend - NVEMAPI_INCLUDE_DIRS: " ${NVEMAPI_INCLUDE_DIRS} )
message( "------------------------------------------------------------------------------------------------")
ADD_SUBDIRECTORY (${PARALIA_BACKEND_SOURCE_DIR}/nvidia-energy-measure)

add_executable (CHLMicrobench2DGrids ${PARALIA_BACKEND_SOURCE_DIR}/microbenchmarks/chl_2D_microbench_grids.cpp ${PARALIA_BACKEND_SOURCE_DIR}/microbenchmarks/helpers.cpp)
add_executable (CHLTestGamalgGEMM ${PARALIA_BACKEND_SOURCE_DIR}/microbenchmarks/chl_test_gamalg_gemm.cpp ${PARALIA_BACKEND_SOURCE_DIR}/microbenchmarks/helpers.cpp)
add_executable (CHLTestMemlocs ${PARALIA_BACKEND_SOURCE_DIR}/microbenchmarks/chl_test_memloc_numbering.cpp ${PARALIA_BACKEND_SOURCE_DIR}/microbenchmarks/helpers.cpp)
add_executable (CHLTest2DThroughput ${PARALIA_BACKEND_SOURCE_DIR}/microbenchmarks/chl_2D_throughput_test.cpp ${PARALIA_BACKEND_SOURCE_DIR}/microbenchmarks/helpers.cpp)

#add_executable (RunCuCuBlasLinkBenchGPUs ${PARALIA_BACKEND_SOURCE_DIR}/microbenchmarks/transfer_link_microbench_GPUs.cpp )

target_link_libraries(CHLMicrobench2DGrids ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)
target_link_libraries(CHLTestGamalgGEMM ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)
target_link_libraries(CHLTestMemlocs ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)
target_link_libraries(CHLTest2DThroughput ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)

#target_link_libraries(RunCuCuBlasLinkBenchGPUs ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)
