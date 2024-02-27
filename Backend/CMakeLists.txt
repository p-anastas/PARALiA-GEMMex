###########################################################################################
#---------------------------------Project-wide Options------------------------------------#
###########################################################################################

cmake_minimum_required (VERSION 3.10)

## Define a name for the backend library wrappers/implementations used in this build, Currently implemented only with CuCuBLAS(==CUDA+CUBLAS)
set(BACKEND "CuCuBLAS")
message("CHANGELINK - BACKEND: " ${BACKEND})
add_definitions(-DBACKEND="${BACKEND}")
## Define the project version (used in logfiles for result testing)

execute_process ( COMMAND bash -c "source ${CMAKE_CURRENT_BINARY_DIR}/../config_karol.sh")

if(NOT DEFINED PROJECT_VERSION)
  set(PROJECT_VERSION "0.1")
endif()
message("CHANGELINK - PROJECT_VERSION: " ${PROJECT_VERSION})
add_definitions(-DVERSION="${PROJECT_VERSION}")

## Define a name for the testbed, to be used for the whole framework setup
set(TESTBED_NAME $ENV{CHANGELINK_SYSTEM})
message("CHANGELINK - TESTBED_NAME: " ${TESTBED_NAME})
add_definitions(-DTESTBED="${TESTBED_NAME}")
set(TESTBED_NAME $ENV{CHANGELINK_SYSTEM} PARENT_SCOPE)

## Define the workspace where execution of the library (benchmarks & deployment) will be invoked from (e.g. to avoid Disk quota problems in /home)
set(CD_LOCATION $ENV{CD_FOLDER})
message("CHANGELINK - CD_LOCATION: " ${CD_LOCATION})

## Define the max GPU memory percentage (%) a problem can use ( not accounting for streams and other control structs; e.g. only the problem data ).
set(PROBLEM_GPU_PERCENTAGE 70)
message("CHANGELINK - PROBLEM_GPU_PERCENTAGE: " ${PROBLEM_GPU_PERCENTAGE})
add_definitions(-DPROBLEM_GPU_PERCENTAGE=${PROBLEM_GPU_PERCENTAGE})

###########################################################################################
#----------------------------------DEBUG/TEST Options-------------------------------------#
###########################################################################################

#Run in debug mode, which includes detailed function calls and info for error-checking
add_definitions(-DCLDEBUG)

#Run in deep debug mode, for command-to-command debugging
#add_definitions(-DDCLDEBUG)

#Run in Link Hop testing mode; transfers using link hops are logged and timed.
#add_definitions(-DCLTEST)

###########################################################################################
#------------------------------------Other options----------------------------------------#
###########################################################################################

## The smart numa alloc method
#add_definitions(-DSMART_NUMALLOC_ROWWISE)
#add_definitions(-DSMART_NUMALLOC_COLWISE)

## The number of parallel workers that can be utilized at each device for comp-comp overlap. 
#add_definitions(-DMAX_BACKEND_L=8)

#add_definitions(-DREDUCE_WORKERS_PERDEV=10)

## Define at which percentage to normalize float values for similar data (e.g. costs, splits etc) to avoid empirical errors leaking into decisions
add_definitions(-DNORMALIZE_NEAR_SPLIT_LIMIT=0.05)

## The workspace size that should be assigned to each cublas handle
## If set to -1, the default ws will be used instead 
## NOTE-FIXME!!! : Currently WS effects the PCIe BW when a cublas kernel is launched. 
## NOTE-FIXME!!! This has a serious overhead if SKs are launched during selection
## NOTE-FIXME!!! Use of SUBKERNELS_FIRE_WHEN_READY is advised for mitigating this
add_definitions(-DWS_SZ=0)

###########################################################################################
#-----------------------------------CUDA definitions--------------------------------------#
###########################################################################################


if(NOT $ENV{CHANGELINK_CUDA_TOOLKIT_PREFIX} STREQUAL "default")
  set (CUDA_PREFIX $ENV{CHANGELINK_CUDA_TOOLKIT_PREFIX})
  set (CMAKE_CUDA_COMPILER ${CUDA_PREFIX}/bin/nvcc)
endif()

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES $ENV{CHANGELINK_CUDA_ARCH})
endif()

set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -w -O3")
message("CHANGELINK - CMAKE_CUDA_FLAGS: " ${CMAKE_CUDA_FLAGS})

set (CUDA_INCLUDE_DIRS ${CUDA_PREFIX}/include ${CUDA_PREFIX}/samples/common/inc)
include_directories(${CUDA_INCLUDE_DIRS})
message("CHANGELINK - CUDA_INCLUDE_DIRS: " ${CUDA_INCLUDE_DIRS})

set (CUDA_LD ${CUDA_PREFIX}/lib64) # or /lib depending on system
set (CUDA_LINK "-L${CUDA_LD} -lcudart -lcublas -lcurand $ENV{CHANGELINK_CUDA_LOAD_COMMAND} -fopenmp")
message("CHANGELINK - CUDA_LINK: " ${CUDA_LINK} )
set (INTERNAL_LINK_LIB "${CUDA_LINK} -lrt -ldl -lnuma")
message("CHANGELINK - INTERNAL_LINK_LIB: " ${INTERNAL_LINK_LIB})

###########################################################################################
#----------------------------------C/C++ definitions--------------------------------------#
###########################################################################################

if(NOT $ENV{CHANGELINK_CXX_PREFIX} STREQUAL "default")
  set (CXX_PREFIX $ENV{CHANGELINK_CXX_PREFIX})
  set (CMAKE_CXX_COMPILER ${CXX_PREFIX}/bin/c++)
  set (CMAKE_C_COMPILER ${CXX_PREFIX}/bin/gcc)
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -lm -Wuninitialized -lgomp -lnuma")
message("CHANGELINK - CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})

#set(CMAKE_C_FLAGS ${CMAKE_CXX_FLAGS})

###########################################################################################
#-------------------------------------Start Project---------------------------------------#
###########################################################################################

project (CHANGELINK VERSION ${PROJECT_VERSION} DESCRIPTION "CHANGELINK interconnect benchmarks" LANGUAGES CUDA CXX)

if($ENV{CHANGELINK_INSTALL_PREFIX} STREQUAL "default")
  message("CHANGELINK - CHANGELINK_INSTALL_PREFIX not provided, project will be installed under build dir.")
  set (CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${TESTBED_NAME}-install)
else()
  set (CMAKE_INSTALL_PREFIX $ENV{CHANGELINK_INSTALL_PREFIX})
endif()
set (CMAKE_INSTALL_PROJECT ${CMAKE_INSTALL_PREFIX})
message("CHANGELINK - CMAKE_INSTALL_PREFIX, CMAKE_INSTALL_PROJECT: " ${CMAKE_INSTALL_PROJECT})

###########################################################################################
#----------------------------------------chl_libs-----------------------------------------#
###########################################################################################


message( "CHANGELINK - Included smart_wrappers DIR: " ${PROJECT_SOURCE_DIR}/smart_wrappers )
include_directories(${PROJECT_SOURCE_DIR}/smart_wrappers)

message( "CHANGELINK - Included grid_amalgamation DIR: " ${PROJECT_SOURCE_DIR}/grid_amalgamation )
include_directories(${PROJECT_SOURCE_DIR}/grid_amalgamation)

set (CHANGELINK_INSTALL ${CMAKE_INSTALL_PROJECT}/CHANGELINK)
file(MAKE_DIRECTORY ${CHANGELINK_INSTALL})

set (CHANGELINK_INSTALL_LIB ${CHANGELINK_INSTALL}/lib)
file(MAKE_DIRECTORY ${CHANGELINK_INSTALL_LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CHANGELINK_INSTALL_LIB})

set (CHANGELINK_INSTALL_INC ${CHANGELINK_INSTALL}/include)
file(MAKE_DIRECTORY ${CHANGELINK_INSTALL_INC})
set(CHANGELINK_INSTALL_INC ${CHANGELINK_INSTALL_INC} PARENT_SCOPE)

add_definitions(-DDEPLOYDB="${CHANGELINK_INSTALL}")
message( "CHANGELINK - DEPLOYDB: " ${CHANGELINK_INSTALL})
file(MAKE_DIRECTORY ${CHANGELINK_INSTALL}/Benchmark-Results)

add_definitions(-DSCONFIG_PATH="$ENV{CHANGELINK_SCONFIG_FILE_PATH}")
message( "CHANGELINK - SCONFIG_PATH: " $ENV{CHANGELINK_SCONFIG_FILE_PATH})

add_library(chl_smart_wrappers SHARED
  ${PROJECT_SOURCE_DIR}/smart_wrappers/general_utils.cpp
  ${PROJECT_SOURCE_DIR}/smart_wrappers/data_helpers.cpp
  ${PROJECT_SOURCE_DIR}/smart_wrappers/numa_awareness.cpp
  ${PROJECT_SOURCE_DIR}/smart_wrappers/${BACKEND}/backend_util_wrappers.cu
  ${PROJECT_SOURCE_DIR}/smart_wrappers/${BACKEND}/backend_class_wrappers.cu
  ${PROJECT_SOURCE_DIR}/smart_wrappers/${BACKEND}/backend_data_wrappers.cu
  ${PROJECT_SOURCE_DIR}/smart_wrappers/${BACKEND}/backend_lib_wrappers.cu
  ${PROJECT_SOURCE_DIR}/smart_wrappers/${BACKEND}/backend_routine_wrappers.cpp
)

target_link_libraries(chl_smart_wrappers ${INTERNAL_LINK_LIB} ${OPENBLAS_LINK})

add_library(chl_grid_amalgamation SHARED
  ${PROJECT_SOURCE_DIR}/grid_amalgamation/chain_transfer.cpp
  ${PROJECT_SOURCE_DIR}/grid_amalgamation/grid_amalgamation.cpp
  ${PROJECT_SOURCE_DIR}/grid_amalgamation/problem_comm_translate.cpp
)

target_link_libraries(chl_grid_amalgamation ${INTERNAL_LINK_LIB} chl_smart_wrappers)

configure_file(${PROJECT_SOURCE_DIR}/smart_wrappers/smart_wrappers.hpp ${CHANGELINK_INSTALL_INC}/chl_smart_wrappers.hpp COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/grid_amalgamation/grid_amalgamation.hpp ${CHANGELINK_INSTALL_INC}/chl_grid_amalgamation.hpp COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/smart_wrappers/${BACKEND}/backend_wrappers.hpp ${CHANGELINK_INSTALL_INC}/backend_wrappers.hpp COPYONLY)

###########################################################################################
#------------------------------------microbenchmarks--------------------------------------#
###########################################################################################

## TODO: Choose between ->
## i) Automatic micro-benchmarking (Population checked within 95% CI) using boost (recommended if boost available).
message("CHANGELINK - Database Builder - Using Boost for microbenchmark statistics")
add_definitions(-DAUTO_BENCH_USE_BOOST)
set(BOOSTO_PREFIX $ENV{CHANGELINK_BOOST_PREFIX})
set (BOOSTO_INCLUDE_DIRS ${BOOSTO_PREFIX}/include)
include_directories(${BOOSTO_INCLUDE_DIRS})
message("CHANGELINK - Database Builder - Included (Boost) DIR: " ${BOOSTO_INCLUDE_DIRS})
## OR
## ii) Benchmark with a predifined itteration number for each problem size.
#add_definitions(-DITER=100)
#message("Database Builder - Performing predefined microbenchmark itterations")

include_directories(${CHANGELINK_INSTALL_INC})
message( "CHANGELINK - Database Builder - Included DIR: " ${CHANGELINK_INSTALL_INC} )

include_directories(${PROJECT_SOURCE_DIR}/microbenchmarks)
message( "CHANGELINK - Included microbenchmarks DIR: " ${PROJECT_SOURCE_DIR}/microbenchmarks)

set(NVEMAPI_PREFIX "${PROJECT_SOURCE_DIR}/nvidia-energy-measure")
set (NVEMAPI_INCLUDE_DIRS ${NVEMAPI_PREFIX}/include) # -L${CUDA_LD}
include_directories(${NVEMAPI_INCLUDE_DIRS} )
message( "CHANGELINK - NVEMAPI_INCLUDE_DIRS: " ${NVEMAPI_INCLUDE_DIRS} )
ADD_SUBDIRECTORY (${PROJECT_SOURCE_DIR}/nvidia-energy-measure)

add_executable (CHLMicrobench2DGrids ${PROJECT_SOURCE_DIR}/microbenchmarks/chl_2D_microbench_grids.cpp ${PROJECT_SOURCE_DIR}/microbenchmarks/helpers.cpp)
add_executable (CHLTestGamalgGEMM ${PROJECT_SOURCE_DIR}/microbenchmarks/chl_test_gamalg_gemm.cpp ${PROJECT_SOURCE_DIR}/microbenchmarks/helpers.cpp)
add_executable (CHLTestMemlocs ${PROJECT_SOURCE_DIR}/microbenchmarks/chl_test_memloc_numbering.cpp ${PROJECT_SOURCE_DIR}/microbenchmarks/helpers.cpp)
add_executable (CHLTest2DThroughput ${PROJECT_SOURCE_DIR}/microbenchmarks/chl_2D_throughput_test.cpp ${PROJECT_SOURCE_DIR}/microbenchmarks/helpers.cpp)

#add_executable (RunCuCuBlasLinkBenchGPUs ${PROJECT_SOURCE_DIR}/microbenchmarks/transfer_link_microbench_GPUs.cpp )

target_link_libraries(CHLMicrobench2DGrids ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)
target_link_libraries(CHLTestGamalgGEMM ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)
target_link_libraries(CHLTestMemlocs ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)
target_link_libraries(CHLTest2DThroughput ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)

#target_link_libraries(RunCuCuBlasLinkBenchGPUs ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)


