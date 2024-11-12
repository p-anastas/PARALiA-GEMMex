###########################################################################################
#---------------------------------Project-wide Options------------------------------------#
###########################################################################################

cmake_minimum_required (VERSION 3.10)

## Define a name for the backend library wrappers/implementations used in this build, Currently implemented only with CuCuBLAS(==CUDA+CUBLAS)
set(BACKEND "CuCuBLAS")
message("PARALiA-GEMMex - BACKEND: " ${BACKEND})
add_definitions(-DBACKEND="${BACKEND}")
## Define the project version (used in logfiles for result testing)

set(PROJECT_VERSION "2.99")
message("PARALiA-GEMMex - PROJECT_VERSION: " ${PROJECT_VERSION})
add_definitions(-DVERSION="${PROJECT_VERSION}")

## Define a name for the testbed, to be used for the whole framework setup
set(TESTBED_NAME $ENV{PARALIA_SYSTEM})
message("PARALiA-GEMMex - TESTBED_NAME: " ${TESTBED_NAME})
add_definitions(-DTESTBED="${TESTBED_NAME}")

add_definitions(-DSYSDB="${CMAKE_CURRENT_SOURCE_DIR}/Deployment_files/${TESTBED_NAME}")
message("PARALiA-GEMMex - SYSDB: " ${CMAKE_CURRENT_SOURCE_DIR}/Deployment_files/${TESTBED_NAME})

## Define the max GPU memory percentage (%) a problem can use (not accounting for streams and other control structs; e.g. only the problem data).
set(PROBLEM_GPU_PERCENTAGE 90)
message("PARALiA-GEMMex - PROBLEM_GPU_PERCENTAGE: " ${PROBLEM_GPU_PERCENTAGE})
add_definitions(-DPROBLEM_GPU_PERCENTAGE=${PROBLEM_GPU_PERCENTAGE})

## Define the (max) number of devices PARALiA will utilize during execution. 
set(NUM_DEVICES $ENV{PARALIA_NUM_DEVICES})
message("PARALiA-GEMMex - NUM_DEVICES: " ${NUM_DEVICES})
add_definitions(-DCHL_WORKERS=${NUM_DEVICES})

## Define the workspace where execution of the library (benchmarks & deployment) will be invoked from (e.g. to avoid Disk quota problems in /home)
set(CD_LOCATION $ENV{CD_FOLDER})
message("PARALiA-GEMMex - CD_LOCATION: " ${CD_LOCATION})

###########################################################################################
#---------------------------Internal resource reuse options-------------------------------#
###########################################################################################

#--------------------------------Tested (production)--------------------------------------#
## Enable Asynchronous execution for Subkernels
add_definitions(-DASYNC_ENABLE)

## Enable GPU backend resourse reuse (streams, handles)
add_definitions(-DQUEUE_REUSE_ENABLE)

## Enable GPU and CPU allocated buffer reuse between different calls
add_definitions(-DBUFFER_REUSE_ENABLE)

## The level of reuse for the autotuning metadata of every routine call by SIMILAR subsequent calls
## NONE -> Rerun full autotuning every time
## MODELS -> Keep the models but rerun the distribution and routing functions
## PROBLEMS -> Keep a cache of size PROBLEM_MD_CACHE with the autotuning MD of DIFFERENT routine calls
##             This enables the full autotuning metadata reuse of up to PROBLEM_MD_CACHE problems. 
#add_definitions(-DMETADATA_REUSE_NONE)
#add_definitions(-DMETADATA_REUSE_MODELS)
add_definitions(-DMETADATA_REUSE_PROBLEMS) 
add_definitions(-DPROBLEM_MD_CACHE=10) # Only relevant for "PROBLEMS" reuse level

###########################################################################################
#---------------------------Hardware/Software Optimizations----------------------------#
###########################################################################################

#--------------------------------Tested (production)--------------------------------------#
## The number of parallel workers that can be utilized at each device for comp-comp overlap. 
set(MAX_BACKEND_L 8)

## The number of sub-buffers FasTCoCoMemcpy2DAsync will use to increase effective BW
## In general effective_trans_time = max_trans_time(route[]) + 1/SBO * sub_trans_time(route_minus_max[])
## Rule of thumb: large SBO -> increased transfer latency, small SBO -> decreased bandwidth
set(STREAMING_BUFFER_OVERLAP 8)

## The output algorithm for obtaining the actual results of each routine. Options:
## - ALGO_WR: The output tiles are processed normally. This entails being fetched at the beggining, 
##            locked while being written on (TODO) and written back after all computations are complete.
## - ALGO_WR_LAZY: The output tiles are fetched from the source lazily. Computations will use
##             local buffers. After all computation on each tile is fired, its source values will be
##             be fetched, reduced locally in the last computation location and written back.
## - ALGO_WREDUCE: The output tiles are not fetched from the source. Instead, computations will use
##             local buffers. After they are complete, tiles are reduced to the output (=source) location.
#set(OUTPUT_ALGO_MODE ALGO_WR)
#set(OUTPUT_ALGO_MODE ALGO_WR_LAZY)
set(OUTPUT_ALGO_MODE ALGO_AUTO)

# Enable overlapping Send & recv Communication
add_definitions(-DENABLE_SEND_RECV_OVERLAP)

## Define if Subkernel operations are fired together with data fetches 
## OR when their data become available
#add_definitions(-DSUBKERNELS_FIRE_LAZY)

#--------------------------------------Experimental----------------------------------------#
## Enable using the CPU as an extra device for part of the total workload, if applicable
#add_definitions(-DENABLE_CPU_WORKLOAD)
#add_definitions(-DOUTPUT_ALGO_MODE="ALGO_WREDUCE")
add_definitions(-DREDUCE_WORKERS_PERDEV=8)

###########################################################################################
#---------------------------General autotuning options------------------------------------#
###########################################################################################

###########################################################################################
#---------------------------------Tiling size selection-----------------------------------#
###########################################################################################

#--------------------------------Tested (production)--------------------------------------#
## Decision regarding perplexing tile selection and workload selection
## by applying estimated overlap sl from tile to device split prediction
add_definitions(-DAPPLY_TILE_SL_TO_WORKLOAD_SPLIT)
## The number of different top "REP_TILE" tile sizes that will be tested.
## Only fruitful for repetitive GEMMs on the same configuration (and with reps >> REP_TILE)
## The maximum allowed value for this is 10 (any higher will be set to 10 instead)
set(REP_TILE 5)

## TODO: For future work could autotune conditions.
## Conditions (with decreasing priority): 
## 1: NO-Imbalance: Spliting creates subkernels that can be distributed equally to devices without WR-sharing
## 2: NO-remainder: split should not create remainder Tiles if possible
## 3: T-min: Tile should be >= than a preset value
add_definitions(-DTILE_MIN=2048)
add_definitions(-DTILE_MIN_SLOWDOWN=0.1)
## 4: SK-num: subkernels per device must be >= MIN_DESIRED_SK_DEV
add_definitions(-DMAX_DESIRED_SK_DEV=256)
add_definitions(-DMAX_DESIRED_SK_DEV_SLOWDOWN=0.1)
## 5: T-max: Tile should be <= than a preset value
add_definitions(-DTILE_MAX=8192)
add_definitions(-DTILE_MAX_SLOWDOWN=0.1)


###########################################################################################
#-----------------------------------Task distributions------------------------------------#
###########################################################################################

#--------------------------------Tested (production)--------------------------------------#
set(DISTRIBUTION "2D-BLOCK-CYCLIC")
# Select which value to use for the largest dimension
# e.g. for D1xD2 = 8 devices, which dimension should be 4 and which 2.
# This impacts the reuse pattern and the actual decomposition so it has potential performance implications
set(ORDER_2DBC D1_lesseq_D2) #Default
#set(ORDER_2DBC "D2_lesseq_D1") 

#--------------------------------------Experimental----------------------------------------#
#add_definitions(-DDISTRIBUTION="SPLIT-CHUNKS-ROBIN")
#add_definitions(-DDISTRIBUTION="ROUND-ROBIN")
#add_definitions(-DDISTRIBUTION="SPLIT-NAIVE")
#add_definitions(-DDISTRIBUTION="SPLIT-CHUNKS-ROBIN-REVERSE")

###########################################################################################
#--------------------------------Fetch routing options--------------------------------#
###########################################################################################
#Select which optimize_tasks() algorithm to use for fetch communication optimization.

#--------------------------------Tested (production)--------------------------------------#
#set (FETCH_ROUTING P2P_FETCH_FROM_INIT)
#set (FETCH_ROUTING P2P_FETCH_FROM_GPU_SERIAL) # This is worse than distance
#set (FETCH_ROUTING P2P_FETCH_FROM_GPU_DISTANCE)
#set (FETCH_ROUTING CHAIN_FETCH_SERIAL) # This is worse than random
#set (FETCH_ROUTING CHAIN_FETCH_RANDOM)
#set (FETCH_ROUTING CHAIN_FETCH_TIME) # This is worse than ETA
set (FETCH_ROUTING CHAIN_FETCH_QUEUE_WORKLOAD)

###########################################################################################
#--------------------------------Writeback routing options--------------------------------#
###########################################################################################
#Select which optimize_tasks() algorithm to use for writeback communication optimization.

#--------------------------------Tested (production)--------------------------------------#
set (WB_ROUTING P2P_TO_INIT)

#--------------------------------------Experimental----------------------------------------#
## Enable the use of transfer hops through memory locations that do not need a tile
## Improves bandwidth, but adds (extra) network load to fast lanes. Currently not integrated in already-chained transfers.
#add_definitions(-DENABLE_TRANSFER_HOPS)

#add_definitions(-DHOP_FETCH_QUEUE_WORKLOAD)
#add_definitions(-DHOP_FETCH_BW_PLUS_ETA)

###########################################################################################
#-----------------------------------Task order options------------------------------------#
###########################################################################################
#Select which optimize_tasks() algorithm to use for selecting the order of firing.

#--------------------------------Tested (production)--------------------------------------#
#set(TASK_ORDER SERIAL)
#set(TASK_ORDER FETCH_MINFETCH)
#set(TASK_ORDER FETCH_MINFETCH_THEN_MINPENDING)
set(TASK_ORDER FETCH_ETA)
#set(TASK_ORDER FETCH_ETA_PLUS_MINPENDING)
#--------------------------------------Experimental----------------------------------------#

###########################################################################################
#---------------------------------Optimization targets------------------------------------#
###########################################################################################

## Define at which percentage to normalize float values for similar data (e.g. costs, splits etc) to avoid empirical errors leaking into decisions
add_definitions(-DNORMALIZE_NEAR_SPLIT_LIMIT=0.05)
#TODO: FIXME !!!! TRIM execution logs automatically, currently by hand for testing

## Define the minimum allowed percentile contribution in total time from each potential unit
## In case of energy, the maximum allowed percentile energy increase by each new unit.
add_definitions(-DMINIMUM_UNIT_CONTRIBUTION=0.00)

## Enable power measuring for benchmarks, and energy-related prediction modes for the Autotuner
add_definitions(-DENABLE_POWA)

## Choose the Unit combination that maximizes total performance (minimum time).
set(PREDICT_OPTIMIZE_TARGET PERF)

## Choose the Unit combination that minimizes consumed energy (J).
## Rarely chooses multi-unit execution, since its almost always less energy-efficient than single-unit.
#set(PREDICT_OPTIMIZE_TARGET ENERGY)

## Choose the Unit combination that maximizes the power-delay product (similar results with energy).
#set(PREDICT_OPTIMIZE_TARGET POWER-DELAY)

## Choose the Unit combination that maximizes the energy-delay product - closer to "PERF" but also slightly considering energy.
#set(PREDICT_OPTIMIZE_TARGET ENERGY-DELAY)

#--------------------------------------Experimental----------------------------------------#

## Choose the Unit combination based on a configurable minimum acceptable percentile performance-improvement-per-J-increase.
## Example for choosing U1(tpred = X, En = J1) vs U2(tpred = Y, En = J2) units with PERPER_LIMIT: if ( X/Y >= PERPER_LIMIT*J2/J1) U2 else U1
#add_definitions(-DPREDICT_OPTIMIZE_TARGET="PERF-PER-J")
## PERPER_LIMIT higher -> "ENERGY", lower -> "PERF"
add_definitions(-DPERPER_LIMIT=0.0)

#add_definitions(-DPREDICT_OPTIMIZE_TARGET="OTHER_TBD")

###########################################################################################
#----------------------------------DEBUG/TEST Options-------------------------------------#
###########################################################################################

#Run in production mode, where checking of correctness for all functions is disabled
#add_definitions(-DPRODUCTION)

#Run in debug mode, which includes detailed function calls and info for error-checking
#add_definitions(-DDEBUG)

#Run in deep debug mode, for command-to-command debugging
#add_definitions(-DDDEBUG)

#Run Grid amalgamantion debug mode
#add_definitions(-DCLDEBUG)

#Run Simple Prediction debug mode (well-presented simple autotuner stuff)
#add_definitions(-DSDEBUG)

#Run Prediction debug mode (autotuner stuff)
#add_definitions(-DPDEBUG)

#Run Prediction deep debug mode (all autotuner stuff)
#add_definitions(-DDPDEBUG)

#Run in caching debug mode, for cache-related debugging
#add_definitions(-DCDEBUG)

#Run in smart_wrappers debug mode
#add_definitions(-DUDEBUG)

#Run in smart_wrappers deep debug mode, for linkmap-related debugging
#add_definitions(-DUDDEBUG)

#Run in testing mode : details about models, decisions and time spend in operations.
#add_definitions(-DTEST)

#Run in Subkernel testing mode : details about time spend in subkernel operations.
#add_definitions(-DSTEST)

#Run in Deep Subkernel testing mode : many details about time spend in subkernel operations.
#add_definitions(-DDSTEST)

#Run in Link Hop testing mode; transfers using link hops are logged and timed.
#add_definitions(-DTTEST)

###########################################################################################
#------------------------------------Other options----------------------------------------#
###########################################################################################

## The number of current production-ready supported datatypes.
add_definitions(-DDTYPE_NUM=3)

## The workspace size that should be assigned to each cublas handle
## If set to -1, the default ws will be used instead 
## NOTE-FIXME!!! : Currently WS effects the PCIe BW when a cublas kernel is launched. 
## NOTE-FIXME!!! This has a serious overhead if SKs are launched during selection
## NOTE-FIXME!!! Use of SUBKERNELS_FIRE_WHEN_READY is advised for mitigating this
add_definitions(-DWS_SZ=-1)

###########################################################################################
#-------------------------------Define library constants----------------------------------#
###########################################################################################

add_definitions(-DMAX_BACKEND_L_IN=${MAX_BACKEND_L})
add_definitions(-DOUTPUT_ALGO_MODE_IN="${OUTPUT_ALGO_MODE}")
add_definitions(-DSTREAMING_BUFFER_OVERLAP_IN=${STREAMING_BUFFER_OVERLAP})
add_definitions(-DDISTRIBUTION_IN="${DISTRIBUTION}")
add_definitions(-DORDER_2DBC_IN="${ORDER_2DBC}")
add_definitions(-DTASK_ORDER_IN="${TASK_ORDER}")
add_definitions(-DFETCH_ROUTING_IN="${FETCH_ROUTING}")
add_definitions(-DWB_ROUTING_IN="${WB_ROUTING}")
add_definitions(-DPREDICT_OPTIMIZE_TARGET_IN="${PREDICT_OPTIMIZE_TARGET}")
add_definitions(-DREP_TILE_IN=${REP_TILE})


###########################################################################################
#-----------------------------------CUDA definitions--------------------------------------#
###########################################################################################

if(NOT $ENV{PARALIA_CUDA_TOOLKIT_PREFIX} STREQUAL "default")
  set (CUDA_PREFIX $ENV{PARALIA_CUDA_TOOLKIT_PREFIX})
  set (CMAKE_CUDA_COMPILER ${CUDA_PREFIX}/bin/nvcc)
endif()

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES $ENV{PARALIA_CUDA_ARCH})
endif()

set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -w -O3 -Wno-pointer-arith -Wno-unused-result")
message("PARALiA-GEMMex - CMAKE_CUDA_FLAGS: " ${CMAKE_CUDA_FLAGS})

set (CUDA_INCLUDE_DIRS ${CUDA_PREFIX}/include ${CUDA_PREFIX}/samples/common/inc)
include_directories(${CUDA_INCLUDE_DIRS})
message("PARALiA-GEMMex - CUDA_INCLUDE_DIRS: " ${CUDA_INCLUDE_DIRS})

set (CUDA_LD ${CUDA_PREFIX}/lib64) # or /lib depending on system
set (CUDA_LINK "-L${CUDA_LD} -lcudart -lcublas -lcurand $ENV{PARALIA_CUDA_LOAD_COMMAND} -fopenmp")
message("PARALiA-GEMMex - CUDA_LINK: " ${CUDA_LINK} )
set (INTERNAL_LINK_LIB "${CUDA_LINK} -lrt -ldl -lnuma")
message("PARALiA-GEMMex - INTERNAL_LINK_LIB: " ${INTERNAL_LINK_LIB})

###########################################################################################
#----------------------------------C/C++ definitions--------------------------------------#
###########################################################################################

if(NOT $ENV{PARALIA_CXX_PREFIX} STREQUAL "default")
  set (CXX_PREFIX $ENV{PARALIA_CXX_PREFIX})
  set (CMAKE_CXX_COMPILER ${CXX_PREFIX}/bin/c++)
  set (CMAKE_C_COMPILER ${CXX_PREFIX}/bin/gcc)
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -lm -lgomp -lnuma -mavx2 -Wno-pointer-arith -Wno-unused-result") 
#  -Wuninitialized 
message("PARALiA-GEMMex - CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})

#set(CMAKE_C_FLAGS ${CMAKE_CXX_FLAGS})

###########################################################################################
#---------------------------------OpenBLAS definitions------------------------------------#
###########################################################################################

if(NOT $ENV{PARALIA_OPENBLAS_PREFIX} STREQUAL "default")
  set(OPENBLAS_PREFIX $ENV{PARALIA_OPENBLAS_PREFIX})
  set (OPENBLAS_INCLUDE_DIRS ${OPENBLAS_PREFIX}/include)
  include_directories(${OPENBLAS_INCLUDE_DIRS})
  message("PARALiA-GEMMex - OPENBLAS_INCLUDE_DIRS: " ${OPENBLAS_INCLUDE_DIRS})
endif()

set (OPENBLAS_LD ${OPENBLAS_PREFIX}/lib)
set (OPENBLAS_LINK  "-L${OPENBLAS_LD} -lopenblas -lgomp" )
message("PARALiA-GEMMex - OPENBLAS_LINK: " ${OPENBLAS_LINK})

###########################################################################################
#-------------------------------------Start Project---------------------------------------#
###########################################################################################

project (PARALiA-GEMMex VERSION ${PROJECT_VERSION} DESCRIPTION "A PARALiA implementation of gemm similar to cublasLt functions, for multi-GPU" LANGUAGES CUDA CXX)

message("PARALiA-GEMMex - Project will be installed under build dir.")
set (CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${TESTBED_NAME}-install)

set (CMAKE_INSTALL_PROJECT ${CMAKE_INSTALL_PREFIX})
message("PARALiA-GEMMex - CMAKE_INSTALL_PREFIX, CMAKE_INSTALL_PROJECT: " ${CMAKE_INSTALL_PROJECT})

set (PARALIA_INSTALL ${CMAKE_INSTALL_PROJECT})
file(MAKE_DIRECTORY ${PARALIA_INSTALL})

set (PARALIA_INSTALL_LIB ${PARALIA_INSTALL}/lib)
file(MAKE_DIRECTORY ${PARALIA_INSTALL_LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PARALIA_INSTALL_LIB})

set (PARALIA_INSTALL_INC ${PARALIA_INSTALL}/include)
file(MAKE_DIRECTORY ${PARALIA_INSTALL_INC})

set (PARALIA_INSTALL_BIN ${PARALIA_INSTALL}/bin)
file(MAKE_DIRECTORY ${PARALIA_INSTALL_BIN})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PARALIA_INSTALL_BIN})

#These should be removed/rechecked at some point...
file(MAKE_DIRECTORY ${PARALIA_INSTALL}/Database)
file(MAKE_DIRECTORY ${PARALIA_INSTALL}/Database/microbenchmarks)

add_definitions(-DSCONFIG_PATH="$ENV{PARALIA_SCONFIG_FILE_PATH}")
message( "PARALiA-GEMMex - Backend - SCONFIG_PATH: " $ENV{PARALIA_SCONFIG_FILE_PATH})

message( "------------------------------------------------------------------------------------------------")
ADD_SUBDIRECTORY (Backend)
ADD_SUBDIRECTORY (Autotuner)
ADD_SUBDIRECTORY (Library)
ADD_SUBDIRECTORY (Benchmarking)

configure_file(${PARALIA_MICROBENCH}/Run_microbenchmarks.in ${PARALIA_INSTALL}/Run_microbenchmarks.sh @ONLY)
configure_file(${BENCH_DIR_BASE}/Targeted_performance_evaluation.in ${PARALIA_INSTALL}/Targeted_performance_evaluation.sh @ONLY)
configure_file(${BENCH_DIR_BASE}/Strong_scaling_evaluation.in ${PARALIA_INSTALL}/Strong_scaling_evaluation.sh @ONLY)

