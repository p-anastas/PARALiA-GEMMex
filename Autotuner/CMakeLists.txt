set ( AUTOTUNE_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
message( "Autotuner - AUTOTUNE_DIR: " ${AUTOTUNE_DIR} )

set (AUTOTUNE_DIR_INC ${AUTOTUNE_DIR}/include)
include_directories( ${AUTOTUNE_DIR_INC})
message( "Autotuner - Included DIR: " ${AUTOTUNE_DIR_INC} )

include_directories(${PARALIA_INSTALL_INC})
message( "Autotuner - Included DIR: " ${PARALIA_INSTALL_INC} )

set (AUTOTUNE_DIR_SRC ${AUTOTUNE_DIR}/src)
message( "Autotuner - Source DIR: " ${AUTOTUNE_DIR_SRC} )

add_library(autotuner SHARED
${AUTOTUNE_DIR_SRC}/TaskDistribution.cpp
${AUTOTUNE_DIR_SRC}/Autotuner.cpp
${AUTOTUNE_DIR_SRC}/Routing.cpp
)

target_link_libraries(autotuner ${INTERNAL_LINK_LIB} chl_grid_amalgamation)
configure_file(${AUTOTUNE_DIR_INC}/Autotuner.hpp ${PARALIA_INSTALL_INC}/Autotuner.hpp COPYONLY)

###########################################################################################
#------------------------------------microbenchmarks--------------------------------------#
###########################################################################################

## TODO: Choose between ->
## i) Automatic micro-benchmarking (Population checked within 95% CI) using boost (recommended if boost available).
message("Autotuner - microbenchmarks - Using Boost for microbenchmark statistics")
add_definitions(-DAUTO_BENCH_USE_BOOST)
set(BOOSTO_PREFIX $ENV{PARALIA_BOOST_PREFIX})
set (BOOSTO_INCLUDE_DIRS ${BOOSTO_PREFIX}/include)
include_directories(${BOOSTO_INCLUDE_DIRS})
message("Autotuner - microbenchmarks - Included (Boost) DIR: " ${BOOSTO_INCLUDE_DIRS})
## OR
## ii) Benchmark with a predifined itteration number for each problem size.
#add_definitions(-DITER=100)
#message("Database Builder - Performing predefined microbenchmark itterations")

include_directories(${PARALIA_INSTALL_INC})
message( "Autotuner - microbenchmarks - Included DIR: " ${PARALIA_INSTALL_INC} )

include_directories(${AUTOTUNE_DIR}/microbenchmarks)
message( "Autotuner - microbenchmarks - Included microbenchmarks DIR: " ${AUTOTUNE_DIR}/microbenchmarks)

set(NVEMAPI_PREFIX "${AUTOTUNE_DIR}/microbenchmarks/nvidia-energy-measure")
set (NVEMAPI_INCLUDE_DIRS ${NVEMAPI_PREFIX}/include) # -L${CUDA_LD}
include_directories(${NVEMAPI_INCLUDE_DIRS} )
message( "Autotuner - microbenchmarks - NVEMAPI_INCLUDE_DIRS: " ${NVEMAPI_INCLUDE_DIRS} )
message( "------------------------------------------------------------------------------------------------")
ADD_SUBDIRECTORY (${NVEMAPI_PREFIX})

## Microbenchmarks for DB
add_executable (bw_bench_broadcast_2D ${AUTOTUNE_DIR}/microbenchmarks/bw_bench_broadcast_2D.cpp ${AUTOTUNE_DIR}/microbenchmarks/helpers.cpp)

target_link_libraries(bw_bench_broadcast_2D ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)

## Testers 
add_executable (CHLTestGamalgGEMM ${AUTOTUNE_DIR}/microbenchmarks/chl_test_gamalg_gemm.cpp ${AUTOTUNE_DIR}/microbenchmarks/helpers.cpp)
add_executable (CHLTestMemlocs ${AUTOTUNE_DIR}/microbenchmarks/chl_test_memloc_numbering.cpp ${AUTOTUNE_DIR}/microbenchmarks/helpers.cpp)

target_link_libraries(CHLTestGamalgGEMM ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)
target_link_libraries(CHLTestMemlocs ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)