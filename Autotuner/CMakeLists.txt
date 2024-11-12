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

## i) Automatic micro-benchmarking (Population checked within 95% CI) using boost.
message("Autotuner - microbenchmarks - Using Boost for microbenchmark statistics")
if(NOT $ENV{PARALIA_BOOST_PREFIX} STREQUAL "default")
    set(BOOSTO_PREFIX $ENV{PARALIA_BOOST_PREFIX})
    set (BOOSTO_INCLUDE_DIRS ${BOOSTO_PREFIX}/include)
    include_directories(${BOOSTO_INCLUDE_DIRS})
    message("Autotuner - microbenchmarks - Included (Boost) DIR: " ${BOOSTO_INCLUDE_DIRS})
endif()

include_directories(${PARALIA_INSTALL_INC})
message( "Autotuner - microbenchmarks - Included DIR: " ${PARALIA_INSTALL_INC} )

set (PARALIA_MICROBENCH ${AUTOTUNE_DIR}/microbenchmarks)
set (PARALIA_MICROBENCH ${PARALIA_MICROBENCH} PARENT_SCOPE)
include_directories(${PARALIA_MICROBENCH})
message( "Autotuner - microbenchmarks - Included microbenchmarks DIR: " ${PARALIA_MICROBENCH})

set(NVEMAPI_PREFIX ${PARALIA_MICROBENCH}/nvidia-energy-measure)
set (NVEMAPI_INCLUDE_DIRS ${NVEMAPI_PREFIX}/include) # -L${CUDA_LD}
include_directories(${NVEMAPI_INCLUDE_DIRS} )
message( "Autotuner - microbenchmarks - NVEMAPI_INCLUDE_DIRS: " ${NVEMAPI_INCLUDE_DIRS} )
message( "------------------------------------------------------------------------------------------------")
ADD_SUBDIRECTORY (${NVEMAPI_PREFIX})

## Microbenchmarks for DB
add_executable (bw_bench_scatter_2D ${PARALIA_MICROBENCH}/bw_bench_scatter_2D.cpp ${PARALIA_MICROBENCH}/helpers.cpp)
add_executable (bw_bench_gather_2D ${PARALIA_MICROBENCH}/bw_bench_gather_2D.cpp ${PARALIA_MICROBENCH}/helpers.cpp)
add_executable (bw_bench_grid_2D ${PARALIA_MICROBENCH}/bw_bench_grid_2D.cpp ${PARALIA_MICROBENCH}/helpers.cpp)

target_link_libraries(bw_bench_scatter_2D ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)
target_link_libraries(bw_bench_gather_2D ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)
target_link_libraries(bw_bench_grid_2D ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)

## Testers 
add_executable (CHLTestGamalgGEMM ${PARALIA_MICROBENCH}/chl_test_gamalg_gemm.cpp ${PARALIA_MICROBENCH}/helpers.cpp)
add_executable (CHLTestMemlocs ${PARALIA_MICROBENCH}/chl_test_memloc_numbering.cpp ${PARALIA_MICROBENCH}/helpers.cpp)

target_link_libraries(CHLTestGamalgGEMM ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)
target_link_libraries(CHLTestMemlocs ${INTERNAL_LINK_LIB} chl_grid_amalgamation nvidia_powa)