set(NVMLAPI_PREFIX ${CMAKE_CURRENT_SOURCE_DIR})
set (NVMLAPI_SOURCE_DIR ${NVMLAPI_PREFIX}/src)
set (NVMLAPI_INCLUDE_DIRS ${NVMLAPI_PREFIX}/include)
include_directories(${NVMLAPI_INCLUDE_DIRS} )
message( "NVML-API - NVMLAPI_INCLUDE_DIRS: " ${NVMLAPI_INCLUDE_DIRS} )

set (NVMLAPI_LINK "-L${CUDA_LD} -lnvidia-ml -L/usr/lib64 -lcuda -lpthread")
message( "NVML-API - NVMLAPI_LINK: " ${NVMLAPI_LINK} )

file(GLOB SOURCE_FILES ${NVMLAPI_SOURCE_DIR}/*.c*)
#message( "NVML-API - SOURCE_FILES: " ${SOURCE_FILES} )
add_library(nvidia_powa ${SOURCE_FILES})
target_link_libraries(nvidia_powa ${INTERNAL_LINK_LIB} ${NVMLAPI_LINK})

add_executable (example ${NVMLAPI_PREFIX}/example_main.cu)
target_link_libraries(example ${INTERNAL_LINK_LIB} nvidia_powa)
message( "------------------------------------------------------------------------------------------------")