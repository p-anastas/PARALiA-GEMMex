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

configure_file(${PARALIA_BACKEND_SOURCE_DIR}/smart_wrappers/chl_smart_wrappers.hpp ${PARALIA_INSTALL_INC}/chl_smart_wrappers.hpp COPYONLY)
configure_file(${PARALIA_BACKEND_SOURCE_DIR}/grid_amalgamation/chl_grid_amalgamation.hpp ${PARALIA_INSTALL_INC}/chl_grid_amalgamation.hpp COPYONLY)
configure_file(${PARALIA_BACKEND_SOURCE_DIR}/smart_wrappers/${BACKEND}/backend_wrappers.hpp ${PARALIA_INSTALL_INC}/backend_wrappers.hpp COPYONLY)

message( "------------------------------------------------------------------------------------------------")
