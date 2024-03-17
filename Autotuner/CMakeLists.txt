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

message( "------------------------------------------------------------------------------------------------")
