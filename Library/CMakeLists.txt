ADD_SUBDIRECTORY (Software-assisted_Memory_Buffer)

set ( LIBSC_DIR ${CMAKE_CURRENT_SOURCE_DIR} )
message( "Library - LIBSC_DIR: " ${LIBSC_DIR} )

set (LIBSC_DIR_INC ${LIBSC_DIR}/include)
include_directories( ${LIBSC_DIR_INC})
message( "Library - Included DIR: " ${LIBSC_DIR_INC} )

set (LIBSC_DIR_INC_BACK ${LIBSC_DIR}/${BACKEND}/)
include_directories( ${LIBSC_DIR_INC_BACK})
message( "Library - Included DIR: " ${LIBSC_DIR_INC_BACK} )

include_directories(${PARALIA_INSTALL_INC})
message( "Library - Included DIR: " ${PARALIA_INSTALL_INC} )

set (LIBSC_DIR_SRC ${LIBSC_DIR}/src)
message( "Library - Source DIR: " ${LIBSC_DIR_SRC} )

set ( BACKEND_LIBSC_DIR ${LIBSC_DIR}/${BACKEND} )
message( "Library - BACKEND_LIBSC_DIR: " ${BACKEND_LIBSC_DIR} )

add_library(paralia SHARED
${LIBSC_DIR_SRC}/extern_definitions.cpp
${LIBSC_DIR_SRC}/DataTile.cpp
${LIBSC_DIR_SRC}/Decomposer.cpp
${LIBSC_DIR_SRC}/Resource_manager.cpp
${LIBSC_DIR_SRC}/Functions/dgemm.cpp
${LIBSC_DIR_SRC}/Functions/sgemm.cpp
${LIBSC_DIR_SRC}/Functions/hgemm.cpp
${LIBSC_DIR}/${BACKEND}/Decom_backend.cpp
)

target_link_libraries(paralia ${INTERNAL_LINK_LIB} autotuner chl_grid_amalgamation caching)
configure_file(${LIBSC_DIR_INC}/PARALiA.hpp ${PARALIA_INSTALL_INC}/PARALiA.hpp COPYONLY)
configure_file(${LIBSC_DIR_INC}/Resource_manager.hpp ${PARALIA_INSTALL_INC}/Resource_manager.hpp COPYONLY)
configure_file(${LIBSC_DIR_INC}/Decomposer.hpp ${PARALIA_INSTALL_INC}/Decomposer.hpp COPYONLY)
configure_file(${LIBSC_DIR_INC}/DataTile.hpp ${PARALIA_INSTALL_INC}/DataTile.hpp COPYONLY)

message( "------------------------------------------------------------------------------------------------")
