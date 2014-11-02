
#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(USolids)

set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_MODULE_PATH})

include_directories(
	${PROJECT_SOURCE_DIR}/include 
)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#

# ??? how to group header files, *.hh to "Header Files" group

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

source_group ("Header Files" FILES ${headers})
source_group ("Source Files" FILES ${sources})

add_library(usolids  ${sources} ${headers})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#

#link_directories (${LINK_DIRECTORIES})

target_link_libraries(usolids ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
###install(TARGETS exampleN01 DESTINATION bin)
