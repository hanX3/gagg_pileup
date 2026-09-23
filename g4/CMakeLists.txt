cmake_minimum_required(VERSION 3.16)
project(GAGGWaveformSim)

find_package(Geant4 REQUIRED ui_all vis_all)
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)

include(${Geant4_USE_FILE})
include(${ROOT_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIRS})

# Create the data directory outside build at CMake configure time.
# The program writes to Constants.hh::DATAPATH, which is "../data" when running from build/.
set(GAGG_DATA_DIR ${PROJECT_SOURCE_DIR}/data)
if(NOT EXISTS ${GAGG_DATA_DIR})
  file(MAKE_DIRECTORY ${GAGG_DATA_DIR})
endif()

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# main.cc is intentionally placed at project root, not inside src/.
# If an old src/main.cc remains, exclude it to avoid duplicate main() definitions.
list(REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/src/main.cc)

add_executable(gagg ${PROJECT_SOURCE_DIR}/main.cc ${sources} ${headers})
target_link_libraries(gagg ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

# CMake reserves the built-in target name "clean". The normal "make clean" still
# removes compiled products. This extra target removes the whole build contents.
add_custom_target(clean_build
  COMMAND ${CMAKE_COMMAND} -E echo "Removing all files under ${CMAKE_BINARY_DIR}"
  COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_BINARY_DIR}/*
)
