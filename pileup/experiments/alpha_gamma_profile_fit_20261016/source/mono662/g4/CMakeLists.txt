cmake_minimum_required(VERSION 3.16)
project(GAGGWaveformSim)

find_package(Geant4 REQUIRED ui_all vis_all)
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)
find_package(nlohmann_json 3.2 QUIET)
if(NOT nlohmann_json_FOUND)
  # ROOT distributions may bundle the header without a separate CMake package.
  find_path(GAGG_JSON_INCLUDE_DIR nlohmann/json.hpp HINTS ${ROOT_INCLUDE_DIRS})
  if(NOT GAGG_JSON_INCLUDE_DIR)
    message(FATAL_ERROR "nlohmann_json >= 3.2 is required (CMake package or ROOT-bundled header)")
  endif()
  file(STRINGS "${GAGG_JSON_INCLUDE_DIR}/nlohmann/json.hpp" GAGG_JSON_VERSION_LINES
    REGEX "^#define NLOHMANN_JSON_VERSION_(MAJOR|MINOR|PATCH) ")
  foreach(component MAJOR MINOR PATCH)
    string(REGEX MATCH "NLOHMANN_JSON_VERSION_${component} ([0-9]+)" version_match "${GAGG_JSON_VERSION_LINES}")
    if(NOT version_match)
      message(FATAL_ERROR "Cannot determine nlohmann_json header version; install its CMake package")
    endif()
    set(GAGG_JSON_${component} "${CMAKE_MATCH_1}")
  endforeach()
  if("${GAGG_JSON_MAJOR}.${GAGG_JSON_MINOR}.${GAGG_JSON_PATCH}" VERSION_LESS "3.2.0")
    message(FATAL_ERROR "nlohmann_json >= 3.2 is required")
  endif()
  add_library(nlohmann_json::nlohmann_json INTERFACE IMPORTED)
  set_target_properties(nlohmann_json::nlohmann_json PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${GAGG_JSON_INCLUDE_DIR}")
  message(STATUS "Using nlohmann_json ${GAGG_JSON_MAJOR}.${GAGG_JSON_MINOR}.${GAGG_JSON_PATCH} headers: ${GAGG_JSON_INCLUDE_DIR}")
endif()

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
target_compile_features(gagg PRIVATE cxx_std_17)
target_link_libraries(gagg ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} nlohmann_json::nlohmann_json)

include(CTest)
if(BUILD_TESTING)
  add_executable(test_sipm_response tests/test_sipm_response.cc src/SiPMResponse.cc)
  target_compile_features(test_sipm_response PRIVATE cxx_std_17)
  target_link_libraries(test_sipm_response PRIVATE nlohmann_json::nlohmann_json)
  add_test(NAME sipm_response COMMAND test_sipm_response ${PROJECT_SOURCE_DIR}/configs/sipm_hamamatsu_s13360_6050pe_v2.json)
endif()

# CMake reserves the built-in target name "clean". The normal "make clean" still
# removes compiled products. This extra target removes the whole build contents.
add_custom_target(clean_build
  COMMAND ${CMAKE_COMMAND} -E echo "Removing all files under ${CMAKE_BINARY_DIR}"
  COMMAND ${CMAKE_COMMAND} -E rm -rf ${CMAKE_BINARY_DIR}/*
)

# Also let the standard "make clean" remove the executable.
set_directory_properties(PROPERTIES
  ADDITIONAL_MAKE_CLEAN_FILES "${CMAKE_BINARY_DIR}/gagg"
)
