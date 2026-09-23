cmake_minimum_required(VERSION 3.16)

project(MakeGAGGStoppingTables_TestEm0Style)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(make_gagg_stopping_tables
  make_gagg_stopping_tables.cc
)

target_link_libraries(make_gagg_stopping_tables
  ${Geant4_LIBRARIES}
)

set_directory_properties(PROPERTIES
  ADDITIONAL_CLEAN_FILES
  "../electron_gagg_stopping.csv;../alpha_gagg_stopping.csv"
)
