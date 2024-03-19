cmake_minimum_required(VERSION 3.16...3.21)
project(123)

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(123 123.cc ${sources} ${headers})
target_link_libraries(123 ${Geant4_LIBRARIES})

set(123_SCRIPTS
  data.out
  123.in
  init_vis.mac
  run.mac
  run1.mac
  vis.mac
  )

foreach(_script ${123_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#add_custom_target(ex2 DEPENDS 123)

install(TARGETS 123 DESTINATION bin)
