cmake_minimum_required(VERSION 3.19)

set(PROJECT_NAME NonLinear)
project(${PROJECT_NAME})

add_subdirectory(core)
add_subdirectory(examples/NSolve)
add_subdirectory(examples/Interpolation)

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Interpolation)
