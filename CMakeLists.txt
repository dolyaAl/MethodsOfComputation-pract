cmake_minimum_required(VERSION 3.19)

set(PROJECT_NAME NonLinear)
project(${PROJECT_NAME})

add_subdirectory(fparser)
add_subdirectory(ev3)
add_subdirectory(core)
#add_subdirectory(examples/NSolve)
#add_subdirectory(examples/Interpolation)
#add_subdirectory(examples/ReverseInterpolation)
#add_subdirectory(examples/NDiff)
add_subdirectory(examples/NInt)
add_subdirectory(examples/NIntMod)
add_subdirectory(examples/ODE)




set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT NDiff)
