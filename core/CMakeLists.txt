cmake_minimum_required(VERSION 3.19)

set(CORE_NAME NonLin_solve)

set(INCLUDES_PUB 
	includes/Functions.h
	includes/FunctionReader.h)

set(SOURCES_PRIV
	source/Functions.cpp
	source/FunctionReader.cpp)

add_library(${CORE_NAME} STATIC 
${INCLUDES_PUB}
${SOURCES_PRIV} $<TARGET_OBJECTS:MAT>) 

add_subdirectory(Matrix )

target_link_libraries(${CORE_NAME} fparser)
target_link_libraries(${CORE_NAME} MAT)
target_link_libraries(${CORE_NAME} ev3)



target_include_directories(${CORE_NAME} PUBLIC includes)
target_compile_features(${CORE_NAME} PUBLIC cxx_std_17)

