cmake_minimum_required(VERSION 3.13)
project (Matrix CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


add_library(MAT OBJECT ${Matrix_SOURCE_DIR}/include/matrix_basic.cpp ${Matrix_SOURCE_DIR}/include/matrix_operations.cpp)

target_include_directories(MAT PUBLIC ${Matrix_SOURCE_DIR}/include)