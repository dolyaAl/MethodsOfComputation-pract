cmake_minimum_required(VERSION 3.19)

set(EXMP2_NAME Interpolation)

add_executable(${EXMP2_NAME}
	source/main.cpp)
target_link_libraries(${EXMP2_NAME} NonLin_solve)
