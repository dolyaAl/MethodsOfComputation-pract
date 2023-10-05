cmake_minimum_required(VERSION 3.19)

set(EXMP_NAME Interpolation)

add_executable(${EXMP_NAME}
	source/main.cpp)
	
target_link_libraries(${EXMP_NAME} NonLin_solve)