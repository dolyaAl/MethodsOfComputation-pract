cmake_minimum_required(VERSION 3.19)

set(EXMP1_NAME NSolve)

add_executable(${EXMP1_NAME}
	source/main.cpp)
target_link_libraries(${EXMP1_NAME} NonLin_solve)
