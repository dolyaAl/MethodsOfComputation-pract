cmake_minimum_required(VERSION 3.19)

set(EXMP4_NAME NInt)

add_executable(${EXMP4_NAME}
	source/main.cpp)
target_link_libraries(${EXMP4_NAME} NonLin_solve)
