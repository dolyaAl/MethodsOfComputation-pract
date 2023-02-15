cmake_minimum_required(VERSION 3.19)

set(INTERFACE_NAME interface)

add_executable(${INTERFACE_NAME}
	source/main.cpp)
target_link_libraries(${INTERFACE_NAME} NonLin_solve)
