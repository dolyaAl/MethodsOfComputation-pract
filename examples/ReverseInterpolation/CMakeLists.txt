cmake_minimum_required(VERSION 3.19)

set(EXMP3_NAME RevInterpolation)

add_executable(${EXMP3_NAME}
	source/main.cpp)
target_link_libraries(${EXMP3_NAME} NonLin_solve)
