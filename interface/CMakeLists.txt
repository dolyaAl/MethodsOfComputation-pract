cmake_minimum_required(VERSION 3.19)

set(INTERFACE_NAME interface)

add_executable(${INTERFACE_NAME}
	source/main.cpp)
target_link_libraries(${INTERFACE_NAME} NonLin_solve)
target_compile_features(${ENGINE_PROJECT_NAME} PUBLIC cxx_std_17)
