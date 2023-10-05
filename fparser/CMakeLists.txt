cmake_minimum_required(VERSION 3.19)

set(FPARSER fparser)

set(INCLUDES_PUB 
	includes/fparser.hh
	fpconfig.hh
	extrasrc/fptypes.hh)

set(SOURCES_PRIV
	fparser.cc 
	fpoptimizer.cc)

add_library(${FPARSER} STATIC 
${INCLUDES_PUB}
${SOURCES_PRIV})

target_include_directories(${FPARSER} PUBLIC includes)
target_compile_features(${FPARSER} PUBLIC cxx_std_17)
