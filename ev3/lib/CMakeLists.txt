
set (headers src/auxiliary.h src/common.h src/exceptions.h  src/expression.h src/operand.h src/parser.h src/tree.h)
set (sources src/parser.cxx src/expression.cxx src/operand.cxx)

add_library (ev3 STATIC ${sources} ${headers})

install (TARGETS ev3
         RUNTIME DESTINATION bin
         LIBRARY DESTINATION lib${LIB_SUFFIX}
         ARCHIVE DESTINATION lib${LIB_SUFFIX})

target_include_directories(ev3 PUBLIC src)

