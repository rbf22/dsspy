find_package(Python REQUIRED COMPONENTS Interpreter Development)
find_package(Boost REQUIRED COMPONENTS python)

add_library(DSSP SHARED dssp-python-plugin.cpp)

target_compile_features(DSSP PUBLIC cxx_std_20)

target_include_directories(DSSP PRIVATE ${Python_INCLUDE_DIRS})
target_link_libraries(DSSP dssp::dssp Boost::python ${Python_LIBRARIES})