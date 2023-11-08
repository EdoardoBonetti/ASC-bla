# cmake_minimum_required(VERSION 3.20)
# project(ASC_bla)
#
# set(CMAKE_CXX_STANDARD 17)
#
# include_directories(src)
#
# find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)
#
# execute_process(
# COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
# OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE PYBIND11_DIR)
# list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_DIR}")
# find_package(pybind11 CONFIG REQUIRED)
#
#
# add_executable(test_vector tests/test_vector.cc)
# add_executable(test_matrix tests/test_matrix.cc)
# add_executable(test_lapack tests/test_lapack.cc)
# include_directories(tests)
#
# target_link_libraries(test_lapack PUBLIC ${LAPACK_LIBRARIES})
#
#
# pybind11_add_module(bla src/bind_bla.cpp)
#
# install(TARGETS bla DESTINATION ASCsoft)
# install(FILES src/vector.h DESTINATION ASCsoft/include)
# install(FILES src/matrix.h DESTINATION ASCsoft/include)

cmake_minimum_required(VERSION 3.20)
project(ASC_bla)

set(CMAKE_CXX_STANDARD 17)

include_directories(src)

find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE PYBIND11_DIR)
list(APPEND CMAKE_PREFIX_PATH "${PYBIND11_DIR}")
find_package(pybind11 CONFIG REQUIRED)

add_executable(test_vector tests/test_vector.cc)

find_package(LAPACK REQUIRED)
add_executable(test_lapack tests/test_lapack.cc)
target_link_libraries(test_lapack PUBLIC ${LAPACK_LIBRARIES})

pybind11_add_module(bla src/bind_bla.cpp)

install(TARGETS bla DESTINATION ASCsoft)
install(FILES src/vector.h DESTINATION ASCsoft/include)
