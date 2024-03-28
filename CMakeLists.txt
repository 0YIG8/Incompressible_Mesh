cmake_minimum_required(VERSION 3.19...3.28)

project(IFD)

set(CMAKE_CXX_COMPILER "clang++")
set(CMAKE_C_COMPILER "clang")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

find_package(Eigen3 REQUIRED)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include_directories(/usr/include/eigen3)

add_executable(main Vector.C Point.C Face.C Cell.C Boundary.C Discretization.C Mesh.C ScalarField.C VectorField.C FileParser.C OutputData.C test_main.C)
target_link_libraries(main PUBLIC Eigen3::Eigen)

execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_BINARY_DIR}/compile_commands.json
                ${CMAKE_SOURCE_DIR}/compile_commands.json)