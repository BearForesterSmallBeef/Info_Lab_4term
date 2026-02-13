cmake_minimum_required(VERSION 3.10)
project(Info_Lab_4term)
set(CMAKE_CXX_STANDARD 17)

find_library(GMSH_LIB gmsh)
find_path(GMSH_INC gmsh.h)

add_executable(torus torus.cpp)
target_include_directories(torus PRIVATE ${GMSH_INC})

target_link_libraries(torus PRIVATE ${GMSH_LIB})

add_executable(cube cube.cpp)
target_include_directories(cube PRIVATE ${GMSH_INC})

target_link_libraries(cube PRIVATE ${GMSH_LIB})

add_executable(circule circule.cpp)
target_include_directories(circule PRIVATE ${GMSH_INC})

target_link_libraries(circule PRIVATE ${GMSH_LIB})

add_executable(cilinder cilinder.cpp)
target_include_directories(cilinder PRIVATE ${GMSH_INC})

target_link_libraries(cilinder PRIVATE ${GMSH_LIB})

add_executable(stl stl.cpp)
target_include_directories(stl PRIVATE ${GMSH_INC})

target_link_libraries(stl PRIVATE ${GMSH_LIB})
