cmake_minimum_required(VERSION 3.10)
project(Info_Lab_4term)
set(CMAKE_CXX_STANDARD 17)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_library(GMSH_LIB gmsh)
find_path(GMSH_INC gmsh.h)

add_executable(cubic2d cubic2d.cpp)
target_link_libraries(cubic2d ${VTK_LIBRARIES})

add_executable(cubic3d cubic3d.cpp)
target_link_libraries(cubic3d ${VTK_LIBRARIES})

add_executable(tetr3d tetr3d.cpp)
target_include_directories(tetr3d PRIVATE ${GMSH_INC})
target_link_libraries(tetr3d PRIVATE ${VTK_LIBRARIES} ${GMSH_LIB})

add_executable(main main.cpp)
target_include_directories(main PRIVATE ${GMSH_INC})
target_link_libraries(main PRIVATE ${VTK_LIBRARIES} ${GMSH_LIB})
