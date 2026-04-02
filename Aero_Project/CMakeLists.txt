cmake_minimum_required(VERSION 3.21)
project(Aero_Project LANGUAGES C CXX)

if(NOT TARGET dolfinx)
  find_package(DOLFINX REQUIRED)
endif()

include(CheckSymbolExists)
set(CMAKE_REQUIRED_INCLUDES ${PETSC_INCLUDE_DIRS})
check_symbol_exists(PETSC_USE_COMPLEX petscsystypes.h PETSC_SCALAR_COMPLEX)
check_symbol_exists(PETSC_USE_REAL_DOUBLE petscsystypes.h PETSC_REAL_DOUBLE)

if(PETSC_SCALAR_COMPLEX EQUAL 1)
  if(PETSC_REAL_DOUBLE EQUAL 1)
    set(SCALAR_TYPE "--scalar_type=complex128")
  else()
    set(SCALAR_TYPE "--scalar_type=complex64")
  endif()
else()
  if(PETSC_REAL_DOUBLE EQUAL 1)
    set(SCALAR_TYPE "--scalar_type=float64")
  else()
    set(SCALAR_TYPE "--scalar_type=float32")
  endif()
endif()

add_custom_command(
  OUTPUT navier_stokes.c navier_stokes.h
  COMMAND ffcx ${CMAKE_CURRENT_SOURCE_DIR}/navier_stokes.py ${SCALAR_TYPE}
  VERBATIM
  DEPENDS navier_stokes.py
  COMMENT "Compile navier_stokes.py using FFCx"
)

add_custom_command(
  OUTPUT navier_stokes_chorin.c navier_stokes_chorin.h
  COMMAND ffcx ${CMAKE_CURRENT_SOURCE_DIR}/navier_stokes_chorin.py ${SCALAR_TYPE}
  VERBATIM
  DEPENDS navier_stokes_chorin.py
  COMMENT "Compile navier_stokes_chorin.py using FFCx"
)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# И обнови цель сборки исполняемого файла:
add_executable(demo_cfd main.cpp navier_stokes_chorin.c)

target_link_libraries(demo_cfd dolfinx)
target_compile_features(demo_cfd PUBLIC cxx_std_20)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wno-comment" HAVE_NO_MULTLINE)
if(HAVE_NO_MULTLINE)
  target_compile_options(demo_cfd PRIVATE "-Wno-comment")
endif()