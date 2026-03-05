cmake_minimum_required(VERSION 3.21)

project(Lab3_FEniCSx LANGUAGES C CXX)

# 1. Поиск библиотеки DOLFINx
if(NOT TARGET dolfinx)
  find_package(DOLFINX REQUIRED)
endif()

# 2. Определение типов данных PETSc (нужно для правильной генерации кода FFCx)
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

# 3. Генерация C-кода из Python-описаний вариационных форм (UFL)

# Для задачи Пуассона
add_custom_command(
  OUTPUT poisson.c poisson.h
  COMMAND ffcx ${CMAKE_CURRENT_SOURCE_DIR}/poisson.py ${SCALAR_TYPE}
  VERBATIM
  DEPENDS poisson.py
  COMMENT "Compile poisson.py using FFCx"
)

# Для задачи упругости (Ам Ням)
add_custom_command(
  OUTPUT elasticity.c elasticity.h
  COMMAND ffcx ${CMAKE_CURRENT_SOURCE_DIR}/elasticity.py ${SCALAR_TYPE}
  VERBATIM
  DEPENDS elasticity.py
  COMMENT "Compile elasticity.py using FFCx"
)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# 4. Создание ПЕРВОГО исполняемого файла (Poisson)
add_executable(demo_poisson main.cpp poisson.c)
target_link_libraries(demo_poisson dolfinx)
target_compile_features(demo_poisson PUBLIC cxx_std_20)

# 5. Создание ВТОРОГО исполняемого файла (Elasticity / Om Nom)
add_executable(demo_elasticity elasticity.cpp elasticity.c)
target_link_libraries(demo_elasticity dolfinx)
target_compile_features(demo_elasticity PUBLIC cxx_std_20)

# 6. Настройки компилятора (предупреждения)
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wno-comment" HAVE_NO_MULTLINE)
if(HAVE_NO_MULTLINE)
  target_compile_options(demo_poisson PRIVATE "-Wno-comment")
  target_compile_options(demo_elasticity PRIVATE "-Wno-comment")
endif()

# 7. Тесты (опционально, если хочешь запускать через ctest)
add_test(NAME test_poisson COMMAND demo_poisson)
add_test(NAME test_elasticity COMMAND demo_elasticity)