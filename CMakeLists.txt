cmake_minimum_required(VERSION 3.0)

project(cfd26)

set(CMAKE_CXX_STANDARD 20)

include(FetchContent)

add_compile_options(
  # $<$<CONFIG:Release>:-Werror> $<$<CONFIG:Debug>:-Wno-error>
  # $<$<CONFIG:RelWithDebInfo>:-Werror>
  -Werror
  -Wall
  -Wextra
  -Wpedantic
  -Wconversion
  -Wunused-const-variable
  -Wduplicated-cond
  -Wnon-virtual-dtor
  -Woverloaded-virtual
  -Wredundant-decls
  -Wsign-promo
  -Wzero-as-null-pointer-constant
  -Wsuggest-override
  -Wshadow
  -Wfloat-equal
  -Wcast-qual
  -Wlogical-op
  -Wcast-align
  # -Wold-style-cast
  -Wmissing-field-initializers)

# Boost search paths
find_package(Boost REQUIRED)
# libformat
find_package(fmt REQUIRED)
# Catch framework
find_package(Catch2 REQUIRED)

# ------ amgcl
FetchContent_Declare(
  amgcl
  GIT_REPOSITORY https://github.com/ddemidov/amgcl
  GIT_TAG 1.4.2)
FetchContent_GetProperties(amgcl)
if(NOT amgcl_POPULATED)
  FetchContent_Populate(amgcl)
endif()
include_directories(SYSTEM ${amgcl_SOURCE_DIR})

# build all binaries to build/bin
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(src)
