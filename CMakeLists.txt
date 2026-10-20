cmake_minimum_required(VERSION 3.24)
project(itchbook VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(ITCHBOOK_BUILD_TESTS "Build unit and property tests" ON)
option(ITCHBOOK_BUILD_BENCH "Build Google Benchmark microbenchmarks" ON)
option(ITCHBOOK_BUILD_FUZZ "Build libFuzzer targets (clang only)" OFF)
option(ITCHBOOK_NATIVE "Tune for the build machine (-march/-mcpu=native)" OFF)
option(ITCHBOOK_WERROR "Treat warnings as errors" OFF)
option(ITCHBOOK_HASH_INDEX "FastBook: hash-table order index instead of the paged direct index (A/B benchmarks)" OFF)
set(ITCHBOOK_SANITIZE "" CACHE STRING "Comma-separated sanitizers, e.g. address,undefined")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(Warnings)
include(Sanitizers)

find_package(ZLIB REQUIRED)

# Core library: ITCH decoding and streaming input. Everything else links this.
add_library(itchbook_core
  src/source.cpp
  src/reference_book.cpp
  src/fast_book.cpp
)
add_library(itchbook::core ALIAS itchbook_core)
target_include_directories(itchbook_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(itchbook_core PUBLIC ZLIB::ZLIB)
itchbook_set_warnings(itchbook_core)
if(ITCHBOOK_HASH_INDEX)
  target_compile_definitions(itchbook_core PUBLIC ITCHBOOK_HASH_INDEX)
endif()

if(ITCHBOOK_NATIVE)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-march=native" HAS_MARCH_NATIVE)
  if(HAS_MARCH_NATIVE)
    target_compile_options(itchbook_core PUBLIC -march=native)
  else()
    target_compile_options(itchbook_core PUBLIC -mcpu=native)
  endif()
endif()

add_executable(itch-stats tools/itch_stats.cpp)
target_link_libraries(itch-stats PRIVATE itchbook::core)
itchbook_set_warnings(itch-stats)

find_package(Threads REQUIRED)
add_executable(itch-bench tools/itch_bench.cpp)
target_link_libraries(itch-bench PRIVATE itchbook::core Threads::Threads)
itchbook_set_warnings(itch-bench)

add_executable(itch-book tools/itch_book.cpp)
target_link_libraries(itch-book PRIVATE itchbook::core)
itchbook_set_warnings(itch-book)

add_executable(itch-queue tools/itch_queue.cpp)
target_link_libraries(itch-queue PRIVATE itchbook::core)
itchbook_set_warnings(itch-queue)

add_executable(itch-lob tools/itch_lob.cpp)
target_link_libraries(itch-lob PRIVATE itchbook::core)
itchbook_set_warnings(itch-lob)

if(ITCHBOOK_BUILD_TESTS OR ITCHBOOK_BUILD_BENCH)
  include(Dependencies)
endif()

if(ITCHBOOK_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(ITCHBOOK_BUILD_BENCH)
  add_subdirectory(bench)
endif()

if(ITCHBOOK_BUILD_FUZZ)
  add_subdirectory(fuzz)
endif()
