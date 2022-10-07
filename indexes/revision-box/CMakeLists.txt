cmake_minimum_required(VERSION 3.16)

project(learnedbench)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


set(CMAKE_CXX_FLAGS "-mbmi2 -O3 -DNDEBUG")
# set(CMAKE_CXX_FLAGS "-mbmi2 -Wall -Wextra")
# set(CMAKE_CXX_FLAGS_DEBUG "-g")
# set(CMAKE_CXX_FLAGS_RELEASE "-O3")

option(PROFILE "whether starting heap profile" OFF)
option(RSMI "whether building RSMI" OFF)


set(BOOST_ROOT "/home/zengyuxiang/boost_1_80_0")
set(Boost_INCLUDE_DIR "/usr/local/include/boost")
set(Boost_LIBRARY_DIR "/usr/local/lib")

SET(EXECUTABLE_OUTPUT_PATH "/home/zengyuxiang/learnedbench/indexes/revision-box/build/bin")

# check dependencies
find_package(TPIE REQUIRED)
find_package(Boost 1.79.0 REQUIRED COMPONENTS program_options system filesystem)
find_package(GEOS REQUIRED)


# the other benchmarks are compiled using cxx17 (required by PGM-Index)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${Boost_INCLUDE_DIR} ${TPIE_INCLUDE_DIRS})

add_executable(main main.cpp)
target_link_libraries(main Boost::program_options GEOS::geos)
  