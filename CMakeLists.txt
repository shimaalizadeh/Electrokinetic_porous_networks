cmake_minimum_required(VERSION 3.20)
project(ek_network LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(LAPACK REQUIRED)

add_executable(ek_solver
  main.cpp
  network_solver.cpp
  network_info.cpp
  table.cpp
  mylapack.cpp
)

target_link_libraries(ek_solver PRIVATE LAPACK::LAPACK)
target_compile_options(ek_solver PRIVATE -Wall -Wextra -Wpedantic)

option(EK_SANITIZE "Build with Address + UB sanitizers" OFF)
if(EK_SANITIZE)
  target_compile_options(ek_solver PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer -g)
  target_link_options(ek_solver PRIVATE -fsanitize=address,undefined)
endif()
