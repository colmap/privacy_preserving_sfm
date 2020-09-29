cmake_minimum_required (VERSION 3.1)
project(re3q3)

if(NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

find_package (Eigen3 REQUIRED NO_MODULE)

add_library(re3q3 INTERFACE)
target_include_directories(re3q3 INTERFACE re3q3)
target_link_libraries(re3q3 INTERFACE Eigen3::Eigen)
target_compile_features(re3q3 INTERFACE cxx_std_14)

add_executable(test_re3q3 test_re3q3.cpp)
target_link_libraries(test_re3q3 re3q3)
target_compile_options(test_re3q3 PRIVATE -Wall -Werror)
