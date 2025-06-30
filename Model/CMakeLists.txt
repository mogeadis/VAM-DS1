# Project Setup
cmake_minimum_required(VERSION 3.22)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(RELEASE_FLAGS /O2 /Ob3 /DNDEBUG /arch:AVX2)
project(VAM-DS1 VERSION 1.0.0)

# Core Library
message(STATUS "Preparing Core Library...")
set(CORE_DIR Core)
set(CORE_TARGET Core)
add_subdirectory(${CORE_DIR} ${CMAKE_BINARY_DIR}/${CORE_TARGET})

# JUCE Plugin
message(STATUS "Preparing JUCE Plugin...")
set(PLUGIN_DIR Plugin)
set(PLUGIN_TARGET Plugin)
add_subdirectory(${PLUGIN_DIR} ${CMAKE_BINARY_DIR}/${PLUGIN_TARGET})

# Main Executable
message(STATUS "Preparing Main Executable...")
set(MAIN_TARGET Main)
add_executable(${MAIN_TARGET} main.cpp)
target_link_libraries(${MAIN_TARGET} PRIVATE ${CORE_TARGET})