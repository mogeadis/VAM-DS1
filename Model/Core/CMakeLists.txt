# Source Files
file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/${CORE_DIR}/*.cpp)
message(STATUS "Collected source files:")
foreach(FILE IN LISTS SOURCES)
    message(STATUS "  ${FILE}")
endforeach()

# Header Files
file(GLOB_RECURSE HEADERS ${CMAKE_SOURCE_DIR}/${CORE_DIR}/*.h ${CMAKE_SOURCE_DIR}/${CORE_DIR}/*.tpp)
message(STATUS "Collected header files:")
foreach(FILE IN LISTS HEADERS)
    message(STATUS "  ${FILE}")
endforeach()

# Header Directories
set(HEADER_DIRECTORIES "")
foreach(HEADER ${HEADERS})
    get_filename_component(HEADER_DIRECTORY ${HEADER} DIRECTORY)
    list(FIND HEADER_DIRECTORIES ${HEADER_DIRECTORY} INDEX)
    if(INDEX EQUAL -1)
        list(APPEND HEADER_DIRECTORIES ${HEADER_DIRECTORY})
    endif()
endforeach()
message(STATUS "Included header directories:")
foreach(DIRECTORY IN LISTS HEADER_DIRECTORIES)
    message(STATUS "  ${DIRECTORY}")
endforeach()

# Fastor Library
set(FASTOR_DIR Fastor)
set(FASTOR_TARGET Fastor)
set(FASTOR_PATH C:/Users/mogeadis/Projects/Cpp/LIBRARIES/${FASTOR_DIR}) # LOCAL PATH
set(BUILD_TESTING OFF CACHE BOOL "")
add_subdirectory(${FASTOR_PATH} ${CMAKE_BINARY_DIR}/${FASTOR_TARGET})

# Core Library
add_library(${CORE_TARGET} STATIC ${SOURCES})
target_include_directories(${CORE_TARGET} PUBLIC ${HEADER_DIRECTORIES})
target_link_libraries(${CORE_TARGET} PUBLIC ${FASTOR_TARGET})
target_compile_options(${CORE_TARGET} PRIVATE $<$<CONFIG:Release>:${RELEASE_FLAGS}>)