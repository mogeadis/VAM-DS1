# Source Files
set(SOURCES_DIR ${PLUGIN_DIR}/src)
file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/${SOURCES_DIR}/*.cpp)
message(STATUS "Collected source files:")
foreach(FILE IN LISTS SOURCES)
    message(STATUS "  ${FILE}")
endforeach()

# Header Files
set(HEADERS_DIR ${PLUGIN_DIR}/include)
file(GLOB_RECURSE HEADERS ${CMAKE_SOURCE_DIR}/${HEADERS_DIR}/*.h)
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

# JUCE Library
set(JUCE_DIR JUCE)
set(JUCE_TARGET JUCE)
set(JUCE_PATH C:/Users/mogeadis/Projects/Cpp/LIBRARIES/${JUCE_DIR}) # LOCAL PATH
add_subdirectory(${JUCE_PATH} ${CMAKE_BINARY_DIR}/${JUCE_TARGET})

# JUCE Plugin
set(JUCE_PLUGIN_TARGET ${PROJECT_NAME})
juce_add_plugin(${JUCE_PLUGIN_TARGET}
    PRODUCT_NAME "VAM DS-1"
    FORMATS Standalone VST3
)
juce_generate_juce_header(${JUCE_PLUGIN_TARGET})
target_include_directories(${JUCE_PLUGIN_TARGET} PRIVATE ${HEADER_DIRECTORIES})
target_sources(${JUCE_PLUGIN_TARGET} PRIVATE ${SOURCES})
target_link_libraries(${JUCE_PLUGIN_TARGET}
    PRIVATE
        ${CORE_TARGET}
        juce::juce_audio_utils
    PUBLIC 
        juce::juce_recommended_config_flags
        juce::juce_recommended_lto_flags
        juce::juce_recommended_warning_flags
)
target_compile_definitions(${JUCE_PLUGIN_TARGET}
    PUBLIC
        JUCE_WEB_BROWSER=0
        JUCE_USE_CURL=0
        JUCE_VST3_CAN_REPLACE_VST2=0
)
target_compile_options(${JUCE_PLUGIN_TARGET} PRIVATE $<$<CONFIG:Release>:${RELEASE_FLAGS}>)