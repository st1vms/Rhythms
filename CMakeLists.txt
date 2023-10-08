cmake_minimum_required(VERSION 3.5)

# Definitions
set(PROJECT_NAME rythms)

# Build directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")

# Default flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c17 -fPIC -Wall -Wextra")

# Set build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Add -g option with Debug build
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g")
endif()

project(${PROJECT_NAME})

# Set global include directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Look for submodules and add their subdirectories
file(GLOB SUBMODULE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/*/CMakeLists.txt")
    foreach(SUBMODULE_PATH ${SUBMODULE_DIRS})
    get_filename_component(SUBMODULE_DIR ${SUBMODULE_PATH} DIRECTORY)
    add_subdirectory(${SUBMODULE_DIR})
endforeach()

file(GLOB ROOT_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.c")
foreach(ROOT_SOURCE ${ROOT_SOURCES})
    get_filename_component(ROOT_SOURCE_NAME ${ROOT_SOURCE} NAME_WE)

    # Add test source for compilation
    add_executable(${ROOT_SOURCE_NAME} ${ROOT_SOURCE})

    # Include headers from submodule folders
    file(GLOB SUBMODULE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/*/CMakeLists.txt")
    foreach(SUBMODULE_PATH ${SUBMODULE_DIRS})
        get_filename_component(SUBMODULE_DIR ${SUBMODULE_PATH} DIRECTORY)
        get_filename_component(LIB_NAME ${SUBMODULE_DIR} NAME_WE)

        # Link submodules libraries
        target_link_libraries(${ROOT_SOURCE_NAME} PRIVATE ${LIB_NAME})
        target_include_directories(${ROOT_SOURCE_NAME} PRIVATE "${SUBMODULE_DIR}/include")

    endforeach()

endforeach()
