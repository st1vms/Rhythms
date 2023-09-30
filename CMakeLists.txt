cmake_minimum_required(VERSION 3.5)

# Definitions
set(PROJECT_NAME rythms)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c17 -fPIC -Wall -Wextra")

project(${PROJECT_NAME})

# Module imports
file(GLOB SUBMODULE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/*/CMakeLists.txt")
foreach(SUBMODULE_DIR ${SUBMODULE_DIRS})
    get_filename_component(SUBMODULE_NAME ${SUBMODULE_DIR} DIRECTORY)
    get_filename_component(SUBMODULE_NAME ${SUBMODULE_NAME} NAME)
    add_subdirectory(${SUBMODULE_NAME})
endforeach()

# Check if the tests directory exists
file(EXISTS tests RESULT_VARIABLE TESTS_DIR_EXISTS)
if(TESTS_DIR_EXISTS)
    # Compile each individual test .c file in the root tests/ project
    file(GLOB ROOT_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.c")
    foreach(ROOT_SOURCE ${ROOT_SOURCES})
        get_filename_component(ROOT_SOURCE_NAME ${ROOT_SOURCE} NAME_WE)
        add_executable(${ROOT_SOURCE_NAME} ${ROOT_SOURCE})

        # Link against submodule libraries
        foreach(SUBMODULE_DIR ${SUBMODULE_DIRS})
            get_filename_component(SUBMODULE_NAME ${SUBMODULE_DIR} DIRECTORY)
            get_filename_component(SUBMODULE_NAME ${SUBMODULE_NAME} NAME)
            target_link_libraries(${ROOT_SOURCE_NAME} PRIVATE  ${SUBMODULE_NAME})
        endforeach()
    endforeach()
endif()

## Each submodule should include this CMakeLists.txt file

# get_filename_component(CURRENT_SOURCE_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
# file(GLOB SOURCE_FILES "src/*.c")
# set(SOURCES ${SOURCE_FILES})
# add_library(${CURRENT_SOURCE_DIR_NAME} STATIC ${SOURCES})
# target_include_directories(${CURRENT_SOURCE_DIR_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
# export(TARGETS ${CURRENT_SOURCE_DIR_NAME} FILE ${CMAKE_CURRENT_SOURCE_DIR}/../build/${CURRENT_SOURCE_DIR_NAME}.cmake)
