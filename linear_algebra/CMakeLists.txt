# Set submodule name
set(PROJECT_NAME linear_algebra)

project(${PROJECT_NAME})

# Find submodule sources
set(SOURCES "")
file(GLOB MODULE_SOURCES "src/*.c")
foreach(MODULE_SOURCE ${MODULE_SOURCES})
    set(SOURCES "${SOURCES} ${MODULE_SOURCE}")
endforeach()
string(STRIP "${SOURCES}" SOURCES)

# Add static library
add_library(${PROJECT_NAME} STATIC ${SOURCES})

# Add submodule include directories
target_include_directories(${PROJECT_NAME} PUBLIC include)