cmake_minimum_required(VERSION 3.12)
cmake_policy(SET CMP0074 NEW) # Find packages using *_ROOT
project(INTERVAL_CCD)

################################################################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

### Configuration
set(INTERVAL_CCD_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(Warnings)
include(UseColors)
include(${PROJECT_NAME}Utils)

#TODO set to off once the main is removed

OPTION(INTERVAL_CCD_WITH_TIMER "Enable profiling timers, for debug"          OFF)

include(${PROJECT_NAME}Dependencies)


# inclusion CCD
add_subdirectory(intervalccd)
prepend_current_path(INTERVAL_SOURCES)
intervalCCD_copy_headers(${INTERVAL_SOURCES})
CCD_set_source_group(${INTERVAL_SOURCES})
add_library(${PROJECT_NAME} ${INTERVAL_SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR}/include)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)
target_link_libraries(${PROJECT_NAME} PRIVATE warnings::all)
find_package(Boost REQUIRED)
target_link_libraries(${PROJECT_NAME} PUBLIC Boost::boost)


#Optional


if(INTERVAL_CCD_WITH_TIMER)
  target_compile_definitions(${PROJECT_NAME} PRIVATE INTERVAL_CCD_USE_TIMER)
endif()

target_compile_definitions(${PROJECT_NAME} PUBLIC NOMINMAX)

