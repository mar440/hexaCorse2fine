project(c2f)
#set(CMAKE_BUILD_TYPE Debug)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
set(CMAKE_BUILD_TYPE Release CACHE STRING "Type" FORCE)
set(VTK_DIR  "$ENV{VTK_ROOT}/build" CACHE PATH "")


cmake_minimum_required(VERSION 2.8)


find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(c2f  ${VTK_LIBRARIES})


