project(c2f)
#set(CMAKE_BUILD_TYPE Debug)



set(CMAKE_BUILD_TYPE Release CACHE STRING "Type" FORCE)
set(VTK_DIR "$ENV{VTK_DIR}/lib/cmake/vtk-6.3" CACHE PATH "")

cmake_minimum_required(VERSION 2.8)
#include(${VTK_USE_FILE})
aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(c2f  )


#cmake_minimum_required(VERSION 2.8)
#project(c2f CXX)
#SET (VTK_DIR "/data_space/Program_files/VTK-build")
#SET (CMAKE_BUILD_TYPE "Debug")
#   find_package(VTK 7.1.0 REQUIRED NO_MODULE)
#include(${VTK_USE_FILE})
#add_executable(u_coarse2fine vtu_coarse2fine.cxx)
#target_link_libraries(u_coarse2fine ${VTK_LIBRARIES} )
