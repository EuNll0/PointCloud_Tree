cmake_minimum_required(VERSION 3.0)
#  Set(CMAKE_C_COMPILER "/usr/bin/clang")
# SET(CMAKE_CXX_COMPILER "/usr/bin/clang++")
project(Aggrate)
include_directories(${PROJECT_SOURCE_DIR}/include/)
set( CMAKE_CXX_STANDARD 17)
# set (CMAKE_BUILD_TYPE "Debug")
set (CMAKE_BUILD_TYPE "Release")
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set( LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

find_package(PCL REQUIRED)    
include_directories(${PCL_INCLUDE_DIRS})               
link_directories(${PCL_LIBRARY_DIRS})                  
add_definitions(${PCL_DEFINITIONS})        

find_package(glog)

set(third_library  ${PCL_LIBRARIES} glog::glog)

#add_subdirectory(src)
add_subdirectory(test)
