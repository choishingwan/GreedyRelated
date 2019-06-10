cmake_minimum_required (VERSION 3.1)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(PROJECT_NAME GreedyRelated)
project(${PROJECT_NAME})
set(CMAKE_CXX_FLAGS "-g -Wall")

option(march "Use --march." OFF)
if(march)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()
set(CMAKE_CXX_STANDARD 11)
# Require (at least) it
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Don't use e.g. GNU extension (like -std=gnu++11) for portability
set(CMAKE_CXX_EXTENSIONS OFF)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
add_library(misc_lib
    ${CMAKE_SOURCE_DIR}/misc.cpp)
add_executable(GreedyRelated ${CMAKE_SOURCE_DIR}/main.cpp)
find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( GreedyRelated PRIVATE ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )
################################
#          Add pthread
################################
find_package (Threads)
target_link_libraries (GreedyRelated PRIVATE ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(GreedyRelated PRIVATE misc_lib)
