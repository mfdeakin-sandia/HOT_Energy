
cmake_minimum_required(VERSION 3.0)

project(hot-error)

find_package(CGAL REQUIRED COMPONENTS Core)
include(${CGAL_USE_FILE})

include_directories(include/hot include/polynomial)

# get xcode project to show include files
file(GLOB_RECURSE INCS "*.h")
file(GLOB_RECURSE INCS "*.hpp")
add_library(hot src/hot/hot.cpp src/hot/WassersteinEdgeEdgeTest.cpp ${INCS})

add_executable(energy src/energy/energy.cpp)

# get xcode project to show include files
source_group("Include Files" FILES ${INCS})

target_compile_options(hot PUBLIC "-g")

set_property(TARGET hot PROPERTY CXX_STANDARD 11)
set_property(TARGET hot PROPERTY CXX_STANDARD_REQUIRED ON)

set_property(TARGET energy PROPERTY CXX_STANDARD 11)
set_property(TARGET energy PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(energy hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

