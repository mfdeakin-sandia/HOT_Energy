
cmake_minimum_required(VERSION 3.0)

project(hot-error)

find_package(CGAL REQUIRED COMPONENTS Core)
include(${CGAL_USE_FILE})

include_directories(include/polynomial)

add_executable(hot src/hot.cpp)
target_compile_options(hot PUBLIC "-g")
set_property(TARGET hot PROPERTY CXX_STANDARD 11)
set_property(TARGET hot PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
