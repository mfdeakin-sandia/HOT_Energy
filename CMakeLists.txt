
cmake_minimum_required(VERSION 3.0)

project(hot-error)

find_package(CGAL REQUIRED COMPONENTS Core)
include(${CGAL_USE_FILE})

include_directories(include/hot include/polynomial)

add_library(hot src/hot/hot.cpp)
add_executable(energy src/energy/energy.cpp)
add_executable(tester src/test/test.cpp)

target_compile_options(hot PUBLIC "-g")

set_property(TARGET hot PROPERTY CXX_STANDARD 11)
set_property(TARGET hot PROPERTY CXX_STANDARD_REQUIRED ON)

set_property(TARGET energy PROPERTY CXX_STANDARD 11)
set_property(TARGET energy PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(energy hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET tester PROPERTY CXX_STANDARD 11)
set_property(TARGET tester PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(tester hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
