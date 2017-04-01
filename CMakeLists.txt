
cmake_minimum_required(VERSION 3.0)

project(hot-error)

find_package(CGAL REQUIRED COMPONENTS Core)
include(${CGAL_USE_FILE})

include_directories(include/hot include/polynomial include/Wasserstein)

# get xcode project to show include files
# ??there should be some way to separate out hot and Wasserstein so they don't appear in each others' lists
file(GLOB_RECURSE INCS "*.h")
file(GLOB_RECURSE INCS "*.hpp")
add_library(hot src/hot/hot.cpp ${INCS})
add_library(wasserstein src/wasserstein/WassersteinEdgeEdgeTest.cpp ${INCS})

add_executable(energy src/energy/energy.cpp)
add_executable(tester src/test/test.cpp)
add_executable(wassertest src/wasserstein/wasserstein.cpp)

# get xcode project to show include files
source_group("Include Files" FILES ${INCS})

target_compile_options(hot PUBLIC "-g")
set_property(TARGET hot PROPERTY CXX_STANDARD 11)
set_property(TARGET hot PROPERTY CXX_STANDARD_REQUIRED ON)

set_property(TARGET energy PROPERTY CXX_STANDARD 11)
set_property(TARGET energy PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(energy hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET tester PROPERTY CXX_STANDARD 11)
set_property(TARGET tester PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(tester hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

target_compile_options(wasserstein PUBLIC "-g")
set_property(TARGET wasserstein PROPERTY CXX_STANDARD 11)
set_property(TARGET wasserstein PROPERTY CXX_STANDARD_REQUIRED ON)

target_compile_options(wassertest PUBLIC "-g")
set_property(TARGET wassertest PROPERTY CXX_STANDARD 11)
set_property(TARGET wassertest PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(wassertest wasserstein ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
