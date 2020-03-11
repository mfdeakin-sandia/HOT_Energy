# CMakeLists.txt for HOT_Energy

cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(hot-error)

find_package(CGAL REQUIRED COMPONENTS Core)
include(${CGAL_USE_FILE})

include_directories(include/hot include/polynomial include/Wasserstein include/optimization include/cgal-kernel)

# get xcode project to show include files
#   ??there should be some way to separate out hot and Wasserstein so they don't appear in each others' lists
#   yes, do two separate GLOBs of subdirectories to reduce these dependencies
file(GLOB_RECURSE INCS "include/hot/*.h" "include/hot/*.hpp")
file(GLOB_RECURSE W_INCS "include/Wasserstein/*.h" "include/Wasserstein/*.hpp")
file(GLOB_RECURSE O_INCS "include/optimization/*.h" "include/optimization/*.hpp")
file(GLOB_RECURSE P_INCS "include/polynomial/*.h" "include/polynomial/*.hpp")

# don't make these libraries, since the triangulation is templated on the underlying field type
# add_library(hot src/hot/hot.cpp src/hot/ply_writer.cpp ${INCS})
# add_library(wasserstein src/wasserstein/WassersteinEdgeEdgeTest.cpp ${INCS})
#
# add_library(build_triangulation src/build_triangulation/build_triangulation.cpp ${INCS})
#

# main executables
add_executable(energy src/energy/energy.cpp ${INCS} ${O_INCS})
add_executable(tester src/test/test.cpp ${INCS} ${P_INCS})
add_executable(wassertest src/wasserstein/wasserstein.cpp src/wasserstein/WassersteinEdgeEdgeTest.cpp ${INCS}  ${W_INCS})

# experiments
add_executable(exp1_constrained_isoscles src/experiments/exp1_constrained_isoscles.cpp ${INCS}) 
add_executable(exp2_vertex_to_fixed_edge src/experiments/exp2_vertex_to_fixed_edge.cpp ${INCS})
add_executable(exp3_collapse_no_large_angles src/experiments/exp3_collapse_no_large_angles.cpp ${INCS})
add_executable(exp4_hex_patch src/experiments/exp4_hex_patch.cpp ${INCS})
add_executable(exp5_horseV src/experiments/exp5_horseV.cpp ${INCS})
add_executable(exp5_horseV_DT src/experiments/exp5_horseV_DT.cpp ${INCS})
add_executable(exp6_rectangle src/experiments/exp6_rectangle.cpp ${INCS})
add_executable(exp7_vertex_to_fixed_edge_correctedformulas src/experiments/exp7_vertex_to_fixed_edge_correctedformulas.cpp ${INCS})

add_executable(lloydsCVT src/energy/lloydsCVT.cpp ${INCS} ${O_INCS})
add_executable(sandbox src/sandbox/sandbox.cpp ${INCS} ${O_INCS})

# NDT vs DT
add_executable(NDTvDT_exp1 src/experiments/NDTvDT/NDTvDT_exp1.cpp ${INCS})
add_executable(NDTvDT_exp2 src/experiments/NDTvDT/NDTvDT_exp2.cpp ${INCS})
add_executable(NDTvDT_exp3 src/experiments/NDTvDT/NDTvDT_exp3.cpp ${INCS})

add_executable(draw_voronoi src/sandbox/draw_voronoi.cpp ${INCS})

# get xcode project to show include files
source_group("Include Files" FILES ${INCS} ${W_INCS} ${O_INCS} ${P_INCS})

####
# target_compile_options(hot PUBLIC "-g")
# set_property(TARGET hot PROPERTY CXX_STANDARD 11)
# set_property(TARGET hot PROPERTY CXX_STANDARD_REQUIRED ON)
####

####
# target_compile_options(hot PUBLIC "-g")
# set_property(TARGET build_triangulation PROPERTY CXX_STANDARD 11)
# set_property(TARGET build_triangulation PROPERTY CXX_STANDARD_REQUIRED ON)
####

set_property(TARGET energy PROPERTY CXX_STANDARD 11)
set_property(TARGET energy PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(energy ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET tester PROPERTY CXX_STANDARD 11)
set_property(TARGET tester PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(tester ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

#
set_property(TARGET exp1_constrained_isoscles PROPERTY CXX_STANDARD 11)
set_property(TARGET exp1_constrained_isoscles PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp1_constrained_isoscles ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET exp2_vertex_to_fixed_edge PROPERTY CXX_STANDARD 11)
set_property(TARGET exp2_vertex_to_fixed_edge PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp2_vertex_to_fixed_edge ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET exp3_collapse_no_large_angles PROPERTY CXX_STANDARD 11)
set_property(TARGET exp3_collapse_no_large_angles PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp3_collapse_no_large_angles ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET exp4_hex_patch PROPERTY CXX_STANDARD 11)
set_property(TARGET exp4_hex_patch PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp4_hex_patch ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET exp5_horseV PROPERTY CXX_STANDARD 11)
set_property(TARGET exp5_horseV PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp5_horseV ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET exp6_rectangle PROPERTY CXX_STANDARD 11)
set_property(TARGET exp6_rectangle PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp6_rectangle ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET lloydsCVT PROPERTY CXX_STANDARD 11)
set_property(TARGET lloydsCVT PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(lloydsCVT ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET exp7_vertex_to_fixed_edge_correctedformulas PROPERTY CXX_STANDARD 11)
set_property(TARGET exp7_vertex_to_fixed_edge_correctedformulas PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp7_vertex_to_fixed_edge_correctedformulas ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET sandbox PROPERTY CXX_STANDARD 11)
set_property(TARGET sandbox PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(sandbox ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET NDTvDT_exp1 PROPERTY CXX_STANDARD 11)
set_property(TARGET NDTvDT_exp1 PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(NDTvDT_exp1 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET NDTvDT_exp2 PROPERTY CXX_STANDARD 11)
set_property(TARGET NDTvDT_exp2 PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(NDTvDT_exp2 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET NDTvDT_exp3 PROPERTY CXX_STANDARD 11)
set_property(TARGET NDTvDT_exp3 PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(NDTvDT_exp3 ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET exp5_horseV_DT PROPERTY CXX_STANDARD 11)
set_property(TARGET exp5_horseV_DT PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp5_horseV_DT ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET draw_voronoi PROPERTY CXX_STANDARD 11)
set_property(TARGET draw_voronoi PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(draw_voronoi ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})


#



#### 
# target_compile_options(wasserstein PUBLIC "-g")
# set_property(TARGET wasserstein PROPERTY CXX_STANDARD 11)
# set_property(TARGET wasserstein PROPERTY CXX_STANDARD_REQUIRED ON)
#### 

target_compile_options(wassertest PUBLIC "-g")
set_property(TARGET wassertest PROPERTY CXX_STANDARD 11)
set_property(TARGET wassertest PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(wassertest ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

