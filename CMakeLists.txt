
cmake_minimum_required(VERSION 3.0)

project(hot-error)

find_package(CGAL REQUIRED COMPONENTS Core)
include(${CGAL_USE_FILE})

include_directories(include/hot include/polynomial include/Wasserstein include/optimization)

# get xcode project to show include files
# ??there should be some way to separate out hot and Wasserstein so they don't appear in each others' lists
file(GLOB_RECURSE INCS "*.h")
file(GLOB_RECURSE INCS "*.hpp")
add_library(hot src/hot/hot.cpp src/hot/ply_writer.cpp ${INCS})
add_library(wasserstein src/wasserstein/WassersteinEdgeEdgeTest.cpp ${INCS})

##
#add_library(build_triangulation src/build_triangulation/build_triangulation.cpp ${INCS})
##

add_executable(energy src/energy/energy.cpp)
add_executable(tester src/test/test.cpp)
add_executable(wassertest src/wasserstein/wasserstein.cpp)

#
add_executable(exp1_constrained_isoscles src/experiments/exp1_constrained_isoscles.cpp) 
add_executable(exp2_vertex_to_fixed_edge src/experiments/exp2_vertex_to_fixed_edge.cpp)
add_executable(exp3_collapse_no_large_angles src/experiments/exp3_collapse_no_large_angles.cpp)
add_executable(exp4_hex_patch src/experiments/exp4_hex_patch.cpp)
add_executable(exp5_horseV src/experiments/exp5_horseV.cpp)
add_executable(exp6_rectangle src/experiments/exp6_rectangle.cpp)
add_executable(lloydsCVT src/energy/lloydsCVT.cpp)
add_executable(exp7_vertex_to_fixed_edge_correctedformulas src/experiments/exp7_vertex_to_fixed_edge_correctedformulas.cpp)
add_executable(sandbox src/sandbox/sandbox.cpp)
add_executable(exp5_horseV_DT src/experiments/exp5_horseV_DT.cpp)

add_executable(NDTvDT_exp1 src/experiments/NDTvDT/NDTvDT_exp1.cpp)
add_executable(NDTvDT_exp2 src/experiments/NDTvDT/NDTvDT_exp2.cpp)
add_executable(NDTvDT_exp3 src/experiments/NDTvDT/NDTvDT_exp3.cpp)

add_executable(draw_voronoi src/sandbox/draw_voronoi.cpp)
#

#add_executable(WassersteinEdgeEdgeTest src/wasserstein/WassersteinEdgeEdgeTest.cpp)

# get xcode project to show include files
source_group("Include Files" FILES ${INCS})

target_compile_options(hot PUBLIC "-g")
set_property(TARGET hot PROPERTY CXX_STANDARD 11)
set_property(TARGET hot PROPERTY CXX_STANDARD_REQUIRED ON)

####
#target_compile_options(hot PUBLIC "-g")
#set_property(TARGET build_triangulation PROPERTY CXX_STANDARD 11)
#set_property(TARGET build_triangulation PROPERTY CXX_STANDARD_REQUIRED ON)

####

set_property(TARGET energy PROPERTY CXX_STANDARD 11)
set_property(TARGET energy PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(energy hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET tester PROPERTY CXX_STANDARD 11)
set_property(TARGET tester PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(tester hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

#
set_property(TARGET exp1_constrained_isoscles PROPERTY CXX_STANDARD 11)
set_property(TARGET exp1_constrained_isoscles PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp1_constrained_isoscles hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET exp2_vertex_to_fixed_edge PROPERTY CXX_STANDARD 11)
set_property(TARGET exp2_vertex_to_fixed_edge PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp2_vertex_to_fixed_edge hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET exp3_collapse_no_large_angles PROPERTY CXX_STANDARD 11)
set_property(TARGET exp3_collapse_no_large_angles PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp3_collapse_no_large_angles hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET exp4_hex_patch PROPERTY CXX_STANDARD 11)
set_property(TARGET exp4_hex_patch PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp4_hex_patch hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET exp5_horseV PROPERTY CXX_STANDARD 11)
set_property(TARGET exp5_horseV PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp5_horseV hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET exp6_rectangle PROPERTY CXX_STANDARD 11)
set_property(TARGET exp6_rectangle PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp6_rectangle hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET lloydsCVT PROPERTY CXX_STANDARD 11)
set_property(TARGET lloydsCVT PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(lloydsCVT hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET exp7_vertex_to_fixed_edge_correctedformulas PROPERTY CXX_STANDARD 11)
set_property(TARGET exp7_vertex_to_fixed_edge_correctedformulas PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp7_vertex_to_fixed_edge_correctedformulas hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET sandbox PROPERTY CXX_STANDARD 11)
set_property(TARGET sandbox PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(sandbox hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET NDTvDT_exp1 PROPERTY CXX_STANDARD 11)
set_property(TARGET NDTvDT_exp1 PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(NDTvDT_exp1 hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET NDTvDT_exp2 PROPERTY CXX_STANDARD 11)
set_property(TARGET NDTvDT_exp2 PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(NDTvDT_exp2 hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET NDTvDT_exp3 PROPERTY CXX_STANDARD 11)
set_property(TARGET NDTvDT_exp3 PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(NDTvDT_exp3 hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET exp5_horseV_DT PROPERTY CXX_STANDARD 11)
set_property(TARGET exp5_horseV_DT PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(exp5_horseV_DT hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})

set_property(TARGET draw_voronoi PROPERTY CXX_STANDARD 11)
set_property(TARGET draw_voronoi PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(draw_voronoi hot ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})


#




target_compile_options(wasserstein PUBLIC "-g")
set_property(TARGET wasserstein PROPERTY CXX_STANDARD 11)
set_property(TARGET wasserstein PROPERTY CXX_STANDARD_REQUIRED ON)

target_compile_options(wassertest PUBLIC "-g")
set_property(TARGET wassertest PROPERTY CXX_STANDARD 11)
set_property(TARGET wassertest PROPERTY CXX_STANDARD_REQUIRED ON)
target_link_libraries(wassertest wasserstein ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})




