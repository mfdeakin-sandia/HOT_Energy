// cgal-kernel.h
// here we define the cgal kernel field once, so that all routines are using the same one and we don't have type conflicts

#ifndef CGAL_KERNEL_HPP
#define CGAL_KERNEL_HPP


// Some choices for kernel K
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian.h>

// caller should define one of these, default is first
// define CGAL_EI
// define CGAL_EE
// define CGAL_CART


#ifdef CGAL_EI
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
#else
#ifdef CGAL_EE
using K = CGAL::Exact_predicates_exact_constructions_kernel;
#else
using K = CGAL::Cartesian<double>;
#endif
#endif

using K_real = K::RT;

using EK=CGAL::Exact_predicates_exact_constructions_kernel; 
using EK_real=EK::RT; 

using CK=CGAL::Exact_predicates_inexact_constructions_kernel;

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/number_utils.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_traits_2.h>
#include <CGAL/Delaunay_triangulation_adaptation_policies_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Arrangement_2.h> // for dealing with half edges 

#include <CGAL/Kernel/global_functions.h>

#include <boost/variant/variant.hpp>

// Newer versions of CGAL define Weighted_point inside of all kernels
// but before that we have to use RegularTriangulationTraits_2 implementations
// Note I'm not certain which version this starts at;
// but it's between 4.9.0 and 4.11.0, so adjust this check as needed
#if CGAL_VERSION_NR > CGAL_VERSION_NUMBER(4, 9, 0)

using DT = CGAL::Delaunay_triangulation_2<K>;
using RegT=CGAL:: Regular_triangulation_2<K>;

#else

#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
using TriTraits = CGAL::Regular_triangulation_euclidean_traits_2<K>;
using DT = CGAL::Delaunay_triangulation_2<TriTraits>;
using RegT=CGAL:: Regular_triangulation_2<TriTraits>;

#endif // CGAL_VERSION_NR

using Regular_triangulation_2 = RegT;

using Face = DT::Face;
using Point = DT::Point;
using Vertex = DT::Vertex;
using midpoint=CGAL::Point_2<K>;
using Face_handle=DT::Face_handle;
using Vertex_handle=DT::Vertex_handle;
using Edge_circulator=DT::Edge_circulator; 
using Vertex_iterator=DT::Vertex_iterator;
using Edge=DT::Edge;
using Wpt=RegT::Weighted_point;


using AT=CGAL::Delaunay_triangulation_adaptation_traits_2<DT>;
using AP=CGAL::Delaunay_triangulation_caching_degeneracy_removal_policy_2<DT>;
using Voronoi_diagram=CGAL::Voronoi_diagram_2<DT,AT,AP>;  

//for weighted triangulations

using weighted_Face_handle=RegT::Face_handle; 


using Triangle = CGAL::Triangle_2<K>;
using Line = CGAL::Line_2<K>;
using Segment = CGAL::Segment_2<K>;
using Vector = CGAL::Vector_2<K>;
using Direction = CGAL::Direction_2<K>;

using Weightpt = RegT::Weighted_point;

using Point_2 = K::Point_2;
using Iso_rectangle_2 = K::Iso_rectangle_2;
using Segment_2 = K::Segment_2;
using Ray_2 = K::Ray_2;
using Line_2 = K::Line_2;

#endif
