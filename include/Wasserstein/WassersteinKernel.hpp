// WassersteinKernel.hpp

#ifndef WASSERSTEINKERNEL_HPP
#define WASSERSTEINKERNEL_HPP

// define basic geometric types, and which CGAL kernel to use

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/number_utils.h>
#include <CGAL/Root_of_traits.h>

#include <boost/variant/variant.hpp>

#include <algorithm>
#include <cmath>
#include <string>

using Real = float;

// exact predicates, but makes debugging and printing floating point in xcode hard
// using K = CGAL::Exact_predicates_exact_constructions_kernel;

// floating point predicates in the plane
using K = CGAL::Cartesian<Real>;

// real type underlying K
using K_real = K::RT;

using DT = CGAL::Delaunay_triangulation_2<K>;
using Point = DT::Point;
using Segment = CGAL::Segment_2<K>;

#endif // WASSERSTEINKERNEL_HPP
