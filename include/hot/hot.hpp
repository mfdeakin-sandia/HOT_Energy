
#ifndef _HOT_HPP_
#define _HOT_HPP_

#include <CGAL/Cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/number_utils.h>
#include <CGAL/Triangulation_2.h>

#include <boost/variant/variant.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional> 

#include <polynomial.hpp>

constexpr const int dims = 2;

using K = CGAL::Cartesian<double>;
// real type underlying K
using K_real = K::RT;

using DT = CGAL::Delaunay_triangulation_2<K>;
using Face = DT::Face;
using Point = DT::Point;
using Vertex = DT::Vertex;
using midpoint=CGAL::Point_2<K>;

using Triangle = CGAL::Triangle_2<K>;
using Line = CGAL::Line_2<K>;
using Segment = CGAL::Segment_2<K>;
using Vector = CGAL::Vector_2<K>;
using Direction = CGAL::Direction_2<K>;

constexpr const int tri_verts = 3;
constexpr const int tri_edges = 3;

void order_points(std::array<Point, tri_verts> &verts);
Triangle face_to_tri(const Face &face);
double signed_dist_circumcenters(const Triangle &tri, int vertex_index);
int sign(double number);

double wasserstein2_edge_edge(const Segment &e0,
                              const Segment &e1);

/* Computes 1 or 2 triangles which can be used for the
 * piecewise integral of the Wasserstein distance with
 * k */
boost::variant<std::array<Triangle, 2>, //boost::variant is like a vector where the entries can be different types
               std::array<Triangle, 1> >
integral_bounds(const Triangle &face);

/* Computes the centroid of the triangle */
Point triangle_centroid(const Triangle &face);

template <typename T, int k>
class triangle_w_helper;

template <typename T>
class triangle_w_helper<T, 2> {
 public:
  K_real operator()(T &bounds, const Point &centroid) {
    K_real integral = 0;
    for(auto area : bounds) {
      /* Start by computing the lines used for boundaries
       * and the width.
       */
      std::array<Point, tri_verts> verts = {  // this is making an array of tri_verts points
          area.vertex(0), area.vertex(1), area.vertex(2)};
      order_points(verts);

      /* Find the index of the non-vertical point
       * It's guaranteed to be either the
       * leftmost or
       * rightmost point
       */
      const int w_idx =
          (verts[0][0] != verts[1][0]) ? 0 : 2;
      // Compute the signed width
      K_real width = verts[w_idx][0] -
                     verts[(w_idx + 1) % tri_verts][0];
      /* Compute the bounding lines
       * x = ax t + bx -> t = (x - bx) /
       * ax
       * y = ay t + by = (ay / ax) x +
       * (-ay bx / ax + by)
       *
       * bx = verts[w_idx][0], by = verts[w_idx][1]
       */
      Vector dir_1 =
          Line(verts[w_idx], verts[(w_idx + 1) % tri_verts])
              .to_vector();
      Vector dir_2 =
          Line(verts[w_idx], verts[(w_idx + 2) % tri_verts])
              .to_vector();
      K_real slope_1 = dir_1[1] / dir_1[0];
      K_real slope_2 = dir_2[1] / dir_2[0];
      if((w_idx == 0 && slope_1 < slope_2) ||
         (w_idx == 2 && slope_1 > slope_2)) {
        std::swap(slope_1, slope_2);
      }
      Numerical::Polynomial<K_real, 1, 1> bound_1;
      bound_1.coeff(1) = slope_1;
      bound_1.coeff(0) =
          -verts[w_idx][0] * slope_1 + verts[w_idx][1];
      Numerical::Polynomial<K_real, 1, 1> bound_2;
      bound_2.coeff(1) = slope_2;
      bound_2.coeff(0) =
          -verts[w_idx][0] * slope_2 + verts[w_idx][1];

      Numerical::Polynomial<K_real, 1, 2> initial_x_root(
          (Tags::Zero_Tag()));
      initial_x_root.coeff(1, 0) = 1;
      initial_x_root.coeff(0, 0) = -centroid[0];

      Numerical::Polynomial<K_real, 1, 2> initial_y_root(
          (Tags::Zero_Tag()));
      initial_y_root.coeff(0, 1) = 1;
      initial_y_root.coeff(0, 0) = -centroid[1];

      Numerical::Polynomial<K_real, 2, 2> initial =
          initial_x_root * initial_x_root +
          initial_y_root * initial_y_root;

      auto y_int = initial.integrate(1);
      Numerical::Polynomial<K_real, 3, 1> upper =
          y_int.var_sub(1, bound_1);
      // 3.55556 y + -1.33333 y ^ 2 + 0.333333 y ^ 3 +
      // -2.66667 x y + x ^ 2 y
      // y = y_1 = -x + 3
      // -1.33333 + -3.55556 x + 4.33333 x ^ 2 + -1 x ^ 3
      Numerical::Polynomial<double, 3, 1> lower =
          y_int.var_sub(1, bound_2);
      auto y_bounded = upper + -lower;

      auto x_int = y_bounded.integrate(0);

      auto left = x_int.slice(0, verts[w_idx][0]).coeff(0);
      auto right =
          x_int.slice(0, verts[(w_idx + 1) % tri_verts][0])
              .coeff(0);
      integral += std::abs(left + -right);
    }
    return integral;
  }	 
};

/* Computes the k Wasserstein distance of the triangular
 * face to it's centroid */
template <int k>
K_real triangle_w(const Triangle &face);

template <>
K_real triangle_w<2>(const Triangle &tri) {
  Point centroid = triangle_centroid(tri);
  boost::variant<std::array<Triangle, 2>,
                 std::array<Triangle, 1> >
      bounds = integral_bounds(tri);
  // Set to NaN until implemented
  K_real distance = 0.0 / 0.0;
  if(bounds.which() == 0) {
    // std::array<Triangle, 2>
    auto area =
        boost::get<std::array<Triangle, 2> >(bounds);
    return triangle_w_helper<std::array<Triangle, 2>, 2>()( //triangle_w_helper gets handed 2 triangles. second instance of 2 is referring to W_2
        area, centroid);
  } else {
    // std::array<Triangle, 1>
    auto area =
        boost::get<std::array<Triangle, 1> >(bounds);
    return triangle_w_helper<std::array<Triangle, 1>, 2>()( //triangle_w_helper gets handed 1 triangle. second instance of 2 is referring to W_2
        area, centroid);
  }
}

/* See "HOT: Hodge-Optimized Triangulations" for details on
 * the energy functional */
template <int k>
K_real hot_energy(const DT &dt) {
  K_real energy = 0;
  for(auto face_itr = dt.finite_faces_begin();
      face_itr != dt.finite_faces_end(); face_itr++) {
    Triangle tri = face_to_tri(*face_itr);
	//std::cout<< "This is the coordinates of a face of triangle: " << tri.vertex(0) << ", " << tri.vertex(1) << ", " << tri.vertex(2) << std::endl; 
    K_real area = tri.area();
    K_real wasserstein = triangle_w<k>(tri);
    energy += wasserstein * area;
  }
  return energy;
}




 // This computes h_(vertex_index) 
double signed_dist_circumcenters(const Triangle &tri, int vertex_index){
	Point point0=tri.vertex(vertex_index);
	Point point1_opp=tri.vertex(vertex_index+1);
	Point point2_opp=tri.vertex(vertex_index+2);
	//std::cout << "point0: " << point0.x() << " " << point0.y() <<std::endl;
	//std::cout << "point1: " << point1_opp.x() << " " << point1_opp.y() <<std::endl;
	//std::cout << "point2: " << point2_opp.x() << " " << point2_opp.y() <<std::endl;
	double cot=((point1_opp -point0)*(point2_opp-point0))/sqrt(squared_distance(point1_opp, point0)*squared_distance(point2_opp, point0)-pow((point1_opp -point0)*(point2_opp-point0),2));
	
	//std::cout<< "cot: " << cot << std::endl;
	//double cotdenom=sqrt(squared_distance(point1_opp, point0)*squared_distance(point2_opp, point0)-pow((point1_opp -point0)*(point2_opp-point0),2));
	//double cotnum=((point1_opp -point0)*(point2_opp-point0));
	//std:: cout<< "cotdenom: " << cotdenom << std::endl;
	//std::cout<< "cotnum: "<< cotnum<<std::endl;
	double length_opp_edge= sqrt(squared_distance(point1_opp, point2_opp)); 
	//std::cout << "legnth_opp_edge: " << length_opp_edge << std::endl; 
	return 0.5*length_opp_edge*cot;	
}



template<const int k, const int star>
double subtri_energy(const Point &xi, const Point &xj, double hk);


template<>
double subtri_energy<2,0>(const Point &xi, const Point &xj, double hk){
	double dij= 0.5*sqrt(pow(xi.x()-xj.x(),2.0)+pow(xi.y()-xj.y(),2.0));
	//std::cout<< "dij: " << dij << std::endl; 
	return pow(dij,3)*hk/2 +dij*pow(hk,3)/6;
}

template<>
double subtri_energy<2,1>(const Point &xi, const Point &xj, double hk){
	double dij= 0.5*sqrt(pow(xi.x()-xj.x(),2.0)+pow(xi.y()-xj.y(),2.0)); // computes distance from xi to midpoint
	return (2.0/3)*(pow(dij,3)*hk+dij*pow(hk,3));
	//int sign=sign(hk+hl);
	//return 2.0/3*(sign*(hk+hl)*pow(dij, 3)+sign*dij*(pow(hk,3)+pow(hl, 3)));
}

template<>
double subtri_energy<2,2>(const Point &xi, const Point &xj, double hk){
	double dij= 0.5*sqrt(pow(xi.x()-xj.x(),2.0)+pow(xi.y()-xj.y(),2.0));
	return pow(dij,3)*hk/6 +dij*pow(hk,3)/2;
}

template<int Wk, int star>
double tri_energy(const Triangle &tri){
	double energy=0;
		for(int k=0; k<3; k++){
			double hk=signed_dist_circumcenters(tri,k);
			//std::cout<< "hk: "<< hk <<std::endl; 
			//std::cout<<std::setw(6) << "subtri eng: " << std::setw(10) << subtri_energy<Wk,star>(tri.vertex(k+1), tri.vertex(k+2), hk)<< std::endl;
		energy+=subtri_energy<Wk,star>(tri.vertex(k+1), tri.vertex(k+2), hk); 
		}
	return energy;
}


template<int Wk, int star>
double hot_energy_density(const DT &dt){
  double energy = 0;
  double mesh_volume=0;

  for(auto face_itr = dt.finite_faces_begin(); face_itr != dt.finite_faces_end(); face_itr++) {
    Triangle tri = face_to_tri(*face_itr);
    energy += tri_energy<Wk,star>(tri);
    mesh_volume+=tri.area();
  }
//	All_edge_iterator edge_itr(t); 
	
//	for(auto edge_itr=dt.All_edges_iterator(); edge_itr!=
//	Segment seg=edge_to_seg(*edge_itr);

	
  //K_real energy_density=energy/mesh_volume;	
  //return energy_density;
  return energy;
}


//template<int Wk, int star>
//double edge_energy(Segment seg, Point xk, Point xl);

//template<>
//double edge_energy(Segment seg, Point xk, Point xl)<2,1>{
//	double hk= sign_dist_circumcenters2(seg, xk);
//	double hl=sign_dist_circumcenters2(seg, xl);
	
//	double energy=0
//	return energy; 
//}

int sign(double number){
	if(number > 0) return 1;
		else return -1;
}

#endif  // _HOT_HPP_
