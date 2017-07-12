#ifndef _ENERGYNOWEIGHTS_HPP_
#define _ENERGYNOWEIGHTS_HPP_

#include <hot.hpp>


#include <CGAL/Cartesian.h>
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

#include <boost/variant/variant.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <list>
#include <vector>
#include <functional> 

#include <polynomial.hpp>

using K = CGAL::Cartesian<double>;
// real type underlying K
using K_real = K::RT;

using EK=CGAL::Exact_predicates_exact_constructions_kernel; // need for point_outside_domain to tell if point is on the boundary
using EK_real=EK::RT; 

//using CK=CGAL::Exact_predicates_inexact_constructions_kernel;

using DT = CGAL::Delaunay_triangulation_2<K>;
using RegT=CGAL::Regular_triangulation_2<K>;


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



int sgn(double number);
Triangle face_to_tri(const Face &face);

template<int Wk, int star>
double Edge_Energy(const Triangle &tri1, int index1, const Triangle &tri2, int index2, bool corrected_formulas);


// This computes h_(vertex_index) 
double signed_dist_circumcenters(const Triangle &tri, int vertex_index){
	Point point0=tri.vertex(vertex_index);
	Point point1_opp=tri.vertex(vertex_index+1);
	Point point2_opp=tri.vertex(vertex_index+2);
	double cot=((point1_opp -point0)*(point2_opp-point0))/sqrt(squared_distance(point1_opp, point0)*squared_distance(point2_opp, point0)-pow((point1_opp -point0)*(point2_opp-point0),2));
	double length_opp_edge= sqrt(squared_distance(point1_opp, point2_opp)); 
	return 0.5*length_opp_edge*cot;	
}

template<const int k, const int star>
double subtri_energy(const Point &xi, const Point &xj, double hk);

template<>
double subtri_energy<2,0>(const Point &xi, const Point &xj, double hk){
	double dij= 0.5*sqrt(pow(xi.x()-xj.x(),2.0)+pow(xi.y()-xj.y(),2.0));
	return pow(dij,3)*hk/2 +dij*pow(hk,3)/6;
}

template<>
double subtri_energy<2,1>(const Point &xi, const Point &xj, double hk){
	double dij= 0.5*sqrt(pow(xi.x()-xj.x(),2.0)+pow(xi.y()-xj.y(),2.0)); // computes distance from xi to midpoint
	return (2.0/3)*(pow(dij,3)*hk+dij*pow(hk,3));

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
			energy+=subtri_energy<Wk,star>(tri.vertex(k+1), tri.vertex(k+2), hk); 
		}
	return energy;
}


//Unless you know triangulation is DT, energy_density_TMethod should not be used. TMethod=Triangle Method. This means energy is computed without regard for neighboring triangles, which in general, is not correct. 

// The TMethod (triangle method) does not handle boundary as in the paper. In the code below, boundary edge can make negative contribution to energy. 
//template<typename T>
template<int Wk, int star>
double energy_density_TMethod(const DT &t){
//double energy_density_TMethod(const T &t, int Wk, int star){
  	double energy = 0;
 
  	for(auto face_itr = t.finite_faces_begin(); face_itr != t.finite_faces_end(); face_itr++) {
    		Triangle tri = face_to_tri(*face_itr);
    		energy += tri_energy<Wk,star>(tri);
 	 }

  return energy;
}





// bool=false means we use the appendix formulas. bool=true means we use our corrected formulas
template<int Wk, int star>
double energy_density_EMethod(const DT &dt, bool corrected_formulas){
  	double energy = 0;
	
	int edgenum=0; 
	for(auto ei=dt.finite_edges_begin();ei!=dt.finite_edges_end(); ei++, edgenum++){

		// check if ei is a boundary edge
		if(dt.is_infinite(ei->first) || dt.is_infinite((ei->first)->neighbor(ei->second))){

			// Two faces can be used to describe each edge. We define finite_tri to be the finite triangle that edge ei bounds. 
			Triangle finite_tri;
			int index_opp_vertex; 

			if(!dt.is_infinite(ei->first)){
				finite_tri=face_to_tri(*(ei->first));
				index_opp_vertex=ei->second;
			}
			else{	Edge mirror_ei=dt.mirror_edge(*ei);
				index_opp_vertex=mirror_ei.second;
				finite_tri=face_to_tri(*(mirror_ei.first));
				
			} 
			double h_index_opp_vertex=signed_dist_circumcenters(finite_tri, index_opp_vertex); 

			//h_index_opp_vertex >0 means that (weighted) circumcenter is inside finite_tri. 
			// if h_index_opp_vertex <=0, then dual edge length 0, no energy contribution

			if(h_index_opp_vertex >0) energy+=subtri_energy<Wk, star>(finite_tri.vertex(index_opp_vertex+1), finite_tri.vertex(index_opp_vertex+2), h_index_opp_vertex); 
			// Author's would stop here. but I think we should add this else clause. Otherwise, when circumcenter is outside boundary triagnle, we don't get a bound for |sigma||*sigma|W(sigma, *sigma). 
			//else energy-=subtri_energy<Wk, star>(finite_tri.vertex(index_opp_vertex+1), finite_tri.vertex(index_opp_vertex+2), h_index_opp_vertex);
			
		}
		
		// if edge is not a boundary edge, we do this else clause
		else{
			Face face1 = *(ei->first); 
    			int index_opp_edge_1 = ei->second;
			
			Edge mirror_ei=dt.mirror_edge(*ei); 
			Face face2=*(mirror_ei.first);
			int index_opp_edge_2=mirror_ei.second;
	
			Triangle tri2=face_to_tri(face2); 
			Triangle tri1=face_to_tri(face1);
			energy+=Edge_Energy<Wk, star>(tri1, index_opp_edge_1, tri2, index_opp_edge_2, corrected_formulas);
		}
		
  	}
	
  return energy;
}

template<int Wk, int star>
double Edge_Energy(const Triangle &tri1, int index1, const Triangle &tri2, int index2, bool corrected_formulas){
	double h_index1=signed_dist_circumcenters(tri1, index1);
	double h_index2=signed_dist_circumcenters(tri2, index2);
	int sign=sgn(h_index1+h_index2); 
	double unsigned_energy=subtri_energy<Wk, star>(tri1.vertex(index1+1), tri1.vertex(index1+2), h_index1)+subtri_energy<Wk, star>(tri2.vertex(index2+1), tri2.vertex(index2+2), h_index2);
	
	if(corrected_formulas){ 
		return sign*unsigned_energy;
	}
	else{ 
		return unsigned_energy; 
	}	

}


int sgn(double number){
	if(number > 0) return 1;
		else return -1;
}

//////////////////////////////////////////////////////////////////////
//////////////////   OTHER ENERGIES /////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

template<int Wk, int star>
double tri_energy_divideArea(const Triangle &tri, int area_pow){
	double tri_area=std::abs(tri.area()); 
	return tri_energy<Wk,star>(tri)/pow(tri_area, area_pow);
}


template<int Wk, int star>
double HOTenergy_divideByTriangleArea(const DT &t, int area_pow){
  	double energy = 0;
 
  	for(auto face_itr = t.finite_faces_begin(); face_itr != t.finite_faces_end(); face_itr++) {
    		Triangle tri = face_to_tri(*face_itr);
    		double triangle_energy= tri_energy<Wk,star>(tri);
		energy+=triangle_energy/pow(std::abs(tri.area()),area_pow); 
 	 }

  return energy;
}
#endif
