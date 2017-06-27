
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
using Face_handle=DT::Face_handle;
using Vertex_handle=DT::Vertex_handle;
using Edge_circulator=DT::Edge_circulator; 
using Vertex_iterator=DT::Vertex_iterator;
using Edge=DT::Edge;


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
int sgn(double number);
void compute_h_deriv(const Point &xi, const Point &xj, const Point &xk, int i, double hk_derv[2]); 

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


template<int Wk, int star>
double Edge_Energy(const Triangle &tri1, int index1, const Triangle &tri2, int index2, bool corrected_formulas){
	double h_index1=signed_dist_circumcenters(tri1, index1);
	double h_index2=signed_dist_circumcenters(tri2, index2);
	int sign=sgn(h_index1+h_index2); 
	//double unsigned_energy_subtri1=subtri_energy<Wk,star>(tri1.vertex(index1+1), tri1.vertex(index1+2), hk);
	//double unsigned_energy_subtri2=subtri_energy<Wk,star>(tri2.vertex(index2+1), tri1.vertex(index2+2), hl);
	//double unsigned_energy=unsigned_energy_subtri1+unsigned_energy_subtri2; 
	double unsigned_energy=subtri_energy<Wk, star>(tri1.vertex(index1+1), tri1.vertex(index1+2), h_index1)+subtri_energy<Wk, star>(tri2.vertex(index2+1), tri1.vertex(index2+2), h_index2);
	
	if(corrected_formulas){ 
		return sign*unsigned_energy;
	}
	else{ 
		return unsigned_energy; 
	}	

}

//Unless you know triangulation is DT, energy_density_TMethod should not be used. TMethod=Triangle Method. This means energy is computed without regard for neighboring triangles, which in general, is not correct. 

// The TMethod (triangle method) does not handle boundary as in the paper. In the code below, boundary edge can make negative contribution to energy. 
template<int Wk, int star>
double energy_density_TMethod(const DT &dt){
  	double energy = 0;
  	double mesh_volume=0;

  	for(auto face_itr = dt.finite_faces_begin(); face_itr != dt.finite_faces_end(); face_itr++) {
    		Triangle tri = face_to_tri(*face_itr);
    		energy += tri_energy<Wk,star>(tri);
    		mesh_volume+=tri.area();
 	 }

  //K_real energy_density=energy/mesh_volume;	
  //return energy_density;
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
	//double mesh_volume=0;
	

//  	for(auto face_itr = dt.finite_faces_begin(); face_itr != dt.finite_faces_end(); face_itr++) {
  // 		Triangle tri = face_to_tri(*face_itr);
   //	 mesh_volume+=tri.area();
  	//}

  //K_real energy_density=energy/mesh_volume;	
  //return energy_density;
  return energy;
}



int sgn(double number){
	if(number > 0) return 1;
		else return -1;
}




//////////////////////////////////////////////////////////////////////////////////
/////////////////////  ENERGY DERIVATIVES /////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

template<typename T>
void energy_gradient(T &triangulation,int Wk, int star, Vertex_handle v, double total_deriv[2], bool corrected_formulas){

//template< > 
//void energy_gradient<2,int star>(const DT &dt, Vertex_handle v, double total_deriv[2], bool //corrected_formulas){
	
	total_deriv[0]=0; 
	total_deriv[1]=0; 


	for(auto ei=triangulation.finite_edges_begin();ei!=triangulation.finite_edges_end(); ei++){
	
		Edge edge=*ei;
		Face edge_face=*(edge.first); 
		int edge_index=edge.second; 
		
		Edge mirror_edge=triangulation.mirror_edge(*ei);
		Face mirror_edge_face=*(mirror_edge.first); 
		int mirror_index=mirror_edge.second; 
		
		Vertex_handle  vk_handle=edge_face.vertex(edge_index); 
		Vertex_handle vl_handle=mirror_edge_face.vertex(mirror_index); 
		Vertex_handle vi_handle=edge_face.vertex(edge_face.cw(edge_index)); 
		Vertex_handle vj_handle=edge_face.vertex(edge_face.ccw(edge_index));
		

		Point pi=vi_handle->point();
		Point pj=(*vj_handle).point(); 
		Point pk=(*vk_handle).point(); 
		Point pl=(*vl_handle).point();


		// for now, we assume equal weights, so dij=dji. 
		double dij=0.5*sqrt(pow(pi.x()-pj.x(),2.0)+pow(pi.y()-pj.y(),2.0));
		double dji=0.5*sqrt(pow(pi.x()-pj.x(),2.0)+pow(pi.y()-pj.y(),2.0)); 

		// intialize derivative arrays
		double dij_derv[2]={0,0};
		double dji_derv[2]={0,0}; 
		double hk_derv[2]={0,0} ; 
		double hl_derv[2]={0,0} ;
		if(v==vl_handle){
			hk_derv[0]=0; 
			hk_derv[1]=0;

			compute_h_deriv(pi,pj,pl, 3, hl_derv);
					
			dij_derv[0]=0;
			dij_derv[1]=0; 

			dji_derv[0]=0;
			dji_derv[1]=0;
		} 

		else if(v==vk_handle){
			compute_h_deriv(pi,pj,pk, 3, hk_derv);

			hl_derv[0]=0; 
			hl_derv[1]=0;

			dij_derv[0]=0;
			dij_derv[1]=0; 

			dij_derv[0]=0;
			dij_derv[1]=0;
		}

		else if(v==vi_handle){
			compute_h_deriv(pi,pj,pk, 1, hk_derv);
			compute_h_deriv(pi,pj,pl, 1, hl_derv);

			dij_derv[0]=-(pj.x()-pi.x())/(2.0*(dij+dji)); 
			dij_derv[1]=-(pj.y()-pi.y())/(2.0*(dij+dji));

			dji_derv[0]=dij_derv[0];  // since for now we assume = weights
			dji_derv[1]=dij_derv[1];  // = weights
		}
				
		else if(v==vj_handle){
			compute_h_deriv(pi,pj,pk, 2, hk_derv);
			compute_h_deriv(pi,pj,pl, 2, hl_derv);
			
			dij_derv[0]=(pj.x()-pi.x())/(2.0*(dij+dji));  // assuming equal weights
			dij_derv[1]=(pj.y()-pi.y())/(2.0*(dij+dji));

			dji_derv[0]=dij_derv[0];
			dji_derv[1]=dij_derv[1];		
		}

		else continue; // this means that the edge engery is unchanged by the movement of v 

	
		double hk=signed_dist_circumcenters(face_to_tri(edge_face), edge_index);
		double hl=signed_dist_circumcenters(face_to_tri(mirror_edge_face), mirror_index); 
		int sign=sgn(hk+hl); 

		// I want to make constants that varry based on what Wk and star are. 
			double constant1; 
			double constant2;
		if(star==0){
			constant1=4.0;
			constant2=12.0; 
		}
		else if(star==1){
			constant1=3.0; 
			constant2=3.0;
		}
		else{
			constant1=12.0; 
			constant2=4.0;

		}
		
		bool boundary_edge= triangulation.is_infinite(edge.first) ||triangulation.is_infinite(mirror_edge.first); 
			
			if(!triangulation.is_infinite(edge.first)){
				if(!boundary_edge || hk>0){
					double edge_deriv[2]={0,0}; 
					for(int coor=0; coor <2 ; coor++){
						edge_deriv[coor]+=(3*pow(dij,2)*hk*dij_derv[coor]+pow(dij,3)*hk_derv[coor])/constant1;
						edge_deriv[coor]+=3*pow(dji,2)*dji_derv[coor]*hk+pow(dji,3)*hk_derv[coor]/constant1;
						edge_deriv[coor]+= dij*3*pow(hk,2)*hk_derv[coor]+dij_derv[coor]*pow(hk,3)/constant2;
						edge_deriv[coor]+=dji*3*pow(hk,2)*hk_derv[coor]+dji_derv[coor]*pow(hk,3)/constant2;
							
						if(star==1){
							if(boundary_edge) total_deriv[coor]+=edge_deriv[coor];
							else if(corrected_formulas) total_deriv[coor]+=sign*edge_deriv[coor]; 						
							else total_deriv[coor]+=edge_deriv[coor]; 
						}
						// otherwise star=0 or 2. for now, we make no special treatment for boundary triangles for these stars.
						else total_deriv[coor]+=edge_deriv[coor]; 
					}
				} 
			}

			if(!triangulation.is_infinite(mirror_edge.first)){
				if(!boundary_edge || hl >0){
					double edge_deriv[2]={0,0}; 
					for(int coor=0; coor <2 ; coor++){
						edge_deriv[coor]+=3*pow(dij,2)*hl*dij_derv[coor]+pow(dij,3)*hl_derv[coor]/constant1;
						edge_deriv[coor]+=3*pow(dji,2)*dji_derv[coor]*hl+pow(dji,3)*hl_derv[coor]/constant1;
						edge_deriv[coor]+= dij*3*pow(hl,2)*hl_derv[coor]+dij_derv[coor]*pow(hl,3)/constant2;
						edge_deriv[coor]+= dji*3*pow(hl,2)*hl_derv[coor]+dji_derv[coor]*pow(hl,3)/constant2;
					
						if(star==1){
							if(boundary_edge) total_deriv[coor]+=edge_deriv[coor];
							else if (corrected_formulas) total_deriv[coor]+=sign*edge_deriv[coor]; 
							else total_deriv[coor]+=edge_deriv[coor];
						}
						// otherwise star=0 or 2. for now, we make no special treatment for boundary triangles for these stars.
						else total_deriv[coor]+=edge_deriv[coor]; 
					}
				} 
			}		
			
		} 
	// rescale to be unit vector
	double grad_length=sqrt(pow(total_deriv[0],2)+ pow(total_deriv[1],2)); 
	std::cout<<"unit vector in direction gradient:" << std::setw(15) << total_deriv[0]/grad_length << std::setw(15) <<total_deriv[1]/grad_length <<std::endl; 
	return; 
}


void compute_h_deriv(const Point &xi, const Point &xj, const Point &xk, int i, double h_derv[2]){

	double xi1=xi.x();  
	double xi2=xi.y(); 
	double xj1=xj.x(); 
	double xj2=xj.y(); 
	double xk1=xk.x(); 
	double xk2=xk.y(); 

	if(i==1){ // diff wrt xi coordinates
		
		h_derv[0]= (2*(pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*(xj1 - xk1)*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2) - 
     ((xi1 - xk1)*(xj1 - xk1) + (xi2 - xk2)*(xj2 - xk2))*(2*(pow(xi1 - xj1,2) + pow(xi2 - xj2,2))*(xj2 - xk2)*(-(xj2*xk1) + xi2*(-xj1 + xk1) + xi1*(xj2 - xk2) + xj1*xk2) + 
        2*(xi1 - xj1)*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2)))/ (4.*pow((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2),1.5));


		h_derv[1]= ((xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2))*((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*(xj2 - xk2)*(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2)) -  ((xi1 - xk1)*(xj1 - xk1) + (xi2 - xk2)*(xj2 - xk2))*((pow(xi1 - xj1,2) + pow(xi2 - xj2,2))*(xj1 - xk1) + (xi2 - xj2)*(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2)))))/ (2.*pow((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2),1.5)); 
	}


	if(i==2){ // diff wrt xj coordinates
		h_derv[0]= ((xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2))*((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*(xi1 - xk1)*(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2)) -  ((xi1 - xk1)*(xj1 - xk1) + (xi2 - xk2)*(xj2 - xk2))*((pow(xi1 - xj1,2) + pow(xi2 - xj2,2))*(xi2 - xk2) + (xi1 - xj1)*(-(xj2*xk1) + xi2*(-xj1 + xk1) + xi1*(xj2 - xk2) + xj1*xk2))))/(2.*pow((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2),1.5));


		h_derv[1]= (2*(pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*(xi2 - xk2)*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2) - 
     ((xi1 - xk1)*(xj1 - xk1) + (xi2 - xk2)*(xj2 - xk2))*(2*(pow(xi1 - xj1,2) + pow(xi2 - xj2,2))*(xi1 - xk1)*(-(xj2*xk1) + xi2*(-xj1 + xk1) + xi1*(xj2 - xk2) + xj1*xk2) - 
        2*(xi2 - xj2)*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2)))/(4.*pow((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2),1.5));
	}

	if(i==3){ //diff wrt xk coordinates
		h_derv[0]= (xj2*pow(xk1,2) + pow(xi1,2)*(xj2 - xk2) + pow(xi2,2)*(xj2 - xk2) + pow(xj1,2)*xk2 + pow(xj2,2)*xk2 - 2*xj1*xk1*xk2 - xj2*pow(xk2,2) + 2*xi1*xk1*(-xj2 + xk2) - 
     xi2*(pow(xj1,2) + pow(xj2,2) - 2*xj1*xk1 + pow(xk1,2) - pow(xk2,2)))/(2.*(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2))*sqrt((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2)));


		h_derv[1]= (-(pow(xj1,2)*xk1) - pow(xj2,2)*xk1 + xj1*pow(xk1,2) + pow(xi1,2)*(-xj1 + xk1) + pow(xi2,2)*(-xj1 + xk1) + 2*xi2*(xj1 - xk1)*xk2 + 2*xj2*xk1*xk2 - xj1*pow(xk2,2) + 
     xi1*(pow(xj1,2) + pow(xj2,2) - pow(xk1,2) - 2*xj2*xk2 + pow(xk2,2)))/(2.*(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2))*sqrt((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2)));
	}

}






#endif  // _HOT_HPP_
