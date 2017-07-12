#ifndef _SB_HPP_
#define _SB_HPP_

#include <cmath>

double perimeter(const Triangle & tri);
Point weighted_circumcenter(const Triangle &tri, double weight[3]); 





//Sb in paper takes powdiff=2 and powarea=1
double triangle_Sb(const Triangle &tri, const Point &wcirc, const double powdist, const double powarea){
	Point bary=centroid(tri); 
	double tri_area=std::abs(tri.area());  
	double dist=sqrt(pow(wcirc.x()-bary.x(),2)+pow(wcirc.y()-bary.y(),2));
	return pow(tri_area, powarea)*pow(dist,powdist); 	
}


double Sb(const RegT &rt,const double powdiff, const double powarea){
  	double energy = 0;
 
  	for(auto face_itr = rt.finite_faces_begin(); face_itr != rt.finite_faces_end(); face_itr++) {
    		Triangle tri = Triangle( face_itr->vertex(0)->point(), face_itr->vertex(1) -> point(), face_itr-> vertex(2) ->point());
		Point wcirc=rt.weighted_circumcenter(face_itr); 
    		double triangle_energy= triangle_Sb(tri,wcirc, powdiff,powarea);
		energy+=triangle_energy; 
 	 }

  return energy;
}


double triangle_Sb_divide_perim4(const Triangle &tri, const Point &wcirc){
	Point bary=centroid(tri); 
	double tri_area=std::abs(tri.area()); 
	double perim =perimeter(tri); 
	double squared_dist=pow(wcirc.x()-bary.x(),2)+pow(wcirc.y()-bary.y(),2);
	return tri_area*squared_dist/pow(perim,4); 	
}

double Sb_divide_perim(const RegT &rt){

	double energy=0; 
	
	for(auto face_itr = rt.finite_faces_begin(); face_itr != rt.finite_faces_end(); face_itr++) {
    		Triangle tri = Triangle( face_itr->vertex(0)->point(), face_itr->vertex(1) -> point(), face_itr-> vertex(2) ->point());
		Point wcirc=rt.weighted_circumcenter(face_itr); 
		energy+=triangle_Sb_divide_perim4(tri,wcirc); 
 	 }

  return energy;
}
	





double perimeter(const Triangle & tri){
	double triangle_perimeter=0;
		triangle_perimeter+= sqrt(CGAL::squared_distance(tri.vertex(0), tri.vertex(1)));
		triangle_perimeter+= sqrt(CGAL::squared_distance(tri.vertex(1), tri.vertex(2)));  
		triangle_perimeter+= sqrt(CGAL::squared_distance(tri.vertex(2), tri.vertex(0)));  
	return triangle_perimeter;
}

Point weighted_circumcenter(const Triangle &tri, double weight[3]){
	double triangle_area=std::abs(tri.area()); 
	
	Point x0_pt=tri.vertex(0); 
	Point x1_pt=tri.vertex(1); 
	Point x2_pt=tri.vertex(2); 
	
	
	double x0[]={(tri.vertex(0)).x(), (tri.vertex(0)).y()}; 
	double x1[]={(tri.vertex(1)).x(), (tri.vertex(1)).y()}; 
	double x2[]={(tri.vertex(2)).x(), (tri.vertex(2)).y()}; 
	double e01[]= {x1[0]-x0[0], x1[1]-x0[1]}; 
	double e02[]= {x2[0]-x0[0], x2[1]-x0[1]};
	
	double e01_perp[2];
	if(CGAL::orientation(x0_pt, x1_pt, x2_pt)==CGAL::LEFT_TURN){
	 	e01_perp[0]=-e01[1];
		e01_perp[1]= e01[0];
	}
	else{
		e01_perp[0]=e01[1];
		e01_perp[1]=-e01[0];
	}


	double e01_perp_length=sqrt(pow(e01_perp[0],2)+pow(e01_perp[1],2));

	double e01_perp_unit[2];
	for(int i=0; i<2; i++){
		e01_perp_unit[i]=e01_perp[i]/e01_perp_length; 
	}

	double e02_perp[2];
		if( CGAL::orientation(x0_pt,x2_pt,x1_pt)==CGAL::LEFT_TURN){
	 		e02_perp[0]=-e02[1];
			e02_perp[1]= e02[0];
		}
		else{
			e02_perp[0]=e02[1];
			e02_perp[1]=-e02[0];
		}
	

	
	double e02_perp_length=sqrt(pow(e02_perp[0],2)+pow(e02_perp[1],2));
	double e02_perp_unit[2];
	for(int i=0; i<2; i++){
		e02_perp_unit[i]=e02_perp[i]/e02_perp_length; 
	}
	
	double scale1=(pow(x0[0]-x1[0],2)+pow(x0[1]-x1[1],2)+weight[0]-weight[1])/(4*triangle_area);
	double scale2=(pow(x0[0]-x2[0],2)+pow(x0[1]-x2[1],2)+weight[0]-weight[2])/(4*triangle_area);

	double wcirc[]={x0[0],x0[1]};
	
	for(int i=0; i<2; i++){
		wcirc[i]+=scale1*e02_perp[i];
		wcirc[i]+=scale2*e01_perp[i]; 
	}
	
	return Point(wcirc[0],wcirc[1]); 
}
#endif
