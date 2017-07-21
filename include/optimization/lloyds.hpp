/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////   LLOYD'S ITERATIONS ////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef _LLOYDS_HPP_
#define _LLOYDS_HPP_


// requires the points to be in order around the perimeter of the polygon. polygon need not be convex 
Point center_mass_polygon( std::vector<Point> points){
	
	double scale_total=0;
	double x_coor=0; 
	double y_coor=0; 
	int size=points.size(); 
	for(int i=0; i <size; i++){
		double scale_i=points[i].x()*points[(i+1)%size].y()-points[(i+1)%size].x()*points[i].y();
		std::cout<<scale_i <<std::endl; 
		scale_total+=scale_i; 
		x_coor+=scale_i*(points[i].x()+points[(i+1)%size].x()); 
		y_coor+=scale_i*(points[i].y()+points[(i+1)%size].y()); 
	}
	x_coor=x_coor/(3*scale_total);
	y_coor=y_coor/(3*scale_total);
	std::cout<<"center of mass: " << x_coor <<", " << y_coor <<std::endl; 
	return Point(x_coor, y_coor);
}

std::vector<Point> lloyds_CVT(std::vector<Point> points, int CVT_iterations, double x_min, double x_max, double y_min, double y_max){
	std::vector<Point> new_sites; 
	while(CVT_iterations >0){
		Voronoi_diagram vd;
		for(auto point:points){
			vd.insert(Wpt(point,0));
		}
		int num_bounded_faces=0;
		std::vector<Point> centroids_Voronoi_regions; 
		for(auto face_itr=vd.bounded_faces_begin(); face_itr!=vd.bounded_faces_end(); face_itr++){
			num_bounded_faces++;
			std::vector<Point> points_in_face;
			
			Voronoi_diagram::Ccb_halfedge_circulator ec_start = (face_itr)->ccb();
			Voronoi_diagram::Ccb_halfedge_circulator ec=ec_start;  
			do{
				points_in_face.push_back(ec ->source() ->point()); 	
			}while(++ec!=ec_start); 
			
			std::cout <<"points in face " << num_bounded_faces <<std::endl;
			for(auto point: points_in_face){std::cout <<point.x() << ", " <<point.y() <<std::endl; 
			}
			
			
		new_sites.push_back(center_mass_polygon(points_in_face)); 
		}	

	
		for(auto unbounded_face_itr=vd.unbounded_faces_begin(); unbounded_face_itr!=vd.unbounded_faces_end(); unbounded_face_itr ++){
			std::vector<Point> points_in_face;
			Voronoi_diagram::Ccb_halfedge_circulator ec_start=(unbounded_face_itr)->ccb(); 
			Voronoi_diagram::Ccb_halfedge_circulator ec=ec_start; 
			do{
				if(ec->is_segment()){
					points_in_face.push_back(ec-> source()->point());
				}
				
			}while(++ec!=ec_start);
		}
	CVT_iterations--;
	std::cout <<"num bouned faces: " << num_bounded_faces <<std::endl; 
	}
	
	
	return new_sites; 
}

RegT build_reg_triangulation(std::vector<Point> points,  int CVT_iterations, double x_min, double x_max, double y_min, double y_max ){

	if(CVT_iterations >0) std::vector<Point> updated_sties=lloyds_CVT(points, CVT_iterations, x_min,  x_max,  y_min,  y_max); 

	RegT rt;

	for(auto point :points){
		rt.insert(Wpt(point, 0));  // intialize with zero weights 
	}
	return rt;
}

#endif // _LLOYD_HPP_
