// exp6_rectangle.cpp

#include <fstream>

#include "hot.hpp"
#include "Sb.hpp"
#include "ply_writer.hpp"

int main(int argc, char **argv) {
	std::ofstream outputFile; 
	double x_coor_initial=.01;
	double y_coor_initial=.01;
	double x_coor=x_coor_initial; 
	double y_coor=y_coor_initial; 
	double step_size=.01; 
	double rec_width=6; 
	double rec_height=1; 
	

	std::vector<Point> points;
	for(int i=0; i<rec_width+1; i++){
		points.push_back(Point(i,0)); 
	}
	points.push_back(Point(rec_width, rec_height)); 
	points.push_back(Point(0,rec_height)); 
	
	
	std::vector<Segment> boundary_segs; 	
	for(int i=0; i<points.size(); i++){
		boundary_segs.push_back(Segment(points[i], points[(i+1)%points.size()])); 
	}
	//star 0
	outputFile.open("exp6a_rectangle_star_0.txt");
	while(x_coor<rec_width-step_size){
		y_coor=.01; 
		while(y_coor<rec_height-step_size){
			Point freept=Point(x_coor, y_coor); 	
			std::vector<Triangle> triangle_array;
				for(int i=0; i<boundary_segs.size(); i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%boundary_segs.size()], freept));
				}
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=tri_energy<2,0>(triangle_array[i]);
				}
				outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		
			y_coor+=step_size;
		}
		x_coor+=step_size;
	}		

	outputFile.close();


	//star 1
	x_coor=x_coor_initial; 
	y_coor=y_coor_initial; 
	outputFile.open("exp6a_rectangle_star_1.txt");
	while(x_coor<rec_width-step_size){
		y_coor=.01; 
		while(y_coor<rec_height-step_size){
			Point freept=Point(x_coor, y_coor); 
			std::vector<Triangle> triangle_array;
				for(int i=0; i<boundary_segs.size(); i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%boundary_segs.size()], freept));
				}
	
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=tri_energy<2,1>(triangle_array[i]);
				}
				outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		
			y_coor+=step_size;
		}
		x_coor+=step_size;
	}		

	outputFile.close();

	//star 2
	x_coor=x_coor_initial; 
	y_coor=y_coor_initial; 
	outputFile.open("exp6a_rectangle_star_2.txt");
	while(x_coor<rec_width-step_size){
		y_coor=.01; 
		while(y_coor<rec_height-step_size){
			Point freept=Point(x_coor, y_coor); 
			std::vector<Triangle> triangle_array;
				for(int i=0; i<boundary_segs.size(); i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%boundary_segs.size()], freept));
				}
	
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=tri_energy<2,2>(triangle_array[i]);
				}
				outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		
			y_coor+=step_size;
		}
		x_coor+=step_size;
	}		

	outputFile.close();



//////////////////////////////////////////////////
////////////////////// Exp 6b ////////////////////
/////////////////////////////////////////////////

	//star 0
	x_coor=x_coor_initial; 
	y_coor=y_coor_initial;
	outputFile.open("exp6b_rectangle_star_0.txt");
	while(x_coor<rec_width-step_size){
		y_coor=.01; 
		while(y_coor<rec_height-step_size){
			Point freept=Point(x_coor, y_coor); 	
			std::vector<Triangle> triangle_array;
				for(int i=0; i<boundary_segs.size(); i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%boundary_segs.size()], freept));
				}
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=tri_energy_divideArea<2,0>(triangle_array[i],2);
				}
				outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		
			y_coor+=step_size;
		}
		x_coor+=step_size;
	}		

	outputFile.close();


	//star 1
	x_coor=x_coor_initial; 
	y_coor=y_coor_initial; 
	outputFile.open("exp6b_rectangle_star_1.txt");
	while(x_coor<rec_width-step_size){
		y_coor=.01; 
		while(y_coor<rec_height-step_size){
			Point freept=Point(x_coor, y_coor); 
			std::vector<Triangle> triangle_array;
				for(int i=0; i<boundary_segs.size(); i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%boundary_segs.size()], freept));
				}
	
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=tri_energy_divideArea<2,1>(triangle_array[i],2);
				}
				outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		
			y_coor+=step_size;
		}
		x_coor+=step_size;
	}		

	outputFile.close();

	//star 2
	x_coor=x_coor_initial; 
	y_coor=y_coor_initial; 
	outputFile.open("exp6b_rectangle_star_2.txt");
	while(x_coor<rec_width-step_size){
		y_coor=.01; 
		while(y_coor<rec_height-step_size){
			Point freept=Point(x_coor, y_coor); 
			std::vector<Triangle> triangle_array;
				for(int i=0; i<boundary_segs.size(); i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%boundary_segs.size()], freept));
				}
	
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=tri_energy_divideArea<2,2>(triangle_array[i],2);
				}
				outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		
			y_coor+=step_size;
		}
		x_coor+=step_size;
	}		

	outputFile.close();


//////////////////////////////////////////////////
////////////////////// Exp 6c ////////////////////
/////////////////////////////////////////////////

	//star 0
	x_coor=x_coor_initial; 
	y_coor=y_coor_initial;
	outputFile.open("exp6c_rectangle_star_0.txt");
	while(x_coor<rec_width-step_size){
		y_coor=.01; 
		while(y_coor<rec_height-step_size){
			Point freept=Point(x_coor, y_coor); 	
			std::vector<Triangle> triangle_array;
				for(int i=0; i<boundary_segs.size(); i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%boundary_segs.size()], freept));
				}
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=tri_energy_divideArea<2,0>(triangle_array[i],1);
				}
				outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		
			y_coor+=step_size;
		}
		x_coor+=step_size;
	}		

	outputFile.close();


	//star 1
	x_coor=x_coor_initial; 
	y_coor=y_coor_initial; 
	outputFile.open("exp6c_rectangle_star_1.txt");
	while(x_coor<rec_width-step_size){
		y_coor=.01; 
		while(y_coor<rec_height-step_size){
			Point freept=Point(x_coor, y_coor); 
			std::vector<Triangle> triangle_array;
				for(int i=0; i<boundary_segs.size(); i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%boundary_segs.size()], freept));
				}
	
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=tri_energy_divideArea<2,1>(triangle_array[i],1);
				}
				outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		
			y_coor+=step_size;
		}
		x_coor+=step_size;
	}		

	outputFile.close();

	//star 2
	x_coor=x_coor_initial; 
	y_coor=y_coor_initial; 
	outputFile.open("exp6c_rectangle_star_2.txt");
	while(x_coor<rec_width-step_size){
		y_coor=.01; 
		while(y_coor<rec_height-step_size){
			Point freept=Point(x_coor, y_coor); 
			std::vector<Triangle> triangle_array;
				for(int i=0; i<boundary_segs.size(); i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%boundary_segs.size()], freept));
				}
	
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=tri_energy_divideArea<2,2>(triangle_array[i],1);
				}
				outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		
			y_coor+=step_size;
		}
		x_coor+=step_size;
	}		

	outputFile.close();

//////////////////////////////////////////////////
////////////////////// Exp 6d ////////////////////
/////////////////////////////////////////////////
	double zero_weights[]={0,0,0}; 

	x_coor=x_coor_initial; 
	y_coor=y_coor_initial;
	outputFile.open("exp6d_rectangle.txt");
	while(x_coor<rec_width-step_size){
		y_coor=.01; 
		while(y_coor<rec_height-step_size){
			Point freept=Point(x_coor, y_coor); 	
			std::vector<Triangle> triangle_array;
				for(int i=0; i<boundary_segs.size(); i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%boundary_segs.size()], freept));
				}
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=triangle_Sb(triangle_array[i],weighted_circumcenter(triangle_array[i],zero_weights),2,1);
				}
				outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		
			y_coor+=step_size;
		}
		x_coor+=step_size;
	}		

	outputFile.close();

//////////////////////////////////////////////////
////////////////////// Exp 6e ////////////////////
/////////////////////////////////////////////////

	x_coor=x_coor_initial; 
	y_coor=y_coor_initial;
	outputFile.open("exp6e_rectangle.txt");
	while(x_coor<rec_width-step_size){
		y_coor=.01; 
		while(y_coor<rec_height-step_size){
			Point freept=Point(x_coor, y_coor); 	
			std::vector<Triangle> triangle_array;
				for(int i=0; i<boundary_segs.size(); i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%boundary_segs.size()], freept));
				}
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=triangle_Sb(triangle_array[i],weighted_circumcenter(triangle_array[i],zero_weights),2,-1);
				}
				outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		
			y_coor+=step_size;
		}
		x_coor+=step_size;
	}		

	outputFile.close();

//////////////////////////////////////////////////
////////////////////// Exp 6f ////////////////////
/////////////////////////////////////////////////

	x_coor=x_coor_initial; 
	y_coor=y_coor_initial;
	outputFile.open("exp6f_rectangle.txt");
	while(x_coor<rec_width-step_size){
		y_coor=.01; 
		while(y_coor<rec_height-step_size){
			Point freept=Point(x_coor, y_coor); 	
			std::vector<Triangle> triangle_array;
				for(int i=0; i<boundary_segs.size(); i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%boundary_segs.size()], freept));
				}
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=triangle_Sb_divide_perim4(triangle_array[i],weighted_circumcenter(triangle_array[i],zero_weights));
				}
				outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		
			y_coor+=step_size;
		}
		x_coor+=step_size;
	}		

	outputFile.close();

	return 0; 
}
