#include <hot.hpp>
#include <Sb.hpp>
#include <ply_writer.hpp>

#include <random>
#include <vector>
#include <fstream>

#include <CGAL/Kernel/global_functions.h>

#include <typeinfo>

#define PI 3.14159265


void hex_data_to_file(double x_min, double x_max, double y_min, double y_max, double step, int viewnum);

// ****USING INCORRECT FORMULAS **********
int main(int argc, char **argv) {

/*
	hex_data_to_file(-1,1,-1,1,.01, 1); 
	
	hex_data_to_file(-1,-.9, -.1, .1, .01, 2);
	 
	hex_data_to_file(-.5,.5,-.5,.5, .01,3);  

	hex_data_to_file(-.5,.5,.8,1, .01, 5); 

	//view at 60 degree corner
	hex_data_to_file(.4,.6, sqrt(3)/2-.1, sqrt(3)/2, .01,6) ; 
*/

//////////////////////////////////////////////
////// Exp 4a. *1 HOT /////////////////////////
//////////////////////////////////////////////////


	std::vector<Point> points={Point(-1,0),  Point(cos(120*PI/180), sin(120*PI/180)),  Point(cos(60*PI/180), sin(60*PI/180)),  Point(1,0),  Point(cos(300*PI/180), sin(300*PI/180)), 				 Point(cos(240*PI/180), sin(240*PI/180))} ; 
	std::vector<Segment> boundary_segs; 
	for(int i=0; i<6; i++){
		boundary_segs.push_back(Segment(points[i], points[(i+1)%6])); 
	}

	std::ofstream outputFile;

	//star 0
	outputFile.open("exp4/exp4a_hex_patch_star_2_0.txt");
	double x_min= -1;
	double y_min=-1; 
	double x_max=1; 
	double y_max=1; 
	double step=.01; 

	double x_coor=x_min; 
	double y_coor=y_min; 
	while(x_coor<x_max){	
		y_coor=y_min; 
		while(y_coor<y_max){		
		Point freept(x_coor, y_coor); 
			int num_intersections=0; 

			//check if freept is contained in the hexagon. if not, don't compute energy for that point. move on. 
			for(int i=0; i<6; i++){		 
				bool intersect=intersection(Segment(Point(2,0), freept), boundary_segs[i]);
				if(intersect)
					num_intersections+=1; 			
			}
			if(!(num_intersections==1)){
				y_coor+=step; 
				 continue; 
			}
			std::vector<Triangle> triangle_array;
			for(int i=0; i<6; i++){
				triangle_array.push_back(Triangle(points[i], points[(i+1)%6], freept));
			}
	
			double total_hex_energy=0;
			for(int i=0; i<triangle_array.size(); i++){
				total_hex_energy+=tri_energy<2,0>(triangle_array[i]);
			}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_hex_energy <<std::endl; 
			
		y_coor+=step; 
		}
		x_coor+=step;
	}
	outputFile.close();

	//star 1
	outputFile.open("exp4/exp4a_hex_patch_star_2_1.txt");

	x_coor=x_min; 
	y_coor=y_min; 
	while(x_coor<x_max){	
		y_coor=y_min; 
		while(y_coor<y_max){		
		Point freept(x_coor, y_coor); 
			int num_intersections=0; 

			//check if freept is contained in the hexagon. if not, don't compute energy for that point. move on. 
			for(int i=0; i<6; i++){		 
				bool intersect=intersection(Segment(Point(2,0), freept), boundary_segs[i]);
				if(intersect)
					num_intersections+=1; 			
			}
			if(!(num_intersections==1)){
				y_coor+=step; 
				 continue; 
			}
			std::vector<Triangle> triangle_array;
			for(int i=0; i<6; i++){
				triangle_array.push_back(Triangle(points[i], points[(i+1)%6], freept));
			}
	
			double total_hex_energy=0;
			for(int i=0; i<triangle_array.size(); i++){
				total_hex_energy+=tri_energy<2,1>(triangle_array[i]);
			}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_hex_energy <<std::endl; 
			
		y_coor+=step; 
		}
		x_coor+=step;
	}
	outputFile.close();

//star 2
	outputFile.open("exp4/exp4a_hex_patch_star_2_2.txt");

	x_coor=x_min; 
	y_coor=y_min; 
	while(x_coor<x_max){	
		y_coor=y_min; 
		while(y_coor<y_max){		
		Point freept(x_coor, y_coor); 
			int num_intersections=0; 

			//check if freept is contained in the hexagon. if not, don't compute energy for that point. move on. 
			for(int i=0; i<6; i++){		 
				bool intersect=intersection(Segment(Point(2,0), freept), boundary_segs[i]);
				if(intersect)
					num_intersections+=1; 			
			}
			if(!(num_intersections==1)){
				y_coor+=step; 
				 continue; 
			}
			std::vector<Triangle> triangle_array;
			for(int i=0; i<6; i++){
				triangle_array.push_back(Triangle(points[i], points[(i+1)%6], freept));
			}
	
			double total_hex_energy=0;
			for(int i=0; i<triangle_array.size(); i++){
				total_hex_energy+=tri_energy<2,2>(triangle_array[i]);
			}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_hex_energy <<std::endl; 
			
		y_coor+=step; 
		}
		x_coor+=step;
	}
	outputFile.close();

//////////////////////////////////////////////
////// Exp 4b. *1 HOT/area^2 /////////////////////////
//////////////////////////////////////////////////

////star 0
	outputFile.open("exp4/exp4b_hex_patch_star_2_0.txt");
	
	x_coor=x_min; 
	y_coor=y_min; 
	while(x_coor<x_max){	
		y_coor=y_min; 
		while(y_coor<y_max){		
		Point freept(x_coor, y_coor); 
			int num_intersections=0; 

			//check if freept is contained in the hexagon. if not, don't compute energy for that point. move on. 
			for(int i=0; i<6; i++){		 
				bool intersect=intersection(Segment(Point(2,0), freept), boundary_segs[i]);
				if(intersect)
					num_intersections+=1; 			
			}
			if(!(num_intersections==1)){
				y_coor+=step; 
				 continue; 
			}
			std::vector<Triangle> triangle_array;
			for(int i=0; i<6; i++){
				triangle_array.push_back(Triangle(points[i], points[(i+1)%6], freept));
			}
	
			double total_hex_energy=0;
			for(int i=0; i<triangle_array.size(); i++){
				total_hex_energy+=tri_energy_divideArea<2,0>(triangle_array[i],2);
			}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_hex_energy <<std::endl; 
			
		y_coor+=step; 
		}
		x_coor+=step;
	}
	outputFile.close();
///star 1
	outputFile.open("exp4/exp4b_hex_patch_star_2_1.txt");
	
	x_coor=x_min; 
	y_coor=y_min; 
	while(x_coor<x_max){	
		y_coor=y_min; 
		while(y_coor<y_max){		
		Point freept(x_coor, y_coor); 
			int num_intersections=0; 

			//check if freept is contained in the hexagon. if not, don't compute energy for that point. move on. 
			for(int i=0; i<6; i++){		 
				bool intersect=intersection(Segment(Point(2,0), freept), boundary_segs[i]);
				if(intersect)
					num_intersections+=1; 			
			}
			if(!(num_intersections==1)){
				y_coor+=step; 
				 continue; 
			}
			std::vector<Triangle> triangle_array;
			for(int i=0; i<6; i++){
				triangle_array.push_back(Triangle(points[i], points[(i+1)%6], freept));
			}
	
			double total_hex_energy=0;
			for(int i=0; i<triangle_array.size(); i++){
				total_hex_energy+=tri_energy_divideArea<2,1>(triangle_array[i],2);
			}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_hex_energy <<std::endl; 
			
		y_coor+=step; 
		}
		x_coor+=step;
	}
	outputFile.close();
//star 2
	outputFile.open("exp4/exp4b_hex_patch_star_2_2.txt");
	
	x_coor=x_min; 
	y_coor=y_min; 
	while(x_coor<x_max){	
		y_coor=y_min; 
		while(y_coor<y_max){		
		Point freept(x_coor, y_coor); 
			int num_intersections=0; 

			//check if freept is contained in the hexagon. if not, don't compute energy for that point. move on. 
			for(int i=0; i<6; i++){		 
				bool intersect=intersection(Segment(Point(2,0), freept), boundary_segs[i]);
				if(intersect)
					num_intersections+=1; 			
			}
			if(!(num_intersections==1)){
				y_coor+=step; 
				 continue; 
			}
			std::vector<Triangle> triangle_array;
			for(int i=0; i<6; i++){
				triangle_array.push_back(Triangle(points[i], points[(i+1)%6], freept));
			}
	
			double total_hex_energy=0;
			for(int i=0; i<triangle_array.size(); i++){
				total_hex_energy+=tri_energy_divideArea<2,2>(triangle_array[i],2);
			}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_hex_energy <<std::endl; 
			
		y_coor+=step; 
		}
		x_coor+=step;
	}
	outputFile.close();


//////////////////////////////////////////////
////// Exp 4c. *1 HOT/area /////////////////////////
//////////////////////////////////////////////////

//star 0
	outputFile.open("exp4/exp4c_hex_patch_star_2_0.txt");
	
	x_coor=x_min; 
	y_coor=y_min; 
	while(x_coor<x_max){	
		y_coor=y_min; 
		while(y_coor<y_max){		
		Point freept(x_coor, y_coor); 
			int num_intersections=0; 

			//check if freept is contained in the hexagon. if not, don't compute energy for that point. move on. 
			for(int i=0; i<6; i++){		 
				bool intersect=intersection(Segment(Point(2,0), freept), boundary_segs[i]);
				if(intersect)
					num_intersections+=1; 			
			}
			if(!(num_intersections==1)){
				y_coor+=step; 
				 continue; 
			}
			std::vector<Triangle> triangle_array;
			for(int i=0; i<6; i++){
				triangle_array.push_back(Triangle(points[i], points[(i+1)%6], freept));
			}
	
			double total_hex_energy=0;
			for(int i=0; i<triangle_array.size(); i++){
				total_hex_energy+=tri_energy_divideArea<2,0>(triangle_array[i],1);
			}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_hex_energy <<std::endl; 
			
		y_coor+=step; 
		}
		x_coor+=step;
	}
	outputFile.close();

// star 1

	outputFile.open("exp4/exp4c_hex_patch_star_2_1.txt");
	
	x_coor=x_min; 
	y_coor=y_min; 
	while(x_coor<x_max){	
		y_coor=y_min; 
		while(y_coor<y_max){		
		Point freept(x_coor, y_coor); 
			int num_intersections=0; 

			//check if freept is contained in the hexagon. if not, don't compute energy for that point. move on. 
			for(int i=0; i<6; i++){		 
				bool intersect=intersection(Segment(Point(2,0), freept), boundary_segs[i]);
				if(intersect)
					num_intersections+=1; 			
			}
			if(!(num_intersections==1)){
				y_coor+=step; 
				 continue; 
			}
			std::vector<Triangle> triangle_array;
			for(int i=0; i<6; i++){
				triangle_array.push_back(Triangle(points[i], points[(i+1)%6], freept));
			}
	
			double total_hex_energy=0;
			for(int i=0; i<triangle_array.size(); i++){
				total_hex_energy+=tri_energy_divideArea<2,1>(triangle_array[i],1);
			}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_hex_energy <<std::endl; 
			
		y_coor+=step; 
		}
		x_coor+=step;
	}
	outputFile.close();

//star 2

	outputFile.open("exp4/exp4c_hex_patch_star_2_2.txt");
	
	x_coor=x_min; 
	y_coor=y_min; 
	while(x_coor<x_max){	
		y_coor=y_min; 
		while(y_coor<y_max){		
		Point freept(x_coor, y_coor); 
			int num_intersections=0; 

			//check if freept is contained in the hexagon. if not, don't compute energy for that point. move on. 
			for(int i=0; i<6; i++){		 
				bool intersect=intersection(Segment(Point(2,0), freept), boundary_segs[i]);
				if(intersect)
					num_intersections+=1; 			
			}
			if(!(num_intersections==1)){
				y_coor+=step; 
				 continue; 
			}
			std::vector<Triangle> triangle_array;
			for(int i=0; i<6; i++){
				triangle_array.push_back(Triangle(points[i], points[(i+1)%6], freept));
			}
	
			double total_hex_energy=0;
			for(int i=0; i<triangle_array.size(); i++){
				total_hex_energy+=tri_energy_divideArea<2,2>(triangle_array[i],1);
			}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_hex_energy <<std::endl; 
			
		y_coor+=step; 
		}
		x_coor+=step;
	}
	outputFile.close();



//////////////////////////////////////////////
////// Exp 4d. Sb /////////////////////////
//////////////////////////////////////////////////

	outputFile.open("exp4/exp4d_hex_patch.txt");
	
	x_coor=x_min; 
	y_coor=y_min;
	double zero_weights[]={0,0,0};  
	while(x_coor<x_max){	
		y_coor=y_min; 
		while(y_coor<y_max){		
		Point freept(x_coor, y_coor); 
			int num_intersections=0; 

			//check if freept is contained in the hexagon. if not, don't compute energy for that point. move on. 
			for(int i=0; i<6; i++){		 
				bool intersect=intersection(Segment(Point(2,0), freept), boundary_segs[i]);
				if(intersect)
					num_intersections+=1; 			
			}
			if(!(num_intersections==1)){
				y_coor+=step; 
				 continue; 
			}
			std::vector<Triangle> triangle_array;
			for(int i=0; i<6; i++){
				triangle_array.push_back(Triangle(points[i], points[(i+1)%6], freept));
			}
	
			double total_hex_energy=0;
			for(int i=0; i<triangle_array.size(); i++){
				total_hex_energy+=triangle_Sb(triangle_array[i],weighted_circumcenter(triangle_array[i], zero_weights),2,1);
			}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_hex_energy <<std::endl; 
			
		y_coor+=step; 
		}
		x_coor+=step;
	}
	outputFile.close();
	
//////////////////////////////////////////////
////// Exp 4e. Sb/area^2 /////////////////////////
//////////////////////////////////////////////////

	outputFile.open("exp4/exp4e_hex_patch.txt");
	
	x_coor=x_min; 
	y_coor=y_min; 
	while(x_coor<x_max){	
		y_coor=y_min; 
		while(y_coor<y_max){		
		Point freept(x_coor, y_coor); 
			int num_intersections=0; 

			//check if freept is contained in the hexagon. if not, don't compute energy for that point. move on. 
			for(int i=0; i<6; i++){		 
				bool intersect=intersection(Segment(Point(2,0), freept), boundary_segs[i]);
				if(intersect)
					num_intersections+=1; 			
			}
			if(!(num_intersections==1)){
				y_coor+=step; 
				 continue; 
			}
			std::vector<Triangle> triangle_array;
			for(int i=0; i<6; i++){
				triangle_array.push_back(Triangle(points[i], points[(i+1)%6], freept));
			}
	
			double total_hex_energy=0;
			for(int i=0; i<triangle_array.size(); i++){
				total_hex_energy+=triangle_Sb(triangle_array[i],weighted_circumcenter(triangle_array[i], zero_weights),2,-1);
			}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_hex_energy <<std::endl; 
			
		y_coor+=step; 
		}
		x_coor+=step;
	}
	outputFile.close();
	
//////////////////////////////////////////////
////// Exp 4f. Sb/perim^4 /////////////////////////
//////////////////////////////////////////////////

	outputFile.open("exp4/exp4f_hex_patch.txt");
	
	x_coor=x_min; 
	y_coor=y_min; 
	while(x_coor<x_max){	
		y_coor=y_min; 
		while(y_coor<y_max){		
		Point freept(x_coor, y_coor); 
			int num_intersections=0; 

			//check if freept is contained in the hexagon. if not, don't compute energy for that point. move on. 
			for(int i=0; i<6; i++){		 
				bool intersect=intersection(Segment(Point(2,0), freept), boundary_segs[i]);
				if(intersect)
					num_intersections+=1; 			
			}
			if(!(num_intersections==1)){
				y_coor+=step; 
				 continue; 
			}
			std::vector<Triangle> triangle_array;
			for(int i=0; i<6; i++){
				triangle_array.push_back(Triangle(points[i], points[(i+1)%6], freept));
			}
	
			double total_hex_energy=0;
			for(int i=0; i<triangle_array.size(); i++){
				total_hex_energy+=triangle_Sb_divide_perim4(triangle_array[i],weighted_circumcenter(triangle_array[i], zero_weights));
			}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_hex_energy <<std::endl; 
			
		y_coor+=step; 
		}
		x_coor+=step;
	}
	outputFile.close();
	
	return 0; 
}


/*void hex_data_to_file(double x_min, double x_max, double y_min, double y_max, double step, int viewnum){

	std::vector<Point> points={Point(-1,0),  Point(cos(120*PI/180), sin(120*PI/180)),  Point(cos(60*PI/180), sin(60*PI/180)),  Point(1,0),  Point(cos(300*PI/180), sin(300*PI/180)), 				 Point(cos(240*PI/180), sin(240*PI/180))} ; 
	std::vector<Segment> boundary_segs; 
	for(int i=0; i<6; i++){
		boundary_segs.push_back(Segment(points[i], points[(i+1)%6])); 
	}

	std::ofstream outputFile;
	outputFile.open("exp4_hex_patch_view"+std::to_string(viewnum)+".txt");
	double x_coor=x_min; 
	double y_coor=y_min; 
	while(x_coor<x_max){	
		y_coor=y_min; 
		while(y_coor<y_max){		
		Point freept(x_coor, y_coor); 
			int num_intersections=0; 

			//check if freept is contained in the hexagon. if not, don't compute energy for that point. move on. 
			for(int i=0; i<6; i++){		 
				bool intersect=intersection(Segment(Point(2,0), freept), boundary_segs[i]);
				if(intersect)
					num_intersections+=1; 			
			}
			if(!(num_intersections==1)){
				y_coor+=step; 
				 continue; 
			}
			std::vector<Triangle> triangle_array;
			for(int i=0; i<6; i++){
				triangle_array.push_back(Triangle(points[i], points[(i+1)%6], freept));
			}
	
			double total_hex_energy=0;
			for(int i=0; i<triangle_array.size(); i++){
				total_hex_energy+=tri_energy<2,1>(triangle_array[i]);
			}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_hex_energy <<std::endl; 
			
		y_coor+=step; 
		}
		x_coor+=step;
	}
	outputFile.close();
	return; 

}
*/
