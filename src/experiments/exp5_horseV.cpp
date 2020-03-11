//exp5_horseV.cpp
#include <fstream>

#include "hot.hpp"
#include "Sb.hpp"
#include "ply_writer.hpp"

int main(int argc, char **argv) {
	std::ofstream outputFile; 


	//EK_real fig_height=8;
	//EK_real fig_width=2;
	//EK_real step_size=.01; 

	//EK_real x_coor= -fig_width/2+.01;
	//EK_real y_coor= .01; 
	 

	// we could do without these doubles if we changed the input type for tri_energy<Wk,star>
	double fig_height=8;
	double fig_width=2; 
	double step_size=.01; 
	
	//double x_coor_d=-.99; 
	//double y_coor_d=.01; 

	double x_coor=-.99; 
	double y_coor=.01; 
	

	//EK::Point_2 exact_points[]={EK::Point_2(-fig_width/2,0), EK::Point_2(0,fig_height/2), EK::Point_2(fig_width/2,0), EK::Point_2(0,fig_height)};
	Point points[]={Point(-fig_width/2,0), Point(0,fig_height/2), Point(fig_width/2,0), Point(0,fig_height)};	

/////////////////////////////////////////////////////////////////
///////////////////// Exp 5a: HOT ///////////////////////////////
/////////////////////////////////////////////////////////////
	//star 0
	outputFile.open("exp5a_horseV_star_2_0.txt");
	while(x_coor< fig_width/2-.01){
		y_coor=.01;
		//y_coor_d=.01; 
		while(y_coor<fig_height){
			Point freept(x_coor, y_coor); 
			//check if freept is contained in non-convex polygon. if not, don't compute energy for that triangulation. move on. 
			/*if(point_outside_domain(Point(x_coor,y_coor), points, points+4, K())){
				y_coor+=step_size;
				//y_coor_d+=step_size_d;
				continue; 
			}
*/
			std::vector<Triangle> triangle_array;
				for(int i=0; i<4; i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%4], freept));
				}
	
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=tri_energy<2,0>(triangle_array[i]);
				}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		y_coor+=step_size; 
		//y_coor_d+=step_size_d; 
		}
		x_coor+=step_size;
		//x_coor_d+=step_size_d;  
	}
	outputFile.close(); 


// star 1

	outputFile.open("exp5a_horseV_star_2_1.txt");
	x_coor= -fig_width/2+.01;
	y_coor= .01;
	while(x_coor< fig_width/2-.01){
		y_coor=.01;
		while(y_coor<fig_height){
			Point freept(x_coor, y_coor); 
			//check if freept is contained in non-convex polygon. if not, don't compute energy for that triangulation. move on. 
			/*if(point_outside_domain(Point(x_coor,y_coor), points, points+4, K())){
				y_coor+=step_size;
				continue; 
			}
*/
			std::vector<Triangle> triangle_array;
				for(int i=0; i<4; i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%4], freept));
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

// star 2

	outputFile.open("exp5a_horseV_star_2_2.txt");
	x_coor= -fig_width/2+.01;
	y_coor= .01;
	while(x_coor< fig_width/2-.01){
		y_coor=.01; 
		while(y_coor<fig_height){
			Point freept(x_coor, y_coor); 
			//check if freept is contained in non-convex polygon. if not, don't compute energy for that triangulation. move on. 
			/*if(point_outside_domain(Point(x_coor,y_coor), points, points+4, K())){
				y_coor+=step_size;
				continue; 
			}
*/
			std::vector<Triangle> triangle_array;
				for(int i=0; i<4; i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%4], freept));
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


////////////////////////////////////////////////////////
//////////////////////  Exp 5b  Hot/area^2 //////////////
//////////////////////////////////////////////////////

	//star 0
	outputFile.open("exp5b_horseV_star_2_0.txt");
	x_coor= -fig_width/2+.01;
	y_coor= .01;
	while(x_coor< fig_width/2-.01){
		y_coor=.01; 
		while(y_coor<fig_height){
			Point freept(x_coor, y_coor); 
			//check if freept is contained in non-convex polygon. if not, don't compute energy for that triangulation. move on. 
			/*if(point_outside_domain(Point(x_coor,y_coor), points, points+4, K())){
				y_coor+=step_size;
				continue; 
			}
*/
			std::vector<Triangle> triangle_array;
				for(int i=0; i<4; i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%4], freept));
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


// star 1

	outputFile.open("exp5b_horseV_star_2_1.txt");
	x_coor= -fig_width/2+.01;
	y_coor= .01;
	while(x_coor< fig_width/2-.01){
		y_coor=.01;
		while(y_coor<fig_height){
			Point freept(x_coor, y_coor); 
			//check if freept is contained in non-convex polygon. if not, don't compute energy for that triangulation. move on. 
			/*if(point_outside_domain(Point(x_coor,y_coor), points, points+4, K())){
				y_coor+=step_size;
				continue; 
			}
*/
			std::vector<Triangle> triangle_array;
				for(int i=0; i<4; i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%4], freept));
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

// star 2

	outputFile.open("exp5b_horseV_star_2_2.txt"); 
	x_coor= -fig_width/2+.01;
	y_coor= .01;
	while(x_coor< fig_width/2-.01){
		y_coor=.01; 
		while(y_coor<fig_height){
			Point freept(x_coor, y_coor); 
			//check if freept is contained in non-convex polygon. if not, don't compute energy for that triangulation. move on. 
			/*if(point_outside_domain(Point(x_coor,y_coor), points, points+4, K())){
				y_coor+=step_size;
				continue; 
			}
*/
			std::vector<Triangle> triangle_array;
				for(int i=0; i<4; i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%4], freept));
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


////////////////////////////////////////////////////////
//////////////////////  Exp 5c  Hot/area //////////////
//////////////////////////////////////////////////////

	//star 0
	outputFile.open("exp5c_horseV_star_2_0.txt"); 
	x_coor= -fig_width/2+.01;
	y_coor= .01;
	while(x_coor< fig_width/2-.01){
		y_coor=.01; 
		while(y_coor<fig_height){
			Point freept(x_coor, y_coor); 
			//check if freept is contained in non-convex polygon. if not, don't compute energy for that triangulation. move on. 
			/*if(point_outside_domain(Point(x_coor,y_coor), points, points+4, K())){
				y_coor+=step_size;
				continue; 
			}
*/
			std::vector<Triangle> triangle_array;
				for(int i=0; i<4; i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%4], freept));
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


// star 1

	outputFile.open("exp5c_horseV_star_2_1.txt"); 
	x_coor= -fig_width/2+.01;
	y_coor= .01;
	while(x_coor< fig_width/2-.01){
		y_coor=.01;
		while(y_coor<fig_height){
			Point freept(x_coor, y_coor); 
			//check if freept is contained in non-convex polygon. if not, don't compute energy for that triangulation. move on. 
			/*if(point_outside_domain(Point(x_coor,y_coor), points, points+4, K())){
				y_coor+=step_size;
				continue; 
			}
*/
			std::vector<Triangle> triangle_array;
				for(int i=0; i<4; i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%4], freept));
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

// star 2

	outputFile.open("exp5c_horseV_star_2_2.txt"); 
	x_coor= -fig_width/2+.01;
	y_coor= .01;
	while(x_coor< fig_width/2-.01){
		y_coor=.01; 
		while(y_coor<fig_height){
			Point freept(x_coor, y_coor); 
			//check if freept is contained in non-convex polygon. if not, don't compute energy for that triangulation. move on. 
			/*if(point_outside_domain(Point(x_coor,y_coor), points, points+4, K())){
				y_coor+=step_size;
				continue; 
			}
*/
			std::vector<Triangle> triangle_array;
				for(int i=0; i<4; i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%4], freept));
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

////////////////////////////////////////////////////////
//////////////////////  Exp 5d  Sb //////////////
//////////////////////////////////////////////////////
	double zero_weights[]={0,0,0};
	
	outputFile.open("exp5d_horseV.txt"); 
	x_coor= -fig_width/2+.01;
	y_coor= .01;
	while(x_coor< fig_width/2-.01){
		y_coor=.01; 
		while(y_coor<fig_height){
			Point freept(x_coor, y_coor); 
			//check if freept is contained in non-convex polygon. if not, don't compute energy for that triangulation. move on. 
			/*if(point_outside_domain(Point(x_coor,y_coor), points, points+4, K())){
				y_coor+=step_size;
				continue; 
			}
*/
			std::vector<Triangle> triangle_array;
				for(int i=0; i<4; i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%4], freept));
				}
	
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=triangle_Sb(triangle_array[i],weighted_circumcenter(triangle_array[i], zero_weights),2,1);
				}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		y_coor+=step_size; 
		}
		x_coor+=step_size;  
	}
	outputFile.close(); 

////////////////////////////////////////////////////////
//////////////////////  Exp 5e  Sb/area^2 //////////////
//////////////////////////////////////////////////////
	
	outputFile.open("exp5e_horseV.txt"); 
	x_coor= -fig_width/2+.01;
	y_coor= .01;
	while(x_coor< fig_width/2-.01){
		y_coor=.01; 
		while(y_coor<fig_height){
			Point freept(x_coor, y_coor); 
			//check if freept is contained in non-convex polygon. if not, don't compute energy for that triangulation. move on. 
			/*if(point_outside_domain(Point(x_coor,y_coor), points, points+4, K())){
				y_coor+=step_size;
				continue; 
			}
*/
			std::vector<Triangle> triangle_array;
				for(int i=0; i<4; i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%4], freept));
				}
	
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=triangle_Sb(triangle_array[i],weighted_circumcenter(triangle_array[i], zero_weights),2,-1);
				}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		y_coor+=step_size; 
		}
		x_coor+=step_size;  
	}
	outputFile.close();

////////////////////////////////////////////////////////
//////////////////////  Exp 5f  Sb/perim^4 //////////////
//////////////////////////////////////////////////////
	
	outputFile.open("exp5f_horseV.txt"); 
	x_coor= -fig_width/2+.01;
	y_coor= .01;
	while(x_coor< fig_width/2-.01){
		y_coor=.01; 
		while(y_coor<fig_height){
			Point freept(x_coor, y_coor); 
			//check if freept is contained in non-convex polygon. if not, don't compute energy for that triangulation. move on. 
			/*if(point_outside_domain(Point(x_coor,y_coor), points, points+4, K())){
				y_coor+=step_size;
				continue; 
			}
*/
			std::vector<Triangle> triangle_array;
				for(int i=0; i<4; i++){
					triangle_array.push_back(Triangle(points[i], points[(i+1)%4], freept));
				}
	
				double total_energy=0;
				for(int i=0; i<triangle_array.size(); i++){
					total_energy+=triangle_Sb_divide_perim4(triangle_array[i],weighted_circumcenter(triangle_array[i], zero_weights));
				}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		y_coor+=step_size; 
		}
		x_coor+=step_size;  
	}
	outputFile.close();

	return 0;
}
