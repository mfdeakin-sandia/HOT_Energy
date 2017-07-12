#include <hot.hpp>
#include <ply_writer.hpp>

#include <random>
#include <vector>
#include <fstream>

#include <CGAL/Kernel/global_functions.h>


#define PI 3.14159265

int main(int argc, char **argv) {
	
	std::ofstream outputFile; 

	Point fixed_pt1(-1,0); 
	Point fixed_pt2(1,0);
	Point fixed_pt3(0,-sqrt(3));

	double x_step=.01; 
	double y_step=.01; 

	double y_coor_start=.01; 
	double free_x_coor=-2;
	double free_y_coor=y_coor_start;

	Triangle bottom_tri(fixed_pt1, fixed_pt2,fixed_pt3);
	Point freept(-.85,.21);
	Triangle tri(freept, fixed_pt1, fixed_pt2);
	//double energy=Edge_Energy<2,1>(Triangle(Point(-.85,.02), fixed_pt1, fixed_pt2), 0, bottom_tri, 2, true); 
	double energy=Edge_Energy<2,1>(tri, 0, bottom_tri, 2, true); 
	std::cout<<energy <<std::endl; 
 /*
	////////star^0
	outputFile.open("exp7_vertex_to_fixed_edge_correctedformulas_2_0.txt"); 
	while(free_x_coor<2){
		free_y_coor=y_coor_start; 
		while(free_y_coor <2){	
			Point freept(free_x_coor, free_y_coor); 
			Triangle tri(freept, fixed_pt1, fixed_pt2); 
			
			double energy=tri_energy<2,0>(tri); 
			outputFile << std::setw(15) << free_x_coor << std::setw(15) << free_y_coor <<  std::setw(15) << energy <<std::endl; 

			free_y_coor+=y_step; 
		}

	free_x_coor+=x_step;
	}

	outputFile.close(); 

	free_x_coor=-2;
	free_y_coor=.1;

	*/
	////////star^1
	outputFile.open("exp7_vertex_to_fixed_edge_corrected.txt"); 
	while(free_x_coor<1.99){
		free_y_coor=y_coor_start; 
		while(free_y_coor <2){	
			Point freept(free_x_coor, free_y_coor); 
			Triangle tri(freept, fixed_pt1, fixed_pt2); 
			double energy=Edge_Energy<2,1>(tri,0,bottom_tri,2,true);
			outputFile << std::setw(15) << free_x_coor << std::setw(15) << free_y_coor <<  std::setw(15) << energy <<std::endl; 

			free_y_coor+=y_step; 
		}

	free_x_coor+=x_step;
	}

	outputFile.close(); 

	free_x_coor=-2;
	free_y_coor=y_coor_start;
	outputFile.open("exp7_vertex_to_fixed_edge_uncorrected.txt"); 
	while(free_x_coor<1.9){
		free_y_coor=y_coor_start; 
		while(free_y_coor <2){	
			Point freept(free_x_coor, free_y_coor); 
			Triangle tri(freept, fixed_pt1, fixed_pt2); 
			double energy=Edge_Energy<2,1>(tri,0,bottom_tri,2,false);
			outputFile << std::setw(15) << free_x_coor << std::setw(15) << free_y_coor <<  std::setw(15) << energy <<std::endl; 

			free_y_coor+=y_step; 
		}

	free_x_coor+=x_step;
	}

	outputFile.close(); 

/*

	///////star^2
	free_x_coor=-2;
	free_y_coor=.1;

	outputFile.open("exp2_vertex_to_fixed_edge_2_2.txt"); 
	while(free_x_coor<2){
		free_y_coor=y_coor_start; 
		while(free_y_coor <2){	
			Point freept(free_x_coor, free_y_coor); 
			Triangle tri(freept, fixed_pt1, fixed_pt2); 
			
			double energy=tri_energy<2,2>(tri); 
			outputFile << std::setw(15) << free_x_coor << std::setw(15) << free_y_coor <<  std::setw(15) << energy <<std::endl; 

			free_y_coor+=y_step; 
		}

	free_x_coor+=x_step;
	}

	outputFile.close(); 
*/

  	return 0;
}
