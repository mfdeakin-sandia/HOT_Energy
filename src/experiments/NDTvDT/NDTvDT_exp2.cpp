#include <hot.hpp>



#include <random>
#include <vector>
#include <fstream>

#include <CGAL/Kernel/global_functions.h>


#define PI 3.14159265

int main(int argc, char **argv) {
	std::ofstream outputFile; 


//star 0	
	outputFile.open("NDTvDT/NDTvDT_exp2_star0.txt"); 
	double height=1;
	while( height>0){
		Triangle DTtri1(Point(-1,0), Point(-1,height), Point(0,-.5));
		Triangle DTtri2(Point(-1,height), Point (1,0), Point(0,-.5)); 
		
		K_real DTtri1_energy=tri_energy<2,0>(DTtri1);
		K_real DTtri2_energy=tri_energy<2,0>(DTtri2);
		K_real DT_energy=DTtri1_energy+DTtri2_energy; 
	
		Triangle NDTtri1(Point(-1,0), Point(1,0), Point(-1,height)); 
		Triangle NDTtri2(Point(-1,0), Point(0,-.5), Point(1,0)); 

		K_real NDTtri1_energy=tri_energy<2,1>(NDTtri1); 
		K_real NDTtri2_energy=tri_energy<2,1>(NDTtri2); 
		K_real NDT_energy=NDTtri1_energy + NDTtri2_energy; 
	
		outputFile<< std::setw(8) << height << std::setw(15) << DT_energy << std::setw(15) << NDT_energy << std::endl;

		height-=.01;

	}
	outputFile.close(); 


	//star 1
	
	outputFile.open("NDTvDT/NDTvDT_exp2_star1.txt"); 
	height=1;	
	while( height>0){
		Triangle DTtri1(Point(-1,0), Point(-1,height), Point(0,-.5));
		Triangle DTtri2(Point(-1,height), Point (1,0), Point(0,-.5)); 
		
		K_real DTtri1_energy=tri_energy<2,1>(DTtri1);
		K_real DTtri2_energy=tri_energy<2,1>(DTtri2);
		K_real DT_energy=DTtri1_energy+DTtri2_energy; 
	
		Triangle NDTtri1(Point(-1,0), Point(1,0), Point(-1,height)); 
		Triangle NDTtri2(Point(-1,0), Point(0,-.5), Point(1,0)); 

		K_real NDTtri1_energy=tri_energy<2,1>(NDTtri1); 
		K_real NDTtri2_energy=tri_energy<2,1>(NDTtri2); 
		K_real NDT_energy=NDTtri1_energy + NDTtri2_energy; 
	
		outputFile<< std::setw(8) << height << std::setw(15) << DT_energy << std::setw(15) << NDT_energy << std::endl;

		height-=.01;

	}
	outputFile.close(); 

	//star 2
	
	outputFile.open("NDTvDT/NDTvDT/NDTvDT_exp2_star2.txt"); 
	height=1;	
	while( height>0){
		Triangle DTtri1(Point(-1,0), Point(-1,height), Point(0,-.5));
		Triangle DTtri2(Point(-1,height), Point (1,0), Point(0,-.5)); 
		
		K_real DTtri1_energy=tri_energy<2,2>(DTtri1);
		K_real DTtri2_energy=tri_energy<2,2>(DTtri2);
		K_real DT_energy=DTtri1_energy+DTtri2_energy; 
		
	
	
		Triangle NDTtri1(Point(-1,0), Point(1,0), Point(-1,height)); 
		Triangle NDTtri2(Point(-1,0), Point(0,-.5), Point(1,0)); 

		K_real NDTtri1_energy=tri_energy<2,2>(NDTtri1); 
		K_real NDTtri2_energy=tri_energy<2,2>(NDTtri2); 
		K_real NDT_energy=NDTtri1_energy + NDTtri2_energy; 
	
		outputFile<< std::setw(8) << height << std::setw(15) << DT_energy << std::setw(15) << NDT_energy << std::endl;

		height-=.01;

	}
	outputFile.close(); 

	return 0;
}
