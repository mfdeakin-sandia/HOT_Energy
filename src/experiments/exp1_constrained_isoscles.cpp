// exp1_constrainted_isoscles.cpp

// Experiment 1: plot or minimize the energy for a triangle that is constrained to be isosceles, with all three vertices on a fixed circle, with two of the vertices moving towards one another and collapsing their shared edge to zero length. Confirm that the energy approaches zero.

#include <fstream>

#include "hot.hpp"
#include "lloyds.hpp"
#include "ply_writer.hpp"
#include "Sb.hpp"

int main(int argc, char **argv) {
	const int Wk=2; 
	Point p_fixed(-1,0); 
	double angle=PI; 
	
	std::ofstream outputFile;
///////////// Exp 1a	
	//star=0
	outputFile.open("exp1a_HOT_constrained_isoscles_star_2_0.txt");
	while(angle >.01){
		Point p_moving1(cos(angle), sin(angle)); 
		Point p_moving2(cos(angle), -sin(angle)); 
		Triangle triangle(p_moving1, p_moving2, p_fixed); 
		double energy= tri_energy<2,0>(triangle);
		outputFile << std::setw(15) << angle << std::setw(15) << energy<<std::endl; 
		angle-=.01;
		//std::vector<Point> tri_pts={p_moving1, p_moving2, p_fixed};
		//std::cout<< center_mass_polygon(tri_pts) << std::endl; 
	}
	outputFile.close(); 

	//star=1
	outputFile.open("exp1a_HOT_constrained_isoscles_star_2_1.txt");
	angle=PI;
	while(angle >.01){
		Point p_moving1(cos(angle), sin(angle)); 
		Point p_moving2(cos(angle), -sin(angle)); 
		Triangle triangle(p_moving1, p_moving2, p_fixed); 
		double energy= tri_energy<2,1>(triangle);
		outputFile << std::setw(15) << angle << std::setw(15) << energy<<std::endl; 
		angle-=.01;
	}
	outputFile.close(); 

	//star=2
	outputFile.open("exp1a_HOT_constrained_isoscles_star_2_2.txt");
	angle=PI; 
	while(angle >.01){
		Point p_moving1(cos(angle), sin(angle)); 
		Point p_moving2(cos(angle), -sin(angle)); 
		Triangle triangle(p_moving1, p_moving2, p_fixed); 
		double energy= tri_energy<2,2>(triangle);
		outputFile << std::setw(15) << angle << std::setw(15) << energy<<std::endl; 
		angle-=.01;
	}
	outputFile.close(); 
	
////////// Exp1b
///// ENERGY: HOT/TRIANGLE_AREA^2 ///////////////////////////////////////////

	//star=0
	outputFile.open("exp1b_HOT_constrained_isoscles_star_2_0.txt");
	angle=PI; 
	while(angle >.01){
		Point p_moving1(cos(angle), sin(angle)); 
		Point p_moving2(cos(angle), -sin(angle)); 
		Triangle triangle(p_moving1, p_moving2, p_fixed); 
		double energy= tri_energy_divideArea<2,0>(triangle,2); // HOT/area^2
		outputFile << std::setw(15) << angle << std::setw(15) << energy<<std::endl; 
		angle-=.01;
	}
	outputFile.close(); 

	//star=1
	outputFile.open("exp1b_HOT_constrained_isoscles_star_2_1.txt");
	angle=PI;
	while(angle >.01){
		Point p_moving1(cos(angle), sin(angle)); 
		Point p_moving2(cos(angle), -sin(angle)); 
		Triangle triangle(p_moving1, p_moving2, p_fixed); 
		double energy= tri_energy_divideArea<2,1>(triangle,2);
		
		//std::cout<< "tri energy: " << tri_energy<2,1>(triangle) << " tri area: " << triangle.area() <<std::endl; 
		outputFile << std::setw(15) << angle << std::setw(15) << energy<<std::endl; 
		angle-=.01;
	}
	outputFile.close(); 

	//star=2
	outputFile.open("exp1b_HOT_constrained_isoscles_star_2_2.txt");
	angle=PI; 
	while(angle >.01){
		Point p_moving1(cos(angle), sin(angle)); 
		Point p_moving2(cos(angle), -sin(angle)); 
		Triangle triangle(p_moving1, p_moving2, p_fixed); 
		double energy= tri_energy_divideArea<2,2>(triangle,2);
		outputFile << std::setw(15) << angle << std::setw(15) << energy<<std::endl; 
		angle-=.01;
	}
	outputFile.close(); 

///////////////Exp 1c
///// ENERGY: HOT/TRIANGLE_AREA ///////////////////////////////////////////

//star=0
	outputFile.open("exp1c_HOT_constrained_isoscles_star_2_0.txt");
	angle=PI; 
	while(angle >.01){
		Point p_moving1(cos(angle), sin(angle)); 
		Point p_moving2(cos(angle), -sin(angle)); 
		Triangle triangle(p_moving1, p_moving2, p_fixed); 
		double energy= tri_energy_divideArea<2,0>(triangle,1); // HOT/area^2
		outputFile << std::setw(15) << angle << std::setw(15) << energy<<std::endl; 
		angle-=.01;
	}
	outputFile.close(); 

	//star=1
	outputFile.open("exp1c_HOT_constrained_isoscles_star_2_1.txt");
	angle=PI;
	while(angle >.01){
		Point p_moving1(cos(angle), sin(angle)); 
		Point p_moving2(cos(angle), -sin(angle)); 
		Triangle triangle(p_moving1, p_moving2, p_fixed); 
		double energy= tri_energy_divideArea<2,1>(triangle,1);
		
		//std::cout<< "tri energy: " << tri_energy<2,1>(triangle) << " tri area: " << triangle.area() <<std::endl; 
		outputFile << std::setw(15) << angle << std::setw(15) << energy<<std::endl; 
		angle-=.01;
	}
	outputFile.close(); 

	//star=2
	outputFile.open("exp1c_HOT_constrained_isoscles_star_2_2.txt");
	angle=PI; 
	while(angle >.01){
		Point p_moving1(cos(angle), sin(angle)); 
		Point p_moving2(cos(angle), -sin(angle)); 
		Triangle triangle(p_moving1, p_moving2, p_fixed); 
		double energy= tri_energy_divideArea<2,2>(triangle,1);
		outputFile << std::setw(15) << angle << std::setw(15) << energy<<std::endl; 
		angle-=.01;
	}
	outputFile.close(); 

//////////Exp1d
////// ENERGY: Sb /////////////////////////////////////////////////////

	Point circumcenter(0,0); 
		//pow distance=2,powarea=-1; i.e. ||b-c||^2 *area
	outputFile.open("exp1d_Sb_constrained_isoscles.txt");
	angle=PI; 
	while(angle >.01){
		Point p_moving1(cos(angle), sin(angle)); 
		Point p_moving2(cos(angle), -sin(angle)); 
		Triangle triangle(p_moving1, p_moving2, p_fixed); 
		double energy= triangle_Sb(triangle, circumcenter,2,1);
		
		outputFile << std::setw(15) << angle << std::setw(15) << energy<<std::endl; 
		angle-=.01;
	}
	outputFile.close(); 

///////// Adjusting powers on terms in Sb ////////////////////////

/////////Exp 1e
//pow distance=2,powarea=-1; i.e. ||b-c||^2 /area
	outputFile.open("exp1e_Sb_constrained_isoscles.txt");
	angle=PI;
	while(angle >.01){
		Point p_moving1(cos(angle), sin(angle)); 
		Point p_moving2(cos(angle), -sin(angle)); 
		Triangle triangle(p_moving1, p_moving2, p_fixed); 
		double energy= triangle_Sb(triangle,circumcenter,2,-1);
		 
		outputFile << std::setw(15) << angle << std::setw(15) << energy<<std::endl; 
		angle-=.01;
	}
	outputFile.close(); 

/////////////////Exp 1f
//////////// Sb/perimeter

outputFile.open("exp1f_Sb_constrained_isoscles.txt");
	angle=PI;
	while(angle >.01){
		Point p_moving1(cos(angle), sin(angle)); 
		Point p_moving2(cos(angle), -sin(angle)); 
		Triangle triangle(p_moving1, p_moving2, p_fixed); 
		double energy= triangle_Sb_divide_perim4(triangle,circumcenter);
		 
		outputFile << std::setw(15) << angle << std::setw(15) << energy<<std::endl; 
		angle-=.01;
	}
	outputFile.close(); 


  	return 0;
}

