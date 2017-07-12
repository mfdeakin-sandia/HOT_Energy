#include <hot.hpp>
#include <Sb.hpp>
#include <ply_writer.hpp>

#include <random>
#include <vector>
#include <fstream>

#include <CGAL/Kernel/global_functions.h>


// ****USING INCORRECT FORMULAS **********
int main(int argc, char **argv) {
	
	//these points stay fixed in the experiment
	Point leftpt(-2,0); 
	Point rightpt(2,0); 
	Point toppt(0,1); 
	Point bottompt(0,-1); 
	double half_length_initial=1.9;


	// This is 1/2 length of horizontal edge we are shrinking
	std::ofstream outputFile; 

///////////////// Exp 3a HOT //////////////////////////////////////////////////////////////

//////////star 0
	double half_length=half_length_initial; 
	outputFile.open("exp3/exp3a_collapse_no_large_angles_2_0.txt"); 
	while(half_length >.01){

		Point inner_right_pt(half_length, 0); 
		Point inner_left_pt(-half_length, 0); 

		Triangle tri1(leftpt, toppt, inner_left_pt); 
		Triangle tri2(inner_left_pt, toppt, inner_right_pt); 
		Triangle tri3(inner_right_pt, toppt, rightpt); 

		Triangle tri4(leftpt, inner_left_pt, bottompt);
		Triangle tri5(inner_left_pt, inner_right_pt, bottompt); 
		Triangle tri6(inner_right_pt, rightpt, bottompt);

		double energy=0; 
		energy+= tri_energy<2,0>(tri1); 
		energy+= tri_energy<2,0>(tri2); 
		energy+= tri_energy<2,0>(tri3); 
		energy+= tri_energy<2,0>(tri4); 
		energy+= tri_energy<2,0>(tri5); 
		energy+= tri_energy<2,0>(tri6); 
		
		outputFile<< std::setw(15) << half_length << std::setw(15) << energy <<std::endl; 
		half_length-= .01; 
	}
	outputFile.close(); 

//////star 1
	half_length=half_length_initial; 
	outputFile.open("exp3/exp3a_collapse_no_large_angles_2_1.txt"); 
	while(half_length >.01){

		Point inner_right_pt(half_length, 0); 
		Point inner_left_pt(-half_length, 0); 

		Triangle tri1(leftpt, toppt, inner_left_pt); 
		Triangle tri2(inner_left_pt, toppt, inner_right_pt); 
		Triangle tri3(inner_right_pt, toppt, rightpt); 

		Triangle tri4(leftpt, inner_left_pt, bottompt);
		Triangle tri5(inner_left_pt, inner_right_pt, bottompt); 
		Triangle tri6(inner_right_pt, rightpt, bottompt);

		double energy=0; 
		energy+= tri_energy<2,1>(tri1); 
		energy+= tri_energy<2,1>(tri2); 
		energy+= tri_energy<2,1>(tri3); 
		energy+= tri_energy<2,1>(tri4); 
		energy+= tri_energy<2,1>(tri5); 
		energy+= tri_energy<2,1>(tri6); 
		
		outputFile<< std::setw(15) << half_length << std::setw(15) << energy <<std::endl; 
		half_length-= .01; 
	}
	outputFile.close();
 
//////star 2
	half_length=1.9; 
	outputFile.open("exp3/exp3a_collapse_no_large_angles_2_2.txt"); 
	while(half_length >.01){

		Point inner_right_pt(half_length, 0); 
		Point inner_left_pt(-half_length, 0); 

		Triangle tri1(leftpt, toppt, inner_left_pt); 
		Triangle tri2(inner_left_pt, toppt, inner_right_pt); 
		Triangle tri3(inner_right_pt, toppt, rightpt); 

		Triangle tri4(leftpt, inner_left_pt, bottompt);
		Triangle tri5(inner_left_pt, inner_right_pt, bottompt); 
		Triangle tri6(inner_right_pt, rightpt, bottompt);

		double energy=0; 
		energy+= tri_energy<2,2>(tri1); 
		energy+= tri_energy<2,2>(tri2); 
		energy+= tri_energy<2,2>(tri3); 
		energy+= tri_energy<2,2>(tri4); 
		energy+= tri_energy<2,2>(tri5); 
		energy+= tri_energy<2,2>(tri6); 
		
		outputFile<< std::setw(15) << half_length << std::setw(15) << energy <<std::endl; 
		half_length-= .01; 
	}
	outputFile.close(); 

//////////////////////////////////////////////////////
////////////////// Exp 3b HOT/area^2	
//////////////////////////////////////////////


//////////star 0
	half_length=half_length_initial; 
	outputFile.open("exp3/exp3b_collapse_no_large_angles_2_0.txt"); 
	while(half_length >.01){

		Point inner_right_pt(half_length, 0); 
		Point inner_left_pt(-half_length, 0); 

		Triangle tri1(leftpt, toppt, inner_left_pt); 
		Triangle tri2(inner_left_pt, toppt, inner_right_pt); 
		Triangle tri3(inner_right_pt, toppt, rightpt); 

		Triangle tri4(leftpt, inner_left_pt, bottompt);
		Triangle tri5(inner_left_pt, inner_right_pt, bottompt); 
		Triangle tri6(inner_right_pt, rightpt, bottompt);

		double energy=0; 
		energy+= tri_energy_divideArea<2,0>(tri1,2); 
		energy+= tri_energy_divideArea<2,0>(tri2,2); 
		energy+= tri_energy_divideArea<2,0>(tri3,2); 
		energy+= tri_energy_divideArea<2,0>(tri4,2); 
		energy+= tri_energy_divideArea<2,0>(tri5,2); 
		energy+= tri_energy_divideArea<2,0>(tri6,2); 
		
		outputFile<< std::setw(15) << half_length << std::setw(15) << energy <<std::endl; 
		half_length-= .01; 
	}
	outputFile.close(); 

//////star 1
	half_length=half_length_initial; 
	outputFile.open("exp3/exp3b_collapse_no_large_angles_2_1.txt"); 
	while(half_length >.01){

		Point inner_right_pt(half_length, 0); 
		Point inner_left_pt(-half_length, 0); 

		Triangle tri1(leftpt, toppt, inner_left_pt); 
		Triangle tri2(inner_left_pt, toppt, inner_right_pt); 
		Triangle tri3(inner_right_pt, toppt, rightpt); 

		Triangle tri4(leftpt, inner_left_pt, bottompt);
		Triangle tri5(inner_left_pt, inner_right_pt, bottompt); 
		Triangle tri6(inner_right_pt, rightpt, bottompt);

		double energy=0; 
		energy+= tri_energy_divideArea<2,1>(tri1,2); 
		energy+= tri_energy_divideArea<2,1>(tri2,2); 
		energy+= tri_energy_divideArea<2,1>(tri3,2); 
		energy+= tri_energy_divideArea<2,1>(tri4,2); 
		energy+= tri_energy_divideArea<2,1>(tri5,2); 
		energy+= tri_energy_divideArea<2,1>(tri6,2); 
		
		outputFile<< std::setw(15) << half_length << std::setw(15) << energy <<std::endl; 
		half_length-= .01; 
	}
	outputFile.close();
 
//////star 2
	half_length=half_length_initial; 
	outputFile.open("exp3/exp3b_collapse_no_large_angles_2_2.txt"); 
	while(half_length >.01){

		Point inner_right_pt(half_length, 0); 
		Point inner_left_pt(-half_length, 0); 

		Triangle tri1(leftpt, toppt, inner_left_pt); 
		Triangle tri2(inner_left_pt, toppt, inner_right_pt); 
		Triangle tri3(inner_right_pt, toppt, rightpt); 

		Triangle tri4(leftpt, inner_left_pt, bottompt);
		Triangle tri5(inner_left_pt, inner_right_pt, bottompt); 
		Triangle tri6(inner_right_pt, rightpt, bottompt);

		double energy=0; 
		energy+= tri_energy_divideArea<2,2>(tri1,2); 
		energy+= tri_energy_divideArea<2,2>(tri2,2); 
		energy+= tri_energy_divideArea<2,2>(tri3,2); 
		energy+= tri_energy_divideArea<2,2>(tri4,2); 
		energy+= tri_energy_divideArea<2,2>(tri5,2); 
		energy+= tri_energy_divideArea<2,2>(tri6,2); 
		
		outputFile<< std::setw(15) << half_length << std::setw(15) << energy <<std::endl; 
		half_length-= .01; 
	}
	outputFile.close(); 

//////////////////////////////////////////////////////
////////////////// Exp 3c HOT/area	
//////////////////////////////////////////////


//////////star 0
	half_length=half_length_initial; 
	outputFile.open("exp3/exp3c_collapse_no_large_angles_2_0.txt"); 
	while(half_length >.01){

		Point inner_right_pt(half_length, 0); 
		Point inner_left_pt(-half_length, 0); 

		Triangle tri1(leftpt, toppt, inner_left_pt); 
		Triangle tri2(inner_left_pt, toppt, inner_right_pt); 
		Triangle tri3(inner_right_pt, toppt, rightpt); 

		Triangle tri4(leftpt, inner_left_pt, bottompt);
		Triangle tri5(inner_left_pt, inner_right_pt, bottompt); 
		Triangle tri6(inner_right_pt, rightpt, bottompt);

		double energy=0; 
		energy+= tri_energy_divideArea<2,0>(tri1,1); 
		energy+= tri_energy_divideArea<2,0>(tri2,1); 
		energy+= tri_energy_divideArea<2,0>(tri3,1); 
		energy+= tri_energy_divideArea<2,0>(tri4,1); 
		energy+= tri_energy_divideArea<2,0>(tri5,1); 
		energy+= tri_energy_divideArea<2,0>(tri6,1); 
		
		outputFile<< std::setw(15) << half_length << std::setw(15) << energy <<std::endl; 
		half_length-= .01; 
	}
	outputFile.close(); 

//////star 1
	half_length=half_length_initial; 
	outputFile.open("exp3/exp3c_collapse_no_large_angles_2_1.txt"); 
	while(half_length >.01){

		Point inner_right_pt(half_length, 0); 
		Point inner_left_pt(-half_length, 0); 

		Triangle tri1(leftpt, toppt, inner_left_pt); 
		Triangle tri2(inner_left_pt, toppt, inner_right_pt); 
		Triangle tri3(inner_right_pt, toppt, rightpt); 

		Triangle tri4(leftpt, inner_left_pt, bottompt);
		Triangle tri5(inner_left_pt, inner_right_pt, bottompt); 
		Triangle tri6(inner_right_pt, rightpt, bottompt);

		double energy=0; 
		energy+= tri_energy_divideArea<2,1>(tri1,1); 
		energy+= tri_energy_divideArea<2,1>(tri2,1); 
		energy+= tri_energy_divideArea<2,1>(tri3,1); 
		energy+= tri_energy_divideArea<2,1>(tri4,1); 
		energy+= tri_energy_divideArea<2,1>(tri5,1); 
		energy+= tri_energy_divideArea<2,1>(tri6,1); 
		
		outputFile<< std::setw(15) << half_length << std::setw(15) << energy <<std::endl; 
		half_length-= .01; 
	}
	outputFile.close();
 
//////star 2
	half_length=half_length_initial; 
	outputFile.open("exp3/exp3c_collapse_no_large_angles_2_2.txt"); 
	while(half_length >.01){

		Point inner_right_pt(half_length, 0); 
		Point inner_left_pt(-half_length, 0); 

		Triangle tri1(leftpt, toppt, inner_left_pt); 
		Triangle tri2(inner_left_pt, toppt, inner_right_pt); 
		Triangle tri3(inner_right_pt, toppt, rightpt); 

		Triangle tri4(leftpt, inner_left_pt, bottompt);
		Triangle tri5(inner_left_pt, inner_right_pt, bottompt); 
		Triangle tri6(inner_right_pt, rightpt, bottompt);

		double energy=0; 
		energy+= tri_energy_divideArea<2,2>(tri1,1); 
		energy+= tri_energy_divideArea<2,2>(tri2,1); 
		energy+= tri_energy_divideArea<2,2>(tri3,1); 
		energy+= tri_energy_divideArea<2,2>(tri4,1); 
		energy+= tri_energy_divideArea<2,2>(tri5,1); 
		energy+= tri_energy_divideArea<2,2>(tri6,1); 
		
		outputFile<< std::setw(15) << half_length << std::setw(15) << energy <<std::endl; 
		half_length-= .01; 
	}
	outputFile.close(); 

//////////////////////////////////////////////////////
////////////////// Exp 3d Energy=Sb	
//////////////////////////////////////////////
	half_length=half_length_initial; 
	outputFile.open("exp3/exp3d_collapse_no_large_angles.txt"); 
	while(half_length >.01){

		Point inner_right_pt(half_length, 0); 
		Point inner_left_pt(-half_length, 0); 

		Triangle tri1(leftpt, toppt, inner_left_pt); 
		Triangle tri2(inner_left_pt, toppt, inner_right_pt); 
		Triangle tri3(inner_right_pt, toppt, rightpt); 

		Triangle tri4(leftpt, inner_left_pt, bottompt);
		Triangle tri5(inner_left_pt, inner_right_pt, bottompt); 
		Triangle tri6(inner_right_pt, rightpt, bottompt);

		double energy=0; 
		double zero_weights[]={0,0,0}; 
		energy+= triangle_Sb(tri1,weighted_circumcenter(tri1,zero_weights),2,1); 
		energy+= triangle_Sb(tri2,weighted_circumcenter(tri2,zero_weights),2,1); 
		energy+= triangle_Sb(tri3,weighted_circumcenter(tri3,zero_weights),2,1); 
		energy+= triangle_Sb(tri4,weighted_circumcenter(tri4,zero_weights),2,1); 
		energy+= triangle_Sb(tri5,weighted_circumcenter(tri5,zero_weights),2,1); 
		energy+= triangle_Sb(tri6,weighted_circumcenter(tri6,zero_weights),2,1); 
		
		outputFile<< std::setw(15) << half_length << std::setw(15) << energy <<std::endl; 
		half_length-= .01; 
	}
	outputFile.close(); 

//////////////////////////////////////////////////////
////////////////// Exp 3e Energy=Sb/area^2
//////////////////////////////////////////////
	half_length=half_length_initial; 
	outputFile.open("exp3/exp3e_collapse_no_large_angles.txt"); 
	while(half_length >.01){

		Point inner_right_pt(half_length, 0); 
		Point inner_left_pt(-half_length, 0); 

		Triangle tri1(leftpt, toppt, inner_left_pt); 
		Triangle tri2(inner_left_pt, toppt, inner_right_pt); 
		Triangle tri3(inner_right_pt, toppt, rightpt); 

		Triangle tri4(leftpt, inner_left_pt, bottompt);
		Triangle tri5(inner_left_pt, inner_right_pt, bottompt); 
		Triangle tri6(inner_right_pt, rightpt, bottompt);

		double energy=0; 
		double zero_weights[]={0,0,0}; 
		energy+= triangle_Sb(tri1,weighted_circumcenter(tri1,zero_weights),2,-1); 
		energy+= triangle_Sb(tri2,weighted_circumcenter(tri2,zero_weights),2,-1); 
		energy+= triangle_Sb(tri3,weighted_circumcenter(tri3,zero_weights),2,-1); 
		energy+= triangle_Sb(tri4,weighted_circumcenter(tri4,zero_weights),2,-1); 
		energy+= triangle_Sb(tri5,weighted_circumcenter(tri5,zero_weights),2,-1); 
		energy+= triangle_Sb(tri6,weighted_circumcenter(tri6,zero_weights),2,-1); 
		
		outputFile<< std::setw(15) << half_length << std::setw(15) << energy <<std::endl; 
		half_length-= .01; 
	}
	outputFile.close(); 

//////////////////////////////////////////////////////
////////////////// Exp 3f Energy=Sb/perim^4
//////////////////////////////////////////////
	half_length=half_length_initial; 
	outputFile.open("exp3/exp3f_collapse_no_large_angles.txt"); 
	while(half_length >.01){

		Point inner_right_pt(half_length, 0); 
		Point inner_left_pt(-half_length, 0); 

		Triangle tri1(leftpt, toppt, inner_left_pt); 
		Triangle tri2(inner_left_pt, toppt, inner_right_pt); 
		Triangle tri3(inner_right_pt, toppt, rightpt); 

		Triangle tri4(leftpt, inner_left_pt, bottompt);
		Triangle tri5(inner_left_pt, inner_right_pt, bottompt); 
		Triangle tri6(inner_right_pt, rightpt, bottompt);

		double energy=0; 
		double zero_weights[]={0,0,0}; 
		energy+= triangle_Sb_divide_perim4(tri1,weighted_circumcenter(tri1,zero_weights)); 
		energy+= triangle_Sb_divide_perim4(tri2,weighted_circumcenter(tri2,zero_weights)); 
		energy+= triangle_Sb_divide_perim4(tri3,weighted_circumcenter(tri3,zero_weights)); 
		energy+= triangle_Sb_divide_perim4(tri4,weighted_circumcenter(tri4,zero_weights)); 
		energy+= triangle_Sb_divide_perim4(tri5,weighted_circumcenter(tri5,zero_weights)); 
		energy+= triangle_Sb_divide_perim4(tri6,weighted_circumcenter(tri6,zero_weights)); 
		
		outputFile<< std::setw(15) << half_length << std::setw(15) << energy <<std::endl; 
		half_length-= .01; 
	}
	outputFile.close(); 


	return 0; 
}
