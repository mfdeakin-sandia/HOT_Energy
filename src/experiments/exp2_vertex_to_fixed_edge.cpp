// exp2_vertex_to_fixed_edge.cpp

#include <fstream>

#include "hot.hpp"
#include "Sb.hpp"
#include "ply_writer.hpp"


int main(int argc, char **argv) {
/*
	std::cout<< std::endl <<"Experiment to see effect of using incorrect appendix formulas" << std::endl;

	double degxi=7*PI/6, degxj=11*PI/6, degxl=3*PI/2; 

	double xkheight=1;


	
	Triangle tri2_stationary( Point(cos(degxi), sin(degxi)), Point(cos(degxj), sin(degxj)), Point(cos(degxl), sin(degxl))); 
	double hl=signed_dist_circumcenters(tri2_stationary, 2);
	//std::cout<< "hl: "<< hl <<std::endl;
	
	double moving_subtri_energy=0, moving_tri_energy=0; 
	std::cout<< std::setw(10) << "xkheight"<< std::setw(15) << "moving_subtri_energy" << std::setw(15) << "moving_tri_energy" << std::endl;
	while(xkheight>-.5){
		Point xk_moving(0, xkheight); 
		Triangle tri1_moving(Point(cos(degxi), sin(degxi)), Point(cos(degxj), sin(degxj)), xk_moving); 
	
		double hk= signed_dist_circumcenters(tri1_moving, 2); 
		//std::cout<<"hk: " << hk<<std::endl;

		moving_subtri_energy+=subtri_energy<2,1>(Point(cos(degxi), sin(degxi)), Point(cos(degxj), sin(degxj)), hk);
		moving_tri_energy=tri_energy<2,1>(tri1_moving);
		//energy+=subtri_energy<2,1>(Point(cos(degxi), sin(degxi)), Point(cos(degxj), sin(degxj)), hl);
		std::cout<< std::setw(10) << xkheight<< std::setw(15) << moving_subtri_energy << std::setw(15) << moving_tri_energy << std::endl;
		xkheight-= .1;
		moving_subtri_energy=0;
		moving_tri_energy=0;

	}
*/

	std::ofstream outputFile; 

	Point fixed_pt1(-1,0); 
	Point fixed_pt2(1,0);

	double x_step=.01; 
	double y_step=.01; 

	double y_coor_start=.01; 
	double free_x_coor=-2;
	double free_y_coor=y_coor_start;

/////////////Exp 2a
	
	////////star^0
	outputFile.open("exp2a_vertex_to_fixed_edge_2_0.txt"); 
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
	free_y_coor=y_coor_start;

	
	////////star^1
	outputFile.open("exp2a_vertex_to_fixed_edge_2_1.txt"); 
	while(free_x_coor<2){
		free_y_coor=y_coor_start; 
		while(free_y_coor <2){	
			Point freept(free_x_coor, free_y_coor); 
			Triangle tri(freept, fixed_pt1, fixed_pt2); 
			
			double energy=tri_energy<2,1>(tri); 
			outputFile << std::setw(15) << free_x_coor << std::setw(15) << free_y_coor <<  std::setw(15) << energy <<std::endl; 

			free_y_coor+=y_step; 
		}

	free_x_coor+=x_step;
	}

	outputFile.close(); 


	///////star^2
	free_x_coor=-2;
	free_y_coor=y_coor_start;

	outputFile.open("exp2a_vertex_to_fixed_edge_2_2.txt"); 
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

////////////////Exp 2b HOT/area^2 //////////////////////////////////////////////////////////////////////////////////
	free_x_coor=-2;
	free_y_coor=y_coor_start;
	////////star^0
	outputFile.open("exp2b_vertex_to_fixed_edge_2_0.txt"); 
	while(free_x_coor<2){
		free_y_coor=y_coor_start; 
		while(free_y_coor <2){	
			Point freept(free_x_coor, free_y_coor); 
			Triangle tri(freept, fixed_pt1, fixed_pt2); 
			
			double energy=tri_energy_divideArea<2,0>(tri,2); 
			outputFile << std::setw(15) << free_x_coor << std::setw(15) << free_y_coor <<  std::setw(15) << energy <<std::endl; 

			free_y_coor+=y_step; 
		}

	free_x_coor+=x_step;
	}

	outputFile.close(); 

	free_x_coor=-2;
	free_y_coor=y_coor_start;

	
	////////star^1
	outputFile.open("exp2b_vertex_to_fixed_edge_2_1.txt"); 
	while(free_x_coor<2){
		free_y_coor=y_coor_start; 
		while(free_y_coor <2){	
			Point freept(free_x_coor, free_y_coor); 
			Triangle tri(freept, fixed_pt1, fixed_pt2); 
			
			double energy=tri_energy_divideArea<2,1>(tri,2); 
			outputFile << std::setw(15) << free_x_coor << std::setw(15) << free_y_coor <<  std::setw(15) << energy <<std::endl; 

			free_y_coor+=y_step; 
		}

	free_x_coor+=x_step;
	}

	outputFile.close(); 


	///////star^2
	free_x_coor=-2;
	free_y_coor=y_coor_start;


	outputFile.open("exp2b_vertex_to_fixed_edge_2_2.txt"); 
	while(free_x_coor<2){
		free_y_coor=y_coor_start; 
		while(free_y_coor <2){	
			Point freept(free_x_coor, free_y_coor); 
			Triangle tri(freept, fixed_pt1, fixed_pt2); 
			
			double energy=tri_energy_divideArea<2,2>(tri,2); 
			outputFile << std::setw(15) << free_x_coor << std::setw(15) << free_y_coor <<  std::setw(15) << energy <<std::endl; 

			free_y_coor+=y_step; 
		}

	free_x_coor+=x_step;
	}

	outputFile.close(); 
/////////////////////////// Exp 2c HOT/area
	free_x_coor=-2;
	free_y_coor=y_coor_start;

	////////star^0
	outputFile.open("exp2c_vertex_to_fixed_edge_2_0.txt"); 
	while(free_x_coor<2){
		free_y_coor=y_coor_start; 
		while(free_y_coor <2){	
			Point freept(free_x_coor, free_y_coor); 
			Triangle tri(freept, fixed_pt1, fixed_pt2); 
			
			double energy=tri_energy_divideArea<2,0>(tri,1); 
			outputFile << std::setw(15) << free_x_coor << std::setw(15) << free_y_coor <<  std::setw(15) << energy <<std::endl; 

			free_y_coor+=y_step; 
		}

	free_x_coor+=x_step;
	}

	outputFile.close(); 

	free_x_coor=-2;
	free_y_coor=y_coor_start;


	
	////////star^1
	outputFile.open("exp2c_vertex_to_fixed_edge_2_1.txt"); 
	while(free_x_coor<2){
		free_y_coor=y_coor_start; 
		while(free_y_coor <2){	
			Point freept(free_x_coor, free_y_coor); 
			Triangle tri(freept, fixed_pt1, fixed_pt2); 
			
			double energy=tri_energy_divideArea<2,1>(tri,1); 
			outputFile << std::setw(15) << free_x_coor << std::setw(15) << free_y_coor <<  std::setw(15) << energy <<std::endl; 

			free_y_coor+=y_step; 
		}

	free_x_coor+=x_step;
	}

	outputFile.close(); 


	///////star^2
	free_x_coor=-2;
	free_y_coor=y_coor_start;


	outputFile.open("exp2c_vertex_to_fixed_edge_2_2.txt"); 
	while(free_x_coor<2){
		free_y_coor=y_coor_start; 
		while(free_y_coor <2){	
			Point freept(free_x_coor, free_y_coor); 
			Triangle tri(freept, fixed_pt1, fixed_pt2); 
			
			double energy=tri_energy_divideArea<2,2>(tri,1); 
			outputFile << std::setw(15) << free_x_coor << std::setw(15) << free_y_coor <<  std::setw(15) << energy <<std::endl; 

			free_y_coor+=y_step; 
		}

	free_x_coor+=x_step;
	}

	outputFile.close(); 

//////////////////////////// Exp 1d ////////////////
/////////////////////////// Energy Sb //////////////////////////
	free_x_coor=-2;
	free_y_coor=y_coor_start;


	outputFile.open("exp2d_vertex_to_fixed_edge.txt"); 
	while(free_x_coor<2){
		free_y_coor=y_coor_start; 
		while(free_y_coor <2){	
			Point freept(free_x_coor, free_y_coor); 
			Triangle tri(freept, fixed_pt1, fixed_pt2); 
			double equal_weights[3]={0,0,0}; 
			Point circ=weighted_circumcenter(tri, equal_weights); 
			double energy=triangle_Sb(tri,circ, 2,1); 
			outputFile << std::setw(15) << free_x_coor << std::setw(15) << free_y_coor <<  std::setw(15) << energy <<std::endl; 

			free_y_coor+=y_step; 
		}

	free_x_coor+=x_step;
	}

	outputFile.close();

//////////////////////////// Exp 1e ////////////////
/////////////////////////// Energy Sb/area^2 //////////////////////////
	free_x_coor=-2;
	free_y_coor=y_coor_start;


	outputFile.open("exp2e_vertex_to_fixed_edge.txt"); 
	while(free_x_coor<2){
		free_y_coor=y_coor_start; 
		while(free_y_coor <2){	
			Point freept(free_x_coor, free_y_coor); 
			Triangle tri(freept, fixed_pt1, fixed_pt2); 
			double equal_weights[3]={0,0,0}; 
			Point circ=weighted_circumcenter(tri, equal_weights); 
			double energy=triangle_Sb(tri,circ, 2,-1); 
			outputFile << std::setw(15) << free_x_coor << std::setw(15) << free_y_coor <<  std::setw(15) << energy <<std::endl; 

			free_y_coor+=y_step; 
		}

	free_x_coor+=x_step;
	}

	outputFile.close();


	free_x_coor= 0; 
	free_y_coor= .2; 
	double step=.01;
		while(free_y_coor >.001){	
			Point freept(free_x_coor, free_y_coor); 
			Triangle tri(freept, fixed_pt1, fixed_pt2); 
			double equal_weights[3]={0,0,0}; 
			Point circ=weighted_circumcenter(tri, equal_weights); 
			std::cout << circ.y() <<std::endl;  

			free_y_coor-=step; 
		}
  
  {
    double equal_weights[3]={0,0,0};
    Triangle sample_tri(Point (0,0), Point(1,0), Point(1,.25));
    Point circ=weighted_circumcenter(sample_tri, equal_weights);
    std::cout <<circ.x() << ", " << circ.y() <<std::endl;
  }


//////////////////////////// Exp 1f ////////////////
/////////////////////////// Energy Sb/perim //////////////////////////
	free_x_coor=-2;
	free_y_coor=y_coor_start;


	outputFile.open("exp2f_vertex_to_fixed_edge.txt"); 
	while(free_x_coor<2){
		free_y_coor=y_coor_start; 
		while(free_y_coor <2){	
			Point freept(free_x_coor, free_y_coor); 
			Triangle tri(freept, fixed_pt1, fixed_pt2); 
			double equal_weights[3]={0,0,0}; 
			Point circ=weighted_circumcenter(tri, equal_weights); 
			double energy=triangle_Sb_divide_perim4(tri,circ); 
			outputFile << std::setw(15) << free_x_coor << std::setw(15) << free_y_coor <<  std::setw(15) << energy <<std::endl; 

			free_y_coor+=y_step; 
		}

	free_x_coor+=x_step;
	}

	outputFile.close();



  	return 0;
}
