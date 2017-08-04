#include <hot.hpp>
#include <Sb.hpp>
#include <ply_writer.hpp>

#include <random>
#include <vector>
#include <fstream>

#include <CGAL/Kernel/global_functions.h>

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

//using EK=CGAL::Exact_predicates_exact_constructions_kernel; 

#define PI 3.14159265

int main(int argc, char **argv) {
	std::ofstream outputFile; 

	double fig_height=8;
	double fig_width=2; 
	double step_size=.1; 

	double x_coor=-.99; 
	double y_coor=.01; 
	

	Point points[]={Point(-fig_width/2,0), Point(0,fig_height/2), Point(fig_width/2,0), Point(0,fig_height)};	

/////////////////////////////////////////////////////////////////
///////////////////// Exp 5a: HOT ///////////////////////////////
/////////////////////////////////////////////////////////////
	//star 0
	outputFile.open("exp5/exp5a_horseV_DT_star_2_1.txt"); 
		y_coor=fig_height/2; 
		while(y_coor<fig_height){
			Point freept(0, y_coor); 
		
			DT dt; 
			dt.insert(points[0]);
			dt.insert(points[1]);
			dt.insert(points[2]);
			dt.insert(points[3]); 
			dt.insert(freept);
			std::cout <<"y= " << y_coor <<std::endl;
			std::cout <<"vertices: " <<std::endl;
			for(auto vertex_itr = dt.finite_vertices_begin(); vertex_itr != dt.finite_vertices_end(); vertex_itr++) {
				Point vertex=vertex_itr->point(); 
				std::cout<< "(" <<vertex.x() <<", " << vertex.y()<<")"<<std::endl; 
			}

			double energy=0; 
			for(auto f_itr = dt.finite_faces_begin(); f_itr != dt.finite_faces_end(); f_itr++) {
				
				Point p0=(f_itr->vertex(0))->point(); 
				Point p1=(f_itr->vertex(1))->point(); 
				Point p2=(f_itr->vertex(2))->point(); 
				std::cout<< "("<< p0.x() <<", " << p0.y() <<"), (" <<p1.x() <<", " <<p1.y()<< "), "  <<"(" <<p2.x() <<", " <<p2.y()<< ")"  <<std::endl; 

				energy+= tri_energy<2,1>(face_to_tri(*f_itr)); 
			
			}
			std::cout<<"mesh energy: " << energy <<std::endl; 
			std::cout <<std::endl; 
			outputFile<< std::setw(15) << y_coor << std::setw(15) << energy <<std::endl; 
			

			
		y_coor+=step_size; 
		}
	
	outputFile.close(); 
//varry horizontal coordinate
		outputFile.open("exp5/exp5a_horseV_DT_star_2_1_hor.txt"); 
		step_size=.01;
		y_coor=6;
		x_coor= -.25; 
		while(x_coor<.25){
			Point freept(x_coor,y_coor); 
		
			DT dt; 
			dt.insert(points[0]);
			dt.insert(points[1]);
			dt.insert(points[2]);
			dt.insert(points[3]); 
			dt.insert(freept);
			//std::cout <<"y= " << y_coor <<std::endl;
			std::cout <<"vertices: " <<std::endl;
			for(auto vertex_itr = dt.finite_vertices_begin(); vertex_itr != dt.finite_vertices_end(); vertex_itr++) {
				Point vertex=vertex_itr->point(); 
				std::cout<< "(" <<vertex.x() <<", " << vertex.y()<<")"<<std::endl; 
			}

			double energy=0; 
			for(auto f_itr = dt.finite_faces_begin(); f_itr != dt.finite_faces_end(); f_itr++) {
				
				Point p0=(f_itr->vertex(0))->point(); 
				Point p1=(f_itr->vertex(1))->point(); 
				Point p2=(f_itr->vertex(2))->point(); 
				//std::cout<< "("<< p0.x() <<", " << p0.y() <<"), (" <<p1.x() <<", " <<p1.y()<< "), "  <<"(" <<p2.x() <<", " <<p2.y()<< ")"  <<std::endl; 

				energy+= tri_energy<2,1>(face_to_tri(*f_itr)); 
			
			}
			//std::cout<<"mesh energy: " << energy <<std::endl; 
			//std::cout <<std::endl; 
			outputFile<< std::setw(15) << x_coor << std::setw(15) << energy <<std::endl; 
			

			
		x_coor+=step_size; 
		}
	
	outputFile.close(); 

// as vertex is moved, update to DT
	outputFile.open("exp5/exp5a_horseV_DTmesh_star_2_1.txt"); 
	Point boundary_pts[]={Point(-fig_width/2,0),Point(fig_width/2,0), Point(0,fig_height)};	
	x_coor= -fig_width/2+.01;
	y_coor= .01;
	while(x_coor< fig_width/2-.01){
		y_coor=.01;
		while(y_coor<fig_height){
			Point freept(x_coor, y_coor); 
			/*if(point_outside_domain(Point(x_coor,y_coor), points, points+4, K())){
				y_coor+=step_size;
				continue; 
			}
			*/
			
			if(point_outside_domain(Point(x_coor,y_coor), boundary_pts, boundary_pts+3, K())){
				y_coor+=step_size;
				continue; 
			}

			DT dt; 
			dt.insert(freept); 
			dt.insert(points[0]);
			dt.insert(points[1]);
			dt.insert(points[2]);
			dt.insert(points[3]); 
			
			double total_energy=0;
			for(auto face_itr=dt.finite_faces_begin(); face_itr!=dt.finite_faces_end();face_itr++){
				total_energy+=tri_energy<2,1>(face_to_tri(*face_itr));
			}
			outputFile<< std::setw(15) <<std::setprecision(5)<<freept.x() << std::setw(15) <<std::setprecision(5)<< freept.y() << std::setw(15) << std::setprecision(5) << total_energy <<std::endl;
		y_coor+=step_size; 
		}
		x_coor+=step_size;
	}
	outputFile.close(); 


	return 0; 
}
