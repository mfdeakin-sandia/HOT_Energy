#include <hot.hpp>
#include <Sb.hpp>
#include <energyWeights.hpp>
#include <lloyds.hpp>
#include <ply_writer.hpp>
//#include <build_triangulation.hpp>

#include <random>
#include <vector>
#include <fstream>

#include <CGAL/Kernel/global_functions.h>


#define PI 3.14159265



int main(int argc, char **argv) {
	//Check center of mass is working 
	
	std::vector<Point> points={Point(0,0), Point(2,0), Point(2,2), Point(0,2)}; 

	Point center=center_mass_polygon(points);
	std::cout<<"(" << center.x() <<", " << center.y() <<")" << std::endl;


	//check out build_triangulation
	//std::vector <Point> intial_points={Point(0,0), Point(2,0), Point(2,2), Point(0,2), Point(1,1)};

	std::vector <Point> intial_points={Point(0,0), Point(2,0), Point(2,2), Point(0,2), Point(1,1),Point(4,0), Point(3,1), Point(4,2)};
	int num_lloyds_iterations=5;  
	double x_min=-1, y_min=-1, y_max=3, x_max=5 ; 
	RegT rt=build_reg_triangulation(intial_points, 10, x_min,y_min,y_max,x_max);
	std::cout<<"Vertices in final trinagulation after Lloyds iterations : " <<std::endl; 
	for(auto v_itr=rt.finite_vertices_begin(); v_itr!=rt.finite_vertices_end(); v_itr++){
		std::cout << "(" << v_itr->point().x()<<", "  << v_itr->point().y() <<")" <<std::endl; 
	}
// because the center of mass of the Voronoi regions of our intial set of points were the initial sites, running Lloyds does not move the points! 



// check Sb energy works 
	RegT random_RT; 
	generate_rand_t(10,random_RT);
	std::cout <<"Sb energy 2,1: " << Sb(random_RT,2,1)<<std::endl; 
	std::cout <<"Sb energy: 2,5: " << Sb(random_RT,2,5)<<std::endl; 
	std::cout <<"*2-HOT_{2,2}/AREA: " << energy_weights_dividebyArea(random_RT,2,2) <<std::endl; 
	std::cout << "*2-HOT_{2,2}:" << energy_weights(random_RT,2,2) <<std::endl; 

// check our weighted_circumcenter function works
	// 1:
	RegT::Face_handle face_handle=random_RT.finite_faces_begin(); 
	Triangle tri((face_handle->vertex(0))->point(), (face_handle->vertex(1))->point(), (face_handle->vertex(2))->point()); 

	double weight[3]={0,0,0};

	std::cout<< "Our code: " << weighted_circumcenter(tri,weight) <<std::endl; 
	std::cout <<"CGAL: " << random_RT.weighted_circumcenter(face_handle) <<std::endl; 

	std::vector<Wpt> weighted_points; 
	for(int i=0; i <20; i++){
		for(int j=0; j <20; j++){
			weighted_points.push_back(Wpt(Point(i,j), 0)); 
		}
	}

	// 2: 
	std::cout<<std::endl; 
	RegT uniform_rt; 
	for(auto wpt: weighted_points){
		uniform_rt.insert(wpt);
	}
	//uniform_rt.insert(weighted_points.begin(), weighted_points.end()); 
	std::cout << uniform_rt.number_of_vertices() <<std::endl; 
	RegT::Face_handle face=uniform_rt.finite_faces_begin(); 
	Triangle tri2((face->vertex(0))->point(), (face->vertex(1))->point(), (face->vertex(2))->point()); 
	
	std::cout<< "Our code: " << weighted_circumcenter(tri2,weight) <<std::endl; 
	std::cout <<"CGAL: " << random_RT.weighted_circumcenter(face) <<std::endl; 

	std::cout<<std::endl; 
	// 3: 	
	Triangle tri3(Point(cos(.75), sin(.75)), Point(cos(1.5), sin(1.5)), Point(cos(2.9), sin(2.9)));
	double tri3_weights[3]={0,0,0};  
	std::cout<< "Our code: " << weighted_circumcenter(tri3, tri3_weights)<<std::endl; 
	std::cout <<"True answer: (0,0)"<<std::endl ; 

// test perimeter function. Checked with wolframalpha 

	std::cout <<"Perimeter tri3: " << perimeter(tri3) <<std::endl; 

double zer_weights[]={0,0,0};

double weights[]={1,3,7};

double xk2=.5;
while(xk2>.01){
	Triangle tri(Point(0,0), Point(0,6), Point(6.578700842374038, xk2));
	Point wcirc=weighted_circumcenter(tri, weights);
	std::cout<< "(" << wcirc.x() <<", " <<wcirc.y() <<")" <<std::endl; 
	xk2-=.01; 
}

xk2=.2;
while(xk2>.001){
	Triangle tri(Point(0,0), Point(0,6), Point(-0.9120341757073718, xk2));
	Point wcirc=weighted_circumcenter(tri, weights);
	std::cout<< "(" << wcirc.x() <<", " <<wcirc.y() <<")" <<std::endl; 
	xk2-=.001; 
}
	
double ex_weights[]={0,0,1.5};

Triangle exampletri(Point(-sqrt(1 - pow((3/20),2)), -3/20), Point(sqrt(1 - pow((3/20),2)), -3/20), Point(0,1)); 

Point ex_tri_wcirc=weighted_circumcenter(exampletri, ex_weights); 
std::cout<< "This is the weighted circ of the example tri:" <<std::endl; 
std::cout<< "(" << ex_tri_wcirc.x() <<", " <<ex_tri_wcirc.y() <<")" <<std::endl; 
	
return 0; 
}
