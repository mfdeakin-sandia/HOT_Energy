#include <hot.hpp>
#include <ply_writer.hpp>

#include <random>
#include <vector>
#include <fstream>

#include <CGAL/Kernel/global_functions.h>


#define CGAL_MESH_2_OPTIMIZER_VERBOSE
//#define CGAL_MESH_2_OPTIMIZERS_DEBUG
//#define CGAL_MESH_2_SIZING_FIELD_USE_BARYCENTRIC_COORDINATES
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/lloyd_optimize_mesh_2.h>
#include <iostream>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_mesh_vertex_base_2<K>                Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K>                  Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>  CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>            Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria>              Mesher;
//typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;


/*
using RNG = std::mt19937_64;
constexpr const float min_pos = 10.0;
constexpr const float max_pos = 20.0;


void generate_rand_dt(int num_points, DT &dt) {
  std::random_device rd;
  RNG engine(rd());
  std::uniform_real_distribution<double> genPos(min_pos,
                                                max_pos);
  for(int i = 0; i < num_points; i++) {
    Point pt(genPos(engine), genPos(engine));
    dt.insert(pt);
  }
}

void generate_rand_cdt(int num_points, CDT &cdt) {
  std::random_device rd;
  RNG engine(rd());
  std::uniform_real_distribution<double> genPos(min_pos,
                                                max_pos);
  for(int i = 0; i < num_points; i++) {
    Point pt(genPos(engine), genPos(engine));
    cdt.insert(pt);
  }
}
*/

int main(int argc, char **argv) {
	DT dt; 
	  CDT cdt;

	const int num_points = 100;
	
  	//generate_rand_dt(num_points, dt);
	generate_rand_t(num_points,cdt);
	
	// insert boundary vertices
	CDT::Vertex_handle va=cdt.insert(Point(10,10));
	CDT::Vertex_handle vb=cdt.insert(Point(20,10));
	CDT::Vertex_handle vc=cdt.insert(Point(20,20));
	CDT::Vertex_handle vd=cdt.insert(Point(10,20));

	cdt.insert(Point(10.5,11.8));
	cdt.insert(Point(17.5,10.9));
	cdt.insert(Point(16.3,14.3));
	cdt.insert(Point(11,12.8));

	cdt.insert_constraint(va,vb);
	cdt.insert_constraint(vb,vc);
	cdt.insert_constraint(vc,vd);
	cdt.insert_constraint(vd,va);	
	
	int num_constrained_edges=0;

	for(auto ei=cdt.finite_edges_begin();ei!=cdt.finite_edges_end(); ei++){
		if(cdt.is_constrained(*ei))  num_constrained_edges++;
	}
	//Observations: when inserting vertices, no constraints are added. it is not assume that the convex hull boundary edges are constrained. 
	//For now we manually insert a boundary box of constrained edges. 
		

	
	std::cout<<"Num constrained edges: " << num_constrained_edges << std::endl; 



  	//CDT::Vertex_handle va = cdt.insert(Point(-2,0));
  	//CDT::Vertex_handle vb = cdt.insert(Point(0,-2));
  	//CDT::Vertex_handle vc = cdt.insert(Point(2,0));
  	//CDT::Vertex_handle vd = cdt.insert(Point(0,1));
  	//cdt.insert(Point(2, 0.6));
  	//cdt.insert_constraint(va, vb);
  	//cdt.insert_constraint(vb, vc);
  //cdt.insert_constraint(vc, vd);
  //cdt.insert_constraint(vd, va);
  //std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
  //std::cout << "Number of faces: " << cdt.number_of_faces() << std::endl;
  std::cout << "Meshing..." << std::endl;
  Mesher mesher(cdt);
  mesher.set_criteria(Criteria(0.125, 0.000005));
  mesher.refine_mesh();
  
  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
  std::cout <<"Number of faces: " <<cdt.number_of_faces()<<std::endl;


int face_num=0; 
	for(auto face_itr = cdt.finite_faces_begin(); face_itr != cdt.finite_faces_end(); face_itr++, face_num++){
		//Triangle face = CDT::face_to_tri(*face_itr);
		//CDT::Face tri=*face_itr;
		Point point0=((*face_itr).vertex(0))->point(); 
		Point point1=((*face_itr).vertex(1))->point();
		Point point2=((*face_itr).vertex(2))->point(); 
	 
		//std::cout<<face_num <<std::endl; 
		//std::cout<<"("<<point0.x() << ", " <<point0.y() <<")" << std::endl; 
		//std::cout<<"("<<point1.x() << ", " <<point1.y() <<")" << std::endl; 
		//std::cout <<"("<<point2.x() << ", " <<point2.y() <<")" << std::endl<<std::endl; 
	}

	std::vector<double> original_x_coor; 
	std::vector<double> original_y_coor;
	for(auto vertex_itr=cdt.finite_vertices_begin(); vertex_itr!=cdt.finite_vertices_end(); vertex_itr++){
		original_x_coor.push_back((vertex_itr->point()).x());
		original_y_coor.push_back((vertex_itr->point()).y());
	}
	
	
	std::cout << "Run Lloyd optimization for CDT...";
	CGAL::lloyd_optimize_mesh_2(cdt, CGAL::parameters::freeze_bound=0,CGAL::parameters::convergence=.000001);
	//  CGAL::parameters::max_iteration_number = 100
	std::cout << " done." << std::endl;
	std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
	std::cout <<"Number of faces: " <<cdt.number_of_faces()<<std::endl;

	int num_unmoved_vertices=0;
	double error=.000000000000000000000000001;
	for(auto vertex_itr=cdt.finite_vertices_begin(); vertex_itr!=cdt.finite_vertices_end(); vertex_itr++){
		for(int i=0; i <original_x_coor.size(); i++){
			//if(((vertex_itr->point()).x()) ==original_x_coor[i] &&((vertex_itr-> point()).y()) ==original_y_coor[i]){
			if(abs((vertex_itr->point()).x() -original_x_coor[i])<error &&abs((vertex_itr-> point()).y() -original_y_coor[i])<error){
				original_x_coor.erase(original_x_coor.begin()+i); 
				original_y_coor.erase(original_y_coor.begin()+i); 
				num_unmoved_vertices++;
				//break;
			}
		
		}

		//std::cout<<"("<< (vertex_itr->point()).x() << ", "<< (vertex_itr->point()).y() << ") " << " Size original_x_coor: " << original_x_coor.size() << std::endl; 

	}
	std::cout<< "Num unmoved vertices: " << num_unmoved_vertices <<std::endl; 
	face_num=0; 
	for(auto face_itr = cdt.finite_faces_begin(); face_itr != cdt.finite_faces_end(); face_itr++, face_num++){
		//Triangle face = CDT::face_to_tri(*face_itr);
		//CDT::Face tri=*face_itr;
		Point point0=((*face_itr).vertex(0))->point(); 
		Point point1=((*face_itr).vertex(1))->point();
		Point point2=((*face_itr).vertex(2))->point(); 
	 
		//std::cout<<face_num <<std::endl; 
		//std::cout<<"("<<point0.x() << ", " <<point0.y() <<")" << std::endl; 
		//std::cout<<"("<<point1.x() << ", " <<point1.y() <<")" << std::endl; 
		//std::cout <<"("<<point2.x() << ", " <<point2.y() <<")" << std::endl<<std::endl; 
	}

	//DT dt; 
	//DT::Vertex_handle va2 = dt.insert(Point(-2,0));
  	//DT::Vertex_handle vb2 = dt.insert(Point(0,-2));
  	//DT::Vertex_handle vc2 = dt.insert(Point(2,0));
  	//DT::Vertex_handle vd2 = dt.insert(Point(0,1));
  	//dt.insert(Point(2, 0.6));
	
//CGAL::lloyd_optimize_mesh_2(dt); 
//	 std::cout << "Number of vertices: " << dt.number_of_vertices() << std::endl;
  //std::cout << "Run Lloyd optimization for DT...";
  //CGAL::lloyd_optimize_mesh_2(cdt,
    //CGAL::parameters::max_iteration_number = 10);
  //std::cout << " done." << std::endl;
  //std::cout << "Number of vertices: " << dt.number_of_vertices() << std::endl;
 	return 0;
}
