
#include <hot.hpp>
#include <ply_writer.hpp>

#include <random>
#include <vector>

#include <CGAL/Kernel/global_functions.h>

void test_tri_w2() {
  Triangle tri(Point(2, 1), Point(3, 1), Point(2, 2));
  boost::variant<std::array<Triangle, 2>,
                 std::array<Triangle, 1> >
      bounds = integral_bounds(tri);
  assert(bounds.which() == 1);
  auto area = boost::get<std::array<Triangle, 1> >(bounds);

  std::cout
      << std::endl
      << "Computing Wasserstein distance to the dual point"
      << std::endl
      << std::endl;
  K_real value = triangle_w<2>(tri);
  std::cout << "Unit Right Triangle expected Wasserstein "
               "Distance: 5/90 (~0.0555555). Calculated "
               "Distance: "
            << value << std::endl
            << std::endl;
}

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

int main(int argc, char **argv) {

  // Test hot energy between triangle and triangle* == point

  DT dt;
std::cout<< "star0-Hot_2,2 Mesh energy density for random DT: " << std::endl;

  const int num_points = 40;

	for(int numiterations=0; numiterations<10; numiterations++){
  		generate_rand_dt(num_points, dt);
		K_real energy=hot_energy_density<2,1>(dt);
		std::cout << std::setw(5) << numiterations <<": " <<  energy << std::endl;
	}


// Make some non-DT triangles, compare to DT

	// Set-up 1 
	const int Wk=2;
	const int star=1; 
	std::cout<< "Set-up 1 "<< std::endl;
	std:: cout<< " Results for Wk=" << Wk << ", star=" <<star << std::endl; 
        std::cout << std:: setw(6) << "" << std::setw(10) << "Tri 1" << std::setw(10) << "Tri 2" << std::setw(10) << "sum" << std::endl; 
	for(double height=1.0; height <10; height++){
		Triangle DTtri1(Point(-1,0), Point(1,0), Point(0,height));
		Triangle DTtri2(Point(-1,0), Point(1,0), Point(0,-height));
		
		K_real DTtri1_energy=tri_energy<Wk,star>(DTtri1);
		K_real DTtri2_energy=tri_energy<Wk,star>(DTtri2);
		K_real DT_energy=DTtri1_energy+DTtri2_energy; 
		
	
	
		Triangle NDTtri1(Point(-1,0), Point(0,height), Point(0,-height));
		Triangle NDTtri2(Point(1,0), Point(0,height), Point(0,-height));

		K_real NDTtri1_energy=tri_energy<Wk,star>(NDTtri1); 
		K_real NDTtri2_energy=tri_energy<Wk,star>(NDTtri2); 
		K_real NDT_energy=NDTtri1_energy + NDTtri2_energy; 


		std:: cout <<std::setw(6) << "DT: "<<  std::setw(10) << DTtri1_energy << std::setw(10) << DTtri2_energy << std::setw(10) << DT_energy << std::endl;
		std:: cout << std::setw(6) << "NDT: " << std::setw(10) << NDTtri1_energy << std::setw(10) << NDTtri2_energy << std::setw(10) << NDT_energy << std::endl<<std::endl;
		
	
	}



// Set-up 2 
	std::cout << std::endl << "Set-up 2: " << std::endl; 
	std:: cout<< " Results for Wk=" << Wk << ", star=" <<star << std::endl; 
	std::cout << std:: setw(6) << "" << std::setw(10) << "Tri 1" << std::setw(10) << "Tri 2" << std::setw(10) << "sum" << std::endl; 
	for(double t=0.01; t >.001; t-=0.001){
		Triangle DTtri1(Point(0,0), Point(t,1), Point(2,0));
		Triangle DTtri2(Point(0,0), Point(t,-1), Point(2,0));
		
		K_real DTtri1_energy=tri_energy<Wk,star>(DTtri1);
		K_real DTtri2_energy=tri_energy<Wk,star>(DTtri2);
		K_real DT_energy=DTtri1_energy+DTtri2_energy; 
	
		Triangle NDTtri1(Point(0,0), Point(t,1), Point(t,-1));
		Triangle NDTtri2(Point(t,1), Point(t,-1), Point(2,0));

		K_real NDTtri1_energy=tri_energy<Wk,star>(NDTtri1); 
		K_real NDTtri2_energy=tri_energy<Wk,star>(NDTtri2); 
		K_real NDT_energy=NDTtri1_energy+NDTtri2_energy; 
		
		//std::cout<< "DTtri1-energy: " << DTtri1_energy << " " << "NDTtri1-energy: " <<NDTtri1_energy <<std::endl;
		//std::cout<< "DTtri2-energy: " << DTtri2_energy << " " << "NDTtri2-energy: " <<NDTtri2_energy <<std::endl;
		//std::cout<< "DT-energy, NDT-energy: " << DT_energy << " , " <<NDT_energy <<std::endl;
		
	//	std::cout<< "t= " << t << std::endl;
	//	std::cout<< "NDTtri1 energy: " << NDTtri1_energy << std::endl;
		std:: cout <<std::setw(6) << "DT: "<<  std::setw(10) << DTtri1_energy << std::setw(10) << DTtri2_energy << std::setw(10) << DT_energy << std::endl;
		std:: cout << std::setw(6) << "NDT: " << std::setw(10) << NDTtri1_energy << std::setw(10) << NDTtri2_energy << std::setw(10) << NDT_energy << std::endl<<std::endl;
	}

	//Point samplepoint1(1.2, 1.4), samplepoint2(1.5, 43.4); 
	//std::cout<< "This is the sample subtriangle energy: " << subtri_energy<2,0>(samplepoint1, samplepoint2 , 5.0) << std::endl; 


	/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

	std::cout<<std::endl<< "hexagon with free vertex experiment: "<< std::endl;
	

	// hexgon perimeter points
	
	std::vector<Point> points={Point(-1,0),  Point(cos(120), sin(120)),  Point(cos(60), sin(60)),  Point(1,0),  Point(cos(300), sin(300)),  Point(cos(240), sin(240))} ; 





	//free point
	Point freept(0.0,0.0); 

	////Create triangles
	//std::vector<Triangle> triangle_array={Triangle(points[0], points[1], freept)}; 


	std::vector<Triangle> triangle_array;

	for(int i=0; i<6; i++){
		triangle_array.push_back(Triangle(points[i], points[(i+1)%6], freept));
	}

	// print triangles vertices 
	std::vector<Triangle>::iterator iter=triangle_array.begin(); // pointer is a pointer to the beginning of the array
	while(iter!=triangle_array.end()){
		std::cout<< *iter << std::endl;
		//std::cout<< *iter.vertex(0) << std::setw(10) << *iter.vertex(1) << std::setw(10) << *iter.vertex(2) << std::endl;
		iter ++;
	}

	double total_hex_energy=0;
	for(int i=0; i<triangle_array.size(); i++){
		total_hex_energy+=tri_energy<Wk,star>(triangle_array[i]);
	}
	std::cout<< "This is the total hexagon energy: " << total_hex_energy << std::endl;

  	write_ply("test.ply", dt);
  	return 0;
}



// std::cout << "Generated the Delaunay Triangulation of "
  //          << num_points << " points" << std::endl;
  //std::cout << "Resulting in " << dt.number_of_faces()
  //          << " faces" << std::endl;
  //int face_idx = 1;
  //for(auto face_itr = dt.finite_faces_begin();
    //  face_itr != dt.finite_faces_end();
     // face_itr++, face_idx++) {
	//std::cout <<"This the circumcetner of the face:" << CGAL::circumcenter(face) << std::endl; 
    
    //Face face = *face_itr;
    //Triangle tri=face_to_tri(face); 
    //std::cout<< "This is the star^0 energy contribution of the current triangle: " << tri_energy<2,0>(tri) << std::endl; 
    //std::cout << "Face " << face_idx << " : ";
    //for(int i = 0; i < tri_verts; i++) {
     // std::cout << " ( " << face.vertex(i)->point()
     //           << " )  ";
   // }
    //std::cout << std::endl;
  //}
 // K_real energy = hot_energy<2>(dt);
  //std::cout << "Mesh energy: " << energy << std::endl;
