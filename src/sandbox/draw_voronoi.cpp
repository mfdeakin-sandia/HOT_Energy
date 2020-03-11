// draw_voronoi.cpp

#include <fstream>
#include <iterator>

#include "hot.hpp"

//A class to recover Voronoi diagram from stream.
//Rays, lines and segments are cropped to a rectangle
//so that only segments are stored
struct Cropped_voronoi_from_delaunay{
  std::list<Segment_2> m_cropped_vd;
  Iso_rectangle_2 m_bbox;
  
  Cropped_voronoi_from_delaunay(const Iso_rectangle_2& bbox):m_bbox(bbox){}
  
  template <class RSL>
  void crop_and_extract_segment(const RSL& rsl){
    CGAL::Object obj = CGAL::intersection(rsl,m_bbox);
    const Segment_2* s=CGAL::object_cast<Segment_2>(&obj);
    if (s) m_cropped_vd.push_back(*s);
  }
  
  void operator<<(const Ray_2& ray)    { crop_and_extract_segment(ray); }
  void operator<<(const Line_2& line)  { crop_and_extract_segment(line); }
  void operator<<(const Segment_2& seg){ crop_and_extract_segment(seg); }
};


int main(){
  //consider some points
  //std::vector<Point_2> points;
	std::vector<Weightpt> wpoints;
  //wpoints.push_back(Weightpt(Point_2(0,0),3));
  //wpoints.push_back(Weightpt(Point_2(1,1),0));
  //wpoints.push_back(Weightpt(Point_2(0,1),1));



  wpoints.push_back(Weightpt(Point_2(0,.5),0)); 

  wpoints.push_back(Weightpt(Point_2(1,3),0)); 

  wpoints.push_back(Weightpt(Point_2(2,4.2),0)); 
//for initial mesh
 // wpoints.push_back(Weightpt(Point_2(4,1.2),0)); 

//for adjusted mesh
  wpoints.push_back(Weightpt(Point_2(3,3.9),0)); 
  
  wpoints.push_back(Weightpt(Point_2(4,4.2),0)); 

  wpoints.push_back(Weightpt(Point_2(1.5,6),0)); 

  wpoints.push_back(Weightpt(Point_2(6,1),0)); 

  wpoints.push_back(Weightpt(Point_2(6,5) ,0)); 


  
  //Delaunay_triangulation_2 dt2;
	Regular_triangulation_2 dt2; 
  //insert points into the triangulation
  dt2.insert(wpoints.begin(),wpoints.end());
  //construct a rectangle
  //Iso_rectangle_2 bbox(-1,-1,2,2);
  Iso_rectangle_2 bbox(-1,-1,8,8);
  Cropped_voronoi_from_delaunay vor(bbox);
  //extract the cropped Voronoi diagram
  dt2.draw_dual(vor);
  //print the cropped Voronoi diagram as segments
 // std::copy(vor.m_cropped_vd.begin(),vor.m_cropped_vd.end(),
   // std::ostream_iterator<Segment_2>(std::cout,"\n"));

//format print out for tikz
	std::cout<< "Voronoi edges: " <<std::endl; 
	for(auto seg_itr=vor.m_cropped_vd.begin(); seg_itr!=vor.m_cropped_vd.end(); seg_itr++){
		std::cout<< " \\draw (" <<(seg_itr->source()).x() <<", " << (seg_itr->source()).y() << ")" << " -- " <<  "(" <<(seg_itr->target()).x() <<", " << (seg_itr->target()).y() << ");"<<std::endl;
	}

std::cout <<"DT edges" <<std::endl;

	for(auto f_itr=dt2.finite_faces_begin(); f_itr!=dt2.finite_faces_end(); f_itr ++){
		std::cout <<" \\draw (" << ((f_itr->vertex(0))->point()).x() <<" , " << ((f_itr->vertex(0))->point()).y() << " ) -- ("<< ((f_itr->vertex(1))->point()).x() <<" , " << ((f_itr->vertex(1))->point()).y() << " ) -- (" <<((f_itr->vertex(2))->point()).x() <<" , " << ((f_itr->vertex(2))->point()).y() << " ) -- (" <<((f_itr->vertex(0))->point()).x() <<" , " << ((f_itr->vertex(0))->point()).y() << " );" << std::endl;
}

	std::vector<Weightpt> wpoints2;




  wpoints2.push_back(Weightpt(Point_2(0,4),0)); 
  wpoints2.push_back(Weightpt(Point_2(5,5),0)); 
wpoints2.push_back(Weightpt(Point_2(9,2),0)); 
wpoints2.push_back(Weightpt(Point_2(7,1),0)); 
wpoints2.push_back(Weightpt(Point_2(3,1),0)); 
wpoints2.push_back(Weightpt(Point_2(1,2),0)); 

//wpoints2.push_back(Weightpt(Point_2(5,1),0)); 
wpoints2.push_back(Weightpt(Point_2(5,1.5),0)); 

 

  Regular_triangulation_2 rt; 
  //insert points into the triangulation
  rt.insert(wpoints2.begin(),wpoints2.end());
  //construct a rectangle

  Iso_rectangle_2 bboxrt(0,0,10,6);
  Cropped_voronoi_from_delaunay vorrt(bboxrt);
  //extract the cropped Voronoi diagram
  rt.draw_dual(vor);

//format print out for tikz
	std::cout<< "Voronoi edges: " <<std::endl; 
	for(auto seg_itr=vorrt.m_cropped_vd.begin(); seg_itr!=vorrt.m_cropped_vd.end(); seg_itr++){
		std::cout<< " \\draw (" <<(seg_itr->source()).x() <<", " << (seg_itr->source()).y() << ")" << " -- " <<  "(" <<(seg_itr->target()).x() <<", " << (seg_itr->target()).y() << ");"<<std::endl;
	}

std::cout <<"RT edges" <<std::endl;

	for(auto f_itr=rt.finite_faces_begin(); f_itr!=rt.finite_faces_end(); f_itr ++){
		std::cout <<" \\draw (" << ((f_itr->vertex(0))->point()).x() <<" , " << ((f_itr->vertex(0))->point()).y() << " ) -- ("<< ((f_itr->vertex(1))->point()).x() <<" , " << ((f_itr->vertex(1))->point()).y() << " ) -- (" <<((f_itr->vertex(2))->point()).x() <<" , " << ((f_itr->vertex(2))->point()).y() << " ) -- (" <<((f_itr->vertex(0))->point()).x() <<" , " << ((f_itr->vertex(0))->point()).y() << " );" << std::endl;
	}

	std::cout <<"Experiment 1" <<std::endl; 
	double freeptx=5; 
	double freepty=4; 
	double energy_U=0; 
	double energy_NU=0; 
	std::ofstream outputFileU;
	std::ofstream outputFileNU;  
	outputFileU.open("patch_energy_barrier/exp1U_patch_energy_barrier.txt"); 
	outputFileNU.open("patch_energy_barrier/exp1NU_patch_energy_barrier.txt"); 
	while(freepty>0){
		energy_U=0; 
		energy_NU=0;
		std::vector<Weightpt> wpoints3;

		//wpoints3.push_back(Weightpt(Point_2(0,4),0)); 
		wpoints3.push_back(Weightpt(Point_2(4,5),0)); 
		wpoints3.push_back(Weightpt(Point_2(9,2),0)); 
		wpoints3.push_back(Weightpt(Point_2(7,1),0)); 
		wpoints3.push_back(Weightpt(Point_2(3,1),0)); 
		//wpoints3.push_back(Weightpt(Point_2(1,2),0)); 
		
		
		for(int i=0; i <4; i++){
			Triangle tri(Point_2(freeptx, freepty), static_cast<Point>(wpoints3[i].point()), static_cast<Point>(wpoints3[(i+1)%4].point()));
			energy_NU+= tri_energy<2,1>(tri); 
		}
		

		Triangle tri1(Point(3,1), Point(7,1), Point(4,-6));
		Triangle tri2(Point(7,1), Point(4,-6), Point(9,2)); 
		energy_NU+=tri_energy<2,1>(tri1); 
		energy_NU+=tri_energy<2,1>(tri2); 

		outputFileNU <<std::setw(10) <<freepty <<std::setw(10) <<energy_NU <<std::endl; 		

		wpoints3.push_back(Weightpt(Point_2(4,-4), 0)); 
		wpoints3.push_back(Weightpt(Point_2(freeptx,freepty),0));
			
		Regular_triangulation_2 rt3; 
		rt3.insert(wpoints3.begin(),wpoints3.end());
		//std::cout <<"RT edges" <<std::endl;
		std::cout << "\\begin{tikzpicture}" <<std::endl; 
		for(auto f_itr=rt3.finite_faces_begin(); f_itr!=rt3.finite_faces_end(); f_itr ++){
	//		std::cout <<" \\draw (" << ((f_itr->vertex(0))->point()).x() <<" , " << ((f_itr->vertex(0))->point()).y() << " ) -- ("<< ((f_itr->vertex(1))->point()).x() <<" , " << ((f_itr->vertex(1))->point()).y() << " ) -- (" <<((f_itr->vertex(2))->point()).x() <<" , " << ((f_itr->vertex(2))->point()).y() << " ) -- (" <<((f_itr->vertex(0))->point()).x() <<" , " << ((f_itr->vertex(0))->point()).y() << " );" << std::endl;
	
			Triangle tri(static_cast<Point>((*f_itr).vertex(0)->point()), static_cast<Point>((*f_itr).vertex(1)->point()), static_cast<Point>((*f_itr).vertex(2)->point())); 
			energy_U+=tri_energy<2,1>(tri); 
		}
		std::cout <<"\\end{tikzpicture}" <<std::endl<<std::endl; 
			
		outputFileU << std::setw(10) << freepty << std::setw(10) << energy_U <<std::endl; 

		freepty-=.5; 	
		
	}
	outputFileU.close();
	outputFileNU.close(); 

}

