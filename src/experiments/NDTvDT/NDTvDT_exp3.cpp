// NDTvDT_exp3.cpp

#include <fstream>

#include "hot.hpp"

double compute_cot(Triangle tri, int vertex_index); 

int main(int argc, char **argv) {
	double x_start=-2; 
	double y_start=.01; 
	double length_start=.1; 
	double angle_start=0; 


	double x=x_start;
	double y=y_start; 
	double x_step=.1; 
	double y_step=.5; 
	double angle=angle_start; 
	double angle_step=.1; 
	double length=length_start; 
	double length_step=.1;

	double total_experiments=0; 
	int num_DTgreaterNDT=0; 
	double num_exp_2_faces=0;
	double DTenergy; 
	double NDTenergy; 

	double cot1=0; 
	double cot2=0; 

	Point fixedpt1(0,0); 
	Point fixedpt2(1,0); 
	
	while(x<2){
		y=y_start; 
		length=length_start;
		angle=angle_start; 
		while(y<4){
			length=length_start; 	
			angle=angle_start; 	
			while(length <2){
				angle=angle_start; 
				while(angle <3.14/2){
					total_experiments++; 
				
					DT dt; 
					dt.insert(fixedpt1);
					dt.insert(fixedpt2);
					dt.insert(Point(x,y)); 
					dt.insert(Point(x+length*cos(angle), y+length*sin(angle))); 
					
					int num_faces=0; 
					for(auto f_itr=dt.finite_faces_begin(); f_itr!=dt.finite_faces_end(); f_itr++){
							num_faces++; 
					}

					if(num_faces!=2){
						angle+=angle_step;
						continue;
					} 

					num_exp_2_faces++;

					double DTenergy=energy_density_TMethod<2,1>(dt);
					double NDTenergy=0; 
					for(auto ei=dt.finite_edges_begin();ei!=dt.finite_edges_end(); ei++)
          {
						// check if ei is a boundary edge
						if(!(dt.is_infinite(ei->first) || dt.is_infinite((ei->first)->neighbor(ei->second))))
            {
              // get both triangles of e1
							Face face1=*(ei->first);
							int index1=ei->second; 
							
							Edge mirror_ei=dt.mirror_edge(*ei);
							Face face2=*(mirror_ei.first);
							int index2=mirror_ei.second;

              // order points
              const int index1p1 = (index1 + 1) % 3;
              const int index1p2 = (index1 + 2) % 3;
              
              auto &p1 = face1.vertex(index1)->point();
              auto &p2 = face1.vertex(index1p1)->point();
              auto &p3 = face1.vertex(index1p2)->point();
              auto &p4 = face2.vertex(index2)->point();
              
              // DT is tri(p1,p2,p3) and tri(p4,p3,p2)
              // p2,p3 is the shared edge
              
              // create triangles with the diagonal flipped
              // tri 124, tri 143
              Triangle tri1(p1, p2, p4);
              Triangle tri2(p1, p4, p3);
              
							NDTenergy=tri_energy<2,1>(tri1)+tri_energy<2,1>(tri2); 	
							
							cot1=compute_cot(tri1,1); 

							cot2=compute_cot(tri2,1);  
						}
					}
					
					
					 
					
					if((DTenergy>NDTenergy) && cot1!=cot2){
						
						//std::cout<< cot1-cot2 <<std::endl;
						//std::cout << "Vertices: " <<std::endl; 
						//for(auto v_itr=dt.finite_vertices_begin(); v_itr!=dt.finite_vertices_end(); v_itr++){
						//	std::cout<< "("<<(v_itr->point()).x()<<", " << (v_itr->point()).y() <<")"<<std::endl; 

						//}

						//std::cout<<"Faces: " <<std::endl; 
						//int num_faces=0; 
						
						//std::cout << "\\begin{tikzpicture}" <<std::endl; 
						//for(auto f_itr=dt.finite_faces_begin(); f_itr!=dt.finite_faces_end(); f_itr++){
						//	num_faces++; 
							//std::cout <<" \\draw (" << ((f_itr->vertex(0))->point()).x() <<" , " << ((f_itr->vertex(0))->point()).y() << " ) -- ("<< ((f_itr->vertex(1))->point()).x() <<" , " << ((f_itr->vertex(1))->point()).y() << " ) -- (" <<((f_itr->vertex(2))->point()).x() <<" , " << ((f_itr->vertex(2))->point()).y() << " ) -- (" <<((f_itr->vertex(0))->point()).x() <<" , " << ((f_itr->vertex(0))->point()).y() << " );" << std::endl;
						//}
						//std::cout <<"\\end{tikzpicture}" <<std::endl<<std::endl; 
						//std::cout <<"DT: " << DTenergy<< std::endl; 
						//std::cout<<"NDT: " << NDTenergy<<std::endl;  

						
						//std::cout<<"DT-NDT: "<< DTenergy-NDTenergy <<std::endl;
						//std::cout <<"DT/NDT: "<<DTenergy/NDTenergy <<std::endl;  
						//std::cout << std::endl; 
						
						num_DTgreaterNDT++; 
						//std::cout<<"cot1:" <<cot1 <<std::endl;
						//std::cout<<"cot2: " << cot2 <<std::endl; 
						std::cout << "\\begin{tikzpicture}" <<std::endl;	
						for(auto f_itr=dt.finite_faces_begin(); f_itr!=dt.finite_faces_end(); f_itr++){
									std::cout <<" \\draw (" << ((f_itr->vertex(0))->point()).x() <<" , " << ((f_itr->vertex(0))->point()).y() << " ) -- ("<< ((f_itr->vertex(1))->point()).x() <<" , " << ((f_itr->vertex(1))->point()).y() << " ) -- (" <<((f_itr->vertex(2))->point()).x() <<" , " << ((f_itr->vertex(2))->point()).y() << " ) -- (" <<((f_itr->vertex(0))->point()).x() <<" , " << ((f_itr->vertex(0))->point()).y() << " );" << std::endl;
							
							//std::cout<<"cot1-cot2:"<< cot1-cot2 <<std::endl;
						}
						std::cout << "\\end{tikzpicture}" <<std::endl;
						std::cout<<std::endl; 
						
					}

					angle+=angle_step; 		
				}
				length+=length_step; 
			}
			y+=y_step; 
		}	
		x+=x_step;		
	}
	std::cout<< "Number of experiments ran: " <<total_experiments <<std::endl; 
	std::cout <<"Number of experiements with DT having 2 faces: " << num_exp_2_faces <<std::endl; 
	std::cout<<"Number of times DT> NDT: " << num_DTgreaterNDT <<std::endl; 
	double percent=num_DTgreaterNDT/total_experiments;
	std::cout<< "(#DT>NDT)/#exp = " << percent<<std::endl; 
	std::cout <<"#DT>NDT, 2 faces/# exp with 2 faces= " << num_DTgreaterNDT/num_exp_2_faces<<std::endl; 
  return 0; 
}

double compute_cot(Triangle tri, int vertex_index){
	Point point0=tri.vertex(vertex_index);
	Point point1_opp=tri.vertex(vertex_index+1);
	Point point2_opp=tri.vertex(vertex_index+2);
	return ((point1_opp -point0)*(point2_opp-point0))/sqrt(squared_distance(point1_opp, point0)*squared_distance(point2_opp, point0)-pow((point1_opp -point0)*(point2_opp-point0),2));


}


