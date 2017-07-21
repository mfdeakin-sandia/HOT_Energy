#ifndef _ANALYTIC_HPP_
#define _ANALYTIC_HPP_

//////////////////////////////////////////////////////////////////////////////////
/////////////////////  ENERGY DERIVATIVES /////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

template<typename T>
void energy_gradient(T &triangulation,int Wk, int star, Vertex_handle v, double total_deriv[2], bool corrected_formulas){

//template< > 
//void energy_gradient<2,int star>(const DT &dt, Vertex_handle v, double total_deriv[2], bool //corrected_formulas){
	
	total_deriv[0]=0; 
	total_deriv[1]=0; 


	for(auto ei=triangulation.finite_edges_begin();ei!=triangulation.finite_edges_end(); ei++){
	
		Edge edge=*ei;
		Face edge_face=*(edge.first); 
		int edge_index=edge.second; 
		
		Edge mirror_edge=triangulation.mirror_edge(*ei);
		Face mirror_edge_face=*(mirror_edge.first); 
		int mirror_index=mirror_edge.second; 
		
		Vertex_handle  vk_handle=edge_face.vertex(edge_index); 
		Vertex_handle vl_handle=mirror_edge_face.vertex(mirror_index); 
		Vertex_handle vi_handle=edge_face.vertex(edge_face.cw(edge_index)); 
		Vertex_handle vj_handle=edge_face.vertex(edge_face.ccw(edge_index));
		

		Point pi=vi_handle->point();
		Point pj=(*vj_handle).point(); 
		Point pk=(*vk_handle).point(); 
		Point pl=(*vl_handle).point();


		// for now, we assume equal weights, so dij=dji. 
		double dij=0.5*sqrt(pow(pi.x()-pj.x(),2.0)+pow(pi.y()-pj.y(),2.0));
		double dji=0.5*sqrt(pow(pi.x()-pj.x(),2.0)+pow(pi.y()-pj.y(),2.0)); 

		// intialize derivative arrays
		double dij_derv[2]={0,0};
		double dji_derv[2]={0,0}; 
		double hk_derv[2]={0,0} ; 
		double hl_derv[2]={0,0} ;
		if(v==vl_handle){
			hk_derv[0]=0; 
			hk_derv[1]=0;

			compute_h_deriv(pi,pj,pl, 3, hl_derv);
					
			dij_derv[0]=0;
			dij_derv[1]=0; 

			dji_derv[0]=0;
			dji_derv[1]=0;
		} 

		else if(v==vk_handle){
			compute_h_deriv(pi,pj,pk, 3, hk_derv);

			hl_derv[0]=0; 
			hl_derv[1]=0;

			dij_derv[0]=0;
			dij_derv[1]=0; 

			dij_derv[0]=0;
			dij_derv[1]=0;
		}

		else if(v==vi_handle){
			compute_h_deriv(pi,pj,pk, 1, hk_derv);
			compute_h_deriv(pi,pj,pl, 1, hl_derv);

			dij_derv[0]=-(pj.x()-pi.x())/(2.0*(dij+dji)); 
			dij_derv[1]=-(pj.y()-pi.y())/(2.0*(dij+dji));

			dji_derv[0]=dij_derv[0];  // since for now we assume = weights
			dji_derv[1]=dij_derv[1];  // = weights
		}
				
		else if(v==vj_handle){
			compute_h_deriv(pi,pj,pk, 2, hk_derv);
			compute_h_deriv(pi,pj,pl, 2, hl_derv);
			
			dij_derv[0]=(pj.x()-pi.x())/(2.0*(dij+dji));  // assuming equal weights
			dij_derv[1]=(pj.y()-pi.y())/(2.0*(dij+dji));

			dji_derv[0]=dij_derv[0];
			dji_derv[1]=dij_derv[1];		
		}

		else continue; // this means that the edge engery is unchanged by the movement of v 

	
		double hk=signed_dist_circumcenters(face_to_tri(edge_face), edge_index);
		double hl=signed_dist_circumcenters(face_to_tri(mirror_edge_face), mirror_index); 
		int sign=sgn(hk+hl); 

		// I want to make constants that varry based on what Wk and star are. 
			double constant1; 
			double constant2;
		if(star==0){
			constant1=4.0;
			constant2=12.0; 
		}
		else if(star==1){
			constant1=3.0; 
			constant2=3.0;
		}
		else{
			constant1=12.0; 
			constant2=4.0;

		}
		
		bool boundary_edge= triangulation.is_infinite(edge.first) ||triangulation.is_infinite(mirror_edge.first); 
			
			if(!triangulation.is_infinite(edge.first)){
				if(!boundary_edge || hk>0){
					double edge_deriv[2]={0,0}; 
					for(int coor=0; coor <2 ; coor++){
						edge_deriv[coor]+=(3*pow(dij,2)*hk*dij_derv[coor]+pow(dij,3)*hk_derv[coor])/constant1;
						edge_deriv[coor]+=3*pow(dji,2)*dji_derv[coor]*hk+pow(dji,3)*hk_derv[coor]/constant1;
						edge_deriv[coor]+= dij*3*pow(hk,2)*hk_derv[coor]+dij_derv[coor]*pow(hk,3)/constant2;
						edge_deriv[coor]+=dji*3*pow(hk,2)*hk_derv[coor]+dji_derv[coor]*pow(hk,3)/constant2;
							
						if(star==1){
							if(boundary_edge) total_deriv[coor]+=edge_deriv[coor];
							else if(corrected_formulas) total_deriv[coor]+=sign*edge_deriv[coor]; 						
							else total_deriv[coor]+=edge_deriv[coor]; 
						}
						// otherwise star=0 or 2. for now, we make no special treatment for boundary triangles for these stars.
						else total_deriv[coor]+=edge_deriv[coor]; 
					}
				} 
			}

			if(!triangulation.is_infinite(mirror_edge.first)){
				if(!boundary_edge || hl >0){
					double edge_deriv[2]={0,0}; 
					for(int coor=0; coor <2 ; coor++){
						edge_deriv[coor]+=3*pow(dij,2)*hl*dij_derv[coor]+pow(dij,3)*hl_derv[coor]/constant1;
						edge_deriv[coor]+=3*pow(dji,2)*dji_derv[coor]*hl+pow(dji,3)*hl_derv[coor]/constant1;
						edge_deriv[coor]+= dij*3*pow(hl,2)*hl_derv[coor]+dij_derv[coor]*pow(hl,3)/constant2;
						edge_deriv[coor]+= dji*3*pow(hl,2)*hl_derv[coor]+dji_derv[coor]*pow(hl,3)/constant2;
					
						if(star==1){
							if(boundary_edge) total_deriv[coor]+=edge_deriv[coor];
							else if (corrected_formulas) total_deriv[coor]+=sign*edge_deriv[coor]; 
							else total_deriv[coor]+=edge_deriv[coor];
						}
						// otherwise star=0 or 2. for now, we make no special treatment for boundary triangles for these stars.
						else total_deriv[coor]+=edge_deriv[coor]; 
					}
				} 
			}		
			
		} 
	// rescale to be unit vector
	double grad_length=sqrt(pow(total_deriv[0],2)+ pow(total_deriv[1],2)); 
	std::cout<<"unit vector in direction gradient:" << std::setw(15) << total_deriv[0]/grad_length << std::setw(15) <<total_deriv[1]/grad_length <<std::endl; 
	return; 
}


void compute_h_deriv(const Point &xi, const Point &xj, const Point &xk, int i, double h_derv[2]){

	double xi1=xi.x();  
	double xi2=xi.y(); 
	double xj1=xj.x(); 
	double xj2=xj.y(); 
	double xk1=xk.x(); 
	double xk2=xk.y(); 

	if(i==1){ // diff wrt xi coordinates
		
		h_derv[0]= (2*(pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*(xj1 - xk1)*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2) - 
     ((xi1 - xk1)*(xj1 - xk1) + (xi2 - xk2)*(xj2 - xk2))*(2*(pow(xi1 - xj1,2) + pow(xi2 - xj2,2))*(xj2 - xk2)*(-(xj2*xk1) + xi2*(-xj1 + xk1) + xi1*(xj2 - xk2) + xj1*xk2) + 
        2*(xi1 - xj1)*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2)))/ (4.*pow((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2),1.5));


		h_derv[1]= ((xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2))*((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*(xj2 - xk2)*(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2)) -  ((xi1 - xk1)*(xj1 - xk1) + (xi2 - xk2)*(xj2 - xk2))*((pow(xi1 - xj1,2) + pow(xi2 - xj2,2))*(xj1 - xk1) + (xi2 - xj2)*(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2)))))/ (2.*pow((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2),1.5)); 
	}


	if(i==2){ // diff wrt xj coordinates
		h_derv[0]= ((xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2))*((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*(xi1 - xk1)*(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2)) -  ((xi1 - xk1)*(xj1 - xk1) + (xi2 - xk2)*(xj2 - xk2))*((pow(xi1 - xj1,2) + pow(xi2 - xj2,2))*(xi2 - xk2) + (xi1 - xj1)*(-(xj2*xk1) + xi2*(-xj1 + xk1) + xi1*(xj2 - xk2) + xj1*xk2))))/(2.*pow((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2),1.5));


		h_derv[1]= (2*(pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*(xi2 - xk2)*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2) - 
     ((xi1 - xk1)*(xj1 - xk1) + (xi2 - xk2)*(xj2 - xk2))*(2*(pow(xi1 - xj1,2) + pow(xi2 - xj2,2))*(xi1 - xk1)*(-(xj2*xk1) + xi2*(-xj1 + xk1) + xi1*(xj2 - xk2) + xj1*xk2) - 
        2*(xi2 - xj2)*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2)))/(4.*pow((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2),1.5));
	}

	if(i==3){ //diff wrt xk coordinates
		h_derv[0]= (xj2*pow(xk1,2) + pow(xi1,2)*(xj2 - xk2) + pow(xi2,2)*(xj2 - xk2) + pow(xj1,2)*xk2 + pow(xj2,2)*xk2 - 2*xj1*xk1*xk2 - xj2*pow(xk2,2) + 2*xi1*xk1*(-xj2 + xk2) - 
     xi2*(pow(xj1,2) + pow(xj2,2) - 2*xj1*xk1 + pow(xk1,2) - pow(xk2,2)))/(2.*(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2))*sqrt((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2)));


		h_derv[1]= (-(pow(xj1,2)*xk1) - pow(xj2,2)*xk1 + xj1*pow(xk1,2) + pow(xi1,2)*(-xj1 + xk1) + pow(xi2,2)*(-xj1 + xk1) + 2*xi2*(xj1 - xk1)*xk2 + 2*xj2*xk1*xk2 - xj1*pow(xk2,2) + 
     xi1*(pow(xj1,2) + pow(xj2,2) - pow(xk1,2) - 2*xj2*xk2 + pow(xk2,2)))/(2.*(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2))*sqrt((pow(xi1,2) + pow(xi2,2) - 2*xi1*xj1 + pow(xj1,2) - 2*xi2*xj2 + pow(xj2,2))*pow(xi2*(xj1 - xk1) + xj2*xk1 - xj1*xk2 + xi1*(-xj2 + xk2),2)));
	}

}




#endif
