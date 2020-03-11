#ifndef _ENERGYWEIGHTS_HPP_
#define _ENERGYWEIGHTS_HPP_

/////////////////////////////////////////////////////////////////////////
//////////////////////////   Energy for weighted triangulations ////////
////////////////////////////////////////////////////////////////

double triangle_energy_weights_dividebyArea(const weighted_Face_handle &face, const Point &wcirc, int Wk, int star);
double triangle_energy_weights(const weighted_Face_handle &face, const Point &wcirc, int Wk, int star);

template<typename T>
double energy_weights(const T &t, int Wk, int star){
  	double energy = 0;

  	for(auto face_itr = t.finite_faces_begin(); face_itr != t.finite_faces_end(); face_itr++) {
		Point wcirc=t.weighted_circumcenter(face_itr);
    		energy+=triangle_energy_weights(face_itr,wcirc, Wk,star);
 	 }

  return energy;
}

template<typename T>
double energy_weights_dividebyArea(const T &t, int Wk, int star){
  	double energy = 0;

  	for(auto face_itr = t.finite_faces_begin(); face_itr != t.finite_faces_end(); face_itr++) {
		Point wcirc=t.weighted_circumcenter(face_itr);
    		energy+=triangle_energy_weights_dividebyArea(face_itr,wcirc,Wk,star);
 	 }
  return energy;
}


double triangle_energy_weights(const weighted_Face_handle &face, const Point &wcirc, int Wk, int star){
	if(Wk!=2){
		std::cout<<"Warning: triangle_energy_weights returning bogus answer because Wk was not 2"<<std::endl;
		return -1;
	}  
	double energy=0; 

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
	
	
	for(int i=0; i<3; i++)
  {
    const auto ip1 = (i+1) % 3;
    const auto ip2 = (i+2) % 3;
    auto wi = face->vertex(i  )->point();
    auto wj = face->vertex(ip1)->point();
    auto wk = face->vertex(ip2)->point();
    
		const double weighti=wi.weight();
		const double weightj=wj.weight();
	
    Point xi(wi);
    Point xj(wj);
    Point xk(wk);

		const double length_eij=sqrt(pow(xi.x()-xj.x(),2.0)+pow(xi.y()-xj.y(),2.0));

		const double dij= (pow(length_eij,2) -weighti+weightj)/(2*length_eij);
		const double dji= (pow(length_eij,2) -weightj+weighti)/(2*length_eij);
		
		const double unsigned_hk=abs((xj.y()-xi.y())*wcirc.x() - (xj.x()-xi.x())*wcirc.y() +xj.x()*xi.y()-xj.y()*xi.x())/length_eij;
    double hk = (CGAL::orientation(xi,xj,xk)==CGAL::orientation(xi,xj,wcirc) ? unsigned_hk : -1.0*unsigned_hk);
		
		energy+=pow(dij,3)*hk/constant1+dij*pow(hk,3)/constant2;
		energy+=pow(dji,3)*hk/constant1+dji*pow(hk,3)/constant2;
	}
	return energy; 
}




double triangle_energy_weights_dividebyArea(const weighted_Face_handle &face, const Point &wcirc, int Wk, int star){
	
	if(Wk!=2){
		std::cout<<"Warning:triangle_energy_weights_dividedbyArea returning bogus answer because Wk was not 2"<<std::endl;
		return -1;
	}  
	
	Triangle tri=Triangle(Point(face->vertex(0)->point()), Point(face->vertex(1)->point()), Point(face->vertex(2)->point()));
	double face_area =std::abs(tri.area()); 	
	
	return triangle_energy_weights(face,wcirc, Wk, star)/face_area; 
}

#endif
