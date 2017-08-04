#include <hot.hpp>
#include <Sb.hpp>
#include <ply_writer.hpp>

#include <random>
#include <vector>
#include <fstream>

#include <CGAL/Kernel/global_functions.h>


int main(int argc, char **argv) {
	double circle radius=1;

	RegT	rt; 

	std::ofstream outputFile; 
	
	outputFile.open("vertices_in_circle.txt");
	outputFile << (v_itr->point()).x() << (v_itr-> point).y() <<std::endl; 
	outputFile.close()

	outputFile.open("triangles_in_circle.txt"); 

	outputFile.close()

	return 0; 
}


