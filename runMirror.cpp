//g++ TestingMirror.cpp mirrorConductor.cpp particleTrajectory.cpp -I/data1/cmswank/BoostCodeSwank/boost_1_64_0 -o testing2 -lboost_thread -lboost_system
#include "particleTrajectory.h"
#include "mirrorConductor.h"
#include "fzerosolve.h"
#include <iostream>
#include <cstdio>
#include <ctime>





//run the parameter file. 
int main(int argc, char *argv[])
{   
   
    std::clock_t start;
    double duration;
    start = std::clock();

    std::string infilename=argv[1];	
	std::string outfilename=argv[2];	
	

	mirrorConductor* mc=new mirrorConductor(infilename,outfilename);


	mc->Conduct();

    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
    
    std::cout<<"duration: "<< duration <<" seconds using "<<mc->cores<<" cores\n";
    delete mc;

	return 0;
}	
	