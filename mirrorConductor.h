#ifndef __MIRRORCONDUCTOR_H__
#define __MIRRORCONDUCTOR_H__

#include <boost/thread.hpp>
#include <fstream>
#include <iostream>
#include <sstream>

class mirrorConductor {
private: 
	int N,timesteps,Nactual; //number of spins/time-steps simulated. Nacutal is the real spins number. 
	int startseed; //random seed. 
	double T; //total time in seconds. 
	double dt; //data save time steps, (does not have to be integer time steps to T.)
	double error; //error per step;
	static constexpr double pi=3.141592653589793;
	double rkdt; //initial guess at the time step of the Runge-Kutta integrator. Should be ~5 times smaller than the collision time of he-3
	std::ofstream  spinfile;
	std::string outfilename;
public:
	std::vector<boost::thread *> bthreads;
	std::vector<particleTrajectory*> tr;
	
	//construction
	mirrorConductor(std::string infilename,std::string ofilename);

	bool FileExists(const std::string& name);
	bool parseParameters(const std::string& name);
	

	//conduct all spins, write to disk once at the end. 
	void Conduct();
	
	//this is the point of contact for the particleTrajectory class
	//particleTrajectory contains all the boundary spin, field information as well as the integrator. ,
	//run a single spin
	void run(particleTrajectory * tr);
	//run a single trajectory but don't do spin integraiton, for some reason. 
	void runNS(particleTrajectory * tr);
	


	void saveSpins();


	double* Position;
	double* Velocity;
	double* Spins;
	//double* BField;
	double*  Phase;
	//double* EField; 
	
	//for parameter file parsing. 
	double ts;
	double nvars;
	double Nsave;

	int geom;
	double diffuse;
	double xwidth;
	double ylength;
	double zheight;
	double vmax;	
	double length;
	double radius;
	double thetam, B0, B0m, phim, eps;
	int cores;
	bool scanbool=false;
	double *scanvar;
	double scanstart;
	double scanstop;
	double scannumber=0;

};
#endif