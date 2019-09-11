#ifndef __PARTICLETRAJECTORY_H__
#define __PARTICLETRAJECTORY_H__
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include "sgn.h"
#include "fzerosolve.h"  //fsolve has boost libraries. 
//#include <boost/math/special_functions.hpp>
//#include <boost/random.hpp>
//#include <boost/random/uniform_real_distribution.hpp>
#define n_gyro -1.83247172e8 //rad/s/T

class particleTrajectory{
private:
	std::string personal ="please don't read this string, it's personal and private." ;

public:	
 		
//boost random engine. 
typedef boost::mt19937 RNG;    // Mersenne Twister
typedef boost::random::uniform_real_distribution<double> DIST_flat;
typedef boost::variate_generator<RNG, DIST_flat> RAND_flat;
RNG rng;
DIST_flat dist_flat;
RAND_flat rand_flat;


static constexpr double m_n=1.67492749804E-27; //neutron mass codata 
static constexpr double pi=3.141592653589793; //ratio of circumference to diameter of circle in euclidian space. 
static constexpr double gravity=9.80665; // m/s^2
double diffuse;
double* oldspin=(double*)new double[3];
double* spin=(double*)new double[3];
double *tempspin=(double*)new double[3];
double tempspin0;
double* oldposition=(double*)new double[3];
double* position=(double*)new double[3];
double* oldvelocity=(double*)new double[3];
double* velocity=(double*)new double[3];
double  oldtime;
double  time;
double  twall; //time until next wall collision. Doesn't need to be complicated with no gas collisions for neutrons. 
double  oldtwall; //save state for presumed complications. 
int whatwall;
int oldwhatwall; //rect:0 is x 1 is y 2 is z. cyl: 0 is circumference 1 is length, sphere is irrelevent. 

//not so virtual doubles...
double xwidth;
double ylength;
double zheight;	
double length;
double radius;

//solver stuff, using sqrts so only want to do it once at the beginning. 
double xmax;
double xguess;
double vmax;
double vmag;
//solver stuff for cylinder only. not really used yet. slower than exact solution. 
double rmax;
double rguess;
double vperp;

//neutron mirror neutron precession parameters (constants). 
double w0;  //This is along z axis only (by definition of Lab frame) (I think). 
double w0m; //assuming mirror n has same n_gyro..
double phim;  //phim and thetam determine the orientation of mirror field with respect to w0. (which is z by def.)                   
double thetam;
double eps;  //coupling strength. 
//derived constants
double theta; //the angle of some nn' spinor basis, 
double Theta; //equation 3
double absw0wm; //used in equation 38. 
double wp; //equation 38
double wm;
double Phip; //equation 37
double Phim;
//time dependent
double mu;
double nuc;
double nus;
double phi;
double phiwall;

//this way seems harder numerically. might be useful if we ever use RK. 
//double Wron;
//double Wron_s;
//double Wron_c;
//double Wron_p;
//double Wron_m;

//necessary for multi threading, used by the conductor for saving this identifies what number neutron it is. 
int nnn;


//Constructor
particleTrajectory(double vmax, double diffuse, int flat_seed) : rng(flat_seed),dist_flat(0.0f,1.0f),rand_flat(rng,dist_flat),vmax(vmax),diffuse(diffuse),time(0.),oldtime(0.){};
	




// move the paritlce to the new position (Object dependent.) 
virtual void moveParticle(double dt){throw std::invalid_argument("this is a place holder, please define a child geometry object\n");};
//start inside of whatever object is defined (object dependent) 
virtual void startInside(){throw std::invalid_argument("this is a place holder, please define a child geometry object\n");};
// find the time to the next boundary collision. 
virtual void nextBoundary(){throw std::invalid_argument("this is a place holder, please define a child geometry object\n");};


//this computes the "integral" from the analytic solution, not really integrating now, was used for rk5. 
void integrateMP();

//set field values for neutron mirror neutron precession and calculate constants
void loadFields(double B0, double B0m, double tempthetam, double tempphim, double eps_temp);

//gets an isotropic velocity distribution after finding the height 
//height is sampled according to the distribution eq 4.7 in Bob G's book. 
double getIsotropic(double floor, double ceiling);



//E.D. Davis neutron mirror-neutron precession solution (at least I hope it is). 
void analyticSolution(double t);

//reset the spin state for adaptive rk (object independent (mostly)) 
//this isn't used right now, but if we go to rk integraiton it will be nice.
void resetState()
{
 	for (int i=0; i<3;i++) 
 	{
 		position[i]=oldposition[i];
 		velocity[i]=oldvelocity[i];
 		spin[i]=oldspin[i];
 	}
 	twall=oldtwall;
 	whatwall=oldwhatwall;
 	time=oldtime;
};

//save the spin state for adaptive rk (object independent)
void storeState()
{
 	for (int i=0; i<3;i++) 
 	{
 		oldposition[i]=position[i];
 		oldvelocity[i]=velocity[i];
 		oldspin[i]=spin[i];
 	}
 	oldtime=time;
 	oldwhatwall=whatwall;
 	oldtwall=twall;
};


};

//chitlins 
class Cylinder: public particleTrajectory {
public:
	Cylinder(double length,double radius, double vmax, double diffuse, int flat_seed): particleTrajectory(vmax,diffuse,flat_seed), length(length), radius(radius){};
	void startInside();
	void nextBoundary();
	void moveParticle(double dt);
    

    double length;
    double radius;
};

class Rectangle: public particleTrajectory {
public:
	Rectangle(double xwidth,double ylength,double zheight, double vmax, double diffuse, int flat_seed): particleTrajectory(vmax,diffuse,flat_seed), xwidth(xwidth), ylength(ylength), zheight(zheight) {};
	void startInside();
	void nextBoundary();
	void moveParticle(double dt);
	
	double xwidth;
	double ylength;
	double zheight;

};

class Sphere: public particleTrajectory {
public:
	Sphere(double radius, double vmax, double diffuse,int flat_seed): particleTrajectory(vmax,diffuse,flat_seed), radius(radius) {};
	void startInside();
	void nextBoundary();
	void moveParticle(double dt);
	
	double radius;
};
#endif

