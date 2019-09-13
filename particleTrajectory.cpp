#include"particleTrajectory.h"
#include "sgn.h"
#include "fzerosolve.h"
#include <boost/math/special_functions.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>


void particleTrajectory::loadFields(double B0, double B0m, double tempthetam, double tempphim, double eps_temp){
		using namespace boost::math;
		w0=n_gyro*B0;  //This is along z axis only (by definition of Lab frame) (I think). 
		//assuming mirror n has same n_gyro..
		w0m=n_gyro*B0m;  //phim and thetam determine the orientation of mirror field with respect to w0. (which is z by def.) 
		phim=tempphim;                  
		thetam=tempthetam;
		
		eps=eps_temp;
		
		//equation 3
		Theta=acos((w0+w0m*cos(thetam))/sqrt(w0*w0+w0m*w0m+2*w0*w0m*cos(thetam)));

		//used to define equation 38. 
		absw0wm=sqrt(pow(w0+w0m*cos(thetam),2)+pow(w0m*sin(thetam)*cos(phim),2)+pow(w0m*sin(thetam)*sin(phim),2));

		//equation 38
		wp = sgn(w0)*sqrt(0.5*(w0*w0+w0m*w0m)+4.*eps*eps+sqrt(0.25*pow(w0m*w0m-w0*w0,2)+4.*eps*eps*absw0wm*absw0wm));
		wm = sgn(w0)*sqrt(0.5*(w0*w0+w0m*w0m)+4.*eps*eps-sqrt(0.25*pow(w0m*w0m-w0*w0,2)+4.*eps*eps*absw0wm*absw0wm));
			
		if (w0m==w0) theta=pi/2;  //this avoids divide by zero errors. c++ doesn't factor them out. 
		else theta=acos(0.5*(w0m*w0m-w0*w0)/sqrt(0.25*pow<2>(w0m*w0m-w0*w0)+4.*eps*eps*absw0wm*absw0wm)); ///should be exact I think. 

		//equation 37
		Phip=asin(w0*w0m*sin(thetam)/(absw0wm*wp));
		Phim=asin(w0*w0m*sin(thetam)/(absw0wm*wm));
}




void particleTrajectory::integrateMP()
{   //coded from E. D. Davis mansuscript (at least I hope it is)
	//ok so its not integrating, could have done it that was, but this way is easier and exacter. 

		//updates our spin and phase vectors assuming wall collision is a measurement. 
	analyticSolution(time-oldtwall);
	

/*
  //Runge-Kutta Cash Karp. put in .h file if that is ever an option. 
  static constexpr double a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875;
  static constexpr double b21 = 0.2;
  static constexpr double b31 = 0.075, b32 = 0.225;
  static constexpr double b41 = 0.3, b42 = -0.9, b43 = 1.2;
  static constexpr double b51 = -2.03703703703703692e-1, b52 = 2.5, b53 = -2.59259259259259256, b54 = 1.29629629629629628;
  static constexpr double b61 = 2.94958043981481469e-2,b62 = 3.41796875e-1,b63 = 4.15943287037037063e-2,b64 = 4.00345413773148140e-1,b65 = 6.1767578125e-2;
  static constexpr double c1 = 9.78835978835978782e-2, c2 = 0.,c3 = 4.02576489533011284e-1,c4 = 2.10437710437710451e-1,c5 = 0.,c6 = 2.89102202145680387e-1;
  static constexpr double d1 = 1.02177372685185189e-1,d2 = 0.,d3 = 3.83907903439153431e-1,d4 = 2.44592737268518517e-1,d5 = 1.93219866071428562e-2,d6 = 0.25;

*/
	return;
}

void particleTrajectory::analyticSolution(double t){
	//t is time from previous wall collision, (at least it should be)
	using namespace boost::math;

	mu=pow<2>(cos(0.5*theta))*cos(0.5*wm*t)+pow<2>(sin(0.5*theta))*cos(0.5*wp*t);
		                                                   
	nuc=pow<2>(cos(0.5*theta))*cos(Theta-Phim)*sin(0.5*wm*t)+pow<2>(sin(0.5*theta))*cos(Theta-Phip)*sin(0.5*wp*t);

	nus=pow<2>(cos(0.5*theta))*sin(Theta-Phim)*sin(0.5*wm*t)+pow<2>(sin(0.5*theta))*sin(Theta-Phip)*sin(0.5*wp*t);

	spin[0]=tempspin0=pow<2>(mu)-pow<2>(nuc)+pow<2>(nus)*cos(2.*phim);

	spin[1]=2.*mu*nuc+pow<2>(nus)*sin(2.*phim);

	spin[2]=2.*nus*(nuc*cos(phim)-mu*sin(phim));	

    //shift by the phase from the last wall collision

	spin[0]=spin[0]*cos(phiwall)-spin[1]*sin(phiwall);

	spin[1]=spin[1]*cos(phiwall)+tempspin0*sin(phiwall);

	phi=atan2(spin[1],spin[0]);
	
	return;
}

double particleTrajectory::getIsotropic(double floor,double ceiling)
{//returns height of particle from Golub's UCN book, distribution from equation 4.7. (althought this is not 100% correct due to parabolic arcs)
 //stores the v^2 distributed from maximum velocity. 
	//not sure if this really works with the sphere, but its probably pretty good.
 	using namespace boost::math;
 	double L = ceiling-floor;
 	//hard v^2 distribution. at the lowest point in the geometry (ergodicity assumed)
	vmag=vmax*cbrt(rand_flat());

	double height;
	if(vmag*vmag>2.*gravity*L){
    	double vfactor=std::sqrt(1.-2.*L*gravity/vmag/vmag);
    	height= 0.5/gravity*(vmag*vmag-cbrt(pow<2>( (vfactor*(pow<3>(vmag)-2*L*gravity*vmag) - pow<3>(vmag) )*rand_flat() + pow<3>(vmag))));      
    }
    else{
    	height=0.5/gravity*vmag*vmag*(1-cbrt(pow<2>(rand_flat())));
    }

    //correct for bad luck and rounding errors... height starts below zero.
    height = (height>0.0) ? height:0.0;

    //correct for bad luck and round errors 2... height starts above ceiling.
  	height = (height+floor>ceiling)? L:height;
     
    //magnitude of the velocity at the height  
    //(it is possible that this populates velocities in places that the neutron cannot replicate from diffuse wall collisions alone, but it should be good enough. )
	double vmagh=std::sqrt(vmag*vmag-2.*gravity*height);
	double thetap=std::acos(2.*rand_flat()-1.);
    double phip= 2.*pi*rand_flat();
    //velocity components ajusted for height. 
	velocity[0]=vmagh*sin(thetap)*cos(phip);
	velocity[1]=vmagh*sin(thetap)*sin(phip);
	velocity[2]=vmagh*cos(thetap);

	  //return in geometry coordinates
	return height+floor;
}



//#####################################
//# implementation of startInside
//#####################################
void Cylinder::startInside()
{//initialize spin state. 
	//position
	do
	{
		position[0]=2.*radius*(rand_flat()-0.5);
		position[1]=2.*radius*(rand_flat()-0.5);
	}
	while(std::sqrt(position[0]*position[0]+position[1]*position[1])>radius);
	
	//fill velocity from vmax and v^2 distribution and get the corrispoinding height sampled from eqn 4.7 in Bob's book. 
	position[2]=getIsotropic(-length/2.,length/2.);
	
	
	for(int i=0; i<3; i++) oldposition[i]=position[i];
	
	//get initial velocity
	
	
	for(int i=0; i<3; i++) oldvelocity[i]=velocity[i];

	//taken out initialization of spin. like, feel free to put it back if you want it 
	spin[0]=1.;
	spin[1]=0.;
	spin[2]=0.;
	phiwall=0.;

	xmax=2.001*std::sqrt(2.*gravity*length/vmag);
	xguess=std::sqrt(2.*gravity*length/vmag);

	vperp=sqrt(velocity[0]*velocity[0]+velocity[1]*velocity[1]);
	

	rmax=vperp/2./radius;
	//do the integral to get a better guess. 
	rguess=vperp/radius;

	for(int i=0; i<3; i++) oldspin[i]=spin[i];
	this->nextBoundary();
	oldtwall=0.;
	return;


}


//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@ nextBoundary
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void Cylinder::nextBoundary(){
	double twalltemp1=(velocity[0]==0 && velocity[1]==0)  ? 1.0E32 :(-velocity[0]*position[0]-velocity[1]*position[1]\
				+std::sqrt(std::pow(radius,2)*((std::pow(velocity[0],2)\
				+std::pow(velocity[1],2)))-std::pow(position[1]*velocity[0]-position[0]*velocity[1],2)))\
		    	/(std::pow(velocity[0],2)+std::pow(velocity[1],2));  //i like rectangles.  
   		
  
  
    //+-zwall=z0+vz*twall-1/2*g*twall^2
   	double twalltop=(velocity[2]-std::sqrt(velocity[2]*velocity[2]+2.0*(position[2]-length/2.0)*gravity))/gravity;
   	double twallbott=(velocity[2]+std::sqrt(velocity[2]*velocity[2]+2.0*(position[2]+length/2.0)*gravity))/gravity;
   	
   	double twalltemp2= twalltop<0. ? twallbott:(twallbott<0. ? twalltop :(twalltop<twallbott ? twalltop:twallbott)); //can't be equal...

   		
		
   		twall=twalltemp2<=twalltemp1 ? time+twalltemp2:time+twalltemp1; //if by some nonsense they are equal, it chooses the top/bottom wall. 
   		whatwall=twalltemp2<twalltemp1 ? 1:0;
   		
   	return;
}





//***********************************************
//*implementation of moveParticle!
//***********************************************
void Cylinder::moveParticle(double dt)
{
      
	if(time+dt<twall)
	{      
		position[0]+=velocity[0]*dt;
    	position[1]+=velocity[1]*dt;
		position[2]+=velocity[2]*dt-0.5*gravity*dt*dt;
		velocity[2]+=-gravity*dt;
		time+=dt;
			//recursive until no wall is hit. 
		return;//the actual exit is here. 
	}

	double dtwall=twall-time;

	analyticSolution(twall-oldtwall);
	phiwall=phi;
	storeState();

	position[0]+=velocity[0]*dtwall;
    position[1]+=velocity[1]*dtwall;
	position[2]+=velocity[2]*dtwall-0.5*gravity*dtwall*dtwall;
	velocity[2]+=-gravity*dtwall;


	
	if(whatwall==1)
	{	
			//might be necessary to clamp to the wall if it hits the corner. 
   			//position[2]=position[2]>0? length/2.0:-length/2.0;
   			
		if(rand_flat()>diffuse)
		{//specular
			
			velocity[2]=-velocity[2];
		}
		else
		{//diffuse bounce on the ends. whats the cdf of cos(theta)^2*sin(theta) again? theta0=arccos((1-rand)^(1/3))
		 //ugg someone from boost should write a diffuse ditrbituion code, (they do have an arcsine dist) this is slow. 
		//maybe template this in a spline for extra speed? lets see how bad it is. 
			//method one.
			double vmag=std::sqrt(velocity[2]*velocity[2]+velocity[0]*velocity[0]+velocity[1]*velocity[1]);
			double temptheta=std::acos(boost::math::cbrt(rand_flat()));
			//don't allow values that can escape from the surface on the cylinder. 
			temptheta=temptheta<(pi/2.0-1E-4)?temptheta:pi/2.0-1E-4;
			double tempphi=2.0*pi*rand_flat();
		    velocity[2]=-vmag*sgn(position[2])*std::cos(temptheta);
		    velocity[0]=vmag*std::sin(temptheta)*std::cos(tempphi);
		    velocity[1]=vmag*std::sin(temptheta)*std::sin(tempphi);
		  
		}

	}
	else
	{
		//might be necessary to clamp to the wall if it hits the corner/edge. 
		//hmm it would be nice to delete this. 
		//double bounce_norm = std::sqrt(position[0]*position[0]+position[1]*position[1]);

		//position[0]=position[0]/bounce_norm*radius;
   		//position[1]=position[1]/bounce_norm*radius;

   		
		if(rand_flat()>diffuse)
		{//specular
			
			double vr=-velocity[0]*position[0]/radius-velocity[1]*position[1]/radius;//sign flipped here 
			double vp=-velocity[0]*position[1]/radius+velocity[1]*position[0]/radius;
			
			
			velocity[0]=-vp*position[1]/radius+vr*position[0]/radius; 
			
			velocity[1]=vr*position[1]/radius+vp*position[0]/radius; 
			
			

			//Normalize(velocity,vmag); I no longer think this is necessary here actually.

			
		}
		else
		{//diffuse bounce on the circumference. 
		double vmag=std::sqrt(velocity[2]*velocity[2]+velocity[0]*velocity[0]+velocity[1]*velocity[1]);
		double temptheta=std::acos(boost::math::cbrt(rand_flat()));
		double tempphi=2.0*pi*rand_flat();
		
		
		double vr=-vmag*std::cos(temptheta); //make negative from reflection.
		double vp=vmag*std::sin(temptheta)*std::sin(tempphi);
		
		velocity[2]=vmag*std::sin(temptheta)*std::cos(tempphi);
        velocity[0]=(-vp*position[1]/radius+vr*position[0]/radius);   //this is all I remember from E&M class. 
		velocity[1]=vr*position[1]/radius+vp*position[0]/radius; 
		
		//Normalize(velocity,vmag); probably not necessary. 
	
		}	
	}

	//should have new velocity[] and position 
	//need to update time
	
	time=twall;
	//maybe we will hit the wall again before dt is over? lets go to infinity and find out!
	this->nextBoundary();
	this->moveParticle(dt-dtwall);
	
	return;//exiting from recursive calls.
}


///:):):):):):):):):):):):):):):):):):):):):):):):)
/// Very happy implementation of startInside. 
///:):):):):):):):):):):):):):):):):):):):):):):):)
void Rectangle::startInside()
{
	position[0]=xwidth*(rand_flat()-0.5);
	position[1]=ylength*(rand_flat()-0.5);
	position[2]=getIsotropic(-0.5*zheight,0.5*zheight);
	
	

	for(int i=0; i<3; i++) oldposition[i]=position[i];
	

	
	
	for(int i=0; i<3; i++) oldvelocity[i]=velocity[i];

	spin[0]=1.;
	spin[1]=0.;
	spin[2]=0.;
	
	phiwall=0.;

	xmax=2.001*std::sqrt(2.*gravity*zheight/vmag);
	xguess=std::sqrt(2.*gravity*zheight/vmag);

	for(int i=0; i<3; i++) oldspin[i]=spin[i];
	this->nextBoundary();
	oldtwall=0.;	
	return;
}


//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//$ Wealthy implementation of nextBoundary
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
void Rectangle::nextBoundary()
{
	//this may look similar to the z walls in the cylinder, but its different.
    //+-zwall=z0+vz*twall-1/2*g*twall^2
   	double twalltop=(velocity[2]-std::sqrt(velocity[2]*velocity[2]+2.0*(position[2]-zheight/2.0)*gravity))/gravity;
   	double twallbott=(velocity[2]+std::sqrt(velocity[2]*velocity[2]+2.0*(position[2]+zheight/2.0)*gravity))/gravity;
   	
   	double twalltemp2= twalltop<0. ? twallbott:(twallbott<0. ? twalltop :(twalltop<twallbott ? twalltop:twallbott)); //can't be equal...

   	double twalltemp0=velocity[0]==0.0 ? 1.0E32 : (xwidth/2.0-sgn(velocity[0])*position[0])/(sgn(velocity[0])*velocity[0]);
	double twalltemp1=velocity[1]==0.0 ? 1.0E32 : (ylength/2.0-sgn(velocity[1])*position[1])/(sgn(velocity[1])*velocity[1]);
   	
   	//goofball inline logic. 		
   	whatwall=twalltemp2<=twalltemp1 ? (twalltemp0<twalltemp2 ? 0:2):(twalltemp0<twalltemp1 ? 0:1);
   	switch (whatwall)
			{	
			case 0: twall=time+twalltemp0;
			break;
			case 1: twall=time+twalltemp1;
			break;
			case 2: twall=time+twalltemp2;
			break;
			default: throw std::invalid_argument("what the fuck?!\n");
			}

   	return;
}

void Rectangle::moveParticle(double dt)
{   
	if(time+dt<twall)
	{      
		position[0]+=velocity[0]*dt;
    	position[1]+=velocity[1]*dt;
		position[2]+=velocity[2]*dt-0.5*gravity*dt*dt;
		velocity[2]+=-gravity*dt;
		time+=dt;
			//recursive until no wall is hit. 
			return;//the actual exit is here at the beginning, Keanu: "whoa". 
	}

	double dtwall=twall-time;

	//interesting must be included here I guess. 
	analyticSolution(twall-oldtwall);
	phiwall=phi;
	storeState();

	position[0]+=velocity[0]*dtwall;
    position[1]+=velocity[1]*dtwall;
	position[2]+=velocity[2]*dtwall-0.5*gravity*dtwall*dtwall;
	velocity[2]+=-gravity*dtwall;


	
	if(whatwall==2)
	{	
			//might be necessary to clamp to the wall if it hits the corner. 
   			//position[2]=position[2]>0? length/2.0:-length/2.0;
   			
		if(rand_flat()>diffuse)
		{//specular
			velocity[2]=-velocity[2];
		}
		else
		{//diffuse bounce on the ends. whats the cdf of cos(theta)^2*sin(theta) again? theta0=arccos((1-rand)^(1/3))
		 //ugg someone from boost should write a diffuse ditrbituion code, (they do have an arcsine dist) this is slow. 
		//maybe template this in a spline for extra speed? lets see how bad it is. 
			//method one.
			double vmag=std::sqrt(velocity[2]*velocity[2]+velocity[0]*velocity[0]+velocity[1]*velocity[1]);
			double temptheta=std::acos(boost::math::cbrt(rand_flat()));
			double tempphi=2.0*pi*rand_flat();
		    velocity[2]=-vmag*sgn(position[2])*std::cos(temptheta);
		    velocity[0]=vmag*std::sin(temptheta)*std::cos(tempphi);
		    velocity[1]=vmag*std::sin(temptheta)*std::sin(tempphi);
		  
		}
	}
	else if(whatwall==1)
	{	
			//might be necessary to clamp to the wall if it hits the corner. 
   			//position[2]=position[2]>0? length/2.0:-length/2.0;
   			
		if(rand_flat()>diffuse)
		{//specular
			
			velocity[1]=-velocity[1];
		}
		else
		{//diffuse bounce on the ends. whats the cdf of cos(theta)^2*sin(theta) again? theta0=arccos((1-rand)^(1/3))
		 //ugg someone from boost should write a diffuse ditrbituion code, (they do have an arcsine dist) this is slow. 
		//maybe template this in a spline for extra speed? lets see how bad it is. 
			//method one.
			double vmag=std::sqrt(velocity[2]*velocity[2]+velocity[0]*velocity[0]+velocity[1]*velocity[1]);
			double temptheta=std::acos(boost::math::cbrt(rand_flat()));
			double tempphi=2.0*pi*rand_flat();
		    velocity[1]=-vmag*sgn(position[1])*std::cos(temptheta);
		    velocity[0]=vmag*std::sin(temptheta)*std::cos(tempphi);
		    velocity[2]=vmag*std::sin(temptheta)*std::sin(tempphi);
		  
		}
	} 
	else
	{	
			//might be necessary to clamp to the wall if it hits the corner. 
   			//position[2]=position[2]>0? length/2.0:-length/2.0;
   			
		if(rand_flat()>diffuse)
		{//specular
			
			velocity[0]=-velocity[0];
		}
		else
		{//diffuse bounce on the ends. whats the cdf of cos(theta)^2*sin(theta) again? theta0=arccos((1-rand)^(1/3))
		 //ugg someone from boost should write a diffuse ditrbituion code, (they do have an arcsine dist) this is slow. 
		//maybe template this in a spline for extra speed? lets see how bad it is. 
			
			double vmag=std::sqrt(velocity[2]*velocity[2]+velocity[0]*velocity[0]+velocity[1]*velocity[1]);
			double temptheta=std::acos(boost::math::cbrt(rand_flat()));
			double tempphi=2.0*pi*rand_flat();
		    velocity[0]=-vmag*sgn(position[0])*std::cos(temptheta);
		    velocity[2]=vmag*std::sin(temptheta)*std::cos(tempphi);
		    velocity[1]=vmag*std::sin(temptheta)*std::sin(tempphi);
		  
		}
	}

//should have new velocity[] and position 
	//need to update time
	
	time=twall;
	//maybe we will hit the wall again before dt is over? lets go to infinity and find out!
	this->nextBoundary();
	this->moveParticle(dt-dtwall);
	
	return;//exiting from recursive calls.

}

//sphere is suck to the max. (exceptions for cylinders along x and y), enter at your own risk.  
void Sphere::startInside()
{
 

   position[2]=getIsotropic(-radius,radius);
  
  
   do 
   {
   position[0]=radius*(2.*rand_flat()-1.);
   position[1]=radius*(2.*rand_flat()-1.);
   } //should be self consistent enough to eventually find a good starting state. presumably there is a better way. 
   while(position[0]*position[0]+position[1]*position[1]+position[2]*position[2]>radius*radius);

	for(int i=0; i<3; i++) oldposition[i]=position[i];
	
	//starting velocity distribution while  accouting for gravity. 
	 //(m/s) 

	for(int i=0; i<3; i++) oldvelocity[i]=velocity[i];

	spin[0]=1.;
	spin[1]=0.;
	spin[2]=0.;
	phiwall=0.;

	xmax=2.000001*std::sqrt(4.*gravity*radius/vmag);
	//this can probably be better with some thought....
	xguess=std::sqrt(4.*gravity*radius/vmag);

	for(int i=0; i<3; i++) oldspin[i]=spin[i];
	
	//taken out collision type so I need to be expliciit. 
	double c4 = 0.25*gravity*gravity;
	double c3 = -velocity[2]*gravity;
	double c2 = (velocity[0]*velocity[0]+velocity[1]*velocity[1]+velocity[2]*velocity[2]-gravity*position[2]);
	double c1 = 2.*velocity[0]*position[0]+2.*velocity[1]*position[1]+2.*velocity[2]*position[2];
	double c  = position[0]*position[0]+position[1]*position[1]+position[2]*position[2]-radius*radius;
	
    
    twall=quartic_root(c,c4,c3,c2,c1,xguess,0.,xmax);
	double x=position[0]+velocity[0]*twall;
    double y=position[1]+velocity[1]*twall;
	double z=position[2]+velocity[2]*twall-0.5*gravity*twall*twall;

	if(isnan(twall) || (x*x+y*y+z*z-radius*radius)<2E-6 || (x*x+y*y+z*z-radius*radius)>-2E-6){// our guess was so bad that everyone died. 
		for(int i = 0; i<101;i++)
		{
			twall=quartic_root(c,c4,c3,c2,c1,(double)i*xmax,0.,xmax);
			//checking to see if it found the correct root. 
			x=position[0]+velocity[0]*twall;
			y=position[1]+velocity[1]*twall;
			z=position[2]+velocity[2]*twall-0.5*gravity*twall*twall;

	    if(!isnan(twall) && (x*x+y*y+z*z-radius*radius)<2E-6 && (x*x+y*y+z*z-radius*radius)>-2E-6) return;
		 	
		}
		throw std::invalid_argument("couldn't find a valid starting state: something is probably wrong with the radius or vmax geometry parameter.\n");
	}
	
	oldtwall=0.;
	return;
}


void Sphere::nextBoundary()
{
	
//need to solve for minimum positive root of:
 // 1/4*g^2*t^4-vz*g*t^3+(vx^2+vy^2+vz^2-g*z)*t^2+(2*vx*x+2*vy*y+2*vz*z)*t +x^2 +y^2 +z^2 -R^2 = 0 
//         c4        c3           c2                    c1                            c
	//on a wall only need to evaluate cube root as one root is zero.  
	double c3 = 0.25*gravity*gravity;
	double c2 = -velocity[2]*gravity;
	double c1 = (velocity[0]*velocity[0]+velocity[1]*velocity[1]+velocity[2]*velocity[2]-gravity*position[2]);
	double c = 2.0*velocity[0]*position[0]+2.0*velocity[1]*position[1]+2.0*velocity[2]*position[2];
   
	//guess is the solution without gravity...
	xguess=mag(2.*(velocity[0]*position[0]+velocity[1]*position[1]+velocity[2]*position[2])/(c1+gravity*position[2]));

	
	twall=cubic_root(c,c3,c2,c1,xguess,0.,xmax)+time;
	
	if(std::isnan(twall)){
    	throw std::invalid_argument("encountered nan value, sorry :.(");
    	
    }
  	
  	//check if it found an imaginary root(this solver is fast, but it gets confused near slightly imaginary roots)
	double x= position[0]+velocity[0]*(twall-time);
    double y= position[1]+velocity[1]*(twall-time);
    double z= position[2]+velocity[2]*(twall-time)-0.5*gravity*(twall-time)*(twall-time);
    
    

    //bad solution if root corrisponds to more than ~1 micron away from the radius . 
    if((x*x+y*y+z*z-radius*radius)>=2E-6 || (x*x+y*y+z*z-radius*radius)<=-2E-6){
    for(int i=0; i<101; i++)
    	{
    	//typically this loop  is run when the neutron reaches the top of its parabolic trajectory, the solver returns this point as the solution,
    	//this loop forces the solver to look from beyond the boundary to the intersection. 
    	xguess=xmax*(1.-(double)i/100.);
    
    	twall=cubic_root(c,c3,c2,c1,xguess,0.0,xmax)+time;
    	
    	x=position[0]+velocity[0]*(twall-time);
    	y= position[1]+velocity[1]*(twall-time);
    	z= position[2]+velocity[2]*(twall-time)-0.5*gravity*(twall-time)*(twall-time);
    	
    	
    	

    	if((x*x+y*y+z*z-radius*radius)<2E-6 && (x*x+y*y+z*z-radius*radius)>-2E-6) return;

    	}
    	//diffuse bounce is not allowed (theta too close to pi/2 (UPDATE,shouldn't happen anymore) and exits or remains on the boundary to numerical precision.). 
    	twall=time;

        //throw std::invalid_argument("invalid boundary collision, sorry :.(");
	}
    
  	
	return;
}

void Sphere::moveParticle(double dt)
{

	if(time+dt<twall)
	{      
		position[0]+=velocity[0]*dt;
    	position[1]+=velocity[1]*dt;
		position[2]+=velocity[2]*dt-0.5*gravity*dt*dt;
		velocity[2]+=-gravity*dt;
		time+=dt;
				//recursive until no wall is hit. 
		return;//the actual exit is here at the beginning, whoa. 
	
	}
		//this bit will become more complicated if we switch to RK. 
		analyticSolution(twall-oldtwall);
		phiwall=phi;
		double dtwall = twall-time;
		time=twall;
		storeState();//twall becomes old twall, 

		
		position[0]+=velocity[0]*dtwall;
    	position[1]+=velocity[1]*dtwall;
		position[2]+=velocity[2]*dtwall-0.5*gravity*dtwall*dtwall;
		velocity[2]=velocity[2]-gravity*dtwall;
		
	
		//need to clamp position to the sphere. it could be up to 1 micron off which will screw up the solver. 
    	double inv_rtemp=1./std::sqrt(position[0]*position[0]+position[1]*position[1]+position[2]*position[2]);
    	for(int i=0; i<3; i++) position[i]=position[i]*inv_rtemp*radius;

    	
        //needed for converting to and from cartesian.
		double rho= std::sqrt(position[0]*position[0]+position[1]*position[1]);
 	if(rand_flat()>diffuse)
 	{//specular. 
 		
 		//radial velocity made negative
 		double vr=-(position[0]*velocity[0]/radius + position[1]*velocity[1]/radius + position[2]*velocity[2]/radius);  
 		
 		double vt=position[2]/radius*position[0]/rho*velocity[0] + position[2]/radius*position[1]/rho*velocity[1]-rho/radius*velocity[2];
 		
 		double vp= - position[1]/rho*velocity[0]+position[0]/rho*velocity[1];

 		velocity[0]=position[0]/radius*vr+position[2]/radius*position[0]/rho*vt-position[1]/rho*vp;
 		velocity[1]=position[1]/radius*vr+position[2]/radius*position[1]/rho*vt+position[0]/rho*vp;
 		velocity[2]=position[2]/radius*vr-rho/radius*vt;
 	}
 	else
 	{//diffuse
 		
		double vmag=std::sqrt(velocity[0]*velocity[0]+velocity[1]*velocity[1]+velocity[2]*velocity[2]);
		double temptheta=std::acos(boost::math::cbrt(rand_flat()));
		//don't allow obnoxious theta values that can escape. (helps with solver stability);
		temptheta=temptheta<(pi/2.0-1E-4)?temptheta:pi/2.0-1E-4;
		
		double tempphi=2.0*pi*rand_flat();   
		double vr=-vmag*std::cos(temptheta);
		double vt=vmag*std::sin(temptheta)*std::cos(tempphi);
		double vp =vmag*std::sin(temptheta)*std::sin(tempphi);

		velocity[0]=position[0]/radius*vr+position[2]/radius*position[0]/rho*vt-position[1]/rho*vp;
 		velocity[1]=position[1]/radius*vr+position[2]/radius*position[1]/rho*vt+position[0]/rho*vp;
 		velocity[2]=position[2]/radius*vr-rho/radius*vt;
 		
 	}

 	
	
	//maybe we will hit the wall again before dt is over? lets go to infinity and find out!
	this->nextBoundary();
	this->moveParticle(dt-dtwall);


	return;//exiting from recursive calls.


}
