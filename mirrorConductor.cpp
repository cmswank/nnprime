//g++ cmsRK.cpp cmsTransport.cpp cmsBoundary.cpp -I/data1/cmswank/BoostCodeSwank/boost_1_64_0 -o Trajtest -lboost_thread -lboost_system
//use above to compile and/or use to generate a makefile. 
#include "particleTrajectory.h"
#include "mirrorConductor.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <boost/thread.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>

mirrorConductor::mirrorConductor(std::string infilename, std::string ofilename)
{

	this->parseParameters(infilename);

  outfilename=ofilename;

	if(this->rkdt>this->dt) this->rkdt=this->dt/5.0; //incase dt is mixed up for some reason
	//this turns out to be impossible in c++. 
	//make dt divisable by rkdt. not good for adaptive solver, but this is exacat solution so no problem. 
	//if(boost::math::fmod(dt,rkdt)!=0.) rkdt=rkdt+boost::math::fmod(dt,rkdt)/((double)((int)(dt/rkdt)));



	
  this->timesteps=(int)(T/dt)+1;

	
	if(N%cores!=0)
	{
		this->Nactual=N-N%cores;
		std::cout<<"Warning: the number of spins has been changed to "<<Nactual<<std::endl;
		std::cout<<"If you want to have exactly N make mod(N,cores)=0"<<std::endl;

	}
	else Nactual=N;
	//prepare file header. (3 doubles that are really ints but matlab and python can't tell the difference)
	ts=(double)this->timesteps;
	nvars=11.0;
	Nsave=(double)Nactual;
	

}



void mirrorConductor::Conduct()
{
	//Memory Heavy, kinda...
	Spins=(double*)malloc(Nactual*timesteps*3*sizeof(double));		
	Position=(double*)malloc(Nactual*timesteps*3*sizeof(double));			
	Velocity=(double*)malloc(Nactual*timesteps*3*sizeof(double));					
	//BField=(double*)malloc(Nactual*timesteps*3*sizeof(double));
	Phase=(double*)malloc(Nactual*timesteps*sizeof(double));
	std::string tempofilename;
  //for loop for scanning parameter, if any. 
  for(int parnum=0;parnum<=scannumber;parnum++){
    //naming for scanning, again if any. 
    if(scanbool) tempofilename="./data/"+outfilename+"_"+std::to_string(parnum);
    else  tempofilename="./data/"+outfilename;
  *scanvar=scanstart+(double)parnum*(scanstop-scanstart)/((double)scannumber);
  this->spinfile.open (tempofilename.c_str(),std::ios::binary | std::ios::out);
  
  spinfile.write(reinterpret_cast <const char*> (&ts),sizeof(double));
  spinfile.write(reinterpret_cast <const char*> (&Nsave),sizeof(double));
  spinfile.write(reinterpret_cast <const char*> (&nvars),sizeof(double));


	
	for(int i = 0; i<N; i++)
		{
			//cast particleTrajectory to its child geometry. based on what we want (or at least what was entered). 
			switch (geom)
			{		
    			case 0:  {Cylinder* cyl = new Cylinder(length,radius,vmax,diffuse,startseed+i);
    					  tr.push_back(dynamic_cast<Cylinder*>(cyl));
    					  tr[i]->nnn=i;} 
        		break;

    			case 1:  {Rectangle* rect= new Rectangle(xwidth,ylength,zheight,vmax,diffuse,startseed+i);
    					  tr.push_back(dynamic_cast<Rectangle*>(rect));
    					  tr[i]->nnn=i;} 
        		
        		break;

    			case 2:  {Sphere* sph= new Sphere(radius,vmax, diffuse, startseed+i);
    					  tr.push_back(dynamic_cast<Sphere*>(sph));
    					  tr[i]->nnn=i;} 
				break;
				//blocked this in the parameter parsing, it wont get here. 
				default: {Cylinder* cyl = new Cylinder(length,radius,vmax, diffuse, startseed+i);
    					  tr.push_back(dynamic_cast<Cylinder*>(cyl));
    					  tr[i]->nnn=i;} 
			}

			

		}
		if (parnum==0) std::cout<<"Run parameters were loaded successfully, running now\n"<<std::flush;
	for(int ii=0; ii<Nactual; ii+=cores)
	{

		for(int i =ii; i<ii+cores; i++)
		{	
				bthreads.push_back(new boost::thread(&mirrorConductor::run,this,tr[i]));
		}

		for(int i = 0; i<cores; i++)
		{	
			
			bthreads[i]->join();

		}


    
		bthreads.erase(bthreads.begin(), bthreads.end());

    if(ii%cores==0) std::cout<<" :) "<<std::flush;
		else if(ii%82==0) std::cout<<" :( "<<std::flush;
		else if(ii%131==0) std::cout<<" (_8(|) "<<std::flush; 
	
	}

	std::cout<<"\n"<<"Saving\n";
	this->saveSpins();
  tr.erase(tr.begin(), tr.end()); 
  }
	std::cout<<"exiting\n";
	return;
}


void mirrorConductor::saveSpins()
{
	double SSX,SSY,SSZ,XX,YY,ZZ,VVX,VVY,VVZ,TT,PPHASE;//,BBX,BBY,BBZ,EEX,EEY,EEZ,
	int numtemp,numtemp1;
	
	for(int i = 0; i<Nactual; i++)
	{
		for(int ii = 0; ii<timesteps; ii++)
		{
		//flattened array coordinate 
		numtemp=i*timesteps*3+3*ii;
    numtemp1=i*timesteps+ii;
		TT=dt*(double)(ii);
		SSX=Spins[numtemp];
		SSY=Spins[numtemp+1];
		SSZ=Spins[numtemp+2];
		PPHASE=Phase[numtemp1];

		XX=Position[numtemp];
		YY=Position[numtemp+1];
		ZZ=Position[numtemp+2];
		

		VVX=Velocity[numtemp];
		VVY=Velocity[numtemp+1];
		VVZ=Velocity[numtemp+2];	
		//BBX=BField[numtemp];
		//BBY=BField[numtemp+1];
		//BBZ=BField[numtemp+2];
		//EEX=EField[numtemp];
		//EEY=EField[numtemp+1];
		//EEZ=EField[numtemp+2];

		spinfile.write(reinterpret_cast <const char*> (&SSX),sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&SSY),sizeof(double)); 		
		spinfile.write(reinterpret_cast <const char*> (&SSZ),sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&PPHASE), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&XX), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&YY), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&ZZ), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&VVX), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&VVY), sizeof(double)); 
		spinfile.write(reinterpret_cast <const char*> (&VVZ), sizeof(double));
		//spinfile.write(reinterpret_cast <const char*> (&BBX), sizeof(double)); 
		//spinfile.write(reinterpret_cast <const char*> (&BBY), sizeof(double)); 
		//spinfile.write(reinterpret_cast <const char*> (&BBZ), sizeof(double)); 
		//spinfile.write(reinterpret_cast <const char*> (&EEX), sizeof(double)); 
		//spinfile.write(reinterpret_cast <const char*> (&EEY), sizeof(double)); 
		//spinfile.write(reinterpret_cast <const char*> (&EEZ), sizeof(double));  
		spinfile.write(reinterpret_cast <const char*> (&TT), sizeof(double));

		}
	}
	spinfile.close();
   	std::cout<<"Saving probably successful.\n";
	
	return;//format file
 	
}


void mirrorConductor::run(particleTrajectory* trp)
{ 
	

  //set up the rest of the trajectory
	trp->startInside();

	//load mirror fields. 
	//loadFields(double B0, double B0m, double tempthetam, double tempphim, double eps_temp)
	trp->loadFields(B0, B0m, thetam, phim, eps);


	//2D array starting place in flattened 3D array
	int n = trp->nnn*timesteps*3;//?
  int n1 =trp->nnn*timesteps;
	double temp_t;
	int numtemp,numtemp1;

	
	    Spins[n]=trp->spin[0]; 
    	Spins[n+1]=trp->spin[1]; 
    	Spins[n+2]=trp->spin[2];
    	Phase[n1]=trp->phi;


    	Position[n]=trp->position[0];
    	Position[n+1]=trp->position[1];
    	Position[n+2]=trp->position[2];
    	


    	Velocity[n]=trp->velocity[0];
    	Velocity[n+1]=trp->velocity[1];
    	Velocity[n+2]=trp->velocity[2];

    	//BField[n]=trp->bfield[0];
    	//BField[n+1]=trp->bfield[1];
    	//BField[n+2]=trp->bfield[2];
    

	for(int i = 1 ; i < (int)(T/dt+1); i++)
    {
    	numtemp=n+3*i;
      numtemp1=n1+i;
    	temp_t=dt*i;
    	
    	//move to the next dt. 
    	for(int ii=0; ii<(int)(dt/rkdt); ii++)
		{
		
    trp->moveParticle(rkdt);
		
      
    	}
     trp->integrateMP();
		//Save point in time
    	Spins[numtemp]=trp->spin[0]; 
    	Spins[numtemp+1]=trp->spin[1]; 
    	Spins[numtemp+2]=trp->spin[2];
      Phase[numtemp1]=trp->phi;
    	
    	Position[numtemp]=trp->position[0];
    	Position[numtemp+1]=trp->position[1];
    	Position[numtemp+2]=trp->position[2];
    	
    	Velocity[numtemp]=trp->velocity[0];
    	Velocity[numtemp+1]=trp->velocity[1];
    	Velocity[numtemp+2]=trp->velocity[2];
    	
          
    	//BField[numtemp]=trp->bfield[0];
    	//BField[numtemp+1]=trp->bfield[1];
    	//BField[numtemp+2]=trp->bfield[2];

    }

   
    return;
}



bool mirrorConductor::parseParameters(const std::string& name) {
    using namespace std;

    int pnum;
    int pcheck=0;

    if(FileExists(name.c_str())) {
      bool scannow=false;
    	cout<<"loading parameter file"<<"\n";
    	ifstream myfile(name.c_str());
    	string paramin;
    	

    	while( getline(myfile,paramin) ){
    	
    	size_t exclude=paramin.find("//");
    	
    	if((int)exclude>=0) paramin=paramin.substr(0,exclude);
   	    
   	    size_t startpos=paramin.find("=");
   	    
   	    size_t endpos=paramin.find(";");
        
        size_t scan=paramin.find(",");
        

   	    if((int)startpos>0 && (int)endpos>(int)startpos)
   	    {   endpos=endpos-startpos;//convert to more logical substring behavior. 
			       
             if((int)scan>-1){
              string scancheck=paramin.substr(scan+1,endpos-1);

              size_t scan2=scancheck.find(",");
              if((int)scan2>-1){
                if(scanbool==true)throw invalid_argument("Sorry, only scanning for one parameter is allowed.\nMay I suggest writing a python or bash script?\n");
                scanbool=true;
                scanstart=atof(paramin.substr(startpos+1,scan-1-startpos).c_str());
                scanstop=atof(paramin.substr(scan+1,scan2).c_str());
                scannumber=atoi(paramin.substr(scan+scan2+2,scan+endpos-2).c_str());
                scannow=true;
                
              }else throw invalid_argument("scan must have a start value, stop value, and number of runs (in that order.)"); 



             }

           //Geometry input, 
   	    	size_t Geometry = paramin.find("geometry");
   	    	if((int)Geometry>-1){
   	    		geom=atoi(paramin.substr(startpos+1,endpos-1).c_str());
   	    		if(geom==0){
   	    			pnum=16;
   	    			pcheck++;
              cout<<"Cylinder Boundary\n";
   	    		}
   	    		else if(geom==1) {pnum=17; pcheck++; cout<<"Rectangle Boundary\n";}
   	    		else if(geom==2){pnum=15; pcheck++; cout<<"Sphere Boundary\n";}
   	    		else{ throw invalid_argument("invalid geometry specification\n"); return false;}
   	    	}

   	    	//radius input. 
   	    	size_t Radius=paramin.find("radius");
   	    	if((int)Radius>-1){
   	    		
   	    	 	if(scannow){ scanvar=&radius;
               cout<<"scanning radius from = "<<scanstart<<" m to "<<scanstop<<" m with "<<scannumber <<" entries\n";
               scannow=false;
            }else {radius=atof(paramin.substr(startpos+1,endpos-1).c_str());
                   cout<<"radius = "<<radius<<" m\n";
            }
            pcheck++;
   	    	}

   	    	//length input. 
   	    	size_t Length=paramin.find("length");
   	    	if((int)Length>-1){

              //y length input. 
          size_t Ylength=paramin.find("ylength");
          if((int)Ylength>-1){
            
            if(scannow){ scanvar=&ylength;
               cout<<"scanning y length from "<<scanstart<<" m to "<<scanstop<<" m with "<<scannumber <<" entries\n";
               scannow=false;
            }else {ylength=atof(paramin.substr(startpos+1,endpos-1).c_str());
                   cout<<"y length = "<<ylength<<" m\n";
            }
            pcheck++;
          }else{
   	    		if(scannow){ scanvar=&length;
               cout<<"scanning length from "<<scanstart<<" m to "<<scanstop<<" m with "<<scannumber <<" entries\n";
               scannow=false;
            }else {length=atof(paramin.substr(startpos+1,endpos-1).c_str());
                   cout<<"length = "<<length<<" m\n";
            }
            pcheck++;
   	    	}
          }

   	    	//z height input. 
   	    	size_t Zheight=paramin.find("zheight");
   	    	if((int)Zheight>-1){
            if(scannow){ scanvar=&zheight;
               cout<<"scanning from z height from "<<scanstart<<" m to "<<scanstop<<" m with "<<scannumber <<" entries\n";
               scannow=false;
            }else {zheight=atof(paramin.substr(startpos+1,endpos-1).c_str());
                   cout<<"z height = "<<zheight<<" m\n";
            }
            pcheck++;
   	    	}

                    //z height input. 
          size_t Vmax=paramin.find("vmax");
          if((int)Vmax>-1){
            if(scannow){ 
               scanvar=&vmax;
               cout<<"scanning from v max from "<<scanstart<<" m/s to "<<scanstop<<" m/s with "<<scannumber <<" entries\n";
               scannow=false;
            }else {vmax=atof(paramin.substr(startpos+1,endpos-1).c_str());
                   cout<<"v max = "<<vmax<<" m/s\n";
            }
            pcheck++;
          }

   	    	size_t Diffuse=paramin.find("diffuse");
   	    	if((int)Diffuse>-1){
            if(scannow){ 
               scanvar=&diffuse;
               cout<<"scanning from diffuse bounce probability from "<<scanstart<<" to "<<scanstop<<" with "<<scannumber <<" entries\n";
               scannow=false;
            }else {diffuse=atof(paramin.substr(startpos+1,endpos-1).c_str());
                   cout<<"diffuse bounce probabililty = "<<diffuse<<"\n";
            } 
            pcheck++;
   	    	}

   	    	//Seed input. 
   	    	size_t Seed=paramin.find("seed");
   	    	if((int)Seed>-1){
   	    		startseed=atoi(paramin.substr(startpos+1,endpos-1).c_str());
   	    		pcheck++;
            cout<<"Starting random number seed = "<<startseed<<"\n";
   	    	} 

   	    	//Numper of particles input. 
   	    	size_t NUM=paramin.find("N");
   	    	if((int)NUM>-1){
   	    		N=atoi(paramin.substr(startpos+1,endpos-1).c_str());
   	    		pcheck++;
            cout<<"number of neutrons simulated "<<N<<"\n";
   	    	}

            //Total Time of sim
   	    	size_t Time=paramin.find("T");
   	    	if((int)Time>-1){
   	    		T=atof(paramin.substr(startpos+1,endpos-1).c_str());
            pcheck++;
   	    	  cout<<"Total time simulated "<<T<<" s\n"; 
          }
   	    	
   	    	//output time steps
   	    	size_t T_step=paramin.find("dt");
   	    	if((int)T_step>-1){
                size_t Xwidth=paramin.find("xwidth");
                size_t S_step=paramin.find("rkdt");
          if((int)S_step>-1){
            rkdt=atof(paramin.substr(startpos+1,endpos-1).c_str());
            pcheck++;
            cout<<"Solver time steps "<<rkdt<<" s\n";
          }else if ((int)Xwidth>-1){
                      //x width input. 
           if(scannow){ scanvar=&xwidth;
               cout<<"scanning x width from "<<scanstart<<" m to "<<scanstop<<" m with "<<scannumber <<" entries\n";
               scannow=false;
            }else {xwidth=atof(paramin.substr(startpos+1,endpos-1).c_str());
                   cout<<"x width = "<<xwidth<<" m\n";
            }
            pcheck++;
          }
            else{
   	    		dt=atof(paramin.substr(startpos+1,endpos-1).c_str());
   	    		pcheck++;
            cout<<"Output time steps "<<dt<<" s\n";
   	    	}
          }

   	    	//Solver time steps
   	  

   	    	//cores
   	    	size_t Cores=paramin.find("cores");
   	    	if((int)Cores>-1){
   	    		cores=atoi(paramin.substr(startpos+1,endpos-1).c_str());
   	    		pcheck++;
   	    	  cout<<"Number of cores "<<cores<<"\n";
          }


   	    	//fields
   	    	size_t b0=paramin.find("B0");
   	    	if((int)b0>-1){
   	    		
          size_t b0m=paramin.find("B0m");
          if((int)b0m>-1){
            //mirror field
           if(scannow){ scanvar=&B0m;
               cout<<"scanning B0m from "<<scanstart<<" T to "<<scanstop<<" T with "<<scannumber <<" entries\n";
               scannow=false;
            }else {B0m=atof(paramin.substr(startpos+1,endpos-1).c_str());
                   cout<<"mirror field strength = "<<B0m<<" T\n";
            }
            pcheck++;
          }
          else
          { //lab field
           if(scannow){ scanvar=&B0;
               cout<<"scanning B0 from "<<scanstart<<" T to "<<scanstop<<" T with "<<scannumber <<" entries\n";
               scannow=false;
            }else {B0=atof(paramin.substr(startpos+1,endpos-1).c_str());
                   cout<<"lab field B0 strength (always along z) = "<<B0<<" T\n";
            }
            pcheck++;
   	    	}
         }
   	    	

			size_t Thetam=paramin.find("thetam");
   	    	if((int)Thetam>-1){
   	    		if(scannow){ scanvar=&thetam;
               cout<<"scanning mirror field theta  from = "<<scanstart<<" radians to "<<scanstop<<" radians with "<<scannumber <<" entries\n";
               scannow=false;
            }else {thetam=atof(paramin.substr(startpos+1,endpos-1).c_str());
                   cout<<"mirror field theta = "<<thetam<<" radians\n";
            }
            pcheck++;
   	    	}		


			size_t Phi_m=paramin.find("phim");
   	    	if((int)Phi_m>-1){
   	    		if(scannow){ 
               scanvar=&phim;
               cout<<"scanning mirror field phi from = "<<scanstart<<" radians to "<<scanstop<<" radians with "<<scannumber <<" entries\n";
               scannow=false;
            }else {phim=atof(paramin.substr(startpos+1,endpos-1).c_str());
                   cout<<"mirror field phi = "<<phim<<" radians\n";
            }
            pcheck++;
   	    	}		

   	    	size_t Eps=paramin.find("eps");
   	    	if((int)Eps>-1){
            if(scannow){ 
              scanvar=&eps;
               cout<<"scanning mirror field coupling from = "<<scanstart<<" radians/s to "<<scanstop<<" radians/s with "<<scannumber <<" entries\n";
               scannow=false;
            }else {eps=atof(paramin.substr(startpos+1,endpos-1).c_str());
                   cout<<"mirror field coupling = "<<eps<<" radians/s\n";
            }
            pcheck++;
   	    	}

          if(scannow==true) throw invalid_argument("scanning not permitted on the requested parameter\n");

		}


		}

    myfile.close();
		cout<<"parameter file closed, aquired "<<pcheck<<" parameters\n";
  


    if(pcheck==pnum) return true;
		else if(pcheck>pnum) {
			cout<<"Warning: unused parameters. \n";
			return true;
		}
		else {
			throw invalid_argument("incorrect parameter file, insufficient number of entries");
			return false;
		}

	}
  else return false;
}

bool mirrorConductor::FileExists(const std::string& name) {
    std::ifstream f(name.c_str());
    if (f.good()) {
        f.close();
        return true;
    } else {
        f.close();
        return false;
    }   
}


