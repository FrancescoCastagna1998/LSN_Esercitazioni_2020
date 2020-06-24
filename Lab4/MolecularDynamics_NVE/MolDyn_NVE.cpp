/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 


  Input(); // the first time the simulation run from a provided configuration fcc
  
  if(restart==0)
  {
	int nconf = 1;
	for(int istep=1; istep <= nstep; ++istep){
     		Move();           //Move particles with Verlet algorithm
  		if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
  		if(istep%10 == 0){
        		Measure();     //Properties measurement
			Equilibration(istep);  // Print the istant values of Temperature       
			//ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
			nconf += 1;
     		}
	}
  
  	ConfFinal();         //Write final configuration to restart
  	ConfFinalOld();      //Write old final configuration before the last one to restart
   }
   if(restart==1)
   {
	InputRestart(); //Restart
  	int nconf = 1;
  	for(int istep=1; istep <= nstep; ++istep){
     			Move();           //Move particles with Verlet algorithm
  			if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
  			if(istep%10 == 0){
        			Measure();     //Properties measurement  
				Equilibration(istep);  // Print the istant values of Temperature       
				//ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
				nconf += 1;
     			}
	}
	ConfFinal();         //Write final configuration to restart
  	ConfFinalOld();      //Write old final configuration before the last one to restart
   }

   if(restart==2)
   {
	InputFinal(); // start the  real measure
	int nconf = 1;
	for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
	{
		Reset(iblk);   //Reset block averages
    		for(int istep=1; istep <= nstep; ++istep){
    			Move();           //Move particles with Verlet algorithm
        		Measure();     //Properties measurement
			Accumulate(); //Update block averages
			if(istep%10 == 0){
				// ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        			nconf += 1;
     			}
  		}
		Averages(iblk);   //Print results for current block
	}
  	ConfFinal();         //Write final configuration to restart
  	ConfFinalOld();      //Write old final configuration before the last one to restart
    }

   return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> restart;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl << endl;

  if(restart==0){
	cout << "The program perform the first necessary run " << endl;
	nstep = nstep*10; // steps number for the first run = 10000;
	cout << "Number of steps = " << nstep << endl << endl; //
  }

  if(restart==1){
	cout << "The program perform the equilibration of the system" << endl;
	cout << "Number of steps to equilibrate the system = " << nstep << endl << endl; // steps number for the equilibration =1000
  }
		
  if(restart==2){
	cout << "The program perform the real measurement" << endl;
	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of steps for each block = " << nstep << endl << endl;
  }

  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature

  n_props = 4; //Number of observables

//measurement of g(r)
  igofr = 4;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;


  if (restart==0)
  {
//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
//Evaluate energy etc. of the initial configuration
   Measure();

//Print initial values for the potential energy 
   cout << "Initial potential energy = " << walker[iv]/(double)npart << endl;
   cout << "Initial kinetic energy = " << walker[ik]/(double)npart << endl;
   cout << "Initial temperature = " << walker[it]/(double)npart << endl;
   cout << "Initial total energy = " << walker[ie]/(double)npart << endl;
   }

   return;
   
}

void InputRestart(void){ 

  //Read initial configuration
  ifstream ReadConf;
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Read old initial configuration
  cout << "Read old initial configuration from file old.0 " << endl << endl;
  ifstream ReadConfOld;
  ReadConfOld.open("old.0");
  for (int i=0; i<npart; ++i){
    ReadConfOld >> xold[i] >> yold[i] >> zold[i];
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  ReadConfOld.close();

  double t=0.; 
  Move(); // compute r(t+dt) with one step of Verlet algorithm 
  for (int i=0; i<npart; ++i){ // compute v(t+dt/2)
  	vx[i] = Pbc(x[i]-xold[i])/delta;
	vy[i] = Pbc(y[i]-yold[i])/delta;
	vz[i] = Pbc(z[i]-zold[i])/delta;
	t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
  }
 
 
  cout << " Trying to rescale the velocities in order to match the temperature Target. " << endl << endl;
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature (t+dt/2)
  cout << "La temperatura stimata in partenza: " << stima_temp << endl;
  cout << "La temperatura attesa: " << temp << endl;
  
  double fs = sqrt(temp/stima_temp); // velocity scale factor 
  
  cout << "The velocity scale factor is: " << fs << endl;
  for (int i=0; i<npart; ++i){
     	vx[i] *= fs;
     	vy[i] *= fs;
     	vz[i] *= fs;

     	xold[i] = Pbc(x[i] - vx[i] * delta);
     	yold[i] = Pbc(y[i] - vy[i] * delta);
     	zold[i] = Pbc(z[i] - vz[i] * delta);
  }
   
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy 
  cout << "Initial potential energy = " << walker[iv]/(double)npart << endl;
  cout << "Initial kinetic energy = " << walker[ik]/(double)npart << endl;
  cout << "Initial temperature = " << walker[it]/(double)npart << endl;
  cout << "Initial total energy = " << walker[ie]/(double)npart << endl;
}

void InputFinal(void) {
  //Read initial configuration
  cout << "Read initial configuration from file config.0" << endl << endl;
  cout << "Read old initial configuration from file old.0 " << endl << endl;
  ifstream ReadLastConf,ReadLastConfOld;
  ReadLastConf.open("config.0");
  ReadLastConfOld.open("old.0");
  for (int i=0; i<npart; ++i){
    ReadLastConf >> x[i] >> y[i] >> z[i];
    ReadLastConfOld >> xold[i] >> yold[i] >> zold[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  Move();
  for (int i=0; i<npart; ++i) {
    vx[i] = Pbc(x[i] - xold[i])/delta;
    vy[i] = Pbc(y[i] - yold[i])/delta;
    vz[i] = Pbc(z[i] - zold[i])/delta;
    xold[i]=Pbc(x[i]-vx[i]*delta);
    yold[i]=Pbc(y[i]-vy[i]*delta);
    zold[i]=Pbc(z[i]-vz[i]*delta);
  }
  return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;

  v = 0.0; //reset observable
  t = 0.0;

//reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;


//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

//update of the histogram of g(r)
     for (bin=igofr; bin<n_props; bin++){
		if (dr>(bin-igofr)*bin_size  && dr<(bin-(igofr-1))*bin_size){
			walker[bin] = walker[bin] + 2;
		}

     }

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    walker[iv]= v; //Potential energy 
    walker[ik] = t; //Kinetic energy 
    walker[it] = (2.0 / 3.0) * t; //Temperature
    walker[ie] = (t+v); //Total energy 

    return;
}

void Equilibration(int istep)
{
   ofstream Temp;
   const int wd=12;

   if(restart==0)
   {	
	Temp.open("output_temp_istant_Argon_gas.dat",ios::app);
	Temp << setw(wd) << istep <<  setw(wd) << walker[it]/(double)npart << endl;
	Temp.close();
   }

   if(restart==1)
   {	
	Temp.open("output_temp_istant_equi_Argon_gas.dat",ios::app);
	Temp << setw(wd) << istep <<  setw(wd) << walker[it]/(double)npart << endl;
	Temp.close();
   }

}

void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
}

void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Print results for current block
{
  double r, gdir;
  ofstream Epot, Ekin, Etot, Temp, Gofr, Gave;
  const int wd=12;

  cout << "Block number " << iblk << endl;

  Epot.open("ave_epot_Argon_solid.out",ios::app);
  Ekin.open("ave_ekin_Argon_solid.out",ios::app);
  Temp.open("ave_temp_Argon_solid.out",ios::app);
  Etot.open("ave_etot_Argon_solid.out",ios::app);
  Gofr.open("output.gofr_solid",ios::app);
  Gave.open("output.gave_solid",ios::app);

  stima_epot = blk_av[iv]/blk_norm/(double)npart; //Potential energy per particle	
  glob_av[iv]  += stima_epot;
  glob_av2[iv] += stima_epot*stima_epot;
  err_epot = Error(glob_av[iv],glob_av2[iv],iblk);
  Epot << setw(wd) << iblk <<  setw(wd) << stima_epot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_epot << endl;

  stima_ekin = blk_av[ik]/blk_norm/(double)npart; //Potential energy per particle	
  glob_av[ik]  += stima_ekin;
  glob_av2[ik] += stima_ekin*stima_ekin;
  err_ekin = Error(glob_av[ik],glob_av2[ik],iblk);
  Ekin << setw(wd) << iblk <<  setw(wd) << stima_ekin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_ekin << endl;

  stima_temp = blk_av[it]/blk_norm/(double)npart; //Potential energy per particle	
  glob_av[it]  += stima_temp;
  glob_av2[it] += stima_temp*stima_temp;
  err_temp = Error(glob_av[it],glob_av2[it],iblk);
  Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
 
  stima_etot = blk_av[ie]/blk_norm/(double)npart; //Potential energy per particle	
  glob_av[ie]  += stima_etot;
  glob_av2[ie] += stima_etot*stima_etot;
  err_etot = Error(glob_av[ie],glob_av2[ie],iblk);
  Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;

  for (int i=igofr; i<n_props; i++){
	r = (i-igofr)*bin_size;
	gdir = blk_av[i]/blk_norm/(rho*npart*(4./3.)*pi*(pow(r+bin_size,3)-pow(r,3))); //Radial distribution function
	glob_av[i] += gdir;
	glob_av2[i] += gdir*gdir;
	err_gdir=Error(glob_av[i],glob_av2[i],iblk);
	Gofr << setw(wd) << iblk <<  setw(wd) << r <<  setw(wd)  << setw(wd) << glob_av[i]/(double)iblk << setw(wd) << err_gdir << endl;
	if (iblk==nblk){
		Gave << setw(wd) << r <<  setw(wd)  << setw(wd) << glob_av[i]/(double)iblk << setw(wd) << err_gdir << endl;
	}
   }


	
  cout << "----------------------------" << endl << endl;

  
  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  Gofr.close();
  Gave.close();
}


void ConfFinal(void){ //Write final configuration and copy it to restart
  ofstream WriteConf;
  ofstream CopyConf;

  cout << "Print final configuration to file config.final and copy it in config.0 " << endl << endl;
  WriteConf.open("config.final");
  CopyConf.open("config.0");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    CopyConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  CopyConf.close();
  return;
}

void ConfFinalOld(void){ //Write old final configuration, which is before the last one, and copy it to restart 
  ofstream WriteConfOld;
  ofstream CopyConfOld;

  cout << "Print old final configuration to file old.final and copy it in old.0 " << endl << endl;
  WriteConfOld.open("old.final");
  CopyConfOld.open("old.0");

  for (int i=0; i<npart; ++i){
    WriteConfOld << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
    CopyConfOld << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConfOld.close();
  CopyConfOld.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk) 
{
    if (iblk==1){ 
	return 0.;
    }
    else{
	return sqrt(fabs((sum2/(double)iblk - pow(sum/(double)iblk,2)))/(double)iblk);
    }
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
