#include <iostream>
#include <fstream>
#include <ostream>
#include <iomanip>
#include <string>
#include <algorithm>

#include "Hydrogen_Atom.h"


using namespace std;

int main()
{ 
  
  Input(); //Inizialization

//Equilibration
  for(int istep=1; istep <= nequistep; ++istep)
  {
	Move();
	Measure();
	if (istep%1000==0) cout << "Number of equilibration time-steps: " << istep << endl;
	Equilibration(istep); //  istant value radial distance r		 
  }
		
//Simulation
   for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
   {
	Reset(iblk);   //Reset block averages
    	for(int istep=1; istep <= nstep; ++istep)
    	{
      		Move();
      		Measure();
      		Accumulate(); //Update block averages
      		if(istep%100 == 0) ConfXYZ(); //Print sampled configurations 
    	}
    	Averages(iblk);   //Print results for current block
    }	

  rnd.SaveSeed();
  return 0;
}


void Input(void)
{
ifstream ReadInput;

  cout << "Hydrogen atom        " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Metropolis algorithm         " << endl << endl;
  cout << "The program uses Bohr radius units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

  
//Read input informations
  ReadInput.open("input.dat");

// Read the starting configuration (starting position in 3D-space)
  ReadInput >> x_start;
  ReadInput >> y_start;
  ReadInput >> z_start;

  Pos.SetX(x_start);
  Pos.SetY(y_start);
  Pos.SetZ(z_start);
  
  cout << "The origin of frame of reference is : O_xyz ( x = " << Ori.GetX() << ", y = " << Ori.GetY() << ", z = " << Ori.GetZ() << " )" << endl << endl; 

  cout << "The starting configuration is : P ( x = " << Pos.GetX() << ", y = " << Pos.GetY() << ", z = " << Pos.GetZ() << " )" << endl << endl; 

//Hydrogen atom wave function
  ReadInput >> wave_function;
  if (wave_function == 0){
	Wave_Fun_HA = new FunOndaHIGS((pow(a0,-(3./2.)))/sqrt(M_PI),1.); //ground state wave function
	cout << "The program samples the squared modulus of ground state wave function: |Psi_1s|^2 " << endl << endl;
  }
  if (wave_function == 1){
	Wave_Fun_HA = new FunOndaHI2P(sqrt(2)*(pow(a0,-(3./2.)))/(8.*sqrt(M_PI)),(1./2.)); //excited level 2P wave function
	cout << "The program samples the squared modulus of excited level 2p wave function: |Psi_2p|^2 " << endl << endl;
  }

// trial transition probability in the metropolis algorithm
  ReadInput >> transition_probability;
  if (transition_probability == 0) cout << "The program perform Metropolis moves with uniform translations for each coordinate" << endl << endl;
  if (transition_probability == 1) cout << "The program perform Metropolis moves with normal multivariate Gaussian translations for each coordinate" << endl << endl;

  if (transition_probability == 0) ReadInput >> delta;
  if (transition_probability == 1) ReadInput >> sigma;

  ReadInput >> nequistep;

  ReadInput >> nblk;

  ReadInput >> nstep;

  if (transition_probability == 0) cout << "Moves parameter = " << delta << endl;
  if (transition_probability == 1) cout << "Sigma parameter = " << sigma << endl;
  cout << "Number of steps for the equilibration = " << nequistep << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;

  ReadInput.close();

//Prepare array for measurement
  ir = 0; //Radial distance

  n_props = 1; //Number of observables


//Evaluate radial distance of the initial configuration
  Measure();

//Print initial values for radial distance
  cout << "Initial radial distance <r> = " << walker[ir] << endl;

}

void Equilibration(int istep)
{
   ofstream Rad_dist;

   if (transition_probability == 0 && wave_function == 0) Rad_dist.open("output.r_equi_GS_Unif",ios::app);
   if (transition_probability == 1 && wave_function == 0) Rad_dist.open("output.r_equi_GS_Gauss",ios::app);
   if (transition_probability == 0 && wave_function == 1) Rad_dist.open("output.r_equi_2P_Unif",ios::app);
   if (transition_probability == 1 && wave_function == 1) Rad_dist.open("output.r_equi_2P_Gauss",ios::app);

   Rad_dist << "   " << istep <<  "   " << walker[ir] << endl;
   Rad_dist.close();
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
   attempted = 0;
   accepted = 0;
}

void Move(void)
{
  double p, old_squared_modulus_wave_function, new_squared_modulus_wave_function; // acceptance probability and old-new |psi|^2
  double xold, yold, zold, xnew, ynew, znew; // old and new configurations

//Old configuraion
  xold = Pos.GetX();
  yold = Pos.GetY();
  zold = Pos.GetZ();

//old density probability
  old_squared_modulus_wave_function = Wave_Fun_HA->Eval_squared_modulus(xold,yold,zold);

//New configuration
  if (transition_probability == 0){ // uniform transition
  	xnew = xold + delta*(rnd.Rannyu() - 0.5);
  	ynew = yold + delta*(rnd.Rannyu() - 0.5);
  	znew = zold + delta*(rnd.Rannyu() - 0.5);
  }
  else { // normal multivariate transition
	xnew = rnd.Gauss(xold,sigma);
  	ynew = rnd.Gauss(yold,sigma);
  	znew = rnd.Gauss(zold,sigma);
  }

//new density probability
  new_squared_modulus_wave_function = Wave_Fun_HA->Eval_squared_modulus(xnew,ynew,znew);


//Metropolis test
  p = min(1., (new_squared_modulus_wave_function)/(old_squared_modulus_wave_function)); // symmetric transition probability

  if(p >= rnd.Rannyu()){
	Pos.SetX(xnew);
	Pos.SetY(ynew);
	Pos.SetZ(znew); 
	accepted++;
  }
  attempted++;
}

void Measure()
{
  double r = 0.0;

//radial distance
  r = Pos.Distanza(Ori);
  walker[ir] = r;
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
    
   ofstream Radial_dist;
   const int wd=12;
    
   cout << "Block number " << iblk << endl;
   cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
   if (transition_probability == 0 && wave_function == 0) Radial_dist.open("output.rAve_GS_Unif",ios::app);
   if (transition_probability == 1 && wave_function == 0) Radial_dist.open("output.rAve_GS_Gauss",ios::app);
   if (transition_probability == 0 && wave_function == 1) Radial_dist.open("output.rAve_2P_Unif",ios::app);
   if (transition_probability == 1 && wave_function == 1) Radial_dist.open("output.rAve_2P_Gauss",ios::app);
    
   stima_radius = blk_av[ir]/blk_norm; //Radial distance
   glob_av[ir] += stima_radius;
   glob_av2[ir] += stima_radius*stima_radius;
   err_radius = Error(glob_av[ir],glob_av2[ir],iblk);
    
    
//Radial distance
   Radial_dist << setw(wd) << iblk <<  setw(wd) << stima_radius << setw(wd) << glob_av[ir]/(double)iblk << setw(wd) << err_radius << endl;

   cout << "----------------------------" << endl << endl;

   Radial_dist.close();
}

void ConfXYZ(void)   //Print current configuration
{     
   ofstream PrintXYZ;

   if (transition_probability == 0 && wave_function == 0) PrintXYZ.open("output.XYZsampled_GS_Unif",ios::app);
   if (transition_probability == 1 && wave_function == 0) PrintXYZ.open("output.XYZsampled_GS_Gauss",ios::app);
   if (transition_probability == 0 && wave_function == 1) PrintXYZ.open("output.XYZsampled_2P_Unif",ios::app);
   if (transition_probability == 1 && wave_function == 1) PrintXYZ.open("output.XYZsampled_2P_Gauss",ios::app);

   PrintXYZ << Pos.GetX() << "  " << Pos.GetY() << "  " << Pos.GetZ() << "  " << endl;
   PrintXYZ.close();
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

