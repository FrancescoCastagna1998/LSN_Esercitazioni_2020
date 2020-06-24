/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <algorithm> //necessary for the min function
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  cout << "This program admits two possibilities to work! " << endl << endl;
  cout << "1) Perform the equilibration of the system before starting the real the simulation. " << endl;
  cout << "2) Perform the real simulation using as the starting configuration the one saved in the file config.final " << endl << endl;
  int answer;
  cout << "Please type the answer number to select the correct option: " << endl;
  cin >> answer;

   
  Input();
  if (answer==1){
	for(int istep=1; istep <= nequistep; ++istep)
	{
		Move(metro);
		Measure();
		if (istep%1000==0) cout << "Number of equilibration time-steps: " << istep << endl;
		Equilibration(istep); //  istant value of internal energy and magnetization 
	}
	ConfFinal(); //Write final configuration
	
	return 0;	
  }
 
  if(answer==2) InputRestart();
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nequistep;

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of steps for the equilibration = " << nequistep << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
  cout << "Initial magnetization = " << walker[im]/(double)nspin << endl;
  cout << "Initial magnetic susceptibility = " << walker[ix]/(double)nspin << endl;
}

void InputRestart(void)
{
  ifstream ReadConfFinal;

  cout << "Read final configuration from file config.final in order to start the real simulation after initial equilibration. " << endl << endl;

  ReadConfFinal.open("config.final");

  for (int i=0; i<nspin; ++i)
  {
	ReadConfFinal >> s[i];
  }

  ReadConfFinal.close();
}
void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {
	sm = s[o];

	energy_old = Boltzmann (sm,o);
	energy_new = Boltzmann (-sm,o);

	p = min(1.,exp(-beta*(energy_new - energy_old))); // Acceptance probability

	double r = rnd.Rannyu(); // Uniform random variable in [0,1)
	if (r < p)
	{
		s[o] = -sm; // Flip the spin
		accepted++;
	}
	else
	{
		s[o] = sm;
	}
     attempted++;		
    }
    else //Gibbs sampling
    {
	sm = 1.;
	energy_up = Boltzmann (sm,o); // sm = +1
	energy_down = Boltzmann (-sm,o); // sm = -1

	p = 1./(1.+exp(-beta*(energy_down - energy_up))); // Prob(s[o] = +1) and Prob(s[o] = -1) = 1-p

	double r = rnd.Rannyu(); // Uniform random variable in [0,1)
	if (r < p)
	{
		s[o] = sm; 
	}
	else
	{
		s[o] = -sm;		
	}
	
	accepted++; // The algorithm acceptes always the move
	attempted++;
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i];
  }
  walker[iu] = u;
  walker[ic] = pow(u,2); // H^2
  walker[im] = m;
  walker[ix] = pow(m,2);
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


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}

void Equilibration(int istep)
{

   ofstream Ene, Mag;
   const int wd=12;

   Ene.open("output.ene_1_T2_h_equi",ios::app);
   Ene << setw(wd) << istep <<  setw(wd) << walker[iu]/(double)nspin << endl;
   Ene.close();

   Mag.open("output.mag_1_T2_h_equi",ios::app);
   Mag << setw(wd) << istep <<  setw(wd) << walker[im]/(double)nspin << endl;
   Mag.close();

}
void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   ofstream Enefin, Heatfin, Magfin, Chifin;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;

    if (h == 0.0) // Termine di accoppiamento nullo
    {
    
	Ene.open("output.ene.1",ios::app);
	stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy	
	glob_av[iu]  += stima_u;
	glob_av2[iu] += stima_u*stima_u;
	err_u = Error(glob_av[iu],glob_av2[iu],iblk);
	Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
	Ene.close();
    	if (iblk==20){
		Enefin.open("output.ene.1_final",ios::app);
		Enefin << setw(wd) << temp <<  setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
		Enefin.close();
    	}

	
	Heat.open("output.heat.1",ios::app);
	stima_c = (pow(beta,2)*((blk_av[ic]/blk_norm) - pow(blk_av[iu]/blk_norm,2)))/(double)nspin; //Heat capacity
	glob_av[ic]  += stima_c;
	glob_av2[ic] += stima_c*stima_c;
	err_c = Error(glob_av[ic],glob_av2[ic],iblk);
	Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
	Heat.close();
	if (iblk==20){
		Heatfin.open("output.heat.1_final",ios::app);
		Heatfin << setw(wd) << temp <<  setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
		Heatfin.close();
    	}

	Chi.open("output.chi.1",ios::app);
	stima_x = beta*(blk_av[ix]/blk_norm)/(double)nspin; // Susceptibility
	glob_av[ix]  += stima_x;
	glob_av2[ix] += stima_x*stima_x;
	err_x = Error(glob_av[ix],glob_av2[ix],iblk);
	Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
	Chi.close();
	if (iblk==20){
		Chifin.open("output.chi.1_final",ios::app);
		Chifin << setw(wd) << temp <<  setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
		Chifin.close();
    	}
    }

    if (h == 0.02) // Termine di accoppiamento non nullo
    {
	Mag.open("output.mag.1",ios::app);
	stima_m = blk_av[im]/blk_norm/(double)nspin; // Magnetization
	glob_av[im]  += stima_m;
	glob_av2[im] += stima_m*stima_m;
	err_m = Error(glob_av[im],glob_av2[im],iblk);
	Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
	Mag.close();
	if (iblk==20){
		Magfin.open("output.mag.1_final",ios::app);
		Magfin << setw(wd) << temp <<  setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
		Magfin.close();
    	}
    }


    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
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
