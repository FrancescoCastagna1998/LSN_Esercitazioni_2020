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
#include <string>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include "Variational_Monte_Carlo.h"

using namespace std;

 
int main (int argc, char *argv[]){


  Input(); //Inizialization

  for(int istep=1; istep <= nequi; ++istep)
  {
	Move();
	Measure();
	if (istep%1000==0) cout << "Number of time-steps: " << istep << endl;
	Equilibration(istep); //  istant values of hamiltonian expextion values
		 
  }


   
   for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
	{
		Reset(iblk);   //Reset block averages
    		for(int istep=1; istep <= nstep; ++istep)
    		{
      			Move();
      			Measure();
      			Accumulate(); //Update block averages
      		}
    		Averages(iblk);   //Print results for current block
  	}	
  
   ConfFinal(); //Write final configuration
   
   return 0;
}

void Input(void){

  ifstream ReadInput;

  //Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();


  cout << "Single quantum particle in 1D      " << endl;
  cout << "Variational Monte Carlo simulation             " << endl << endl;
  cout << "Trial Wave Function psi_T(x;mu,sigma) = exp() " << endl << endl; // completa
  cout << "External potential V(x) = x^4 - 5/2 * x^2" << endl << endl;
  cout << "The program uses units with hbar = 1 and m = 1  " << endl;

  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> mu;
  cout << "The first variational parameter of the trial wave function = " << mu << endl;  
  ReadInput >> sigma;
  cout << "The second variational parameter of the trial wave function = " << sigma << endl;

// Trial wave function
  Psi_trial = new TrialWaveFunction(mu,sigma);

// External potential
  Ext_Pot = new Potential((-5)/2.);

// starting position
  ReadInput >> start_position;
  x_pos = start_position;

  cout << "The initial position is x = " << x_pos << endl << endl;

  ReadInput >> delta;

  ReadInput >> nequi;

  ReadInput >> nblk;

  ReadInput >> nstep;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of steps for equilibration = " << nequi << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;

//Prepare arrays for measurements
  ih = 0; //Hamiltonian expectation values
  n_props = 1; //Number of observables

  ReadInput >> nbins;				
  ReadInput >> histogram_start;				
  ReadInput >> histogram_end;

// Histogram with the sampled configurations of |psi_trial(x)|^2
  ihist = 1;
  n_props = n_props + nbins;
  bin_size = (histogram_end - histogram_start)/(double)nbins;

  cout << "The histogram with the sample configurations has " << nbins << " bins " << " with bins dimension = " << bin_size << endl;
  cout << "The histogram starts from x = " << histogram_start << " , and ends to x = " << histogram_end << endl << endl;

  ReadInput.close();

//Evaluate hamiltonian expectation value of the initial position
  Measure();

//Print initial value for the hamiltonian expectation value
  cout << "Initial hamiltonian expectation value = " << walker[ih] << endl << endl;

}
void Move(void)
{
  double p;
  double xold,xnew;

//Old
  xold = x_pos;
//New
  xnew = x_pos + delta*(rnd.Rannyu() - 0.5);

//Metropolis test
  p = min(1., (pow(Psi_trial->Eval(xnew),2))/(pow(Psi_trial->Eval(xold),2))); 
  if(p>=rnd.Rannyu()){
	x_pos = xnew;
	accepted++;
  }
  else{
	x_pos = xold;	
  }
  attempted++;
}


void Measure()
{
  int bin;

//reset the hystogram of |psi_trial(x)|^2
  for (int k=ihist; k<ihist+nbins; ++k) walker[k]=0.0;

//update of the histogram of |psi_trial(x)|^2
  for (bin=ihist; bin<n_props; bin++){
	if ( (x_pos > (histogram_start + (bin-ihist)*bin_size))  && (x_pos < (histogram_start + bin*bin_size))){
		walker[bin] = walker[bin] + 1;
	}
  }
  walker[ih] = ((-1/2.)*Psi_trial->EvalDerSec(x_pos) + Ext_Pot->Eval(x_pos)*Psi_trial->Eval(x_pos))/(Psi_trial->Eval(x_pos));
}

void Equilibration(int istep)
{
   ofstream Ham;
   const int wd=12;
    
   Ham.open("output.expection_values_hamiltonian_equi",ios::app);
   Ham << setw(wd) << istep <<  setw(wd) << walker[ih] << endl;
   Ham.close();
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

void Averages(int iblk) //Print results for current block
{
   double x_current, PsiT2;
   double norma = 0.0; // normalization histogram
   ofstream Ham, PsiT_2;
   const int wd=12;
    
   cout << "Block number " << iblk << endl;
   cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
   Ham.open("output.expection_values_hamiltonian",ios::app);
   PsiT_2.open("ProbabilityGS.dat",ios::app);

    
   stima_ham = blk_av[ih]/blk_norm; // hamiltonian expectation value
   glob_av[ih] += stima_ham;
   glob_av2[ih] += stima_ham*stima_ham;
   err_ham = Error(glob_av[ih],glob_av2[ih],iblk);
   Ham << setw(wd) << iblk <<  setw(wd) << stima_ham << setw(wd) << glob_av[ih]/(double)iblk << setw(wd) << err_ham << endl;


 for (int i=ihist; i<n_props; i++){
	x_current = histogram_start + (i-ihist+0.5)*bin_size;
	PsiT2 = blk_av[i]/blk_norm; // |psi_trial(x)|^2
	glob_av[i] += PsiT2;
	glob_av2[i] += PsiT2*PsiT2;
	err_PsiT2=Error(glob_av[i],glob_av2[i],iblk);
	if (iblk==nblk){
		for (int j=ihist; j<n_props; j++){
			norma += (glob_av[j]/(double)iblk);
		}
		norma *= bin_size;
		PsiT_2 << setw(wd) << x_current <<  setw(wd) << glob_av[i]/(double)iblk/norma << setw(wd) << err_PsiT2/norma << endl;
	}
    }
    
   cout << "----------------------------" << endl << endl;

   Ham.close();
   PsiT_2.close();
   
}

void ConfFinal(void)
{
  ofstream WritePos;

  cout << "Print final position to file pos.final " << endl << endl;
  WritePos.open("pos.final");

  WritePos << " x = " << x_pos << endl;

  WritePos.close();

  rnd.SaveSeed();
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
