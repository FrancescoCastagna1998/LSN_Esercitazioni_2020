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
#include <cmath>
#include <algorithm>
#include <ostream>
#include <iomanip>
#include "European_Option.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Input(); //Inizialization

   for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
	{
		Reset(iblk);   //Reset block averages
    		for(int istep=1; istep <= nstep; ++istep)
    		{
      			Move1();
			S_T2 = starting_asset_price;
			for (int iint=1; iint<=ninterval; iint++)
			{ 
				Move2();
			}
      			Measure();
      			Accumulate();   //Update block averages
    		}
    		Averages(iblk);   //Print results for current block
  	}


   rnd.SaveSeed();
   return 0;
}

void Input(void)
{
   ifstream ReadInput;

//Inizialize random generator
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
   ReadInput.open("Input.dat");

   //partiamo campionando direttamente il prezzo del bene al tempo di consegna T
   //sfrutteremo il modello del moto browniano geometrico per la descrizione dell'asset price in termini di equazioni differenziali stocastiche
   cout << "European option pricing!" << endl;
   cout << "Monte Carlo simulations" << endl;
   cout << "The time evolution of asset price exhibits a Geometric Brownian Motion GBM(r,sigma^2) " << endl << endl;
   
   ReadInput >> starting_asset_price;
   cout << "Asset price for t=0, S(0): " << starting_asset_price << endl;

   ReadInput >> delivery_time;
   cout << "Deleivery time: " << delivery_time << endl;

   ReadInput >> strike_price;
   cout << "Strike price: " << strike_price << endl;

   ReadInput >> r;
   cout << "Risk-free interest rate: " << r << endl;

   ReadInput >> volatility;
   cout << "Volatility: " << volatility << endl;

   ReadInput >> ninterval;
   deltaT = delivery_time/double(ninterval);
   cout << "Number of intervals which divides [0,T]: " << ninterval << ", with amplitude: " <<  deltaT << endl << endl;

   cout << "Perform the Monte carlo simulation for call and put-option price: " << endl;
   cout << "1) Sampling directly the final price asset S(t)" << endl;
   cout << "2) Sampling the discretized path of the asset price dividing [0,T] in " << ninterval << " time intervals "  << endl << endl;

   ReadInput >> nblk;
   ReadInput >> nstep;

   cout << "Number of blocks = " << nblk << endl;
   cout << "Number of steps in one block = " << nstep << endl << endl;

   ReadInput.close();

//Prepare arrays for measurements
   // first Monte Carlo Simulation
   iput1 = 0; 
   icall1 = 1;
   // second Monte Carlo Simulation
   iput2 = 2; 
   icall2 = 3;

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
void Move1(void)
{
   S_T1 = starting_asset_price*exp((r-pow(volatility,2)/2)*delivery_time+volatility*rnd.Gauss(0,delivery_time));
}
void Move2(void)
{
   S_T2 = S_T2*exp((r-pow(volatility,2)/2)*deltaT+volatility*rnd.Gauss(0,1)*sqrt(deltaT));
}
void Measure(void)
{
   double  Put_1, Call_1, Put_2, Call_2;

   Call_1 = exp(-r*delivery_time)*max(0.,S_T1-strike_price);
   Put_1 = exp(-r*delivery_time)*max(0.,strike_price-S_T1);

   Option_Price[icall1] = Call_1;
   Option_Price[iput1] = Put_1; 

   Call_2 = exp(-r*delivery_time)*max(0.,S_T2-strike_price);
   Put_2 = exp(-r*delivery_time)*max(0.,strike_price-S_T2);

   Option_Price[icall2] = Call_2;
   Option_Price[iput2] = Put_2; 
		
}
void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + Option_Price[i];
   }
   blk_norm = blk_norm + 1.0;
}
void Averages(int iblk) //Print results for current block
{
    
    ofstream call_1, put_1, call_2, put_2;
    const int wd=12;

    
    call_1.open("Call_sampling_final_asset_price.txt",ios::app);
    put_1.open("Put_sampling_final_asset_price.txt",ios::app);
    call_2.open("Call_sampling_dis.txt",ios::app);
    put_2.open("Put_sampling_dis.txt",ios::app);
    
    stima_call1 = blk_av[icall1]/blk_norm; 
    glob_av[icall1] += stima_call1;
    glob_av2[icall1] += stima_call1*stima_call1;
    err_call1=Error(glob_av[icall1],glob_av2[icall1],iblk);

    stima_put1 = blk_av[iput1]/blk_norm; 
    glob_av[iput1] += stima_put1;
    glob_av2[iput1] += stima_put1*stima_put1;
    err_put1=Error(glob_av[iput1],glob_av2[iput1],iblk);

    stima_call2 = blk_av[icall2]/blk_norm; 
    glob_av[icall2] += stima_call2;
    glob_av2[icall2] += stima_call2*stima_call2;
    err_call2=Error(glob_av[icall2],glob_av2[icall2],iblk);

    stima_put2 = blk_av[iput2]/blk_norm; 
    glob_av[iput2] += stima_put2;
    glob_av2[iput2] += stima_put2*stima_put2;
    err_put2=Error(glob_av[iput2],glob_av2[iput2],iblk);


//Call_1
    call_1 << setw(wd) << iblk <<  setw(wd) << stima_call1 << setw(wd) << glob_av[icall1]/(double)iblk << setw(wd) << err_call1 << endl;
//Put_1
    put_1 << setw(wd) << iblk <<  setw(wd) << stima_put1 << setw(wd) << glob_av[iput1]/(double)iblk << setw(wd) << err_put1 << endl;
//Call_2
    call_2 << setw(wd) << iblk <<  setw(wd) << stima_call2 << setw(wd) << glob_av[icall2]/(double)iblk << setw(wd) << err_call2 << endl;
//Put_2
    put_2 << setw(wd) << iblk <<  setw(wd) << stima_put2 << setw(wd) << glob_av[iput2]/(double)iblk << setw(wd) << err_put2 << endl;


    if (iblk%10 == 0){
	cout << "Block number " << iblk << endl;
	cout << "----------------------------" << endl << endl;
    }

    call_1.close();
    put_1.close();
    call_2.close();
    put_2.close();
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
