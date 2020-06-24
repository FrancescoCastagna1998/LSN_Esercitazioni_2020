#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <ostream>
#include <iomanip>
#include "Buffon.h"

using namespace std;

int main (int argc, char *argv[]){

   Input(); //Inizialization

   for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
	{
		Reset(iblk);   //Reset block averages
    		for(int ithrow=1; ithrow <= nthrow; ++ithrow)
    		{
      			Move();
      			Measure();
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
   ReadInput.open("input.dat");

   cout << "Buffon experiment" << endl << endl;
   cout << "Monte Carlo simulation" << endl << endl;
   
   ReadInput >> L;
   cout << "The lenght of the needle:  " << L << endl;
   ReadInput >> d;
   cout << "The distance between the straight lines on the horizontal plane :  " << d << endl << endl;

   ReadInput >> nblk;
   ReadInput >> nthrow;

   cout << "Number of blocks = " << nblk << endl;
   cout << "Number of steps in one block = " << nthrow << endl << endl;

   ReadInput.close();

}
void Reset(int iblk) //Reset block averages
{  
   if(iblk == 1)
   {
	glob_av = 0;
	glob_av2 = 0;
   }
   N_hit = 0;
}
void Move(void)
{
// expoloit the simmetry of the system
// generate the position of the needle's length centre in [-d/2; d/2)
   y_centre = rnd.Rannyu(-d/2.,d/2.); 
// generate Theta uniformly in the interval [0,pigreco) without using pigreco
   Theta = rnd.Unif_Theta();
}
void Measure(void)
{
   double y1, y2;

   y1 = y_centre + (L/2.)*sin(Theta);
   y2 = y_centre - (L/2.)*sin(Theta);
// check if the needle intersects the horizontal line
   if((y1>=0 && y2<=0) || (y1<=0 && y2>=0)){  
	N_hit++;
   }
}
void Averages(int iblk) //Print results for current block
{  
    ofstream Pi;
    const int wd=12;
    
    Pi.open("pigreco.txt",ios::app);
    
    stima_pigreco = (2.*L*nthrow)/(N_hit*d); 
    glob_av += stima_pigreco;
    glob_av2 += stima_pigreco*stima_pigreco;
    err_pigreco=Error(glob_av,glob_av2,iblk);

//Pigeco
    Pi << setw(wd) << iblk <<  setw(wd) << stima_pigreco << setw(wd) << glob_av/(double)iblk << setw(wd) << err_pigreco << endl;

    if (iblk%10 == 0){
	cout << "Block number " << iblk << endl;
	cout << "----------------------------" << endl << endl;
    }

    Pi.close();
}
double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}
