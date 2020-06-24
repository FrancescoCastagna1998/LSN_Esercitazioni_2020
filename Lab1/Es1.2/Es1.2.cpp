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
#include "random.h"
#include "Vettore.h"
#include "funzioni.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
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

   
	
   int M = 10000; //numero di realizzazioni per la variabile somma

   cout << "This program checks the Central Limit Theorem and the Levy Theorem " << endl;

//media su una singola variabile casuale (N=1)

   Vettore st1(M); //contenitore per numeri generati uniformemente tra [0.5,6.5)
   Vettore exp1(M);//contenitore per numeri generati secondo una pdf esponenziale
   Vettore cl1(M);//contenitore per numeri generati secondo una pdf di Cauchy-Lorentz

   for(int i=0; i<st1.GetN(); i++){
   	st1.SC(i,rnd.Rannyu(0.5,6.5));
   }
   Stampa ("st1.txt",st1);

   for(int i=0; i<exp1.GetN(); i++){
   	exp1.SC(i,rnd.Exp(1.));
   }
   Stampa ("exp1.txt",exp1);

   for(int i=0; i<cl1.GetN(); i++){
   	cl1.SC(i,rnd.CL(0.,1.));
   }
   Stampa ("cl1.txt",cl1); 

//media su 2 variabili casuali (N=2)

   Vettore st2(M);
   Vettore exp2(M);
   Vettore cl2(M);
   int N = 2;

   for(int i=0; i<st2.GetN(); i++){
	double k=0.;
	for (int j=0; j<N; j++){
   		k+=rnd.Rannyu(0.5,6.5);
	}
	st2.SC(i,(k/N));
   }
   Stampa ("st2.txt",st2);

   for(int i=0; i<exp2.GetN(); i++){
	double k=0.;
	for (int j=0; j<N; j++){
   		k+=rnd.Exp(1.);
	}
	exp2.SC(i,(k/N));
   }
   Stampa ("exp2.txt",exp2);

   for(int i=0; i<cl2.GetN(); i++){
	double k=0.;
	for (int j=0; j<N; j++){
   		k+=rnd.CL(0.,1.);
	}
	cl2.SC(i,(k/N));
   }
   Stampa ("cl2.txt",cl2);

//media su 10 variabili casuali (N=10)

   Vettore st10(M);
   Vettore exp10(M);
   Vettore cl10(M);
   N = 10;

   for(int i=0; i<st10.GetN(); i++){
	double k=0.;
	for (int j=0; j<N; j++){
   		k+=rnd.Rannyu(0.5,6.5);
	}
	st10.SC(i,(k/N));
   }
   Stampa ("st10.txt",st10);

   for(int i=0; i<exp10.GetN(); i++){
	double k=0.;
	for (int j=0; j<N; j++){
   		k+=rnd.Exp(1.);
	}
	exp10.SC(i,(k/N));
   }
   Stampa ("exp10.txt",exp10);

   for(int i=0; i<cl10.GetN(); i++){
	double k=0.;
	for (int j=0; j<N; j++){
   		k+=rnd.CL(0.,1.);
	}
	cl10.SC(i,(k/N));
   }
   Stampa ("cl10.txt",cl10);

//media su 100 variabili casuali (N=100)

   Vettore st100(M);
   Vettore exp100(M);
   Vettore cl100(M);
   N = 100;

   for(int i=0; i<st100.GetN(); i++){
	double k=0.;
	for (int j=0; j<N; j++){
   		k+=rnd.Rannyu(0.5,6.5);
	}
	st100.SC(i,(k/N));
   }
   Stampa ("st100.txt",st100);

   for(int i=0; i<exp100.GetN(); i++){
	double k=0.;
	for (int j=0; j<N; j++){
   		k+=rnd.Exp(1.);
	}
	exp100.SC(i,(k/N));
   }
   Stampa ("exp100.txt",exp100);

   for(int i=0; i<cl100.GetN(); i++){
	double k=0.;
	for (int j=0; j<N; j++){
   		k+=rnd.CL(0.,1.);
	}
	cl100.SC(i,(k/N));
   }
   Stampa ("cl100.txt",cl100);
 

   rnd.SaveSeed();
   return 0;
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
