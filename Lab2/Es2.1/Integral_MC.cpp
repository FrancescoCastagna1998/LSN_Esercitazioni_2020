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
#include "FunzioneBase.h"
#include "Coseno.h"
#include "CosenoIS.h"
#include "funzioni.h"
#include <cmath>


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

   cout << "This program provides the Monte Carlo evaluation of the same integral with two different algorithm: " << endl;
   cout << "1) sampling random variables uniformly distributed in [0,1) " << endl;
   cout << "2) using importance sampling " << endl; 

   int M = 1000000; //MC steps
   int N = 100; //Blocks
   int L = M/N; //MC steps for each block


   FunzioneBase* cos = new Coseno(M_PI/2,M_PI/2); // integranda

   Vettore INT(N); // integral estimation
   Vettore INT2(N); // squared integral extimation
   Vettore sumprog(N);
   Vettore su2prog(N);
   Vettore errprog(N);


   
   for(int i=0; i<N; i++){
	double temp=0.;
	for(int j=0; j<L; j++){
		temp+=cos->Eval(rnd.Rannyu());
	}
	INT.SC(i,temp/L);
	INT2.SC(i,pow(INT.GC(i),2));
   }

// blocking method
   for(int i=0; i<N; i++){
	double sum1=0;
	double sum2=0;
   	for(int j=0; j<i+1; j++){
		sumprog.SC(i,sum1 += INT.GC(j));
		su2prog.SC(i,sum2 += INT2.GC(j));
	}
	sumprog.SC(i,(sumprog.GC(i))/(i+1));//valore medio cumulativo
	su2prog.SC(i,(su2prog.GC(i))/(i+1));//valore medio quadro cumulativo
	errprog.SC(i,calcolaErrore(sumprog,su2prog,i));//incertezza statistica
   }
   
   Stampa("int.txt",sumprog ,errprog,1);

//ora mediante la tecnica del campionamento rilevante prendo come nuova distribuzione di probabilità che meglio approssima l'integranda
//la retta passante per (pigreco/2,0) e (1,0), la normalizzo e uso il metodo dell'inversa della funzione cumulativa per estrarre la 
//variabile casuale ditribuita in questo modo

   FunzioneBase* cosIS = new CosenoIS(M_PI/4,M_PI/2); // new integrand after the application of the importance sampling

   for(int i=0; i<N; i++){
	double temp1=0.;
	for(int j=0; j<L; j++){
		temp1+=cosIS->Eval(rnd.Retta());
	}
	INT.SC(i,temp1/L);
	INT2.SC(i,pow(INT.GC(i),2));
   }

// blocking method
   for(int i=0; i<N; i++){
	double sum1=0;
	double sum2=0;
   	for(int j=0; j<i+1; j++){
		sumprog.SC(i,sum1 += INT.GC(j));
		su2prog.SC(i,sum2 += INT2.GC(j));
	}
	sumprog.SC(i,(sumprog.GC(i))/(i+1));//valore medio cumulativo
	su2prog.SC(i,(su2prog.GC(i))/(i+1));//valore medio quadro cumulativo
	errprog.SC(i,calcolaErrore(sumprog,su2prog,i));//incertezza statistica
   }
   
   Stampa("int_IS.txt",sumprog ,errprog,1);//stampo i valori dell'integrale calcolato sfruttando la tecninca dell'importance sampling


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
