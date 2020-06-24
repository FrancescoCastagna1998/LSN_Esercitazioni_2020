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
#include "Posizione.h"
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

   Posizione ORI;

   int M = 10000; //numero totale delle simulazioni del Random Walk 3D
   int N = 100; // numero di blocchi in cui suddividere le M simulazioni e numero di passi massimo nel RW
   int L = M/N;
   double a = 1.; //passo reticolare

   Vettore RW3D(N);//vale sia per il caso discreto che per il continuo
   Vettore DIST(N); 
   Vettore DIST2(N);
   Vettore temp(N);
   Vettore errprog(N);


   for(int i=0; i<N; i++){ //numero dei blocchi
	for(int j=0; j<L; j++){ //numero di RWs per blocco
		Posizione p;
		for(int l=0; l<N; l++){//numero dei passi
			p.RWdiscreto(rnd.Axis(),rnd.Verso(),a);
			temp.SC(l,pow(p.Distanza(ORI),2));
			RW3D.SC(l,RW3D.GC(l)+temp.GC(l));
			}
		}
	for(int k=0; k<N; k++){
		RW3D.SC(k,sqrt(RW3D.GC(k)/L));
		DIST.SC(k,DIST.GC(k)+RW3D.GC(k));//la somma delle radici di moduli quadri mediati sui L RWs al k-esimo espeimento
		DIST2.SC(k,DIST2.GC(k)+pow(RW3D.GC(k),2));
		RW3D.SC(k,0.);//riazzero per il successivo blocco di RWs
	}
    }
	
   for(int i=0; i<N; i++){
	DIST.SC(i,DIST.GC(i)/N);
	DIST2.SC(i,DIST2.GC(i)/N);
	errprog.SC(i,calcolaErrore(DIST,DIST2,i)/sqrt(N-1));//incertezza statistica
   }
 
   Stampa("RW3Ddiscreto.txt",DIST,errprog,1);

   //passiamo ora al caso continuo 

   for(int i=0; i<N; i++){//riazzero i valori per riutilizzarli
	DIST.SC(i,0);
	DIST2.SC(i,0);
   }	

   for(int i=0; i<N; i++){ //numero dei blocchi
	for(int j=0; j<L; j++){//numero di RWs per blocco
		Posizione p;
		for(int l=0; l<N; l++){//numero dei passi
			p.RWcontinuo(rnd.Longitudine(),rnd.Latitudine(),a);
			temp.SC(l,pow(p.Distanza(ORI),2));
			RW3D.SC(l,RW3D.GC(l)+temp.GC(l));
			}
		}
	for(int k=0; k<N; k++){
		RW3D.SC(k,sqrt(RW3D.GC(k)/L));
		DIST.SC(k,DIST.GC(k)+RW3D.GC(k));//la somma delle radici di moduli quadri mediati sui L RWs al k-esimo espeimento
		DIST2.SC(k,DIST2.GC(k)+pow(RW3D.GC(k),2));
		RW3D.SC(k,0.);//riazzero per il successivo blocco di RWs
	}
    }
	
   for(int i=0; i<N; i++){
	DIST.SC(i,DIST.GC(i)/N);
	DIST2.SC(i,DIST2.GC(i)/N);
	errprog.SC(i,calcolaErrore(DIST,DIST2,i)/sqrt(N-1));//incertezza statistica
   }
 
   Stampa("RW3Dcontinuo.txt",DIST,errprog,1);



   

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
