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

   int M = 1000000; //numero totale di variabili casuali generate 
   int N = 100; //numero di blocchi(o esperimenti)
   int L = M/N;//numero di variabili casuali in ciascun blocco
   Vettore r(M);
   Vettore sumprog(N); //vettore che contiene i valori medi cumulativi progressivi
   Vettore su2prog(N); //vettore che contiene il quadrato dei valori medi cumulativi progressivi
   Vettore errprog(N); //vettore che contiene le incertezze statistiche progressive
   Vettore ave(N); //valori medi per ogni blocco
   Vettore av2(N);// quadrato dei valori medi per ogni blocco

   cout << "This program check the properties of Pseudo-Random Number generator with 3 tests: " << endl;
   cout << "1) Average value for a Random variable uniformly ditributed in the interval [0,1) " << endl;
   cout << "2) Variance for a Random variable uniformly ditributed in the interval [0,1) " << endl;
   cout << "3) Chi quadro " << endl;
 
//riempio il vettore r con il generatore di numeri casuali uniformemente distribuiti in [0,1)	
   for(int i=0; i<r.GetN(); i++){
   	r.SC(i,rnd.Rannyu());
   }

   for(int i=0; i<N; i++){
	double sum=0;
	for(int j=0; j<L; j++){
		int k = j+i*L;
		sum += r.GC(k);
   	}
	ave.SC(i,sum/L);
	av2.SC(i,pow(ave.GC(i),2));
   }

   for(int i=0; i<N; i++){
	double sum1=0;
	double sum2=0;
   	for(int j=0; j<i+1; j++){
		sumprog.SC(i,sum1 += ave.GC(j)); // somma progressiva dei valori medi
		su2prog.SC(i,sum2 += av2.GC(j));// somma progressiva del quadrato dei valori medi
	}
	sumprog.SC(i,(sumprog.GC(i))/(i+1));//valore medio cumulativo
	su2prog.SC(i,(su2prog.GC(i))/(i+1));//valore medio quadro cumulativo
	errprog.SC(i,calcolaErrore(sumprog,su2prog,i));//incertezza statistica
   }
   
  
   Stampa("ave.txt",sumprog ,errprog,1);
   


//Ripetiamo lo stesso procedimento per la varianza

   for(int i=0; i<N; i++){
	double sum=0;
	for(int j=0; j<L; j++){
		int k = j+i*L;
		sum += pow(r.GC(k)-(0.5),2);
   	}
	ave.SC(i,sum/L);
	av2.SC(i,pow(ave.GC(i),2));
   }

   for(int i=0; i<N; i++){
	double sum1=0;
	double sum2=0;
   	for(int j=0; j<i+1; j++){
		sumprog.SC(i,sum1 += ave.GC(j));
		su2prog.SC(i,sum2 += av2.GC(j));
	}
	sumprog.SC(i,(sumprog.GC(i))/(i+1));
	su2prog.SC(i,(su2prog.GC(i))/(i+1));
	errprog.SC(i,calcolaErrore(sumprog,su2prog,i));//incertezza statistica della varianza
   }

   Stampa("variance.txt",sumprog ,errprog,1);

//Test del chi^2

   M=100; //numero di sotto-intervalli dell'intervallo [0,1)
   N=10000; // numeri casuali distribuiti uniformemente nell'intervallo [0,1) per la stima del j-esimo valore del Chi^2
   int R = M*M;
   Vettore conteggio(R);//conteggio cumulativo dei numeri generati che cade nell'M-esimo intervallo per ognuno dei 100 test
   Vettore chi_2(M);
   double temp;

   	for (int i=0; i<M; i++ ){ // ciclo sugli M tests
		for (int l=0; l<M; l++ ){ // ciclo sugli M sottointervalli
			int count = 0;
			for (int j=0; j<N; j++ ){ // genero N numeri casuali
				temp = rnd.Rannyu();
				int k=l+i*M;
 				if((temp > double (l)/M)&&(temp < double(l+1)/M )){
					count++;
					conteggio.SC(k,count);
				}
			}
		}
	}
   	for (int i=0; i<chi_2.GetN(); i++ ){
		for (int l=0; l<M; l++ ){
			double sum3=0;
			for (int j=0; j<M; j++){
				int k = j+i*M;
				chi_2.SC(i,sum3+=pow((conteggio.GC(k)-N/M), 2)/(double(N/M)));
			}
   		}
	}
  
   Stampa("chi_2.txt",chi_2,1);


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
