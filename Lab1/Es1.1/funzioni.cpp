#include "funzioni.h"
#include "Vettore.h"

double calcolaErrore( const Vettore & v,const Vettore & w, int n ) {
	if (n==0){
		return 0;
	}
	return sqrt(fabs((w.GC(n)-pow(v.GC(n),2))/n));
}

void Stampa ( const Vettore & v){
	for (int i=0; i<v.GetN(); i++) 
	cout << v.GC(i) << endl;
}

void Stampa ( const char* Filename ,const Vettore & v,const Vettore & w,int n){
	ofstream fout (Filename);
	for (int i=0; i<v.GetN(); i++){ 
		fout << n*(i+1) << " " << setprecision(10) << v.GC(i) << " " << setprecision(10) << w.GC(i)  << endl;
	}
	fout.close();
}
void Stampa ( const char* Filename ,const Vettore & v, int n){
	ofstream fout (Filename);
	for (int i=1; i<v.GetN(); i++){ 
		fout << n*i << " " << setprecision(10) << v.GC(i) << endl;
	}
	fout.close();
}
