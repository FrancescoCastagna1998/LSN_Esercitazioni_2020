#include "funzioni.h"
#include "Vettore.h"

double calcolaErrore( const Vettore & v,const Vettore & w, int n ) {
	if (n==0){
		return 0;
	}
	return sqrt(fabs((w.GC(n)-pow(v.GC(n),2))));
}
void Stampa ( const Vettore & v){
	for (int i=0; i<v.GetN(); i++) 
	cout << v.GC(i) << endl;
}

void Stampa ( const char* Filename ,const Vettore & v,const Vettore & w,int n){
	ofstream fout (Filename);
	for (int i=0; i<v.GetN(); i++){ 
		fout << setw(10) << n*(i+1) << setw(15) << setprecision(10) << v.GC(i) << setw(20) << setprecision(10) << w.GC(i)  << endl;
	}
	fout.close();
}
void Stampa ( const char* Filename ,const Vettore & v){
	ofstream fout (Filename);
	for (int i=0; i<v.GetN(); i++){ 
		fout << i << "  " << setprecision(10) << v.GC(i) << endl;
	}
	fout.close();
}
