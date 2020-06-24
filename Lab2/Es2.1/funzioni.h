#include <iostream>
#include <cstdlib>
#include <fstream>
#include "Vettore.h"
#include <cmath>
#include <iomanip>

double calcolaErrore( const Vettore &,const Vettore &, int ) ; //calcolo per l'incertezza statistica

void Stampa ( const char* Filename ,const Vettore & v,const Vettore & w, int n);
void Stampa ( const char* Filename ,const Vettore & v);
void Stampa ( const Vettore &);
