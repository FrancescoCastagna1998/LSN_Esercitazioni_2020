#ifndef __Buffon__
#define __Buffon__

//Random numbers
#include "random.h"

int seed[4];
Random rnd;

//Buffon experiment parameters
double L, d, y_centre, Theta; 
int N_hit;

//parameters, observables
double pigreco;

//averages
double glob_av,glob_av2;
double stima_pigreco,err_pigreco;


//simulation
int nblk, nthrow;

//functions
void Input(void);
void Reset(int);
void Move(void);
void Measure(void);
void Averages(int);
double Error(double,double,int);

#endif
