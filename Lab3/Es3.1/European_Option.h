#ifndef __EO__
#define __EO__

//Random numbers
#include "random.h"

int seed[4];
Random rnd;

//parameters, observables
const int n_props = 4;
int iput1, iput2 , icall1, icall2;
double S_T1, S_T2;
double Option_Price[n_props];

// averages
double blk_av[n_props],blk_norm;
double glob_av[n_props],glob_av2[n_props];
double stima_call1,stima_put1,stima_call2,stima_put2,err_call1,err_put1,err_call2,err_put2;

//European Option parameters
double delivery_time, strike_price, r, volatility, starting_asset_price,deltaT;
int ninterval;
//simulation
int nblk, nstep;



//functions
void Input(void);
void Reset(int);
void Move1(void);
void Move2(void);
void Measure(void);
void Accumulate(void);
void Averages(int);
double Error(double,double,int);

#endif
