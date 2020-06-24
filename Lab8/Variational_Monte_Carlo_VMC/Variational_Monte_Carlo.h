#ifndef __VMC__
#define __VMC__


#include "random.h"
#include "FunzioneBase.h"
#include "TrialWaveFunction.h" 
#include "Potential.h" 

//Random numbers
int seed[4];
Random rnd;

//variational parameters
double mu,sigma;

//Trial wave function
FunzioneBase* Psi_trial; 
//External potential
FunzioneBase* Ext_Pot;

//position
double x_pos;

//parameters, observables
const int m_props=1000;
int n_props, ih, ihist; //hamiltonian and histogram indexes
double bin_size,nbins, histogram_start, histogram_end;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_ham,err_ham, err_PsiT2;

// simulation
int nstep, nblk, nequi;
double delta;
double start_position;

//equilibration
int answer;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void Measure(void);
void Equilibration(int);
void ConfFinal(void);
double Error(double,double,int);

#endif

