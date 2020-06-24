#ifndef __NVT__
#define __NVT__

#include <cmath>

#include "Posizione.h"
#include "FunzioneBase.h"
#include "FunOndaHIGS.h" 
#include "FunOndaHI2P.h" 
//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//constants with Bohr Radius

double epsilon0 = 8.854187817*pow(10,-12); // dielectric constant
double e = -1.602176634*pow(10,-19); // electron charge
double hbar = 1.054571817*pow(10,-34); // h bar
double me = 9.1093837015*pow(10,-31); // electron mass
double a0 = (4*M_PI*epsilon0*pow(hbar,2))/(me*pow(e,2)); // Bohr radius

// Hydrogen atom wave functions : ground state and excited level 2p
FunzioneBase* Wave_Fun_HA;

//"bool" variables
int wave_function; //which wave function
int transition_probability; //which probability transition


//parameters, observables
const int m_props=100; 
int ir, n_props;
double walker[m_props];


// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_radius,err_radius;

//configuration ando origin of frame of reference
Posizione Pos;
Posizione Ori; // default constructor: x_Ori = 0 , y_Ori = 0 , z_Ori = 0 

// simulation
int nequistep, nstep, nblk;
double delta, sigma, x_start, y_start, z_start;


//functions
void Input(void);
void Reset(int);
void Equilibration(int);
void Move(void);
void Measure(void);
void Accumulate(void);
void Averages(int);
void ConfXYZ(void);
double Error(double,double,int);

#endif

