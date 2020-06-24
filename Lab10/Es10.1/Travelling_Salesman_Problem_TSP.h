#ifndef __TSP__
#define __TSP__

#include <cmath>
#include <algorithm>
#include <iomanip>
#include "random.h"
#include "city.h"
#include "Individual.h"



//Random numbers
int seed[4];
Random* rnd = new Random;

//parameters
int ncity,nsteps,igen,iblk;
double starting_temp, end_temp, evolution_temp,beta;
//probability
double prob_gen_mut1,prob_gen_mut2, prob_gen_mut3;
double accepted,attempted;


//position of the first city randomly placed on a circumfernce
double xstart_city;
double ystart_city;

//cities geometry
int circle_or_square;

//Individual
Individual starting_individual;
Individual Evolution_Individual;




//functions
void Input(void);
void Move_Mutation(void);
void PrintResults(void);
void PrintBest(void);
void Restart(void);


#endif

