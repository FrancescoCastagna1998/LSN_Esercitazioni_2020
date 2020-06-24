
#include "Travelling_Salesman_Problem_TSP.h"
#include "funzioni.h"


using namespace std;

int main (int argc, char *argv[]){

   Input(); //initialitation
   
   while(evolution_temp >= end_temp){

	accepted = 0.0;
	attempted = 0.0;
	for(int istep=1; istep<=nsteps; istep++)
	{
		Move_Mutation(); //Move the individual with mutations
		if (!Evolution_Individual.Check()) cerr << "PROBLEM: The population doesn't fulfil the bonds of the problem" << endl << endl; //check
		PrintResults(); 	
   	}
	Restart();
   }

   PrintBest();
			  
   rnd->SaveSeed();
   return 0;
}


void Input(void){

  ifstream ReadInput;

  //Read seed for random numbers
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
            rnd->SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   cout << "This program solve the Travelling salesman Problem with a Simulated Annealing algorithm" << endl << endl;
   cout << "The bonds are: the salesman must visit one and only one every city and must be back to the fisrt city at the end of the path!" << endl << endl;
   cout << "Monte Carlo simulation             " << endl << endl;
   cout << "The Metropolis algorithm proposes a genetic mutation as move " << endl;
   cout << "Boltzmann weight exp(- beta * (L(x') - L(x))), beta = 1/T " << endl << endl;
   cout << "The program uses k_b = 1 unit " << endl << endl;


  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> ncity;
  cout << "The path (individual) is made of " << ncity << " cities" << endl;

  ReadInput >> circle_or_square;

  ReadInput >> starting_temp;

  if (circle_or_square==1) beta = 1.0/starting_temp;
  if (circle_or_square==0) beta = 1.0/starting_temp;

  cout << "The intial fictitious temperature: " << starting_temp << " and its inverse beta : "  << beta << endl;
 
  ReadInput >> end_temp; 	

  cout << "The final fictitious temperature in annealing schedule: " << end_temp << " and its inverse beta : "  << 1.0/end_temp << endl;

  ReadInput >> nsteps; //Monte Carlo steps
  cout << "The generations (MC-steps) produced by tha algorithm for each temperature : " << nsteps << endl << endl;

  ReadInput >> xstart_city; 
  ReadInput >> ystart_city;

//The starting and returning city of each paths
  if (circle_or_square==1)  //city randomly placed on a circumference
  {
	starting_individual.NewCity(xstart_city,ystart_city,1);               
	double r = 1.; // circumference's radius
	double theta;

	cout << "The cities are randomly placed on a circumference of radius " << r << endl;
	cout << "The position of the first city on the plane is: x = " << xstart_city << " y = " << ystart_city  << endl << endl;

// Generate the starting individual
	for(int icity=2; icity<=ncity; icity++) {
		double x=0.;
		double y=0.;	
		theta = rnd->Rannyu(0., 2.*M_PI);
		x = r*cos(theta);
		y = r*sin(theta);
		starting_individual.NewCity(x, y, icity);
   	}
   }
   else //city randomly placed inside a square
   {
	double l = 2.; // square's side
	cout << "The cities are randomly placed inside a square of side " << l << endl;
// Generate the starting individual
	for(int icity=1; icity<=ncity; icity++) {
		double x=0.;
		double y=0.;	
		x = rnd->Rannyu(-l/2.,l/2.);
		y = rnd->Rannyu(-l/2.,l/2.);
		starting_individual.NewCity(x, y, icity);
   	}
   }


   ReadInput.close();


//shuffle the order of the cities of the starting individual
   starting_individual.Shuffle();

// Check if the starting individual fulfils the bonds
   if (starting_individual.Check()) cout << "Correct initialization" << endl << endl;

// Print the starting individual
   starting_individual.PrintInd("Starting_Individual.dat");

   if (circle_or_square==1)  //city randomly placed on a circumference
   {
   	ofstream Length ("Individual_Evolution_circle.dat",ios::app);
   	int wd=12;
   	Length << setw(wd) << 0 << setw(wd) << starting_individual.CostFun() << endl;
   	Length.close();
   }
   if (circle_or_square==0) //city randomly placed inside a square
   {
   	ofstream Length ("Individual_Evolution_square.dat",ios::app);
   	int wd=12;
   	Length << setw(wd) << 0 << setw(wd) << starting_individual.CostFun() << endl;
   	Length.close();
   }


//copy the starting_individual in Evolution_Individual
   for(int icity=0; icity<ncity; icity++){
   	Evolution_Individual.NewCity(starting_individual.GetCity(icity)); 
   }

   evolution_temp = starting_temp;
   beta = 1.0/evolution_temp;
   
   igen = 0;
   iblk = 0;
}
void Move_Mutation(void)
{
   Individual Old_Individual;
   Individual New_Individual;
   double p, Length_Old, Length_New;
   double r;

   for(int icity=0; icity<ncity; icity++){
	Old_Individual.NewCity(Evolution_Individual.GetCity(icity));
   	New_Individual.NewCity(Evolution_Individual.GetCity(icity)); 
   }
   r = rnd->Rannyu();
//first mutation
   if (r <= 1./3.)
   {
	New_Individual.Gen_Mut_1((int)(ncity*(rnd->Rannyu())),(int)(ncity*(rnd->Rannyu())));
	New_Individual.Gen_Mut_1((int)(ncity*(rnd->Rannyu())),(int)(ncity*(rnd->Rannyu())));
   }	
//third mutation
   if ((r > 1./3.) && (r <= 2./3.))
   {
	New_Individual.Gen_Mut_2((int)(ncity*(rnd->Rannyu())),rnd,(int)((ncity/2-1)*(rnd->Rannyu())));
	New_Individual.Gen_Mut_2((int)(ncity*(rnd->Rannyu())),rnd,(int)((ncity/2-1)*(rnd->Rannyu())));
   }	
//forth mutation	
   if (r > 2./3.)
   {
	New_Individual.Gen_Mut_3((int)(ncity*(rnd->Rannyu())),(int)((ncity-1)*(rnd->Rannyu())));
	New_Individual.Gen_Mut_3((int)(ncity*(rnd->Rannyu())),(int)((ncity-1)*(rnd->Rannyu())));
   }

//New and old "energy of the system"
   Length_Old = Old_Individual.CostFun();
   Length_New = New_Individual.CostFun();

//Metropolis test
   p = min(1.,exp(beta*(Length_Old-Length_New)));

   if(p >= rnd->Rannyu()){
	for(int icity=0; icity<ncity; icity++){
		Evolution_Individual.SetCity(icity,New_Individual.GetCity(icity));
	} 
	accepted = accepted + 1.0;
   }
   else{
	for(int icity=0; icity<ncity; icity++){
		Evolution_Individual.SetCity(icity,Old_Individual.GetCity(icity));
	}
   }  
   attempted = attempted + 1.0;
}
void PrintResults(void)
{
//Print the length of the best path
   if (igen%50 == 0)
   {
	if (circle_or_square==1)  //city randomly placed on a circumference
	{
		ofstream Length ("Individual_Evolution_circle.dat",ios::app);
		int wd=12;
		Length << setw(wd) << igen + 1 << setw(wd) << Evolution_Individual.CostFun() << endl;
		Length.close();
	}
	if (circle_or_square==0) //city randomly placed inside a square
	{
		ofstream Length ("Individual_Evolution_square.dat",ios::app);
		int wd=12;
		Length << setw(wd) << igen + 1 << setw(wd) << Evolution_Individual.CostFun() << endl;
		Length.close();
	}
   }

   igen = igen + 1.0; 
}
void Restart(void)
{
    cout << "Block number " << iblk << endl;
    cout << "Temperature T = " << evolution_temp << endl;
    cout << "Beta = " << beta << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;

//Annealing Schedule
    if (circle_or_square==1)  //city randomly placed on a circumference
    {
	if (evolution_temp >= 0.5) evolution_temp = evolution_temp - 0.05; // Low the temperature
	if ((evolution_temp < 0.5) && (evolution_temp >= 0.1)) evolution_temp = evolution_temp - 0.02;
	if ((evolution_temp < 0.1) && (evolution_temp >= 0.01)) evolution_temp = evolution_temp - 0.01;
	if ((evolution_temp < 0.01) && (evolution_temp >= end_temp)) evolution_temp = evolution_temp - 0.002;
    }
    else //city randomly placed inside a square
    {
	if (evolution_temp >= 0.8) evolution_temp = evolution_temp - 0.05; // Low the temperature
	if ((evolution_temp < 0.8) && (evolution_temp >= 0.2)) evolution_temp = evolution_temp - 0.02;
	if ((evolution_temp < 0.2) && (evolution_temp >= 0.01)) evolution_temp = evolution_temp - 0.01;
	if ((evolution_temp < 0.01) && (evolution_temp >= end_temp)) evolution_temp = evolution_temp - 0.002;
    }

    beta = 1.0/evolution_temp;

    cout << "----------------------------" << endl << endl;

    iblk = iblk + 1.0;
}

void PrintBest(void){
//city randomly placed on a circumference
	if (circle_or_square==1) Evolution_Individual.PrintInd("Best_Individual_circle.dat");
//city randomly placed inside a square  
	if (circle_or_square==0) Evolution_Individual.PrintInd("Best_Individual_square.dat");
}



