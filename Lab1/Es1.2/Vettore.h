#pragma once //#ifndef__Vettore_h    #define__Vettore_h 

#include <iostream>
#include <cstdlib>

using namespace std;

class Vettore {
	public:
		Vettore();		// Costruttore di default
		Vettore(int N);		// Costruttore con dimensione del vettore

		~Vettore();		//Distruttore

		Vettore(const Vettore&); // Copy Constructor
		Vettore& operator=(const Vettore&); // overloading operatore di assegnazione
		double& operator[](int) const;

		int GetN() const { return m_N;}; // accede ala dimensione del vettore
		void SC (int, double); // modifica la componente i-esima
		double GC (int) const; // accede ala componente i-esima

		void Scambia( int, int );
		
		double Media() const;
  		double DevSt() const;


	protected:

		int m_N;
		double* m_v;
};
//#endif //__Vettore_h 



