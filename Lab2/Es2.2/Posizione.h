#ifndef __Posizione_h__
#define __Posizione_h__

#include <iostream>
using namespace std;

class Posizione {
	public:	
		Posizione(); //costruttore di default
		Posizione(double,double,double);

		double GetX() const;//coordinate cartesiane
		double GetY() const;
		double GetZ() const;

		double GetR() const;//coordinate sferiche 
		double GetTheta() const;
		double GetPhi() const;

		double GetRho() const; //coordinate cilindriche
		double Distanza(const Posizione&) const;// distanza da un punto
		void RWdiscreto(int asse, int verso, double a); // mi indica dove avanzare randomicamente nel reticolo dicreto di passo reticolare a
		void RWcontinuo(double Longitudine, double Latitudine, double a); //mi indica come cambia la posizione randomicamente uin na direzione casuale uniformemnte distribuita sull'angolo solido per un passo a
	private:
		double m_x, m_y, m_z;
};
#endif
