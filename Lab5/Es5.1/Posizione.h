#ifndef __Posizione_h__
#define __Posizione_h__

#include <iostream>
#include "random.h"
#include "FunzioneBase.h"
using namespace std;

class Posizione {
	public:	
		Posizione(); //costruttore di default
		Posizione(double,double,double);

		double GetX() const;//coordinate cartesiane
		double GetY() const;
		double GetZ() const;

		void SetX(double);
		void SetY(double);
		void SetZ(double);

		double GetR() const;//coordinate sferiche 
		double GetTheta() const;
		double GetPhi() const;

		double GetRho() const; //coordinate cilindriche
		double Distanza(const Posizione&) const;// distanza da un punto

	private:
		double m_x, m_y, m_z;
};
#endif
