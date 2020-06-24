#include "FunzioneBase.h"
#include <cmath>

class Coseno: public FunzioneBase {

	public:

		Coseno() { m_a=0; m_b=0;};
		Coseno(double a, double b) {m_a=a; m_b=b;};
		
		~Coseno(){};

		double Geta() const {return m_a; };
		double Getb() const {return m_b; };

		void Seta(double a) {m_a=a;};
		void Setb(double b) {m_b=b;};
		
		double Eval (double x) const {return m_a*cos(m_b*x);};	


	private: 
		double m_a, m_b;
};
