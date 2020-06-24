#include "FunzioneBase.h"
#include <cmath>

class CosenoIS: public FunzioneBase {

	public:

		CosenoIS() { m_a=0; m_b=0;};
		CosenoIS(double a, double b) {m_a=a; m_b=b;};
		
		~CosenoIS(){};

		double Geta() const {return m_a; };
		double Getb() const {return m_b; };

		void Seta(double a) {m_a=a;};
		void Setb(double b) {m_b=b;};
		
		double Eval (double x) const {return (m_a*cos(m_b*x))/(1.-x);};	


	private: 
		double m_a, m_b;
};
