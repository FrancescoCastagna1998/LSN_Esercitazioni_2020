#include "FunzioneBase.h"
#include <cmath>

class Potential: public FunzioneBase {

	public:

		Potential() {m_a=0;};
		Potential(double a) {m_a=a;};
		
		~Potential(){};

		double Get() const {return m_a;};

		void Set(double a) {m_a=a;};
		
		virtual double Eval (double x) const {return pow(x,4) + m_a*pow(x,2) ;};

		virtual double EvalDerSec (double x) const {return 12*pow(x,2) + 2*m_a ;}; // second derivative of the external Potential	


	private: 
		double m_a; // x^2 coefficient
};
