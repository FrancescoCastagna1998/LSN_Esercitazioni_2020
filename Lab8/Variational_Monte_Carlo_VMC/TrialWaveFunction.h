#include "FunzioneBase.h"
#include <cmath>

class TrialWaveFunction: public FunzioneBase {

	public:

		TrialWaveFunction() { m_mu=0; m_sigma=0;}; //constructor
		TrialWaveFunction(double mu, double sigma) {m_mu=mu; m_sigma=sigma;};
		
		~TrialWaveFunction(){}; //distructor

		double Getmu() const {return m_mu; };
		double Getsigma() const {return m_sigma; };

		void Setmu(double mu) {m_mu=mu;};
		void Setsigma(double sigma) {m_sigma=sigma;};
		
		virtual double Eval (double x) const; // wave function

		virtual double EvalDerSec (double x) const; // second derivative of the wave function	


	private: 
		// the two variational parameters
		double m_mu; // distance from the origin of the Gaussians average values
		double m_sigma; // the width of the Gaussians under square root
};
