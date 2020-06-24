#include "FunzioneBase.h"

class FunOndaHI2P: public FunzioneBase {

	public:

		FunOndaHI2P() { m_a=0; m_b=0;};
		FunOndaHI2P(double a, double b) {m_a=a; m_b=b;};
		
		~FunOndaHI2P(){};

		double Geta() const {return m_a; };
		double Getb() const {return m_b; };

		void Seta(double a) {m_a=a;};
		void Setb(double b) {m_b=b;};
		
		virtual double Eval (double,double,double) const;
		virtual double Eval_squared_modulus (double,double,double) const; 	


	private: 
		double m_a, m_b; //normalization and esponent constants
};
