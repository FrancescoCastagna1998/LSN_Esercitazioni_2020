#include "TrialWaveFunction.h"
#include <cmath>

double TrialWaveFunction::Eval (double x) const {
	//esponenti della Gaussiana
	double esp_sum = (-1.)*pow((x + m_mu),2)/(2*pow(m_sigma,2));
	double esp_diff = (-1.)*pow((x - m_mu),2)/(2*pow(m_sigma,2)); 

	return exp(esp_diff) + exp(esp_sum);
}

double TrialWaveFunction::EvalDerSec (double x) const {
	// esponenti della Gaussiana
	double esp_sum = (-1.)*pow((x + m_mu),2)/(pow(m_sigma,2));
	double esp_diff = (-1.)*pow((x - m_mu),2)/(pow(m_sigma,2));
	double a = 1/(pow(m_sigma,2));

	return exp(esp_diff/2.)*((-1)*esp_diff*a-a) + exp(esp_sum/2.)*((-1.)*esp_sum*a-a);
}
