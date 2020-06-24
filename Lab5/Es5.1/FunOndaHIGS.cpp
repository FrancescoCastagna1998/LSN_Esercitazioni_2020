#include "FunOndaHIGS.h"
#include <cmath>

double FunOndaHIGS::Eval (double x, double y, double z) const {

   return m_a*exp((-1.)*m_b*sqrt(x*x+y*y+z*z));
}

double FunOndaHIGS::Eval_squared_modulus (double x,double y,double z) const {

   return pow(m_a*exp((-1.)*m_b*sqrt(x*x+y*y+z*z)),2);
}

