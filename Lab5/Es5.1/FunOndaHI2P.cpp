#include "FunOndaHI2P.h"
#include <cmath>

double FunOndaHI2P::Eval (double x, double y, double z) const {

   return m_a*z*exp((-1.)*m_b*sqrt(x*x+y*y+z*z));
}

double FunOndaHI2P::Eval_squared_modulus (double x,double y,double z) const {

   return pow(m_a*z*exp((-1.)*m_b*sqrt(x*x+y*y+z*z)),2);
}

