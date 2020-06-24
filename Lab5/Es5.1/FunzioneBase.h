#ifndef __FunzioneBase_h__
#define __FunzioneBase_h__

#include<iostream>

class FunzioneBase{ // classe astratta per generica funzione
	public:
		virtual double Eval (double x,double y, double z) const = 0;
		virtual double Eval_squared_modulus (double x,double y, double z) const = 0;
};

#endif
