#ifndef __FunzioneBase_h__
#define __FunzioneBase_h__

#include<iostream>

class FunzioneBase{ // classe astratta per generica funzione
	public:
		virtual double Eval (double x) const = 0;
		virtual double EvalDerSec (double x) const = 0;
};

#endif
