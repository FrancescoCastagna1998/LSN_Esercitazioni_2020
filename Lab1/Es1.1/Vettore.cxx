#include "Vettore.h"
//#define NDEBUG
#include <assert.h>
#include <cmath>

Vettore::Vettore() {
	m_N = 0;
	m_v = NULL;
}
Vettore::Vettore(int N) {
	assert((N >=0 ) && "Errore : la dimensione deve essere positiva");	
	if( N < 0 ) {
		cerr << "Errore: la dimensione del vettore deve essere positiva" << endl;
		exit(1);
	} else {
		m_N = N;
		m_v = new double[N];
		for (int i =0; i<N; i++)
			m_v[i] = 0;	
	}
}
Vettore::~Vettore() {
	delete[] m_v;
}
void Vettore::SC (int i, double a){
	assert(( m_N > i ) && "Errore: l'indice e' troppo grande");
	if( i<m_N ) {
		m_v[i] = a;
	}else {
		cerr << "Errore: indice" << i << " , dimensione" << m_N << endl;
		exit(2);
	}
}
double Vettore::GC (int i) const{
	assert(( m_N > i ) && "Errore: l'indice e' troppo grande");
	if( m_N > i) {
		return m_v[i];
	}else {
		cerr << "Errore: indice" << i << " , dimensione" << m_N << endl;
		exit(3);
	}
}
void Vettore::Scambia (int i, int j){
	double temp = GC(i);
	SC(i,GC(j));
	SC(j,temp);
}
Vettore::Vettore(const Vettore& V){
	m_N = V.GetN();
	m_v = new double [m_N];
	for (int i =0; i < m_N; i++)
			m_v[i] = V.GC(i);	
}
Vettore& Vettore::operator=(const Vettore& V){
	m_N = V.GetN();
	if ( m_v ) 
		delete[] m_v;
	m_v = new double [m_N];
	for (int i =0; i < m_N; i++)
			m_v[i] = V.GC(i);
	return *this; // utile per restituire una copia dell'oggetto
}
double& Vettore::operator[](int i) const {
	if ( i < m_N ) {
		return m_v[i];
	} else {
		cerr << "Errore: indice " << i << " e dimensione " << m_N << endl;
		exit(4);
	}
}
double Vettore::Media() const {
	double s = 0;
	for (int i=0; i<m_N; i++) s += m_v[i];
	return s/m_N;
}

double Vettore::DevSt() const {
	double s = 0;
	for (int i = 0; i<m_N; i++) s += pow(m_v[i]-Media(), 2);
	return pow( s/(m_N-1), 0.5);
}
