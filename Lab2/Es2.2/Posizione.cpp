#include "Posizione.h"
#include <cmath>

Posizione::Posizione(){
	m_x = 0.;
	m_y = 0.;
	m_z = 0.;
}
Posizione::Posizione(double x, double y, double z){
	m_x = x;
	m_y = y;
	m_z = z;
}
double Posizione::GetX() const{
	return m_x;
}
double Posizione::GetY() const{
	return m_y;
}
double Posizione::GetZ() const{
	return m_z;
}
double Posizione::GetR() const{
	return sqrt(pow(m_x,2) + pow(m_y,2) + pow(m_z,2));
}
double Posizione::GetTheta() const{
	return acos(m_z/GetR());
}
double Posizione::GetPhi() const{
	return atan2(m_y,m_x);
}
double Posizione::GetRho() const{
	return sqrt(pow(m_x,2) + pow(m_y,2));
}
double Posizione::Distanza(const Posizione& p) const{
	return sqrt (pow(GetX()-p.GetX(),2) + pow(GetY()-p.GetY(),2) + pow(GetZ()-p.GetZ(),2));
}
void Posizione::RWdiscreto(int asse, int verso, double a){
	if(asse==1){//asse x
		m_x+=a*verso;
	}
	if(asse==2){// asse y
		m_y+=a*verso;
	}
	if(asse==3){// asse z
		m_z+=a*verso;
	}
}
void Posizione::RWcontinuo(double Longitudine, double Latitudine, double a) {
  m_x+=a*sin(Latitudine)*cos(Longitudine);
  m_y+=a*sin(Latitudine)*sin(Longitudine);
  m_z+=a*cos(Latitudine);
}




