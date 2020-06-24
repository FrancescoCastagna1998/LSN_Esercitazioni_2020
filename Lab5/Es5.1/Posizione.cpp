#include "Posizione.h"
#include <cmath>
#include <algorithm>

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
void Posizione::SetX(double x){
	m_x = x;
}
void Posizione::SetY(double y){
	m_y = y;
}
void Posizione::SetZ(double z){
	m_z = z;
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
  	
	





