// 3-step strong order 2 stochastic Runge-Kutta solver

#ifndef INCLUDED_SRK32_hpp_
#define INCLUDED_SRK32_hpp_

#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <time.h>

using namespace std;

class SRKintegrater
{
protected:
  double t,t0, // time t and its initial value t0
    A0[3][3],A1[3][3],B0[3][3],B1[3][3],Alpha[3],Beta1[3],Beta2[3],C0[3],C1[3]; // params. for stochastic Runge-Kutta
  vector<double> dW; // dW[alpha], alpha : noise d.o.f.
  vector< vector<double> > xx,xxi,aa[3],H0[3]; // xx[xp][I], xp : phi or pi, I : fields d.o.f.
  // xx, xxi : phase-space value and its initial condition
  //aa, H0 : for SRK
  vector< vector<double> > aIs,vIs; // aIs[I][alpha] only for Idim == noisedim
  // for calc. of adiabatic direction
  vector< vector< vector<double> > > bb[3],Hk[3]; //bb[3][alpha][xp][I]
  // bb, Hk : for SRK
  int xpdim,Idim,noisedim; // xpdim = 1 or 2, Idim : fields d.o.f., noisedim : noise d.o.f.

public:
  SRKintegrater(){}
  SRKintegrater(vector< vector<double> > &XPi, double T0, int NoiseDim); // XPi : initial phase space value, T0 : initial time, NoiseDim : noise d.o.f.
  void SRK2(double dt); // proceed SRK 1 step
  void coeff(double dt, int step); // set coefficients of SDE
  double e1(vector<double> &X, vector<double> &P); // first SR param. -\dot{H}/H^2
  double return_H(); // return current Hubble 
  double return_V(); // return_current potential
  double return_t(); // return current time
  double return_xp(int xp, int I); // return current phase-space value
  double return_e1(); // return current -\dot{H}/H^2
  double vielbein(vector< vector<double> > &XP, int I, int alpha); // projection from field coordinate I to noise frame alpha
  double eIsigma(vector< vector<double> > &XP, int I); // projection to adiabatic direction
  double eIs(vector< vector<double> > &XP, int I, int alpha); // projection to entropic direction
  
  virtual double H(vector<double> &X, vector<double> &P);
  virtual double V(vector<double> &X);
  virtual double VI(vector<double> &X, int I);
  virtual double metric(vector<double> &X, int I, int J);
  virtual double inversemetric(vector<double> &X, int I, int J);
  virtual double affine(vector<double> &X, int I, int J, int K);
  virtual double derGamma(vector<double> &X, int I, int J, int K, int L);

  // ---------------------------------------------------
  // solve dxI = DI dN + gIa dWa
  // difusion term DIJ = gIa gJa
  virtual double DI(int xp, int I, vector< vector<double> > &psv);
  virtual double DIJ(int xpI, int I, int xpJ, int J, vector< vector<double> > &psv);
  virtual double gIa(int xp, int I, int alpha, vector< vector<double> > &psv);
  // ---------------------------------------------------
};

double Uniform();
double rand_normal(double mu, double sigma);

#endif
