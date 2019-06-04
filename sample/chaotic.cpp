#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "chaotic" // model name

// ---------- box size & step h ------------
#define XMIN 0
#define XMAX 15
#define PMIN (-1e-4)
#define PMAX 0
#define HX 0.01
#define HPMIN (1e-6)
#define HPOV (1./20)
// -----------------------------------------

// ---------- for PDE ----------------------
#define MAXSTEP 10000 // max recursion
#define TOL (1e-10) // tolerance
// -----------------------------------------

// ---------- potential parameter ----------
#define MPHI (1e-5)
// -----------------------------------------

#define RHOC (MPHI*MPHI) // end of inflation

// ---------- for SDE ----------------------
#define RECURSION 100 // recursion for power spectrum
#define PHIIN 13 // i.c. for phi
#define PIN -(5e-6) // i.c. for pi
#define TIMESTEP (1e-2) // time step : delta N
// -----------------------------------------

// ---------- for power spectrum -----------
#define DELTAN 0.1 // calc. PS every DELTAN e-folds
#define NMAX 40 // calc. PS for 0--NMAX e-folds
// -----------------------------------------


int main(int argc, char** argv)
{
  // ---------- start stop watch ----------
  struct timeval tv;
  struct timezone tz;
  double before, after;
  
  gettimeofday(&tv, &tz);
  before = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  // --------------------------------------


  // ---------- set box ---------------
  double h = HX, sitev = XMIN;
  vector<double> site;
  vector< vector<double> > xpsite;
  vector< vector< vector<double> > > sitepack;
  while(sitev <= XMAX) {
    site.push_back(sitev);
    sitev += h;
  }
  xpsite.push_back(site);
  sitepack.push_back(xpsite);
  site.clear();
  xpsite.clear();

  sitev = PMIN;
  while (sitev <= PMAX) {
    h = max(fabs(sitev)*HPOV,HPMIN);

    site.push_back(sitev);
    sitev += h;
  }
  xpsite.push_back(site);
  sitepack.push_back(xpsite);
  site.clear();
  xpsite.clear();
  // ----------------------------------

  vector<double> params = {MAXSTEP,TOL,2,RHOC,(double)sitepack[0].size(),TIMESTEP,NMAX,DELTAN,RECURSION}; // set parameters

  vector< vector<double> > xpi = {{PHIIN},{PIN}}; // set i.c. for inflationary trajectories

  StocDeltaN sdn(MODEL,sitepack,xpi,0,params); // declare the system

  sdn.sample(); // obtain 1 sample path
  //sdn.sample_plot(); // plot obtained sample path
  
  //sdn.solve(); // solve PDE & SDE to obtain power spectrum
  //sdn.f_logplot(0); // show plot of <N>
  //sdn.f_logplot(1); // show plot of <delta N^2>
  //sdn.calP_plot(); // show plot of power spectrum of zeta
  

  // ---------- stop stop watch ----------
  gettimeofday(&tv, &tz);
  after = (double)tv.tv_sec + (double)tv.tv_usec * 1.e-6;
  cout << after - before << " sec." << endl;
  // -------------------------------------
}


// ---------- Lagrangian params. X[0]=phi ----------

double StocDeltaN::V(vector<double> &X)
{
  return 1./2*MPHI*MPHI*X[0]*X[0];
}

double StocDeltaN::VI(vector<double> &X, int I) // \partial_I V
{
  return MPHI*MPHI*X[0];
}

double StocDeltaN::metric(vector<double> &X, int I, int J) // G_IJ
{
  return 1;
}

double StocDeltaN::inversemetric(vector<double> &X, int I, int J) // G^IJ
{
  return 1;
}

double StocDeltaN::affine(vector<double> &X, int I, int J, int K) // \Gamma^I_JK
{
  return 0;
}

double StocDeltaN::derGamma(vector<double> &X, int I, int J, int K, int L) // \partial_L \Gamma^I_JK
{
  return 0;
}

