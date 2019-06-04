#include "../source/StocDeltaN.hpp"
#include <sys/time.h>

#define MODEL "hilltop_SR" // model name

// ---------- box size & step h ------------
#define XMIN -21
#define XMAX 21
#define HMIN (1e-8)
#define HOV (1./20)
// -----------------------------------------

// ---------- for PDE ----------------------
#define MAXSTEP 100000000 // max recursion
#define TOL 1e-10 // tolerance
// -----------------------------------------

// ---------- potential parameter ----------
#define LAMBDA 0.01
#define MU 20
// -----------------------------------------

#define RHOC (LAMBDA*LAMBDA*LAMBDA*LAMBDA*(sqrt(1+2*MU*MU)-1)/MU/MU) // end of inflation

// ---------- for SDE ----------------------
#define RECURSION 100 // recursion for power spectrum
#define XIN (1e-6) // i.c. for phi
#define TIMESTEP (1e-1) // time step : delta N
// -----------------------------------------

// ---------- for power spectrum -----------
#define DELTAN 10 // calc. PS every DELTAN e-folds
#define NMAX 2800 // calc. PS for 0--NMAX e-folds
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
  double h, sitev = XMIN;
  vector<double> site;
  vector< vector<double> > xsite;
  vector< vector< vector<double> > > sitepack;
  while (sitev <= XMAX) {
    h = max(fabs(sitev)*HOV,HMIN);
    
    site.push_back(sitev);
    sitev += h;
  }
  xsite.push_back(site);
  sitepack.push_back(xsite);
  // ----------------------------------

  vector<double> params = {MAXSTEP,TOL,2,RHOC,(double)sitepack[0].size(),TIMESTEP,NMAX,DELTAN,RECURSION}; // set parameters

  vector< vector<double> > xpi = {{XIN}}; // set i.c. for inflationary trajectories

  StocDeltaN sdn(MODEL,sitepack,xpi,0,params); // declare the system
  
  sdn.sample(); // obtain 1 sample path
  //sdn.sample_logplot(); // plot obtained sample path
  
  //sdn.solve(); // solve PDE & SDE to obtain power spectrum
  //sdn.f_loglinearplot(0); // show plot of <N>
  //sdn.f_loglinearplot(1); // show plot of <delta N^2>
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
  return LAMBDA*LAMBDA*LAMBDA*LAMBDA*(1-X[0]*X[0]/MU/MU);
}

double StocDeltaN::VI(vector<double> &X, int I) // \partial_I V
{
  return -2*X[0]*LAMBDA*LAMBDA*LAMBDA*LAMBDA/MU/MU;
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
