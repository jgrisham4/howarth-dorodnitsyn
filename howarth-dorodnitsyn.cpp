#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <initializer_list>
#include <cmath>
#include "gss.h"

// Time step
#define DT 0.1

// Global variables
std::vector<std::array<double,4> > state_hist;

//-------------------------------------------------------------------
// First-order form of the equation
// f''' + ff'' = 0
//-------------------------------------------------------------------

void deriv(const std::vector<double>& f, std::vector<double>& dfdeta, double eta) {

  dfdeta[0] = f[1];
  dfdeta[1] = f[2];
  dfdeta[2] = -f[0]*f[2];

}

//-------------------------------------------------------------------
// Function for storing ODE solution
//-------------------------------------------------------------------

void output(const std::vector<double>& f, const double eta) {
  state_hist.push_back({eta,f[0],f[1],f[2]});
}

//-------------------------------------------------------------------
// Function for solving the ODE using boost odeint
//-------------------------------------------------------------------

double solve_ode(const double fpp0) {

  // Initial conditions
  if (state_hist.size()!=0) {
    state_hist.clear();
  }
  std::vector<double> state({0.0,0.0,fpp0});

  // Setting up for RK4
  boost::numeric::odeint::runge_kutta4<std::vector<double> > stepper;

  // Integrating ODE
  size_t steps = boost::numeric::odeint::integrate_const(stepper,deriv,state,0.0,4.5,DT,output);

  // Returning objective function
  return pow(state[1]-1.0,2);

}

//-------------------------------------------------------------------
// Function for integrating using the trapezoidal rule
//-------------------------------------------------------------------

double integrate(const int i, const double rhoe, const double nue, const double Ue, const double x, const std::vector<double>& eta, const std::vector<double>& rho) {
  
  double sum=0.0;
  double Yim1,Yi;

  for (int j=0; j<i; ++j) {
    Yim1 = eta[j]*sqrt(2.0*nue*x/Ue);
    Yi   = eta[j+1]*sqrt(2.0*nue*x/Ue);
    sum += 0.5*(rhoe/rho[j+1] + rhoe/rho[j])*fabs(Yi - Yim1);
  }

  return sum;

}


//-------------------------------------------------------------------
// Main
//-------------------------------------------------------------------

int main() {

  // Inputs
  int size;
  double Me = 2.64;
  double gamma = 1.4;
  double R = 286.9;       // [J/(kg K)]
  double Te = 273.0;      // [K]
  double pe = 107.574;    // [Pa]
  double x = 1.0;         // [m] -- streamwise location at which information will be determined
  double nue = 13.30e-6;  // [m^2/s]
  double Tw, Ue, cp, rhoe;
  std::vector<double> T, rho, y, eta_vec;

  // Finding freestream properties
  Ue = Me*sqrt(gamma*R*Te);
  rhoe = pe/(R*Te);
  std::cout << "Ue = " << Ue << " m/s" << std::endl;

  // Finding adiabatic wall temp
  cp = gamma*R/(gamma - 1.0);
  Tw = Te + Ue*Ue/(2.0*cp);

  // Calling Golden Section Search code to determine the initial condition which satisfies the equations
  double fpp0 = gss(1.0,0.01,2.0,1e-8,&solve_ode);

  // Now have Blasius solution in x,Y.  Need to transform back to x,y by integrating.
  // First must determine temperature profile
  size = state_hist.size();
  T.resize(size);
  rho.resize(size);
  eta_vec.resize(size);
  double fp;
  for (int i=0; i<size; ++i) {
    eta_vec[i] = state_hist[i][0];
    fp = state_hist[i][2];
    T[i] = Te*(1.0 + (gamma-1.0)/2.0*Me*Me*(1.0 - fp*fp) + (Tw/Te - 1.0 - (gamma-1.0)/2.0*Me*Me)*(1-fp));
    rho[i] = pe/(R*T[i]);
  }

  // Integrating
  y.resize(size);
  for (int i=0; i<size; ++i) {
    y[i] = integrate(i,rhoe,nue,Ue,x,eta_vec,rho);
  }

  // Writing data
  std::ofstream outfile("results.dat");
  outfile.setf(std::ios_base::scientific);
  outfile.precision(10);
  outfile << "# y, y/2 sqrt(Ue/(nue x)), u/Ue, T [K]\n";
  for (int i=0; i<size; ++i) {
    outfile << y[i] << " " << 0.5*y[i]*sqrt(Ue/(nue*x)) << " " << state_hist[i][2] << " " << T[i] << "\n";
  }
  outfile.close();

  return 0;

}
