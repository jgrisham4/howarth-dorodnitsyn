#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <initializer_list>
#include <cmath>
#include "viscosity.h"
#include "nlopt.hpp"
#include "utils.h"

// Time step
#define DT 0.001

// Global variables
std::vector<std::array<double,6> > state_hist;
input_data<double> inp = read_input_file<double>("input");


//-------------------------------------------------------------------
// First-order form of the equation
// f''' + ff'' = 0
//-------------------------------------------------------------------

void deriv(const std::vector<double>& f, std::vector<double>& dfdeta, double eta) {

  dfdeta[0] = f[1];
  dfdeta[1] = f[2];
  dfdeta[2] = -f[0]*f[2];
  dfdeta[3] = f[4];
  dfdeta[4] = -f[0]*f[4] - (inp.gamma-1.0)*inp.Me*inp.Me*inp.Te*f[2]*f[2];

}

//-------------------------------------------------------------------
// Function for storing ODE solution
//-------------------------------------------------------------------

void output(const std::vector<double>& f, const double eta) {
  state_hist.push_back({eta,f[0],f[1],f[2],f[3],f[4]});
}

//-------------------------------------------------------------------
// Function for solving the ODE using boost odeint
//-------------------------------------------------------------------

double solve_ode(const std::vector<double>& X, std::vector<double>& gradX, void* f_data) {

  // Initial conditions
  if (state_hist.size()!=0) {
    state_hist.clear();
  }
  double fpp0 = X[0];
  double Tprime = X[1];
  double Tw = *(double *) f_data;
  std::vector<double> state({0.0,0.0,fpp0,Tw,Tprime});

  // Setting up for RK4
  boost::numeric::odeint::runge_kutta4<std::vector<double> > stepper;

  // Integrating ODE
  size_t steps = boost::numeric::odeint::integrate_const(stepper,deriv,state,0.0,4.5,DT,output);

  // Returning objective function
  return pow(state[1]-1.0,2) + pow(state[3]-inp.Te,2);

}

//-------------------------------------------------------------------
// Function for integrating using the trapezoidal rule
//-------------------------------------------------------------------

double integrate(const int i, const double rhoe, const double nue, const double Ue, const double x, const std::vector<double>& eta, const std::vector<double>& T) {

  double sum=0.0;
  double Yim1,Yi;

  for (int j=0; j<i; ++j) {
    Yim1 = eta[j]*sqrt(2.0*nue*x/Ue);
    Yi   = eta[j+1]*sqrt(2.0*nue*x/Ue);
    sum += 0.5*(T[j+1]/inp.Te + T[j]/inp.Te)*(Yi - Yim1);
  }

  return sum;

}

//-------------------------------------------------------------------
// Main
//-------------------------------------------------------------------

int main() {

  // Inputs
  int size;
  double R = 287.0;       // [J/(kg K)]
  double x = 2.0;         // [m] -- streamwise location at which information will be determined
  double nue;
  double Ue, cp, rhoe, pe;
  std::vector<double> T, rho, y, eta_vec;

  // Finding freestream properties
  Ue = inp.Me*sqrt(inp.gamma*R*inp.Te);

  // Computing the state using the Reynolds number
  rhoe = mu<double>(inp.Te)*inp.ReL/(Ue*inp.L);
  nue = mu<double>(inp.Te)/rhoe;
  pe = rhoe*R*inp.Te;
  std::cout << "pe = " << pe << " Pa" << std::endl;

  // Finding adiabatic wall temp
  cp = inp.gamma*R/(inp.gamma - 1.0);
  double* Tw;
  double tmp = inp.Te + Ue*Ue/(2.0*cp);
  Tw = &tmp;

  // Outputting some info
  std::cout << "Te = " <<  inp.Te << " K\n";
  std::cout << "Tw = " << *Tw << " K\n";
  std::cout << "Ue = " <<  Ue << " m/s\n";

  // Calling nlopt code to determine the initial condition which satisfies the equations
  double Fopt;
  std::vector<double> Xopt({1.0,10.0});
  nlopt::opt opt(nlopt::LN_NELDERMEAD,2);
  std::vector<double> lb({0.01,0.01});
  std::vector<double> ub({5.0,1000.0});
  opt.set_min_objective(solve_ode,Tw);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  opt.optimize(Xopt,Fopt);

  // Now have Blasius solution in x,Y.  Need to transform back to x,y by integrating.
  // First must determine temperature profile
  size = state_hist.size();
  T.resize(size);
  eta_vec.resize(size);
  for (int i=0; i<size; ++i) {
    eta_vec[i] = state_hist[i][0];
    T[i] = state_hist[i][4];
  }

  // Integrating
  y.resize(size);
  for (int i=0; i<size; ++i) {
    y[i] = integrate(i,rhoe,nue,Ue,x,eta_vec,T);
  }

  // Writing data
  std::ofstream outfile("results.dat");
  outfile.setf(std::ios_base::scientific);
  outfile.precision(10);
  outfile << "# y, y sqrt(Ue/(nue x)), u/Ue, T/Te, eta_incomp\n";
  for (int i=0; i<size; ++i) {
    //outfile << y[i] << " " << y[i]*sqrt(Ue/(2.0*nue*x)) << " "
    outfile << y[i] << " " << y[i]*sqrt(Ue/(nue*x)) << " "
      << state_hist[i][2] << " " << T[i]/inp.Te << " " << eta_vec[i] << "\n";
  }
  outfile.close();

  return 0;

}
