#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <initializer_list>
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
// Main
//-------------------------------------------------------------------

int main() {

  // Calling Golden Section Search code to determine the initial condition which satisfies the equations
  double fpp0 = gss(1.0,0.01,2.0,1e-8,&solve_ode);

  // Writing state
  std::ofstream outfile("results.dat");
  outfile.setf(std::ios_base::scientific);
  outfile.precision(10);
  outfile << "# eta f(eta) f'(eta) f''(eta) \n";
  for (std::array<double,4> current_state : state_hist) {
    for (double value : current_state) {
      outfile << value << " ";
    }
    outfile << "\n";
  }
  outfile.close();

  return 0;

}
