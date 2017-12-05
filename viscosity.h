#ifndef VISCOSITYHEADER
#define VISCOSITYHEADER

// Function for Sutherland's law
template <typename U>
U mu(U T) {
  U T0,S,mu0;
  T0 = (U) 273.1;
  S = (U) 110.6;
  mu0 = (U) 1.716e-5;
  return (U) mu0*(pow(T/T0,1.5)*(T0+S)/(T+S));
}

#endif
