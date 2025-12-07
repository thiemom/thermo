#ifndef TRANSPORT_H
#define TRANSPORT_H

#include <vector>

// Collision integral functions
double linear_interp(double x, const std::vector<double>& x_values, const std::vector<double>& y_values);
double omega22(double T, double well_depth);

// Transport properties
double viscosity(double T, double P, const std::vector<double>& X);
double thermal_conductivity(double T, double P, const std::vector<double>& X);
double prandtl(double T, double P, const std::vector<double>& X);
double kinematic_viscosity(double T, double P, const std::vector<double>& X);
double thermal_diffusivity(double T, double P, const std::vector<double>& X);  // α = k/(ρ·cp) [m²/s]
double reynolds(double T, double P, const std::vector<double>& X, double V, double L);  // Re = ρVL/μ [-]
double peclet(double T, double P, const std::vector<double>& X, double V, double L);    // Pe = VL/α [-]

#endif // TRANSPORT_H
