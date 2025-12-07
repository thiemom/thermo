#ifndef COMPRESSIBLE_H
#define COMPRESSIBLE_H

#include "state.h"
#include <cstddef>

// -------------------------------------------------------------
// Compressible isentropic flow for ideal gas with variable cp(T)
// -------------------------------------------------------------
// Uses ideal gas equation of state (PV = nRT) with temperature-dependent
// thermodynamic properties (cp, h, s) from NASA polynomial fits.
// This is sometimes called "thermally perfect" gas modeling.
//
// Isentropic relations are solved numerically since gamma varies with T.

// -------------------------------------------------------------
// Compressible flow solution
// -------------------------------------------------------------

// Result of a compressible flow calculation
struct CompressibleFlowSolution {
    State stagnation;    // Stagnation state (T0, P0, h0, s0, ...)
    State outlet;        // Outlet static state (T, P, rho, h, s, a)
    double v = 0.0;      // Outlet velocity [m/s]
    double M = 0.0;      // Outlet Mach number [-]
    double mdot = 0.0;   // Mass flow rate [kg/s]
    bool choked = false; // True if flow is choked (M = 1 at throat)
};

// -------------------------------------------------------------
// Forward problem: given geometry and pressures, find mass flow
// -------------------------------------------------------------

// Isentropic nozzle flow: given stagnation conditions, back pressure, and
// effective flow area (A * Cd), compute mass flow and outlet state.
//
// T0      : stagnation temperature [K]
// P0      : stagnation pressure [Pa]
// P_back  : back pressure [Pa]
// A_eff   : effective flow area [m²] (= A * Cd)
// X       : mole fractions
// tol     : convergence tolerance (default: 1e-8)
// max_iter: maximum iterations (default: 50)
//
// Returns CompressibleFlowSolution with outlet state at throat/exit.
// If P_back <= P_critical, flow is choked and outlet is at sonic conditions.
CompressibleFlowSolution nozzle_flow(
    double T0, double P0, double P_back, double A_eff,
    const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

// -------------------------------------------------------------
// Inverse problems: solve for unknown given mass flow
// -------------------------------------------------------------

// Find effective area for given mass flow rate.
// Throws if mdot_target exceeds choked mass flow.
double solve_A_eff_from_mdot(
    double T0, double P0, double P_back, double mdot_target,
    const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

// Find back pressure for given mass flow rate.
// Returns the subsonic solution (P_back > P_critical).
// Throws if mdot_target exceeds choked mass flow.
double solve_P_back_from_mdot(
    double T0, double P0, double A_eff, double mdot_target,
    const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

// Find stagnation pressure for given mass flow rate.
// Throws if no solution exists.
double solve_P0_from_mdot(
    double T0, double P_back, double A_eff, double mdot_target,
    const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

// -------------------------------------------------------------
// Utility functions
// -------------------------------------------------------------

// Compute critical (sonic) pressure ratio P*/P0 for given stagnation state.
// Uses isentropic expansion to find where M = 1.
double critical_pressure_ratio(
    double T0, double P0, const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

// Compute Mach number from pressure ratio P/P0 (isentropic).
double mach_from_pressure_ratio(
    double T0, double P0, double P,
    const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

// Compute mass flux G = rho * v [kg/(m²·s)] at given pressure for isentropic
// expansion from stagnation conditions.
double mass_flux_isentropic(
    double T0, double P0, double P,
    const std::vector<double>& X,
    double tol = 1e-8, std::size_t max_iter = 50);

#endif // COMPRESSIBLE_H
