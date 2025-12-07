// Compressible isentropic flow for ideal gas with variable cp(T).
// Solves isentropic relations numerically using Newton iteration.

#include "compressible.h"
#include "thermo.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

using combaero::thermo::R_GAS;

namespace {

// Reference pressure for entropy calculations
constexpr double P_REF = 101325.0;

// Solve for temperature at pressure P such that s(T, P) = s0 (isentropic)
// Uses Newton iteration with ds/dT = cp/T
double solve_T_isentropic(double P, double s0, double T_init,
                          const std::vector<double>& X,
                          double tol, std::size_t max_iter) {
    double T = T_init;
    
    // Clamp temperature bounds to valid thermo data range
    constexpr double T_MIN = 200.0;
    constexpr double T_MAX = 6000.0;
    
    for (std::size_t it = 0; it < max_iter; ++it) {
        double s_curr = s(T, X, P, P_REF);
        double F = s_curr - s0;
        
        // ds/dT at constant P = cp/T
        double dF = cp(T, X) / T;
        
        if (std::abs(dF) < 1e-30) break;
        
        double dT = -F / dF;
        
        // Limit step size to avoid overshooting
        double max_step = 0.5 * T;
        if (std::abs(dT) > max_step) {
            dT = (dT > 0) ? max_step : -max_step;
        }
        
        T += dT;
        
        // Clamp temperature to reasonable bounds
        T = std::max(T_MIN, std::min(T, T_MAX));
        
        if (std::abs(dT) < tol * (1.0 + std::abs(T))) {
            return T;
        }
    }
    
    return T;
}

// Compute mass flux G = rho * v for isentropic expansion from (T0, P0) to P
// Returns G in kg/(m²·s)
double compute_mass_flux(double T0, double P0, double P, double s0, double h0,
                         const std::vector<double>& X,
                         double tol, std::size_t max_iter) {
    if (P >= P0) return 0.0;
    
    double T = solve_T_isentropic(P, s0, T0, X, tol, max_iter);
    
    double h_curr = h(T, X);
    double dh = h0 - h_curr;
    
    if (dh <= 0.0) return 0.0;
    
    double v = std::sqrt(2.0 * dh);
    double rho_curr = density(T, P, X);
    
    return rho_curr * v;
}

// Golden section search to find pressure that maximizes mass flux (critical point)
double find_critical_pressure(double T0, double P0, double s0, double h0,
                              const std::vector<double>& X,
                              double tol, std::size_t max_iter) {
    const double phi = 0.5 * (3.0 - std::sqrt(5.0));  // Golden ratio conjugate
    
    // For ideal gas, P*/P0 ~ 0.528 for gamma=1.4
    // Search in a reasonable range around this
    double a = 0.3 * P0;   // Lower bound (below typical critical ratio)
    double b = 0.99 * P0;  // Upper bound (just below stagnation)
    
    double x1 = b - phi * (b - a);
    double x2 = a + phi * (b - a);
    
    double f1 = compute_mass_flux(T0, P0, x1, s0, h0, X, tol, max_iter);
    double f2 = compute_mass_flux(T0, P0, x2, s0, h0, X, tol, max_iter);
    
    for (std::size_t it = 0; it < max_iter; ++it) {
        if (std::abs(b - a) < tol * (1.0 + std::abs(a) + std::abs(b))) break;
        
        if (f1 < f2) {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + phi * (b - a);
            f2 = compute_mass_flux(T0, P0, x2, s0, h0, X, tol, max_iter);
        } else {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - phi * (b - a);
            f1 = compute_mass_flux(T0, P0, x1, s0, h0, X, tol, max_iter);
        }
    }
    
    return (f1 > f2) ? x1 : x2;
}

}  // anonymous namespace

// -------------------------------------------------------------
// Forward problem: nozzle flow
// -------------------------------------------------------------

CompressibleFlowSolution nozzle_flow(
    double T0, double P0, double P_back, double A_eff,
    const std::vector<double>& X,
    double tol, std::size_t max_iter) {
    
    if (T0 <= 0.0) {
        throw std::invalid_argument("nozzle_flow: T0 must be positive");
    }
    if (P0 <= 0.0) {
        throw std::invalid_argument("nozzle_flow: P0 must be positive");
    }
    if (P_back <= 0.0) {
        throw std::invalid_argument("nozzle_flow: P_back must be positive");
    }
    if (A_eff <= 0.0) {
        throw std::invalid_argument("nozzle_flow: A_eff must be positive");
    }
    
    CompressibleFlowSolution sol;
    
    // Set up stagnation state
    sol.stagnation.T = T0;
    sol.stagnation.P = P0;
    sol.stagnation.X = X;
    
    double h0 = sol.stagnation.h();
    double s0 = sol.stagnation.s();
    
    // Find critical pressure (where mass flux is maximum)
    double P_crit = find_critical_pressure(T0, P0, s0, h0, X, tol, max_iter);
    double G_crit = compute_mass_flux(T0, P0, P_crit, s0, h0, X, tol, max_iter);
    
    double P_outlet;
    if (P_back <= P_crit) {
        // Choked flow: outlet is at critical conditions
        P_outlet = P_crit;
        sol.choked = true;
    } else {
        // Subsonic flow: outlet is at back pressure
        P_outlet = P_back;
        sol.choked = false;
    }
    
    // Compute outlet state
    double T_outlet = solve_T_isentropic(P_outlet, s0, T0, X, tol, max_iter);
    sol.outlet.T = T_outlet;
    sol.outlet.P = P_outlet;
    sol.outlet.X = X;
    
    double h_outlet = sol.outlet.h();
    double dh = h0 - h_outlet;
    
    sol.v = (dh > 0.0) ? std::sqrt(2.0 * dh) : 0.0;
    
    double a_outlet = sol.outlet.a();
    sol.M = (a_outlet > 0.0) ? sol.v / a_outlet : 0.0;
    
    double rho_outlet = sol.outlet.rho();
    double G_outlet = rho_outlet * sol.v;
    
    sol.mdot = A_eff * G_outlet;
    
    return sol;
}

// -------------------------------------------------------------
// Inverse problems
// -------------------------------------------------------------

double solve_A_eff_from_mdot(
    double T0, double P0, double P_back, double mdot_target,
    const std::vector<double>& X,
    double tol, std::size_t max_iter) {
    
    if (mdot_target <= 0.0) {
        throw std::invalid_argument("solve_A_eff_from_mdot: mdot_target must be positive");
    }
    
    // Compute mass flux at operating conditions
    State stag;
    stag.T = T0;
    stag.P = P0;
    stag.X = X;
    double h0 = stag.h();
    double s0 = stag.s();
    
    // Find critical conditions
    double P_crit = find_critical_pressure(T0, P0, s0, h0, X, tol, max_iter);
    double G_crit = compute_mass_flux(T0, P0, P_crit, s0, h0, X, tol, max_iter);
    
    double P_outlet = (P_back <= P_crit) ? P_crit : P_back;
    double G_outlet = compute_mass_flux(T0, P0, P_outlet, s0, h0, X, tol, max_iter);
    
    if (G_outlet <= 0.0) {
        throw std::runtime_error("solve_A_eff_from_mdot: zero mass flux at operating conditions");
    }
    
    return mdot_target / G_outlet;
}

double solve_P_back_from_mdot(
    double T0, double P0, double A_eff, double mdot_target,
    const std::vector<double>& X,
    double tol, std::size_t max_iter) {
    
    if (mdot_target <= 0.0) {
        throw std::invalid_argument("solve_P_back_from_mdot: mdot_target must be positive");
    }
    if (A_eff <= 0.0) {
        throw std::invalid_argument("solve_P_back_from_mdot: A_eff must be positive");
    }
    
    State stag;
    stag.T = T0;
    stag.P = P0;
    stag.X = X;
    double h0 = stag.h();
    double s0 = stag.s();
    
    // Find critical conditions
    double P_crit = find_critical_pressure(T0, P0, s0, h0, X, tol, max_iter);
    double G_crit = compute_mass_flux(T0, P0, P_crit, s0, h0, X, tol, max_iter);
    double mdot_choked = A_eff * G_crit;
    
    if (mdot_target > mdot_choked * (1.0 + tol)) {
        throw std::runtime_error("solve_P_back_from_mdot: mdot_target exceeds choked mass flow");
    }
    
    // If essentially choked, return critical pressure
    if (mdot_target >= mdot_choked * (1.0 - tol)) {
        return P_crit;
    }
    
    // Bisection search for P_back in [P_crit, P0]
    double G_target = mdot_target / A_eff;
    double P_lo = P_crit;
    double P_hi = P0;
    
    for (std::size_t it = 0; it < max_iter; ++it) {
        double P_mid = 0.5 * (P_lo + P_hi);
        double G_mid = compute_mass_flux(T0, P0, P_mid, s0, h0, X, tol, max_iter);
        
        if (std::abs(G_mid - G_target) < tol * G_target) {
            return P_mid;
        }
        
        // Mass flux decreases as P increases toward P0
        if (G_mid > G_target) {
            P_lo = P_mid;
        } else {
            P_hi = P_mid;
        }
        
        if ((P_hi - P_lo) < tol * P_lo) {
            return P_mid;
        }
    }
    
    return 0.5 * (P_lo + P_hi);
}

double solve_P0_from_mdot(
    double T0, double P_back, double A_eff, double mdot_target,
    const std::vector<double>& X,
    double tol, std::size_t max_iter) {
    
    if (mdot_target <= 0.0) {
        throw std::invalid_argument("solve_P0_from_mdot: mdot_target must be positive");
    }
    if (A_eff <= 0.0) {
        throw std::invalid_argument("solve_P0_from_mdot: A_eff must be positive");
    }
    if (P_back <= 0.0) {
        throw std::invalid_argument("solve_P0_from_mdot: P_back must be positive");
    }
    
    // Bisection search for P0
    // P0 must be > P_back for any flow
    double P0_lo = P_back * 1.001;
    double P0_hi = P_back * 1000.0;  // Start with large upper bound
    
    // Expand upper bound if needed
    for (int expand = 0; expand < 10; ++expand) {
        auto sol = nozzle_flow(T0, P0_hi, P_back, A_eff, X, tol, max_iter);
        if (sol.mdot >= mdot_target) break;
        P0_hi *= 10.0;
    }
    
    // Bisection
    for (std::size_t it = 0; it < max_iter; ++it) {
        double P0_mid = 0.5 * (P0_lo + P0_hi);
        auto sol = nozzle_flow(T0, P0_mid, P_back, A_eff, X, tol, max_iter);
        
        if (std::abs(sol.mdot - mdot_target) < tol * mdot_target) {
            return P0_mid;
        }
        
        if (sol.mdot < mdot_target) {
            P0_lo = P0_mid;
        } else {
            P0_hi = P0_mid;
        }
        
        if ((P0_hi - P0_lo) < tol * P0_lo) {
            return P0_mid;
        }
    }
    
    return 0.5 * (P0_lo + P0_hi);
}

// -------------------------------------------------------------
// Utility functions
// -------------------------------------------------------------

double critical_pressure_ratio(
    double T0, double P0, const std::vector<double>& X,
    double tol, std::size_t max_iter) {
    
    State stag;
    stag.T = T0;
    stag.P = P0;
    stag.X = X;
    double h0 = stag.h();
    double s0 = stag.s();
    
    double P_crit = find_critical_pressure(T0, P0, s0, h0, X, tol, max_iter);
    return P_crit / P0;
}

double mach_from_pressure_ratio(
    double T0, double P0, double P,
    const std::vector<double>& X,
    double tol, std::size_t max_iter) {
    
    if (P >= P0) return 0.0;
    
    State stag;
    stag.T = T0;
    stag.P = P0;
    stag.X = X;
    double h0 = stag.h();
    double s0 = stag.s();
    
    double T = solve_T_isentropic(P, s0, T0, X, tol, max_iter);
    
    double h_curr = h(T, X);
    double dh = h0 - h_curr;
    
    if (dh <= 0.0) return 0.0;
    
    double v = std::sqrt(2.0 * dh);
    double a_curr = speed_of_sound(T, X);
    
    return (a_curr > 0.0) ? v / a_curr : 0.0;
}

double mass_flux_isentropic(
    double T0, double P0, double P,
    const std::vector<double>& X,
    double tol, std::size_t max_iter) {
    
    State stag;
    stag.T = T0;
    stag.P = P0;
    stag.X = X;
    double h0 = stag.h();
    double s0 = stag.s();
    
    return compute_mass_flux(T0, P0, P, s0, h0, X, tol, max_iter);
}
