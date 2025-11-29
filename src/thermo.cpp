#include "../include/thermo_transport.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iostream>

// Core thermodynamic utilities factored out of thermo_transport.cpp.
// Public API remains declared in thermo_transport.h / thermo.h.

// Conversion factor from J/mol to J/kg
double J_per_mol_to_J_per_kg(double value, double molar_mass) {
    return value * 1000.0 / molar_mass; // molar_mass is in g/mol, convert to kg/mol
}

// Species name lookup functions
std::string species_name(std::size_t species_index) {
    if (species_index >= species_names.size()) {
        throw std::out_of_range("Species index out of range");
    }
    return species_names[species_index];
}

std::size_t species_index_from_name(const std::string& name) {
    auto it = species_index.find(name);
    if (it == species_index.end()) {
        throw std::out_of_range("Species name not found: " + name);
    }
    return static_cast<std::size_t>(it->second);
}

// Calculate mixture molecular weight [g/mol]
double mwmix(const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }

    double sum = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        sum += X[i] * molar_masses[i];
    }
    return sum;
}

// NASA polynomial evaluation for Cp/R
double cp_R(std::size_t species_idx, double T) {
    const NASA_Coeffs& nasa = nasa_coeffs[species_idx];

    const std::vector<double>* coeffs;
    if (T < nasa.T_mid) {
        if (T < nasa.T_low) {
            std::cerr << "Warning: Temperature " << T << " K is below valid range (" << nasa.T_low
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the lowest valid temperature coefficients
            T = nasa.T_low;
        }
        coeffs = &nasa.low_coeffs;
    } else {
        if (T > nasa.T_high) {
            std::cerr << "Warning: Temperature " << T << " K is above valid range (" << nasa.T_high
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the highest valid temperature coefficients
            T = nasa.T_high;
        }
        coeffs = &nasa.high_coeffs;
    }

    // NASA polynomial for Cp/R: a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^4
    return (*coeffs)[0] + (*coeffs)[1] * T + (*coeffs)[2] * T * T +
           (*coeffs)[3] * T * T * T + (*coeffs)[4] * T * T * T * T;
}

// NASA polynomial evaluation for H/RT
double h_RT(std::size_t species_idx, double T) {
    const NASA_Coeffs& nasa = nasa_coeffs[species_idx];

    const std::vector<double>* coeffs;
    if (T < nasa.T_mid) {
        if (T < nasa.T_low) {
            std::cerr << "Warning: Temperature " << T << " K is below valid range (" << nasa.T_low
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the lowest valid temperature coefficients
            T = nasa.T_low;
        }
        coeffs = &nasa.low_coeffs;
    } else {
        if (T > nasa.T_high) {
            std::cerr << "Warning: Temperature " << T << " K is above valid range (" << nasa.T_high
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the highest valid temperature coefficients
            T = nasa.T_high;
        }
        coeffs = &nasa.high_coeffs;
    }

    // NASA polynomial for H/RT: a1 + a2*T/2 + a3*T^2/3 + a4*T^3/4 + a5*T^4/5 + a6/T
    return (*coeffs)[0] + (*coeffs)[1] * T / 2.0 + (*coeffs)[2] * T * T / 3.0 +
           (*coeffs)[3] * T * T * T / 4.0 + (*coeffs)[4] * T * T * T * T / 5.0 +
           (*coeffs)[5] / T;
}

// NASA polynomial evaluation for S/R
double s_R(std::size_t species_idx, double T) {
    const NASA_Coeffs& nasa = nasa_coeffs[species_idx];

    const std::vector<double>* coeffs;
    if (T < nasa.T_mid) {
        if (T < nasa.T_low) {
            std::cerr << "Warning: Temperature " << T << " K is below valid range (" << nasa.T_low
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the lowest valid temperature coefficients
            T = nasa.T_low;
        }
        coeffs = &nasa.low_coeffs;
    } else {
        if (T > nasa.T_high) {
            std::cerr << "Warning: Temperature " << T << " K is above valid range (" << nasa.T_high
                      << " K) for species " << species_names[species_idx] << std::endl;
            std::cerr << "Results may be less accurate but represent the best available estimate." << std::endl;
            // Use the highest valid temperature coefficients
            T = nasa.T_high;
        }
        coeffs = &nasa.high_coeffs;
    }

    // NASA polynomial for S/R: a1*ln(T) + a2*T + a3*T^2/2 + a4*T^3/3 + a5*T^4/4 + a7
    return (*coeffs)[0] * std::log(T) + (*coeffs)[1] * T + (*coeffs)[2] * T * T / 2.0 +
           (*coeffs)[3] * T * T * T / 3.0 + (*coeffs)[4] * T * T * T * T / 4.0 +
           (*coeffs)[6];
}

// Dimensionless Gibbs free energy G/(R*T) = H/(R*T) - S/R
double g_over_RT(std::size_t species_idx, double T) {
    return h_RT(species_idx, T) - s_R(species_idx, T);
}

// Number of thermo species used in the internal tables
std::size_t num_species() {
    return species_names.size();
}

// Molar mass lookup [g/mol]
double species_molar_mass(std::size_t species_index) {
    if (species_index >= molar_masses.size()) {
        throw std::out_of_range("species_molar_mass: species_index out of range");
    }
    return molar_masses[species_index];
}

double species_molar_mass_from_name(const std::string& name) {
    const std::size_t idx = species_index_from_name(name);
    return species_molar_mass(idx);
}

// Calculate heat capacity at constant pressure [J/(mol·K)]
double cp(double T, const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }

    double cp_mix = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        cp_mix += X[i] * cp_R(i, T) * R_GAS;
    }
    return cp_mix;
}

// Calculate enthalpy [J/mol]
double h(double T, const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }

    double h_mix = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        h_mix += X[i] * h_RT(i, T) * R_GAS * T;
    }
    return h_mix;
}

// Calculate entropy [J/(mol·K)]
double s(double T, const std::vector<double>& X, double P, double P_ref) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }

    double s_mix = 0.0;

    // Calculate pure species entropies and mixing term
    for (std::size_t i = 0; i < X.size(); ++i) {
        if (X[i] > 0.0) {
            // Pure species entropy at reference pressure
            double s_i = s_R(i, T) * R_GAS;

            // Pressure correction: -R * ln(P/P_ref)
            s_i -= R_GAS * std::log(P / P_ref);

            // Mixing entropy: -R * ln(X_i)
            s_i -= R_GAS * std::log(X[i]);

            s_mix += X[i] * s_i;
        }
    }

    return s_mix;
}

// Calculate heat capacity at constant volume [J/(mol·K)]
double cv(double T, const std::vector<double>& X) {
    // Cv = Cp - R for ideal gases
    return cp(T, X) - R_GAS;
}

// Calculate derivative of enthalpy with respect to temperature [J/(mol·K)]
// For ideal gas, this is the same as Cp
double dh_dT(double T, const std::vector<double>& X) {
    return cp(T, X);
}

// Calculate derivative of entropy with respect to temperature [J/(mol·K^2)]
double ds_dT(double T, const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }

    double ds_dT_mix = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        if (X[i] > 0.0) {
            // Derivative of NASA polynomial for S/R with respect to T
            const NASA_Coeffs& nasa = nasa_coeffs[i];

            const std::vector<double>* coeffs;
            if (T < nasa.T_mid) {
                if (T < nasa.T_low) {
                    throw std::out_of_range("Temperature below valid range for species " + species_names[i]);
                }
                coeffs = &nasa.low_coeffs;
            } else {
                if (T > nasa.T_high) {
                    throw std::out_of_range("Temperature above valid range for species " + species_names[i]);
                }
                coeffs = &nasa.high_coeffs;
            }

            // d(S/R)/dT = a1/T + a2 + a3*T + a4*T^2/2 + a5*T^3/3
            double ds_dT_i = (*coeffs)[0] / T + (*coeffs)[1] + (*coeffs)[2] * T +
                             (*coeffs)[3] * T * T / 2.0 + (*coeffs)[4] * T * T * T / 3.0;

            ds_dT_mix += X[i] * ds_dT_i * R_GAS;
        }
    }

    return ds_dT_mix;
}

// Calculate derivative of heat capacity with respect to temperature [J/(mol·K^2)]
double dcp_dT(double T, const std::vector<double>& X) {
    if (X.size() != species_names.size()) {
        throw std::invalid_argument("Mole fraction vector size does not match number of species");
    }

    double dcp_dT_mix = 0.0;
    for (std::size_t i = 0; i < X.size(); ++i) {
        const NASA_Coeffs& nasa = nasa_coeffs[i];

        const std::vector<double>* coeffs;
        if (T < nasa.T_mid) {
            if (T < nasa.T_low) {
                throw std::out_of_range("Temperature below valid range for species " + species_names[i]);
            }
            coeffs = &nasa.low_coeffs;
        } else {
            if (T > nasa.T_high) {
                throw std::out_of_range("Temperature above valid range for species " + species_names[i]);
            }
            coeffs = &nasa.high_coeffs;
        }

        // d(Cp/R)/dT = a2 + 2*a3*T + 3*a4*T^2 + 4*a5*T^3
        double dcp_dT_i = (*coeffs)[1] + 2.0 * (*coeffs)[2] * T +
                          3.0 * (*coeffs)[3] * T * T + 4.0 * (*coeffs)[4] * T * T * T;

        dcp_dT_mix += X[i] * dcp_dT_i * R_GAS;
    }

    return dcp_dT_mix;
}
