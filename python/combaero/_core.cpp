#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>

#include "thermo_transport.h"
#include "equilibrium.h"

namespace py = pybind11;

// Build WgsConfig using species indices looked up by name
static combaero::equilibrium::WgsConfig make_wgs_cfg()
{
    combaero::equilibrium::WgsConfig cfg{};
    cfg.i_CO  = species_index_from_name("CO");
    cfg.i_H2O = species_index_from_name("H2O");
    cfg.i_CO2 = species_index_from_name("CO2");
    cfg.i_H2  = species_index_from_name("H2");
    return cfg;
}

static std::vector<double> to_vec(
    py::array_t<double, py::array::c_style | py::array::forcecast> arr)
{
    auto buf = arr.unchecked<1>();
    std::vector<double> v(buf.shape(0));
    for (ssize_t i = 0; i < buf.shape(0); ++i)
        v[static_cast<std::size_t>(i)] = buf(i);
    return v;
}

PYBIND11_MODULE(_core, m)
{
    m.doc() = "Python bindings for combaero core";

    // Species metadata helpers
    m.def("num_species", &num_species, "Number of thermo species in the internal tables.");

    m.def(
        "species_name",
        &species_name,
        py::arg("index"),
        "Return the canonical species name for a given index."
    );

    m.def(
        "species_index_from_name",
        &species_index_from_name,
        py::arg("name"),
        "Return the species index for a given canonical name."
    );

    m.def(
        "species_molar_mass",
        &species_molar_mass,
        py::arg("index"),
        "Return the molar mass [g/mol] for a given species index."
    );

    m.def(
        "species_molar_mass_from_name",
        &species_molar_mass_from_name,
        py::arg("name"),
        "Return the molar mass [g/mol] for a given species name."
    );

    m.def(
        "mixture_h",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return h(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mixture enthalpy at temperature T and mole fractions X."
    );

    // Thermodynamic mixture properties
    m.def(
        "cp",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return cp(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mixture isobaric heat capacity Cp(T, X) [J/(mol·K)]."
    );

    m.def(
        "h",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return h(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mixture enthalpy h(T, X) [J/mol]."
    );

    m.def(
        "s",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr,
           double P,
           double P_ref)
        {
            auto X = to_vec(X_arr);
            return s(T, X, P, P_ref);
        },
        py::arg("T"),
        py::arg("X"),
        py::arg("P"),
        py::arg("P_ref") = 101325.0,
        "Mixture entropy s(T, X, P, P_ref) [J/(mol·K)]."
    );

    m.def(
        "cv",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return cv(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Mixture isochoric heat capacity Cv(T, X) [J/(mol·K)]."
    );

    m.def(
        "density",
        [](double T,
           double P,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return density(T, P, X);
        },
        py::arg("T"),
        py::arg("P"),
        py::arg("X"),
        "Mixture density rho(T, P, X) [kg/m^3]."
    );

    m.def(
        "specific_gas_constant",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return specific_gas_constant(X);
        },
        py::arg("X"),
        "Mixture specific gas constant R_s(X) [J/(kg·K)]."
    );

    m.def(
        "isentropic_expansion_coefficient",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return isentropic_expansion_coefficient(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Isentropic expansion coefficient gamma(T, X) = Cp/Cv."
    );

    m.def(
        "speed_of_sound",
        [](double T,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return speed_of_sound(T, X);
        },
        py::arg("T"),
        py::arg("X"),
        "Speed of sound c(T, X) [m/s]."
    );

    // Composition conversion helpers
    m.def(
        "mole_to_mass",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_arr)
        {
            auto X = to_vec(X_arr);
            return mole_to_mass(X);
        },
        py::arg("X"),
        "Convert mole fractions X to mass fractions Y.")
    ;

    m.def(
        "mass_to_mole",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> Y_arr)
        {
            auto Y = to_vec(Y_arr);
            return mass_to_mole(Y);
        },
        py::arg("Y"),
        "Convert mass fractions Y to mole fractions X.")
    ;

    // Equivalence ratio (mole basis)
    m.def(
        "equivalence_ratio_mole",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> X_mix_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_ox_arr)
        {
            auto X_mix  = to_vec(X_mix_arr);
            auto X_fuel = to_vec(X_fuel_arr);
            auto X_ox   = to_vec(X_ox_arr);
            return equivalence_ratio_mole(X_mix, X_fuel, X_ox);
        },
        py::arg("X_mix"),
        py::arg("X_fuel"),
        py::arg("X_ox"),
        "Equivalence ratio phi (mole basis) for mixture X_mix formed from fuel X_fuel and oxidizer X_ox.")
    ;

    m.def(
        "set_equivalence_ratio_mole",
        [](double phi,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_ox_arr)
        {
            auto X_fuel = to_vec(X_fuel_arr);
            auto X_ox   = to_vec(X_ox_arr);
            return set_equivalence_ratio_mole(phi, X_fuel, X_ox);
        },
        py::arg("phi"),
        py::arg("X_fuel"),
        py::arg("X_ox"),
        "Construct mole-fraction mixture X_mix with target phi from fuel X_fuel and oxidizer X_ox.")
    ;

    // Equivalence ratio (mass basis)
    m.def(
        "equivalence_ratio_mass",
        [](py::array_t<double, py::array::c_style | py::array::forcecast> Y_mix_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_ox_arr)
        {
            auto Y_mix  = to_vec(Y_mix_arr);
            auto Y_fuel = to_vec(Y_fuel_arr);
            auto Y_ox   = to_vec(Y_ox_arr);
            return equivalence_ratio_mass(Y_mix, Y_fuel, Y_ox);
        },
        py::arg("Y_mix"),
        py::arg("Y_fuel"),
        py::arg("Y_ox"),
        "Equivalence ratio phi (mass basis) for mixture Y_mix formed from fuel Y_fuel and oxidizer Y_ox.")
    ;

    m.def(
        "set_equivalence_ratio_mass",
        [](double phi,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_fuel_arr,
           py::array_t<double, py::array::c_style | py::array::forcecast> Y_ox_arr)
        {
            auto Y_fuel = to_vec(Y_fuel_arr);
            auto Y_ox   = to_vec(Y_ox_arr);
            return set_equivalence_ratio_mass(phi, Y_fuel, Y_ox);
        },
        py::arg("phi"),
        py::arg("Y_fuel"),
        py::arg("Y_ox"),
        "Construct mass-fraction mixture Y_mix with target phi from fuel Y_fuel and oxidizer Y_ox.")
    ;

    m.def(
        "adiabatic_T_wgs",
        [](double T_in,
           py::array_t<double, py::array::c_style | py::array::forcecast> X_in_arr,
           py::object T_guess_obj)
        {
            auto X_in = to_vec(X_in_arr);
            if (X_in.empty())
                throw std::runtime_error("X_in must be non-empty");

            double T_guess = T_in + 500.0;
            if (!T_guess_obj.is_none())
                T_guess = T_guess_obj.cast<double>();

            // 1 mol basis
            std::vector<double> n_in = X_in;

            double H_target = h(T_in, X_in);
            auto cfg = make_wgs_cfg();

            double T_eq = combaero::equilibrium::solve_adiabatic_T_wgs(
                n_in, H_target, T_guess, cfg
            );
            return T_eq;
        },
        py::arg("T_in"),
        py::arg("X_in"),
        py::arg("T_guess") = py::none(),
        "Adiabatic flame temperature with WGS limited equilibrium.\n\n"
        "T_in   : inlet temperature [K]\n"
        "X_in   : inlet mole fractions (1D NumPy array)\n"
        "T_guess: optional initial guess for T [K]"
    );
}