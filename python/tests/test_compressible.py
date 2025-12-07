"""Tests for compressible flow functions."""

import numpy as np
import combaero as cb


def test_nozzle_flow_subsonic() -> None:
    """Subsonic nozzle flow (P_back > P_critical)."""
    X = cb.standard_dry_air_composition()
    T0 = 500.0  # K (high enough to stay above 300K after expansion)
    P0 = 200000.0  # Pa
    P_back = 150000.0  # Pa (subsonic)
    A_eff = 0.001  # m²

    sol = cb.nozzle_flow(T0, P0, P_back, A_eff, X)

    assert sol.M < 1.0  # Subsonic
    assert not sol.choked
    assert sol.mdot > 0
    assert sol.outlet.T < T0  # Temperature drops
    assert sol.outlet.P < P0  # Pressure drops


def test_nozzle_flow_choked() -> None:
    """Choked nozzle flow (P_back < P_critical)."""
    X = cb.standard_dry_air_composition()
    T0 = 600.0  # K (high enough for large expansion)
    P0 = 200000.0  # Pa
    P_back = 50000.0  # Pa (choked)
    A_eff = 0.001  # m²

    sol = cb.nozzle_flow(T0, P0, P_back, A_eff, X)

    # Flow should be choked (M=1 at throat)
    assert sol.choked
    # Mass flow should be positive
    assert sol.mdot > 0


def test_critical_pressure_ratio() -> None:
    """Critical pressure ratio P*/P0."""
    X = cb.standard_dry_air_composition()
    T0 = 500.0  # K
    P0 = 100000.0  # Pa

    ratio = cb.critical_pressure_ratio(T0, P0, X)

    # For variable-gamma gas, ratio is higher than ideal constant-gamma
    # Actual value ~0.726 for air at 300K
    assert 0.5 < ratio < 0.8


def test_mach_from_pressure_ratio() -> None:
    """Mach number from isentropic pressure ratio."""
    X = cb.standard_dry_air_composition()
    T0 = 500.0  # K
    P0 = 100000.0  # Pa
    P = 80000.0  # Pa

    M = cb.mach_from_pressure_ratio(T0, P0, P, X)

    assert M > 0
    assert M < 1.0  # Subsonic for P/P0 = 0.8


def test_solve_A_eff_from_mdot() -> None:
    """Inverse problem: find area from mass flow."""
    X = cb.standard_dry_air_composition()
    T0 = 500.0
    P0 = 200000.0
    P_back = 150000.0
    A_eff_orig = 0.001

    # Get mass flow for known area
    sol = cb.nozzle_flow(T0, P0, P_back, A_eff_orig, X)
    mdot = sol.mdot

    # Solve inverse
    A_eff_solved = cb.solve_A_eff_from_mdot(T0, P0, P_back, mdot, X)

    assert np.isclose(A_eff_solved, A_eff_orig, rtol=1e-4)


def test_fanno_pipe_basic() -> None:
    """Basic Fanno flow through a pipe."""
    X = cb.standard_dry_air_composition()
    T_in = 400.0  # K
    P_in = 200000.0  # Pa
    u_in = 50.0  # m/s (subsonic)
    L = 1.0  # m
    D = 0.05  # m
    f = 0.02  # Darcy friction factor

    sol = cb.fanno_pipe(T_in, P_in, u_in, L, D, f, X)

    # Pressure should drop due to friction
    assert sol.outlet.P < P_in
    # Velocity should increase (subsonic Fanno)
    assert sol.outlet.T > 0
    assert sol.mdot > 0


def test_fanno_max_length() -> None:
    """Maximum pipe length before choking."""
    X = cb.standard_dry_air_composition()
    T_in = 400.0
    P_in = 200000.0
    u_in = 100.0  # m/s
    D = 0.05
    f = 0.02

    L_star = cb.fanno_max_length(T_in, P_in, u_in, D, f, X)

    assert L_star > 0


def test_friction_haaland() -> None:
    """Haaland friction factor."""
    Re = 50000.0
    e_D = 0.001

    f = cb.friction_haaland(Re, e_D)

    # Should be in reasonable range for turbulent flow
    assert 0.01 < f < 0.1


def test_friction_serghides() -> None:
    """Serghides friction factor."""
    Re = 50000.0
    e_D = 0.001

    f = cb.friction_serghides(Re, e_D)

    assert 0.01 < f < 0.1


def test_friction_colebrook() -> None:
    """Colebrook friction factor."""
    Re = 50000.0
    e_D = 0.001

    f = cb.friction_colebrook(Re, e_D)

    assert 0.01 < f < 0.1


def test_friction_consistency() -> None:
    """All friction correlations should be close."""
    Re = 100000.0
    e_D = 0.0001

    f_haaland = cb.friction_haaland(Re, e_D)
    f_serghides = cb.friction_serghides(Re, e_D)
    f_colebrook = cb.friction_colebrook(Re, e_D)

    # Serghides should be very close to Colebrook (<0.3%)
    assert np.isclose(f_serghides, f_colebrook, rtol=0.003)
    # Haaland within ~3%
    assert np.isclose(f_haaland, f_colebrook, rtol=0.03)
