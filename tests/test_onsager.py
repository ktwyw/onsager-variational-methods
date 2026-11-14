"""Tests of the onsager module: each function against a result known in closed form."""
import math

import numpy as np
import pytest
import sympy as sp

import onsager as on


def test_single_particle():
    x, v, zeta, k = sp.symbols("x v zeta k", positive=True)
    eq = on.kinetic_equations(k * x**2 / 2, zeta * v**2 / 2, [x], [v])
    assert sp.simplify(eq[v] + k * x / zeta) == 0


def test_rayleighian_and_friction_matrix():
    x1, x2, v1, v2, a, b, c = sp.symbols("x1 x2 v1 v2 a b c")
    Phi = (a * v1**2 + 2 * b * v1 * v2 + c * v2**2) / 2
    assert on.friction_matrix(Phi, [v1, v2]) == sp.Matrix([[a, b], [b, c]])
    R = on.rayleighian(x1 * x2, Phi, [x1, x2], [v1, v2])
    assert sp.expand(R - Phi - x2 * v1 - x1 * v2) == 0


def test_energy_balance_identity():
    """At the minimum, dA/dt = -2 Phi."""
    x1, x2, v1, v2 = sp.symbols("x1 x2 v1 v2")
    A = x1**2 + x1 * x2 + 3 * x2**2
    Phi = (2 * v1**2 + 2 * v1 * v2 + 5 * v2**2) / 2
    eq = on.kinetic_equations(A, Phi, [x1, x2], [v1, v2])
    Adot = sum(sp.diff(A, xi) * eq[vi] for xi, vi in ((x1, v1), (x2, v2)))
    assert sp.simplify(Adot + 2 * Phi.subs(eq)) == 0


def test_rod_drift_angle():
    """A slender rod (zeta_perp = 2 zeta_par) drifts sideways by at most asin(1/3)."""
    X, Z, vx, vz, th, z0, F = sp.symbols("X Z v_x v_z theta zeta_0 F", positive=True)
    u = sp.Matrix([sp.sin(th), sp.cos(th)])
    zeta = z0 * (u * u.T) + 2 * z0 * (sp.eye(2) - u * u.T)
    v = sp.Matrix([vx, vz])
    Phi = (v.T * zeta * v)[0] / 2
    eq = on.kinetic_equations(F * Z, Phi, [X, Z], [vx, vz])
    angle = sp.lambdify(th, sp.atan(-eq[vx] / eq[vz]))
    t = np.linspace(0, np.pi / 2, 20001)
    assert abs(np.max(np.abs(angle(t))) - math.asin(1 / 3)) < 1e-7


def test_constraint_gives_hydrostatic_pressure():
    h1, h2, q = sp.symbols("h1 h2 q")
    v1, v2, Q, p1, p2 = sp.symbols("v1 v2 Q p1 p2")
    S1, S2, rg, Rh = sp.symbols("S1 S2 rho_g R_h", positive=True)
    A = rg * (S1 * h1**2 + S2 * h2**2) / 2
    eq = on.kinetic_equations(A, Rh * Q**2 / 2, [h1, h2, q], [v1, v2, Q],
                              constraints=[S1 * v1 + Q, S2 * v2 - Q], multipliers=[p1, p2])
    assert sp.simplify(eq[p1] - rg * h1) == 0
    assert sp.simplify(eq[p2] - rg * h2) == 0
    assert sp.simplify(eq[Q] - rg * (h1 - h2) / Rh) == 0


def test_nonquadratic_dissipation_is_refused():
    x, v = sp.symbols("x v")
    with pytest.raises(ValueError):
        on.kinetic_equations(x**2, v**4, [x], [v])


def test_rate_function():
    x, v, k = sp.symbols("x v k")
    eq = on.kinetic_equations(k * x**2 / 2, v**2 / 2, [x], [v])
    f = on.rate_function(eq, [x], [v], {k: 3.0})
    assert f(0.0, [2.0])[0] == pytest.approx(-6.0)


def test_friction_of_a_gaussian_profile():
    """For a Gaussian of N particles and width s, v = x sdot/s and zeta_ss = zeta N."""
    N, s = 2.5, 0.7
    x = np.linspace(-12 * s, 12 * s, 40001)
    c = N / (math.sqrt(2 * math.pi) * s) * np.exp(-x**2 / (2 * s**2))
    dc_ds = c * (x**2 / s**3 - 1 / s)
    Z = on.friction_matrix_1d(x, c, [dc_ds], zeta=1.3)
    assert Z[0, 0] == pytest.approx(1.3 * N, rel=1e-6)


def test_friction_of_a_linear_profile():
    """c = 1 - x/l on 0 < x < l with rate ldot: zeta_ll = (11/48) l."""
    ell = 2.0
    x = np.linspace(0, ell, 20001)
    Z = on.friction_matrix_1d(x, 1 - x / ell, [x / ell**2])
    assert Z[0, 0] == pytest.approx(11 / 48 * ell, rel=1e-5)


def test_friction_of_power_and_exponential_profiles():
    """(1 - x)**3: I(3) = 21/480;  exp(-x): 5  (Chapter 5)."""
    x = np.linspace(0, 1, 20001)
    Z = on.friction_matrix_1d(x, (1 - x) ** 3, [3 * (1 - x) ** 2 * x])
    assert Z[0, 0] == pytest.approx(21 / 480, rel=1e-6)
    x = np.linspace(0, 40, 400001)
    Z = on.friction_matrix_1d(x, np.exp(-x), [x * np.exp(-x)])
    assert Z[0, 0] == pytest.approx(5.0, rel=1e-6)


def test_friction_matrix_is_symmetric_and_positive():
    x = np.linspace(0, 1, 4001)
    s = 1 - x
    c = s**2.5
    dc_dl = 2.5 * s**1.5 * x
    dc_dn = np.where(s > 0, c * np.log(np.maximum(s, 1e-300)), 0.0)
    Z = on.friction_matrix_1d(x, c, [dc_dl, dc_dn])
    assert Z[0, 1] == Z[1, 0]
    assert np.all(np.linalg.eigvalsh(Z) > 0)


def test_flux_from_rate():
    x = np.linspace(0, 1, 1001)
    assert on.flux_from_rate(x, np.ones_like(x))[0] == pytest.approx(1.0)
    assert on.flux_from_rate(x, np.ones_like(x), closed="left")[-1] == pytest.approx(-1.0)
    with pytest.raises(ValueError):
        on.flux_from_rate(x, x, closed="both")


def test_finite_volume_diffusion():
    """Free diffusion of a Gaussian: conserved, and equal to the exact solution."""
    edges = np.linspace(-10, 10, 801)
    xc = 0.5 * (edges[1:] + edges[:-1])
    gauss = lambda x, var: np.exp(-x**2 / (2 * var)) / np.sqrt(2 * np.pi * var)
    t = np.array([0.0, 0.5, 2.0])
    _, c = on.smoluchowski_fv(gauss(xc, 0.25), edges, t)
    dx = edges[1] - edges[0]
    assert np.allclose(c.sum(axis=1) * dx, 1.0, atol=1e-10)
    assert on.sup_error(c[-1], gauss(xc, 0.25 + 2 * 2.0)) < 2e-4


def test_finite_volume_equilibrium_in_a_trap():
    """In U = x**2/2 the profile relaxes to the Boltzmann distribution."""
    edges = np.linspace(-8, 8, 401)
    xc = 0.5 * (edges[1:] + edges[:-1])
    c0 = np.exp(-(xc - 2) ** 2 / 0.5) / np.sqrt(0.5 * np.pi)
    _, c = on.smoluchowski_fv(c0, edges, np.array([0.0, 20.0]), U=lambda x: x**2 / 2)
    assert on.sup_error(c[-1], np.exp(-xc**2 / 2) / np.sqrt(2 * np.pi)) < 1e-4


def test_duct():
    assert on.duct_flow_rate() == pytest.approx(0.5623078, abs=2e-6)
    # a very wide duct tends to the plane channel, Q = (2/3) per unit width
    assert on.duct_flow_rate(b=50.0) / 100 == pytest.approx(2 / 3, rel=0.02)
    # flow rate from the velocity field
    g, w = np.polynomial.legendre.leggauss(60)
    X, Y = np.meshgrid(g, g)
    Q = np.sum(np.outer(w, w) * on.duct_velocity(X, Y))
    assert Q == pytest.approx(on.duct_flow_rate(), rel=1e-6)
    assert on.duct_velocity(1.0, 0.3) == pytest.approx(0.0, abs=1e-12)


def test_compare():
    eps = np.logspace(-3, -1, 9)
    slope, pref = on.loglog_slope(eps, 3 * eps**2)
    assert slope == pytest.approx(2.0) and pref == pytest.approx(3.0)
    assert on.sup_error([1, 2, 3], [1, 2.5, 3]) == 0.5
    assert on.relative_error(1.1, 1.0) == pytest.approx(0.1)
    with pytest.raises(ValueError):
        on.loglog_slope([1.0], [1.0])
