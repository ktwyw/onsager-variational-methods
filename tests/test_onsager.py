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


# ---------------------------------------------------------------- thin films

def _cap(n=400, length=2.4, b=0.01):
    edges = np.linspace(0.0, length, n + 1)
    r = 0.5 * (edges[1:] + edges[:-1])
    return edges, r, np.maximum(1 - r**2, 0.0) + b


def test_thin_film_conserves_volume_and_spreads():
    edges, r, h0 = _cap()
    t = np.array([0.0, 1.0, 10.0])
    r, h = on.thin_film_fv(h0, edges, t)
    volume = 2 * np.pi * np.sum(h * r, axis=1) * (edges[1] - edges[0])
    assert np.ptp(volume) < 1e-10
    assert h.min() > 0
    radii = [on.fit_cap(r, hk, b=0.01)[1] for hk in h]
    assert radii[0] == pytest.approx(1.0, abs=2e-3)
    assert radii[0] < radii[1] < radii[2]
    # Tanner's law with its logarithm, dR/dt = 4 / (l R**9) with l = ln(h_c / b) - 1:
    # the mean rate between t = 1 and t = 10 lies between the rates at the two ends
    rate = (radii[2] - radii[1]) / 9.0
    tanner = [4 / ((np.log(1 / (R**2 * 0.01)) - 1) * R**9) for R in radii[1:]]
    assert 0.8 * tanner[1] < rate < 1.2 * tanner[0]


def test_thin_film_partial_wetting_reaches_the_equilibrium_angle():
    edges, r, h0 = _cap(n=500, length=2.0, b=0.005)
    r, h = on.thin_film_fv(h0, edges, np.array([0.0, 300.0]), theta_e=0.6, b=0.005)
    hc, R = on.fit_cap(r, h[-1], b=0.005)
    assert 2 * hc / R == pytest.approx(0.6, rel=0.06)
    assert R == pytest.approx((2 / 0.6) ** (1 / 3), rel=0.03)


def test_thin_film_planar_and_bad_arguments():
    edges, r, h0 = _cap(n=300, length=2.0)
    r, h = on.thin_film_fv(h0, edges, np.array([0.0, 1.0]), geometry="planar")
    assert np.ptp(np.sum(h, axis=1)) < 1e-10
    with pytest.raises(ValueError):
        on.thin_film_fv(h0, edges, np.array([0.0, 1.0]), geometry="spherical")
    with pytest.raises(ValueError):
        on.thin_film_fv(h0, edges, np.array([0.0, 1.0]), theta_e=0.5)


def test_fit_cap_recovers_a_parabola():
    r = np.linspace(0.005, 3, 300)
    h = np.maximum(0.7 * (1 - r**2 / 1.6**2), 0) + 0.02
    hc, R = on.fit_cap(r, h, b=0.02)
    assert hc == pytest.approx(0.7) and R == pytest.approx(1.6)


# ---------------------------------------------------------------- coalescence

def test_hopper_solution_limits_and_area():
    assert on.hopper_time(1.0) == 0.0
    assert on.hopper_neck(1.0) == 0.0 and on.hopper_neck(0.0) == pytest.approx(np.sqrt(2))
    assert on.hopper_length(1.0) == pytest.approx(2.0)          # two touching circles of radius 1
    for m in (0.9, 0.5, 0.1):
        X, Y = on.hopper_shape(m, n=20001)
        area = 0.5 * abs(np.sum(X[:-1] * Y[1:] - X[1:] * Y[:-1]))
        assert area == pytest.approx(2 * np.pi, rel=1e-5)          # area of two unit circles
        assert np.max(X) == pytest.approx(on.hopper_length(m), rel=1e-6)
        assert X[0] == pytest.approx(0, abs=1e-12) and Y[0] == pytest.approx(on.hopper_neck(m), rel=1e-12)
        polygon = np.sum(np.hypot(np.diff(X), np.diff(Y)))
        assert on.hopper_perimeter(m) == pytest.approx(polygon, rel=1e-5)
    assert on.hopper_perimeter(0.0) == pytest.approx(2 * np.sqrt(2) * np.pi)


def test_hopper_time_asymptotes_and_inverse():
    # late times: m ~ exp(-gamma t / eta R_f), R_f = sqrt(2)
    t1, t2 = on.hopper_time(1e-3), on.hopper_time(1e-4)
    assert (t2 - t1) / np.log(10) == pytest.approx(np.sqrt(2), rel=1e-3)
    # early times: t = pi x / |ln(x / 8)| (Eggers, Lister and Stone 1999), to 10 %
    m = 1 - 1e-6
    x = on.hopper_neck(m)
    assert on.hopper_time(m) == pytest.approx(np.pi * x / abs(np.log(x / 8)), rel=0.1)
    for t in (0.01, 0.5, 3.0):
        assert on.hopper_time(on.hopper_parameter(t)) == pytest.approx(t, rel=1e-9)


# ---------------------------------------------------------------- drying drops

def test_drying_film_conserves_solute_and_dries():
    b = 0.005
    edges = np.linspace(0.0, 1.5, 301)
    r = 0.5 * (edges[1:] + edges[:-1])
    cap = np.maximum(1 - r**2, 0.0)
    out = on.drying_film_fv(cap + b, cap, edges, evaporation=1.0, theta_e=2.0, b=b, stop_volume=0.2, n_out=20)
    assert out["volume"][-1] == pytest.approx(0.2 * out["volume"][0], rel=1e-3)
    assert np.ptp(out["solute"]) < 1e-9
    assert np.all(np.diff(out["volume"]) < 0)
    # the contact line recedes: less liquid beyond r = 0.9 at the end
    assert np.all(out["h"][-1][r > 0.9] < 2 * b)


# ---------------------------------------------------------------- sedimentation

def test_sedimentation_matches_mason_weaver_and_conserves():
    edges = np.linspace(0.0, 1.0, 401)
    z = 0.5 * (edges[1:] + edges[:-1])
    t = np.array([0.0, 0.05, 0.5, 5.0])
    z, phi = on.sedimentation_fv(np.ones_like(z), edges, t, alpha=5.0)
    assert np.ptp(phi.sum(axis=1)) * (z[1] - z[0]) < 1e-12
    for k in (1, 2, 3):
        assert np.max(np.abs(phi[k] - on.mason_weaver(z, t[k], 5.0))) < 2e-4
    # equilibrium: the barometric profile
    assert np.max(np.abs(phi[-1] - 5 * np.exp(-5 * z) / (1 - np.exp(-5)))) < 1e-4


def test_mason_weaver_limits():
    z = np.linspace(0, 1, 11)
    assert np.allclose(on.mason_weaver(z, 50.0, 3.0), 3 * np.exp(-3 * z) / (1 - np.exp(-3)))
    zz = np.linspace(0, 1, 2001)
    assert np.trapezoid(on.mason_weaver(zz, 0.1, 3.0), zz) == pytest.approx(1.0, rel=1e-6)


def test_sedimentation_front_moves_at_the_hindered_velocity():
    edges = np.linspace(0.0, 1.0, 801)
    z = 0.5 * (edges[1:] + edges[:-1])
    K = lambda p: (1 - p) ** 4.65
    cs = lambda p: p * (1 + p + p**2 - p**3) / (1 - p) ** 3
    z, phi = on.sedimentation_fv(np.full_like(z, 0.1), edges, np.array([0.0, 0.5]), 60.0, cs, K)
    front = z[np.where(phi[-1] > 0.05)[0].max()]
    assert front == pytest.approx(1 - 0.5 * K(0.1), abs=0.01)


# ------------------------------------------------------------------ gels

def test_gel_sphere_roots_limits():
    assert np.allclose(on.gel_sphere_roots(0.0, 4), np.pi * np.arange(1, 5), atol=1e-10)
    # shear-dominated gel: the slowest rate tends to 15 K / f a^2 (uniform dilation)
    alpha = on.gel_sphere_roots(1e4, 1)[0]
    assert (1 + 4e4 / 3) * alpha**2 == pytest.approx(15.0, rel=1e-3)


def test_gel_sphere_series_is_tanaka_fillmore_for_zero_shear():
    t = 0.03
    n = np.arange(1, 400)
    expected = 1 - np.sum(6 / (n * np.pi) ** 2 * np.exp(-(n * np.pi) ** 2 * t))   # swelling of the volume
    r = np.linspace(0, 1, 4001)
    u = on.gel_sphere_series(r, t, 0.0)
    dilation = np.gradient(u, r) + 2 * np.divide(u, r, out=np.ones_like(u), where=r > 0)
    assert np.trapezoid(dilation * r**2, r) == pytest.approx(expected, abs=2e-4)   # mean of div u over 3 = (u(1) for volume)
    assert u[-1] == pytest.approx(expected, abs=1e-6)


def test_gel_finite_elements_match_the_series():
    t = np.array([0.01, 0.1])
    x, u = on.gel_fe(300, 0.6, t)
    for k in range(2):
        assert np.max(np.abs(u[k] - on.gel_sphere_series(x, t[k], 0.6))) < 2e-5
    x, u = on.gel_fe(300, 0.6, t, geometry="sheet")
    for k in range(2):
        assert np.max(np.abs(u[k] - on.gel_sheet_series(x, t[k], 0.6))) < 2e-5


# ------------------------------------------------------------------ phase separation

def _fh(chi):
    f = lambda p: p * np.log(p) + (1 - p) * np.log(1 - p) + chi * p * (1 - p)
    df = lambda p: np.log(np.clip(p, 1e-12, 1) / np.clip(1 - p, 1e-12, 1)) + chi * (1 - 2 * p)
    return f, df


def test_cahn_hilliard_linear_growth_rates():
    chi = 2.5
    _, df = _fh(chi)
    n, box = 256, 2 * np.pi / 0.25
    rng = np.random.default_rng(0)
    t = np.linspace(0, 20, 11)
    phi = on.cahn_hilliard_spectral(0.5 + 1e-6 * rng.standard_normal(n), box, t, df, 0.01)
    amp = np.abs(np.fft.rfft(phi, axis=1))
    k = 2 * np.pi * np.fft.rfftfreq(n, d=box / n)
    for m in (1, 2, 3):
        rate = np.polyfit(t[3:], np.log(amp[3:, m]), 1)[0]
        assert rate == pytest.approx(-k[m] ** 2 * (4 - 2 * chi + k[m] ** 2), rel=1e-3)
    assert np.ptp(phi.mean(axis=1)) < 1e-12


def test_cahn_hilliard_interface_tension():
    chi = 3.0
    f, df = _fh(chi)
    from scipy.optimize import brentq
    from scipy.integrate import quad
    pe = brentq(lambda p: df(p), 1e-9, 0.4999)
    gamma = quad(lambda p: np.sqrt(max(2 * (f(p) - f(pe)), 0.0)), pe, 1 - pe, limit=200)[0]
    n, box = 512, 100.0
    x = np.arange(n) * box / n
    phi = on.cahn_hilliard_spectral(np.where((x > 25) & (x < 75), 1 - pe, pe), box, [0.0, 100.0], df, 0.01)[-1]
    k = 2 * np.pi * np.fft.fftfreq(n, d=box / n)
    grad = np.real(np.fft.ifft(1j * k * np.fft.fft(phi)))
    excess = np.sum(f(phi) - f(pe) + 0.5 * grad**2) * box / n / 2
    assert excess == pytest.approx(gamma, rel=1e-4)


def test_cahn_hilliard_free_energy_decreases_in_2d():
    chi = 2.5
    f, df = _fh(chi)
    rng = np.random.default_rng(2)
    n, box = 64, 64.0
    t = np.array([0.0, 20.0, 40.0, 80.0, 160.0])
    snaps = on.cahn_hilliard_spectral(0.5 + 0.01 * rng.standard_normal((n, n)), box, t, df, 0.5)
    k = 2 * np.pi * np.fft.fftfreq(n, d=box / n)
    kx, ky = np.meshgrid(k, k, indexing="ij")
    def energy(p):
        h = np.fft.fft2(p)
        gx, gy = np.real(np.fft.ifft2(1j * kx * h)), np.real(np.fft.ifft2(1j * ky * h))
        return np.mean(f(p) + 0.5 * (gx**2 + gy**2))
    e = [energy(p) for p in snaps]
    assert all(b < a for a, b in zip(e, e[1:]))
    assert np.ptp(snaps.mean(axis=(1, 2))) < 1e-12


# ------------------------------------------------------------------ polymer chains

def test_rouse_brownian_dynamics_equilibrium_and_decay():
    rng = np.random.default_rng(5)
    t, ends, middle = on.rouse_brownian_dynamics(8, 3000, 0.005, 6.0, 200, rng)
    size = np.mean(np.sum(ends**2, axis=2), axis=1)
    assert np.all(np.abs(size / 7.0 - 1) < 0.06)                      # <R^2> = (N - 1) b^2 throughout
    corr = np.mean(np.sum(ends * ends[0], axis=2), axis=1) / 7.0
    lam1 = 12 * np.sin(np.pi / 16) ** 2                                  # slowest Rouse rate, N = 8
    assert corr[-1] == pytest.approx(np.exp(-lam1 * t[-1]), abs=0.05)


def test_fene_chain_steady_state():
    N, Q = 20, 10.0
    lam1 = 12 * np.sin(np.pi / (2 * N)) ** 2
    assert on.fene_chain_steady(N, 0.95 * lam1, Q) is None                # below the threshold: no stretch
    x = on.fene_chain_steady(N, 3 * lam1, Q)
    assert np.allclose(x, -x[::-1], atol=1e-8)                           # symmetric about the stagnation point
    s = np.diff(x)
    T = 3.0 * s / (1 - (s / Q) ** 2)
    assert np.allclose(np.diff(T), -3 * lam1 * x[1:-1], rtol=1e-6, atol=1e-8)  # force balance on the inner beads


# ------------------------------------------------------------------ nematic liquid crystals

def test_freedericksz_twist_angle():
    assert on.freedericksz_twist_angle(0.99) == 0.0
    h = 1 + 1e-6
    assert on.freedericksz_twist_angle(h) == pytest.approx(2 * np.sqrt(h - 1), rel=1e-5)
    assert on.freedericksz_twist_angle(20.0) == pytest.approx(np.pi / 2, abs=1e-6)


def test_director_fd_matches_the_exact_twist_and_the_linear_rate():
    n = 199
    z = np.arange(1, n + 1) / (n + 1)
    h = 2.0
    torque = lambda th: h * h * np.sin(th) * np.cos(th)
    final = on.director_fd(0.01 * np.sin(np.pi * z), np.array([0.0, 60.0]), torque, 1 / np.pi**2)[-1]
    assert final.max() == pytest.approx(on.freedericksz_twist_angle(h), abs=1e-4)
    t = np.linspace(0.0, 1.0, 6)
    early = on.director_fd(1e-7 * np.sin(np.pi * z), t, torque, 1 / np.pi**2)
    assert np.polyfit(t, np.log(early[:, n // 2]), 1)[0] == pytest.approx(h * h - 1, rel=1e-3)
