"""Reference solutions against which the variational results are checked."""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags

__all__ = ["smoluchowski_fv", "duct_flow_rate", "duct_velocity", "thin_film_fv",
           "hopper_time", "hopper_parameter", "hopper_neck", "hopper_length", "hopper_shape",
           "hopper_perimeter", "drying_film_fv", "sedimentation_fv", "mason_weaver",
           "gel_sphere_roots", "gel_sphere_series", "gel_sheet_series", "gel_fe",
           "cahn_hilliard_spectral", "rouse_brownian_dynamics", "fene_chain_steady",
           "director_fd", "freedericksz_twist_angle"]


def smoluchowski_fv(c0, edges, t_eval, mu_excess=None, U=None, kT=1.0, zeta=1.0,
                    rtol=1e-8, atol=1e-11):
    """Finite-volume solution of one-dimensional diffusion in conservative form.

        dc/dt = -dJ/dx ,   J = -(c/zeta) d/dx [ kT ln c + mu_excess(c) + U(x) ] ,

    with J = 0 at both ends.  The ideal part of the flux is -(kT/zeta) dc/dx; the
    rest is the mobility c/zeta times the gradient of mu_excess(c) + U(x).

    c0 : initial cell averages (length n);  edges : cell edges (length n + 1);
    mu_excess : function of c, the non-ideal part of the chemical potential;
    U : function of x, an external potential.  Fluxes are evaluated on the cell
    faces, so the number of particles is conserved to rounding error.

    Returns (centres, c) with c of shape (len(t_eval), n).
    """
    edges = np.asarray(edges, dtype=float)
    xc = 0.5 * (edges[1:] + edges[:-1])
    dx = np.diff(edges)
    dxc = np.diff(xc)
    Uc = np.zeros_like(xc) if U is None else U(xc)

    def rhs(t, c):
        extra = Uc if mu_excess is None else Uc + mu_excess(c)
        J = -(kT * np.diff(c) + 0.5 * (c[1:] + c[:-1]) * np.diff(extra)) / (zeta * dxc)
        J = np.concatenate(([0.0], J, [0.0]))
        return -np.diff(J) / dx

    n = len(xc)
    pattern = diags([np.ones(n - 1), np.ones(n), np.ones(n - 1)], [-1, 0, 1])   # tridiagonal Jacobian
    sol = solve_ivp(rhs, (0.0, t_eval[-1]), np.asarray(c0, dtype=float), t_eval=t_eval,
                    method="BDF", rtol=rtol, atol=atol, jac_sparsity=pattern)
    if not sol.success:   # pragma: no cover
        raise RuntimeError(sol.message)
    return xc, sol.y.T


def duct_flow_rate(b=1.0, terms=200):
    """Exact flow rate in a rectangular duct, -1 < x < 1, -b < y < b.

    Lengths are in units of the half-width and the result is Q * eta / (G a**4),
    where G is the pressure gradient.  For the square duct (b = 1) it is 0.56231.
    """
    s = sum(math.tanh(n * math.pi * b / 2) / n**5 for n in range(1, 2 * terms, 2))
    return 4 * b / 3 * (1 - 192 / (math.pi**5 * b) * s)


def duct_velocity(x, y, b=1.0, terms=200):
    """Exact velocity w * eta / (G a**2) in the rectangular duct of duct_flow_rate."""
    x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
    w = np.zeros(x.shape)
    for k, n in enumerate(range(1, 2 * terms, 2)):
        arg = n * math.pi / 2
        # cosh(arg*y)/cosh(arg*b), written so that nothing overflows
        ratio = np.exp(arg * (np.abs(y) - b)) * (1 + np.exp(-2 * arg * np.abs(y))) / (1 + math.exp(-2 * arg * b))
        w += (-1) ** k / n**3 * (1 - ratio) * np.cos(arg * x)
    return 16 / math.pi**3 * w


def thin_film_fv(h0, edges, t_eval, geometry="axisymmetric", bond=0.0, theta_e=0.0, b=None,
                 exponents=(3, 2), rtol=1e-6, atol=1e-9):
    """Finite-volume solution of the thin-film (lubrication) equation in scaled form,

        dh/dt = -div( h**3 grad( lap h - bond * h + Pi(h) ) ) ,

    for a film of thickness h(r, t) on a flat solid, with no flux through either end
    of the domain.  Lengths along the solid are in units of L, thicknesses in units
    of H and times in units of 3 eta L**4 / (gamma H**3); bond = rho g L**2 / gamma.

    geometry : "axisymmetric" (r is the distance from the axis, edges[0] = 0) or
               "planar" (r is a Cartesian coordinate).
    theta_e, b : the solid is covered by a precursor film of thickness b, which
               should be part of h0.  With theta_e = 0 the liquid wets the solid
               completely and Pi = 0.  With theta_e > 0 (the equilibrium slope, in
               units of H / L) the disjoining pressure

                   Pi = B [ (b/h)**n - (b/h)**m ],  B = (n-1)(m-1) theta_e**2 / (2 b (n-m)),

               makes b the equilibrium film and theta_e the equilibrium contact
               angle (Schwartz and Eley, 1998); (n, m) = exponents.

    Fluxes are evaluated on the cell faces, so the volume is conserved to rounding
    error; the mobility on a face is 2 h1**2 h2**2 / (h1 + h2), which keeps h positive.

    Returns (centres, h) with h of shape (len(t_eval), n).
    """
    edges = np.asarray(edges, dtype=float)
    rc = 0.5 * (edges[1:] + edges[:-1])
    dr = np.diff(edges)
    drc = np.diff(rc)
    if geometry == "axisymmetric":
        wf, wc = edges, rc                       # measure: r dr
    elif geometry == "planar":
        wf, wc = np.ones_like(edges), np.ones_like(rc)
    else:
        raise ValueError("geometry must be 'axisymmetric' or 'planar'")
    n_, m_ = exponents
    if theta_e > 0:
        if b is None:
            raise ValueError("a precursor thickness b is needed when theta_e > 0")
        B = (n_ - 1) * (m_ - 1) * theta_e**2 / (2 * b * (n_ - m_))

    def rhs(t, h):
        slope = np.concatenate(([0.0], wf[1:-1] * np.diff(h) / drc, [0.0]))
        p = -np.diff(slope) / (wc * dr) + bond * h          # pressure in units of gamma H / L**2
        if theta_e > 0:
            x = b / h
            p = p - B * (x**n_ - x**m_)
        mobility = 2 * h[1:] ** 2 * h[:-1] ** 2 / (h[1:] + h[:-1])
        q = np.concatenate(([0.0], -wf[1:-1] * mobility * np.diff(p) / drc, [0.0]))
        return -np.diff(q) / (wc * dr)

    n = len(rc)
    pattern = diags([np.ones(n - abs(k)) for k in range(-2, 3)], list(range(-2, 3)))   # five bands
    sol = solve_ivp(rhs, (0.0, t_eval[-1]), np.asarray(h0, dtype=float), t_eval=t_eval,
                    method="BDF", rtol=rtol, atol=atol, jac_sparsity=pattern)
    if not sol.success:   # pragma: no cover
        raise RuntimeError(sol.message)
    return rc, sol.y.T


# ------------------------------------------------- coalescence of two viscous cylinders

def hopper_time(m):
    """Time at which Hopper's exact solution has the parameter m (0 <= m <= 1).

    Two equal viscous cylinders of radius a0 touch at t = 0 (m = 1) and coalesce
    under surface tension into one cylinder of radius sqrt(2) a0 (m = 0), with
    no outer fluid and no inertia (Hopper, 1984, 1990).  Time is in units of
    eta a0 / gamma:

        t = (pi sqrt(2) / 4) * integral_{m^2}^{1} ds / (s sqrt(1 + s) K(s)),

    with K the complete elliptic integral of the first kind, K(s) = ellipk(s).
    """
    from scipy.integrate import quad
    from scipy.special import ellipk

    def one(mm):
        if mm >= 1.0:
            return 0.0
        if mm <= 0.0:
            return np.inf
        f = lambda s: 1.0 / (s * np.sqrt(1.0 + s) * ellipk(s))
        return np.pi * np.sqrt(2.0) / 4.0 * quad(f, mm * mm, 1.0, limit=200)[0]
    m = np.asarray(m, dtype=float)
    return np.vectorize(one)(m) if m.ndim else one(float(m))


def hopper_parameter(t):
    """Inverse of hopper_time: the parameter m at time t (units of eta a0 / gamma)."""
    from scipy.optimize import brentq

    def one(tt):
        if tt <= 0:
            return 1.0
        # for large t, m ~ exp(-t / sqrt 2): bracket generously
        lo = max(np.exp(-tt / np.sqrt(2.0)) * 1e-6, 1e-300)
        return brentq(lambda mm: hopper_time(mm) - tt, lo, 1.0, xtol=1e-15, rtol=1e-13)
    t = np.asarray(t, dtype=float)
    return np.vectorize(one)(t) if t.ndim else one(float(t))


def hopper_neck(m):
    """Half-width of the neck, in units of a0: sqrt(2) (1 - m) / sqrt(1 + m^2)."""
    m = np.asarray(m, dtype=float)
    return np.sqrt(2.0) * (1 - m) / np.sqrt(1 + m * m)


def hopper_length(m):
    """Half-length of the pair along the line of centres, in units of a0."""
    m = np.asarray(m, dtype=float)
    return np.sqrt(2.0) * (1 + m) / np.sqrt(1 + m * m)


def hopper_shape(m, n=721):
    """Outline (X, Y) of Hopper's shape with parameter m, in units of a0.

    The line of centres is the X axis and the neck lies on the Y axis.  The
    outline is the image of the unit circle under the conformal map
    z = C zeta / (1 + m zeta^2), C = sqrt(2) (1 - m^2) / sqrt(1 + m^2), rotated
    by 90 degrees.
    """
    zeta = np.exp(2j * np.pi * np.arange(n) / (n - 1))
    C = np.sqrt(2.0) * (1 - m * m) / np.sqrt(1 + m * m)
    z = 1j * C * zeta / (1 + m * zeta * zeta)
    return z.real, z.imag


def hopper_perimeter(m, n=1 << 14):
    """Perimeter of Hopper's shape with parameter m < 1, in units of a0.

    Computed as the integral of |dz/dzeta| around the unit circle; the trapezoidal
    rule converges geometrically because the integrand is smooth and periodic.
    """
    zeta = np.exp(2j * np.pi * np.arange(n) / n)
    C = np.sqrt(2.0) * (1 - m * m) / np.sqrt(1 + m * m)
    return 2 * np.pi * float(np.mean(np.abs(C * (1 - m * zeta**2) / (1 + m * zeta**2) ** 2)))


def drying_film_fv(h0, s0, edges, evaporation, theta_e, b, stop_volume=0.03, t_max=1e4,
                   n_out=200, rtol=1e-7, atol=1e-10):
    """An evaporating, partially wetting drop that carries a dilute solute.

    Solves, in the scaled units of thin_film_fv (axisymmetric),

        dh/dt = -(1/r) d/dr(r q) - J max(0, 1 - b/h),   q = -h^3 dp/dr,
        p = -(1/r) d/dr(r dh/dr) - B [(b/h)^3 - (b/h)^2],   B = theta_e^2 / b,

    on a precursor film of thickness b that neither flows nor evaporates, and
    then carries the solute (amount per unit area s = h * concentration) with the
    liquid flux, ds/dt = -(1/r) d/dr(r q s/h), by a first-order upwind scheme.
    The solute does not act back on the liquid.  Solute left in the precursor
    film cannot move: it is the deposit.

    Integration stops when the liquid volume above the precursor film has fallen
    to stop_volume times its initial value.

    Returns a dict with r (cell centres), t (n_out times), h and s (arrays of
    shape (n_out, n)), volume and solute (totals at those times).
    """
    from scipy.sparse import diags

    edges = np.asarray(edges, dtype=float)
    rc = 0.5 * (edges[1:] + edges[:-1])
    dr = np.diff(edges)
    drc = np.diff(rc)
    B = theta_e**2 / b

    def flux(h):
        slope = np.concatenate(([0.0], edges[1:-1] * np.diff(h) / drc, [0.0]))
        p = -np.diff(slope) / (rc * dr) - B * ((b / h) ** 3 - (b / h) ** 2)
        mobility = 2 * h[1:] ** 2 * h[:-1] ** 2 / (h[1:] + h[:-1])
        return -mobility * np.diff(p) / drc              # liquid flux on the interior faces

    def divergence(face_values):
        F = np.concatenate(([0.0], edges[1:-1] * face_values, [0.0]))
        return np.diff(F) / (rc * dr)

    def rhs(t, h):
        return -divergence(flux(h)) - evaporation * np.clip(1 - b / h, 0, None)

    def volume(h):
        return 2 * np.pi * np.sum((h - b) * rc * dr)

    v0 = volume(np.asarray(h0, dtype=float))
    stop = lambda t, h: volume(h) - stop_volume * v0
    stop.terminal = True
    n = len(rc)
    pattern = diags([np.ones(n - abs(k)) for k in range(-2, 3)], list(range(-2, 3)))
    film = solve_ivp(rhs, (0.0, t_max), np.asarray(h0, dtype=float), method="BDF", rtol=rtol, atol=atol,
                     jac_sparsity=pattern, dense_output=True, events=stop)
    if not film.success:   # pragma: no cover
        raise RuntimeError(film.message)
    t_out = np.linspace(0.0, film.t[-1], n_out)

    s = np.asarray(s0, dtype=float).copy()
    S = np.empty((n_out, n)); S[0] = s
    t, k = 0.0, 1
    while k < n_out:
        h = film.sol(t)
        q = flux(h)
        speed = np.max(np.abs(q) / np.minimum(h[1:], h[:-1]))
        dt = min(0.4 * dr.min() / max(speed, 1e-12), t_out[k] - t)
        concentration = s / h
        upwind = np.where(q > 0, concentration[:-1], concentration[1:])
        s = s - dt * divergence(q * upwind)
        t += dt
        if t >= t_out[k] - 1e-14:
            S[k] = s
            k += 1
    H = film.sol(t_out).T
    return dict(r=rc, t=t_out, h=H, s=S, volume=np.array([volume(x) for x in H]),
                solute=2 * np.pi * np.sum(S * rc * dr, axis=1))


# --------------------------------------------------------------- sedimentation

def sedimentation_fv(phi0, edges, t_eval, alpha, osmotic=None, hindrance=None, rtol=1e-8, atol=1e-11):
    """Sedimentation and diffusion of a suspension in a closed vertical cell.

    Solves, with z upwards in units of the height L of the cell and t in units of
    L / U0 (U0 the sedimentation velocity of a single particle),

        dphi/dt = -dj/dz ,   j = -K(phi) [ (1/alpha) d Pi(phi)/dz + phi ] ,

    with j = 0 at the bottom and the top.  alpha = L U0 / D0 = L / l_g is the
    height of the cell in units of the gravitational length l_g = kT / m_b g.
    Pi = osmotic(phi) is the osmotic pressure in units of kT / v_c (default:
    the ideal value phi) and K = hindrance(phi) the hindered settling function
    (default 1).  phi0 holds the initial cell averages on the cells between
    `edges`.  Returns (centres, phi) with phi of shape (len(t_eval), n).
    """
    edges = np.asarray(edges, dtype=float)
    zc = 0.5 * (edges[1:] + edges[:-1])
    dz = np.diff(edges)
    dzc = np.diff(zc)
    P = (lambda p: p) if osmotic is None else osmotic
    K = (lambda p: np.ones_like(p)) if hindrance is None else hindrance

    def rhs(t, phi):
        face = 0.5 * (phi[1:] + phi[:-1])
        j = -K(face) * (np.diff(P(phi)) / (alpha * dzc) + face)
        return -np.diff(np.concatenate(([0.0], j, [0.0]))) / dz

    n = len(zc)
    pattern = diags([np.ones(n - 1), np.ones(n), np.ones(n - 1)], [-1, 0, 1])
    sol = solve_ivp(rhs, (0.0, t_eval[-1]), np.asarray(phi0, dtype=float), t_eval=t_eval,
                    method="BDF", rtol=rtol, atol=atol, jac_sparsity=pattern)
    if not sol.success:   # pragma: no cover
        raise RuntimeError(sol.message)
    return zc, sol.y.T


def mason_weaver(z, t, alpha, n_terms=300, n_quad=600):
    """Exact solution of the dilute problem of sedimentation_fv (Mason and Weaver, 1924).

    The suspension starts uniform, phi = 1, at t = 0.  The solution is the
    barometric profile alpha exp(-alpha z) / (1 - exp(-alpha)) plus the modes

        exp(-alpha z / 2) [cos(n pi z) - alpha / (2 n pi) sin(n pi z)] exp(-lambda_n t),
        lambda_n = n^2 pi^2 / alpha + alpha / 4 ,

    which are orthogonal with the weight exp(alpha z).  The series converges
    slowly at short times when alpha is large (lambda_n grows only as n^2/alpha);
    alpha should not exceed a few hundred.
    """
    z = np.asarray(z, dtype=float)
    x, w = np.polynomial.legendre.leggauss(n_quad)
    zq, wq = 0.5 * (x + 1), 0.5 * w
    norm = -np.expm1(-alpha)
    out = alpha * np.exp(-alpha * z) / norm
    deviation = (1 - alpha * np.exp(-alpha * zq) / norm) * np.exp(alpha * zq)
    for n in range(1, n_terms + 1):
        k = n * np.pi
        mode_q = np.exp(-alpha * zq / 2) * (np.cos(k * zq) - alpha / (2 * k) * np.sin(k * zq))
        a = np.sum(wq * deviation * mode_q) / np.sum(wq * mode_q**2 * np.exp(alpha * zq))
        out = out + a * np.exp(-(k * k / alpha + alpha / 4) * t) * np.exp(-alpha * z / 2) * (
            np.cos(k * z) - alpha / (2 * k) * np.sin(k * z))
    return out


# ------------------------------------------------------------------ gels

def gel_sphere_roots(G_over_K, n):
    """First n positive roots alpha of  G/M = alpha^2 / (4 - 4 alpha cot alpha),  M = K + 4G/3.

    These fix the relaxation rates  (M/f a^2) alpha_n^2  of a swelling gel sphere of radius a
    with a stress-free surface (Peters and Candau; Li and Tanaka).  For G = 0 they are n pi.
    """
    from scipy.optimize import brentq
    ratio = G_over_K / (1 + 4 * G_over_K / 3)
    def h(a):
        return 4 * ratio * (math.sin(a) - a * math.cos(a)) - a * a * math.sin(a)
    grid = np.linspace(1e-6, (n + 2) * math.pi, 64 * (n + 2))
    vals = np.array([h(a) for a in grid])
    roots = []
    for k in range(len(grid) - 1):
        if vals[k] == 0.0 or vals[k] * vals[k + 1] < 0:
            roots.append(brentq(h, grid[k], grid[k + 1], xtol=1e-14))
            if len(roots) == n:
                break
    return np.array(roots)


def _j1(x):
    x = np.asarray(x, dtype=float)
    small = np.abs(x) < 1e-3
    xs = np.where(small, 1.0, x)
    return np.where(small, x / 3 - x**3 / 30, np.sin(xs) / xs**2 - np.cos(xs) / xs)


def gel_sphere_series(r, t, G_over_K, n_terms=200, n_quad=800):
    """Exact radial displacement u(r, t) of a gel sphere swelling from u = 0 to u = r.

    Linear theory in units of the radius a, the bulk modulus K and the time f a^2 / K, with f
    the friction between network and solvent.  The deviation from the final state is a sum of
    modes j1(alpha_n r) exp(-lambda_n t), lambda_n = (1 + 4G/3K) alpha_n^2, with alpha_n from
    `gel_sphere_roots`; the modes are orthogonal with weight r^2.  G = 0 is the case solved by
    Tanaka and Fillmore.  For an array of times the result has shape (len(t), len(r)).
    """
    alphas = gel_sphere_roots(G_over_K, n_terms)
    rates = (1 + 4 * G_over_K / 3) * alphas**2
    x, w = np.polynomial.legendre.leggauss(n_quad)
    x, w = 0.5 * (x + 1), 0.5 * w
    phi = _j1(alphas[:, None] * x[None, :])
    coeff = -(phi @ (w * x**3)) / (phi**2 @ (w * x**2))
    r = np.asarray(r, dtype=float)
    modes = _j1(alphas[:, None] * np.atleast_1d(r)[None, :])          # (n_terms, len(r))
    times = np.atleast_1d(np.asarray(t, dtype=float))
    u = np.atleast_1d(r)[None, :] + (coeff[None, :] * np.exp(-np.outer(times, rates))) @ modes
    if np.ndim(t) == 0:
        return u[0].reshape(r.shape)
    return u.reshape((len(times),) + r.shape)


def gel_sheet_series(z, t, G_over_K, n_terms=400):
    """Exact displacement u(z, t) of a gel sheet on a rigid, impermeable substrate (z = 0)
    swelling from u = 0 to u = z through its free surface z = 1.

    Units: thickness h, bulk modulus K, time f h^2 / K.  The sheet swells uniaxially, so the
    modulus is M = K + 4G/3 and u obeys du/dt = (M/K) d^2u/dz^2.  The same function describes a
    sheet squeezed by a porous piston, with the sign of u reversed.
    """
    z = np.asarray(z, dtype=float)
    u = z.copy()
    for k in range(n_terms):
        m = 2 * k + 1
        q = m * math.pi / 2
        u = u - (-1) ** k * 8 / (m * math.pi) ** 2 * np.sin(q * z) * math.exp(-(1 + 4 * G_over_K / 3) * q * q * t)
    return u


def gel_fe(n_elements, G_over_K, t_eval, geometry="sphere"):
    """Finite-element solution of the linear swelling of a gel sphere or a gel sheet.

    The Rayleighian is discretised with piecewise-linear displacements: the free energy
    gives the stiffness matrix, the friction f u'^2/2 the mass matrix, and the evolution is
    then exact in time (generalised eigenvectors).  Same units and initial state as
    `gel_sphere_series` (geometry "sphere") and `gel_sheet_series` (geometry "sheet").
    Returns (nodes, u) with u of shape (len(t_eval), n_elements + 1).
    """
    from scipy.linalg import eigh
    x = np.linspace(0.0, 1.0, n_elements + 1)
    n = n_elements + 1
    stiff, mass = np.zeros((n, n)), np.zeros((n, n))
    gq, gw = np.polynomial.legendre.leggauss(4)
    modulus = 1 + 4 * G_over_K / 3
    for e in range(n_elements):
        x0, length = x[e], x[e + 1] - x[e]
        idx = np.ix_([e, e + 1], [e, e + 1])
        for q, wq in zip(gq, gw):
            s = 0.5 * (q + 1)
            r, w = x0 + s * length, 0.5 * wq * length
            phi, dphi = np.array([1 - s, s]), np.array([-1.0, 1.0]) / length
            if geometry == "sphere":
                dil, shear = dphi + 2 * phi / r, dphi - phi / r
                k = np.outer(dil, dil) + (4 / 3) * G_over_K * np.outer(shear, shear)
                weight = r * r
            elif geometry == "sheet":
                k, weight = modulus * np.outer(dphi, dphi), 1.0
            else:
                raise ValueError("geometry must be 'sphere' or 'sheet'")
            stiff[idx] += w * weight * k
            mass[idx] += w * weight * np.outer(phi, phi)
    rates, modes = eigh(stiff[1:, 1:], mass[1:, 1:])
    coeff = modes.T @ mass[1:, 1:] @ (-x[1:])
    u = np.array([x[1:] + modes @ (coeff * np.exp(-rates * t)) for t in np.atleast_1d(t_eval)])
    return x, np.concatenate((np.zeros((u.shape[0], 1)), u), axis=1)


# ------------------------------------------------------------------ phase separation

def cahn_hilliard_spectral(phi0, box, t_eval, dfdphi, dt, stabilisation=2.0):
    """Pseudo-spectral solution of the Cahn-Hilliard equation on a periodic interval or square,

        dphi/dt = laplacian( f'(phi) - laplacian(phi) ),

    in units in which the gradient coefficient and the mobility are 1.  phi0 is a 1D array or a
    square 2D array on a uniform grid of side `box`; dfdphi(phi) is f'(phi).  The time stepping
    is second order (semi-implicit BDF2): the fourth-order term is implicit, f' is extrapolated,
    and a stabilising term stabilisation * laplacian(phi^{n+1} - phi_extrapolated) allows large
    steps.  The mean of phi is conserved exactly.  The times in t_eval must be multiples of dt.

    Returns an array of snapshots, shape (len(t_eval),) + phi0.shape.
    """
    phi = np.array(phi0, dtype=float)
    n = phi.shape[0]
    k = 2 * np.pi * np.fft.fftfreq(n, d=box / n)
    if phi.ndim == 1:
        q2 = k**2
    elif phi.ndim == 2 and phi.shape[1] == n:
        kx, ky = np.meshgrid(k, k, indexing="ij")
        q2 = kx**2 + ky**2
    else:
        raise ValueError("phi0 must be 1D or a square 2D array")
    t_eval = np.asarray(t_eval, dtype=float)
    targets = np.rint(t_eval / dt).astype(int)
    if np.any(np.abs(targets * dt - t_eval) > 1e-9 * max(1.0, t_eval.max())):
        raise ValueError("the times in t_eval must be multiples of dt")
    S = stabilisation
    hat, nonlin = np.fft.fftn(phi), np.fft.fftn(dfdphi(phi))
    hat_old = nonlin_old = None
    out, idx = [], 0
    while idx < len(targets) and targets[idx] == 0:
        out.append(phi.copy())
        idx += 1
    for step in range(1, targets[-1] + 1):
        if hat_old is None:                                    # first step: semi-implicit Euler
            new = (hat - dt * q2 * nonlin + dt * S * q2 * hat) / (1 + dt * q2**2 + dt * S * q2)
        else:
            guess = 2 * hat - hat_old
            new = (4 * hat - hat_old - 2 * dt * q2 * (2 * nonlin - nonlin_old) + 2 * dt * S * q2 * guess) \
                / (3 + 2 * dt * q2**2 + 2 * dt * S * q2)
        hat_old, nonlin_old = hat, nonlin
        hat = new
        phi = np.real(np.fft.ifftn(hat))
        nonlin = np.fft.fftn(dfdphi(phi))
        while idx < len(targets) and step == targets[idx]:
            out.append(phi.copy())
            idx += 1
    return np.array(out)


# ------------------------------------------------------------------ polymer chains

def rouse_brownian_dynamics(n_beads, n_chains, dt, t_end, stride, rng, spring=3.0):
    """Brownian dynamics of free Rouse chains, started from equilibrium.

    Units: kT = 1, bead friction 1, spring constant `spring` (3 for a mean-square bond length 1).
    Each bead obeys dR = F dt + sqrt(2 dt) xi (Euler-Maruyama), with the spring forces F and a
    random force whose strength is fixed by the fluctuation-dissipation theorem.

    Returns (times, end_to_end, middle): the end-to-end vectors, shape (n_times, n_chains, 3), and
    the positions of the middle bead, same shape, recorded every `stride` steps.
    """
    bonds = rng.standard_normal((n_chains, n_beads - 1, 3)) / math.sqrt(spring)
    R = np.concatenate((np.zeros((n_chains, 1, 3)), np.cumsum(bonds, axis=1)), axis=1)
    n_steps = int(round(t_end / dt))
    noise = math.sqrt(2 * dt)
    times, ends, middle = [], [], []
    for step in range(n_steps + 1):
        if step % stride == 0:
            times.append(step * dt)
            ends.append(R[:, -1] - R[:, 0])
            middle.append(R[:, n_beads // 2].copy())
        tension = spring * np.diff(R, axis=1)
        F = np.zeros_like(R)
        F[:, :-1] += tension
        F[:, 1:] -= tension
        R += F * dt + noise * rng.standard_normal(R.shape)
    return np.array(times), np.array(ends), np.array(middle)


def fene_chain_steady(n_beads, rate, q_max, spring=3.0):
    """Stretched steady state of a bead-spring chain with FENE springs in an elongational flow.

    Friction 1 per bead, strain rate `rate` along x, springs f(s) = spring * s / (1 - s^2/q_max^2).
    The mean positions satisfy rate * x_n + T_n - T_{n-1} = 0, with T_n the tension of bond n and
    no tension beyond the ends; the tension is marched from one end and the end position is found
    so that the chain is centred on the stagnation point.  Returns the bead positions x_n, or
    None if only the unstretched state exists.
    """
    from scipy.optimize import brentq
    def stretch(T):
        if abs(T) < 1e-14:
            return T / spring
        return (-spring + math.sqrt(spring**2 + 4 * T * T / q_max**2)) / (2 * T / q_max**2)
    def march(x0):
        x, T = [x0], 0.0
        for _ in range(n_beads - 1):
            T -= rate * x[-1]
            x.append(x[-1] + stretch(T))
        return np.array(x)
    def mismatch(x0):
        return march(x0)[-1] + x0
    grid = -np.linspace(1e-3, 0.4999 * (n_beads - 1) * q_max, 4000)[::-1]
    values = [mismatch(v) for v in grid]
    for i in range(len(grid) - 1):
        if values[i] * values[i + 1] < 0:
            return march(brentq(mismatch, grid[i], grid[i + 1], xtol=1e-13))
    return None


# ------------------------------------------------------------------ nematic liquid crystals

def director_fd(theta0, t_eval, torque, diffusivity=1.0, rtol=1e-9, atol=1e-12):
    """Finite-difference solution of  d(theta)/dt = diffusivity * theta_zz + torque(theta)  on 0 < z < 1,
    with theta = 0 at both walls (strong anchoring).

    theta0 holds the values at the n interior points z_i = i/(n + 1); torque(theta) is the
    field or flow torque divided by the rotational viscosity.  Returns the snapshots at t_eval,
    shape (len(t_eval), n).
    """
    from scipy.sparse import diags
    theta0 = np.asarray(theta0, dtype=float)
    n = len(theta0)
    dz2 = (1.0 / (n + 1)) ** 2
    def rhs(t, th):
        lap = np.empty_like(th)
        lap[1:-1] = th[2:] - 2 * th[1:-1] + th[:-2]
        lap[0] = th[1] - 2 * th[0]
        lap[-1] = th[-2] - 2 * th[-1]
        return diffusivity * lap / dz2 + torque(th)
    pattern = diags([np.ones(n - 1), np.ones(n), np.ones(n - 1)], [-1, 0, 1])
    sol = solve_ivp(rhs, (t_eval[0], t_eval[-1]), theta0, t_eval=t_eval, method="BDF",
                    rtol=rtol, atol=atol, jac_sparsity=pattern)
    return sol.y.T


def freedericksz_twist_angle(h):
    """Exact twist angle at the mid-plane of a twist cell above the Freedericksz threshold.

    h is the field over its threshold value; the angle phi_m solves (2/pi) K(sin phi_m) = h, with K
    the complete elliptic integral of the first kind.  Returns 0 for h <= 1.
    """
    from scipy.optimize import brentq
    from scipy.special import ellipkm1
    if h <= 1:
        return 0.0
    # solve for p = 1 - m = cos^2(phi_m) on a log scale: for strong fields p is far below 1e-16
    log_p = brentq(lambda q: 2 / math.pi * ellipkm1(math.exp(q)) - h, -700.0, 0.0, xtol=1e-14)
    return math.acos(math.sqrt(math.exp(log_p)))
