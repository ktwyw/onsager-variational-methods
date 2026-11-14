"""Reference solutions against which the variational results are checked."""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags

__all__ = ["smoluchowski_fv", "duct_flow_rate", "duct_velocity", "thin_film_fv"]


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
