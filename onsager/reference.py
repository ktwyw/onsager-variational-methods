"""Reference solutions against which the variational results are checked."""

from __future__ import annotations

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags

__all__ = ["smoluchowski_fv", "duct_flow_rate", "duct_velocity"]


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
