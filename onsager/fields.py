"""Reducing a one-dimensional diffusion problem to a few parameters.

Particles of number density c(x, t), each with friction constant zeta against the
solvent, move with velocity v(x, t); the flux is J = c v and the dissipation
function is

    Phi = (zeta/2) * integral of J**2 / c dx .

If the profile is assumed to depend on time only through parameters a_i(t), then
dc/dt = sum_i (dc/da_i) adot_i, conservation gives J = sum_i J_i adot_i, and
Phi = (1/2) sum_ij zeta_ij adot_i adot_j with

    zeta_ij = zeta * integral of J_i J_j / c dx .
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid

__all__ = ["flux_from_rate", "friction_matrix_1d"]


def flux_from_rate(x, dcdt, closed="right"):
    """Flux J(x) that satisfies dc/dt = -dJ/dx on the grid x.

    closed="right": J vanishes at the last grid point, J(x) = integral from x to
    the right end of dc/dt.  closed="left": J vanishes at the first grid point,
    J(x) = -integral from the left end to x of dc/dt.
    """
    dcdt = np.asarray(dcdt, dtype=float)
    cum = cumulative_trapezoid(dcdt, x, initial=0.0)
    if closed == "right":
        return cum[-1] - cum
    if closed == "left":
        return -cum
    raise ValueError("closed must be 'left' or 'right'")


def friction_matrix_1d(x, c, dc_da, zeta=1.0, closed="right", floor=1e-12):
    """Friction matrix of the parameters of a trial profile.

    x : grid;  c : profile on the grid;  dc_da : list of arrays dc/da_i on the grid.
    Points where c is below floor times its maximum carry no particles to speak of
    and are left out of the integral (there J/c would be rounding error divided by
    rounding error).
    """
    c = np.asarray(c, dtype=float)
    fluxes = [flux_from_rate(x, d, closed) for d in dc_da]
    occupied = c > floor * c.max()
    n = len(fluxes)
    Z = np.empty((n, n))
    for i in range(n):
        for j in range(i, n):
            integrand = np.zeros_like(c)
            integrand[occupied] = fluxes[i][occupied] * fluxes[j][occupied] / c[occupied]
            Z[i, j] = Z[j, i] = zeta * np.trapezoid(integrand, x)
    return Z
