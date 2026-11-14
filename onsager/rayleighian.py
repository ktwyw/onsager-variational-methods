"""The Onsager principle for a finite number of variables, in SymPy.

A system is described by slow variables x = (x1, ..., xn), a free energy A(x)
and a dissipation function Phi(x, xdot) that is quadratic in the rates,

    Phi = (1/2) sum_ij zeta_ij(x) xdot_i xdot_j .

The Rayleighian is  R = Phi + sum_i (dA/dx_i) xdot_i .  The rates that the system
chooses are those that minimise R; if the rates are tied by constraints that are
linear in them, R is minimised subject to the constraints, with one Lagrange
multiplier for each.
"""

from __future__ import annotations

import sympy as sp
from sympy.solvers.solveset import NonlinearError

__all__ = ["rayleighian", "friction_matrix", "kinetic_equations", "rate_function"]


def rayleighian(A, Phi, x, xdot):
    """R = Phi + dA/dt, with dA/dt = sum_i (dA/dx_i) xdot_i."""
    return Phi + sum(sp.diff(A, xi) * vi for xi, vi in zip(x, xdot))


def friction_matrix(Phi, xdot):
    """The matrix zeta_ij = d^2 Phi / d xdot_i d xdot_j (symmetric by construction)."""
    return sp.hessian(Phi, list(xdot))


def kinetic_equations(A, Phi, x, xdot, constraints=(), multipliers=()):
    """Minimise the Rayleighian with respect to the rates.

    Parameters
    ----------
    A, Phi : SymPy expressions for the free energy and the dissipation function.
    x, xdot : sequences of symbols, the variables and their rates.  A rate that
        has no variable of its own (a flow rate, say) can be listed in xdot with
        a dummy entry in x on which A does not depend.
    constraints : expressions g_k(x, xdot), linear in the rates, that must vanish.
    multipliers : one symbol per constraint.  The quantity made stationary is
        R - sum_k multiplier_k * g_k.

    Returns
    -------
    dict mapping each rate (and each multiplier) to its expression in x.

    The function does not check that Phi is positive definite.  If it is not, the
    stationary point returned is not a minimum and the model is unphysical.
    """
    x, xdot = list(x), list(xdot)
    constraints, multipliers = list(constraints), list(multipliers)
    if len(constraints) != len(multipliers):
        raise ValueError("one multiplier is needed for each constraint")
    R = rayleighian(A, Phi, x, xdot) - sum(m * g for m, g in zip(multipliers, constraints))
    unknowns = xdot + multipliers
    equations = [sp.diff(R, v) for v in xdot] + constraints
    try:
        M, rhs = sp.linear_eq_to_matrix(equations, unknowns)
    except NonlinearError as exc:
        raise ValueError("Phi must be quadratic, and the constraints linear, in the rates") from exc
    solution = M.LUsolve(rhs)
    return {v: sp.simplify(s) for v, s in zip(unknowns, solution)}


def rate_function(equations, x, xdot, parameters=None):
    """Turn the result of kinetic_equations into f(t, x) for scipy's solve_ivp.

    parameters is a dict {symbol: number} for the constants in the equations.
    """
    parameters = parameters or {}
    exprs = [equations[v].subs(parameters) for v in xdot]
    f = sp.lambdify(list(x), exprs, "numpy")
    return lambda t, y: f(*y)
