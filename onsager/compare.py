"""Comparing a variational result with a reference solution."""

from __future__ import annotations

import numpy as np

__all__ = ["sup_error", "relative_error", "loglog_slope"]


def sup_error(approx, exact):
    """Largest absolute difference between two arrays: the uniform error on a grid."""
    return float(np.max(np.abs(np.asarray(approx) - np.asarray(exact))))


def relative_error(approx, exact):
    """(approx - exact) / exact, for numbers or arrays."""
    return np.asarray(approx, dtype=float) / np.asarray(exact, dtype=float) - 1.0


def loglog_slope(x, y, fit_range=None):
    """Least-squares slope of log(y) against log(x).

    Returns (slope, prefactor): y is approximately prefactor * x**slope over the
    fitted range (xmin, xmax), or over all points if fit_range is None.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = (y > 0) & (x > 0) & np.isfinite(y)
    if fit_range is not None:
        keep &= (x >= fit_range[0]) & (x <= fit_range[1])
    if keep.sum() < 2:
        raise ValueError("need at least two positive values in the fitted range")
    slope, intercept = np.polyfit(np.log(x[keep]), np.log(y[keep]), 1)
    return slope, float(np.exp(intercept))
