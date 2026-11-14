"""onsager: small, readable tools for the Onsager variational principle.

The module accompanies the notebooks of onsager-variational-methods.  It provides

* rayleighian:  the Rayleighian of a system with a few variables and the kinetic
                equations that minimise it, with constraints, in SymPy;
* fields:       the friction matrix of a trial concentration profile;
* reference:    reference solutions used for the checks (conservative finite-volume
                solvers for diffusion and for the thin-film equation; flow in a
                rectangular duct);
* compare:      uniform and relative errors, slopes on log-log plots, and the
                apparent radius of a drop;

plus the plotting style used by every figure.
"""

from .compare import fit_cap, loglog_slope, relative_error, sup_error
from .fields import flux_from_rate, friction_matrix_1d
from .plotting import panel_label, use_style
from .rayleighian import friction_matrix, kinetic_equations, rate_function, rayleighian
from .reference import duct_flow_rate, duct_velocity, smoluchowski_fv, thin_film_fv

__version__ = "0.1.0.dev0"

__all__ = [
    "rayleighian", "friction_matrix", "kinetic_equations", "rate_function",
    "flux_from_rate", "friction_matrix_1d",
    "smoluchowski_fv", "duct_flow_rate", "duct_velocity", "thin_film_fv",
    "sup_error", "relative_error", "loglog_slope", "fit_cap",
    "use_style", "panel_label",
]
