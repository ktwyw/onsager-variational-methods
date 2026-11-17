"""Animation at the top of README.md: a spreading drop (Chapter 7).

A drop that wets the solid completely spreads on a precursor film.  The profile is
computed from the thin-film equation (on.thin_film_fv); the variational result,
Eq. (7.12), uses only the radius of a spherical cap.  Left: the drop (thickness
exaggerated) with the fitted cap.  Right: the radius against time, thin-film equation
(dots) and Eq. (7.12) (line), approaching Tanner's law R ~ t^(1/10).
Output: docs/images/spreading.gif.

Run from the repository root:  python scripts/make_readme_drop_animation.py   (needs Pillow)
"""
import io
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.integrate import solve_ivp

try:
    import onsager as on
except ImportError:          # running from a clone without `pip install -e .`
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    import onsager as on

on.use_style()

b = 3e-3
edges = np.linspace(0.0, 3.2, 1601)
r = 0.5 * (edges[1:] + edges[:-1])
times = np.concatenate(([0.0], np.logspace(-2, 3, 46)))
r, h = on.thin_film_fv(np.maximum(1 - r**2, 0.0) + b, edges, times)
fits = np.array([on.fit_cap(r, hk, b=b) for hk in h])
start = 21                                       # t = 1.7: start the kinetic equation from the computed radius

def rate(t, R):
    theta = 2 / R[0] ** 3
    return [theta**3 / (2 * (math.log(1 / (R[0] ** 2 * b)) - 1))]

t_var = np.logspace(math.log10(times[start]), 3, 200); t_var[0], t_var[-1] = times[start], 1e3
R_var = solve_ivp(rate, (times[start], 1e3), [fits[start, 1]], t_eval=t_var, rtol=1e-10).y[0]

frames = []
for k in range(1, len(times)):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.7), gridspec_kw=dict(width_ratios=[1.25, 1]))
    rr = np.concatenate((-r[::-1], r))
    hh = np.concatenate((h[k][::-1], h[k]))
    ax1.fill_between(rr, 0, hh, color="C0", alpha=0.3)
    ax1.plot(rr, hh, "C0", lw=2, label="thin-film equation")
    hc, R = fits[k]
    xc = np.linspace(-R, R, 200)
    ax1.plot(xc, hc * (1 - xc**2 / R**2) + b, "k--", lw=1, label="spherical cap")
    ax1.plot([-3.2, 3.2], [0, 0], "k", lw=2)
    ax1.set_xlim(-3.0, 3.0); ax1.set_ylim(-0.02, 1.05)
    ax1.set_xlabel("$r/R_0$"); ax1.set_ylabel("$h/H_0$ (exaggerated)")
    label = f"{times[k]:.0f}" if times[k] >= 10 else f"{times[k]:.2g}"
    ax1.text(-2.85, 0.92, f"t = {label}", fontsize=11)
    ax1.legend(loc="upper right", fontsize=8.5)

    ax2.loglog(times[1:k + 1], fits[1:k + 1, 1], "o", color="C0", ms=4, label="thin-film equation")
    ax2.loglog(t_var, R_var, "C3", lw=1.6, label="one variable, Eq. (7.12)")
    tt = np.array([30.0, 1000.0])
    ax2.loglog(tt, 1.15 * tt**0.1, "k:", lw=1.2); ax2.text(80, 1.45, r"slope 1/10", fontsize=10)
    ax2.set_xlim(0.008, 1300); ax2.set_ylim(0.95, 3.2)
    ax2.set_yticks([1, 2, 3]); ax2.set_yticklabels(["1", "2", "3"]); ax2.minorticks_off()
    ax2.set_xlabel(r"time $t\,\gamma H_0^3/3\eta R_0^4$"); ax2.set_ylabel("radius $R/R_0$")
    ax2.legend(loc="upper left", fontsize=8.5)
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=80)
    plt.close(fig)
    buffer.seek(0)
    frames.append(Image.open(buffer).convert("P", palette=Image.ADAPTIVE, colors=64))

frames += [frames[-1]] * 12                      # hold the last frame
out = os.path.join(os.path.dirname(__file__), "..", "docs", "images", "spreading.gif")
frames[0].save(out, save_all=True, append_images=frames[1:], duration=110, loop=0, optimize=True)
print(f"{len(frames)} frames -> {os.path.normpath(out)} ({os.path.getsize(out) / 1e3:.0f} kB)")
