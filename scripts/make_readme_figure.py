"""Overview figure of README.md: four results from Part I.  Output: docs/images/overview.png.

Run from the repository root:  python scripts/make_readme_figure.py
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erfc

try:
    import onsager as on
except ImportError:          # running from a clone without `pip install -e .`
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    import onsager as on

on.use_style()
fig, axes = plt.subplots(2, 2, figsize=(9.4, 7.6))

# (a) the Rayleighian of a particle (Chapter 1)
ax = axes[0, 0]
dU = lambda x: x**3 - x - 0.2
v = np.linspace(-1.2, 1.2, 200)
for xp, c in ((-0.3, "C0"), (0.4, "C1"), (0.9, "C3")):
    ax.plot(v, 0.5 * v**2 + dU(xp) * v, color=c, label=f"$x = {xp}$")
    ax.plot(-dU(xp), -0.5 * dU(xp) ** 2, "o", color=c)
ax.axhline(0, color="k", lw=0.6); ax.set_ylim(-0.25, 0.8)
ax.set_xlabel(r"velocity $\dot x$"); ax.set_ylabel(r"Rayleighian $R = \Phi + \dot A$")
ax.legend(loc="upper center"); on.panel_label(ax, "(a)")

# (b) flow in a square duct: exact and one trial function (Chapter 1)
ax = axes[0, 1]
gr = np.linspace(-1, 1, 161)
GX, GY = np.meshgrid(gr, gr)
levels = np.linspace(0.03, 0.27, 5)
ax.contour(GX, GY, on.duct_velocity(GX, GY), levels=levels, colors="k", linewidths=1.3)
ax.contour(GX, GY, 5 / 16 * (1 - GX**2) * (1 - GY**2), levels=levels, colors="C3", linewidths=1.3, linestyles="--")
ax.plot([], [], "k-", label="exact"); ax.plot([], [], "--", color="C3", label="one trial function")
ax.set_aspect("equal"); ax.set_xlabel("$x/a$"); ax.set_ylabel("$y/a$")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.13), ncol=2); on.panel_label(ax, "(b)", x=0.46, y=0.54)

# (c) drift of a settling rod (Chapter 2)
ax = axes[1, 0]
th = np.linspace(0, np.pi / 2, 400)
for ratio, c, lab in ((2.0, "C3", r"$\zeta_\perp/\zeta_\parallel = 2$ (slender rod)"), (1.5, "C0", r"$\zeta_\perp/\zeta_\parallel = 1.5$")):
    ax.plot(np.degrees(th), np.degrees(np.arctan((ratio - 1) * np.sin(th) * np.cos(th) / (ratio * np.cos(th) ** 2 + np.sin(th) ** 2))),
            color=c, label=lab)
ax.axhline(np.degrees(math.asin(1 / 3)), color="k", ls=":", lw=1)
ax.text(2, 20.1, r"$\arcsin\frac{1}{3} = 19.5^\circ$", fontsize=9.5)
ax.set_xlim(0, 90); ax.set_ylim(0, 24)
ax.set_xlabel(r"inclination of the rod, $\theta$ (degrees)"); ax.set_ylabel(r"drift angle $\alpha$ (degrees)")
ax.legend(loc="lower center"); on.panel_label(ax, "(c)", x=0.9)

# (d) diffusion into a half-space: trial profiles (Chapter 5)
ax = axes[1, 1]
I = lambda n: (5 * n + 6) / ((n + 1) ** 2 * (n + 2) * (n + 3))
g = lambda n: (2 * n + 1) / (n + 1) ** 2
xp = np.linspace(0, 6, 600)
ax.plot(xp, erfc(xp / 2), "k-", lw=2, label="exact, erfc")
for n, c in ((1, "C0"), (2, "C2"), (3, "C3")):
    lam = math.sqrt(2 * g(n) / I(n))
    ax.plot(xp, np.where(xp < lam, (1 - np.minimum(xp / lam, 1)) ** n, 0), "--", color=c, lw=1.4,
            label=rf"$n = {n}$: uptake {100 * (lam / (n + 1) * math.sqrt(math.pi) / 2 - 1):+.1f} %")
ax.set_xlabel(r"$x/\sqrt{Dt}$"); ax.set_ylabel("$c/c_0$")
ax.legend(loc="upper right"); on.panel_label(ax, "(d)", x=0.04, y=0.12)

out = os.path.join(os.path.dirname(__file__), "..", "docs", "images", "overview.png")
fig.savefig(out, dpi=160)
print("wrote", os.path.normpath(out))
