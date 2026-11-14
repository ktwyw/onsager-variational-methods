"""Animation at the top of README.md: a trial solution against the exact one.

Diffusion into a half-space from a wall held at concentration c0 (Chapter 5).  The
exact profile is c0 erfc(x / 2 sqrt(D t)).  The trial profiles are c0 (1 - x/l)**n,
with the penetration length l(t) given by the Onsager principle: l = lambda sqrt(D t),
lambda**2 = 2 g(n) / I(n).  Left: the profiles as time goes on.  Right: the amount
taken up, accumulated in time.  Output: docs/images/penetration.gif.

Run from the repository root:  python scripts/make_readme_animation.py   (needs Pillow)
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
from scipy.special import erfc

try:
    import onsager as on
except ImportError:          # running from a clone without `pip install -e .`
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    import onsager as on

on.use_style()

I = lambda n: (5 * n + 6) / ((n + 1) ** 2 * (n + 2) * (n + 3))
g = lambda n: (2 * n + 1) / (n + 1) ** 2
lam = lambda n: math.sqrt(2 * g(n) / I(n))
trials = ((1, "C0", "straight line, $n = 1$"), (3, "C3", "$n = 3$"))

x = np.linspace(0, 6, 500)
times = np.linspace(0.05, 2.0, 48) ** 2
frames = []
for k, t in enumerate(times):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.9))
    ax1.fill_between(x, erfc(x / (2 * math.sqrt(t))), color="0.88")
    ax1.plot(x, erfc(x / (2 * math.sqrt(t))), "k-", lw=2, label="exact")
    for n, c, lab in trials:
        ell = lam(n) * math.sqrt(t)
        ax1.plot(x, np.where(x < ell, (1 - np.minimum(x / ell, 1)) ** n, 0), "--", color=c, label=lab)
    ax1.set_xlim(0, 6); ax1.set_ylim(0, 1.05)
    ax1.set_xlabel("distance from the wall, $x$"); ax1.set_ylabel("concentration $c/c_0$")
    ax1.legend(loc="upper right"); ax1.set_title(f"$Dt = {t:4.2f}$", loc="left")

    tt = times[: k + 1]
    ax2.plot(tt, 2 / math.sqrt(math.pi) * np.sqrt(tt), "k-", lw=2, label=r"exact, $1.128\sqrt{Dt}$")
    for n, c, lab in trials:
        ax2.plot(tt, lam(n) / (n + 1) * np.sqrt(tt), "--", color=c, label=rf"${lam(n) / (n + 1):.3f}\sqrt{{Dt}}$")
    ax2.set_xlim(0, 4); ax2.set_ylim(0, 2.7)
    ax2.set_xlabel("time $Dt$"); ax2.set_ylabel("amount taken up, $M/c_0$")
    ax2.legend(loc="lower right"); ax2.set_title("one parameter, $\\ell(t)$, from the Onsager principle", loc="left")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=90)
    plt.close(fig)
    buf.seek(0)
    frames.append(Image.open(buf).convert("P", palette=Image.ADAPTIVE, colors=128))

out = os.path.join(os.path.dirname(__file__), "..", "docs", "images", "penetration.gif")
frames[0].save(out, save_all=True, append_images=frames[1:] + [frames[-1]] * 14, duration=110, loop=0, optimize=True)
print("wrote", os.path.normpath(out), f"({os.path.getsize(out) / 1e6:.2f} MB)")
