from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from dmhd.spectral import diffspec

Ns, errs = [], []
for N in range(4, 50, 2):
    h = 2 * np.pi / N
    x = h * np.arange(N)
    v = np.exp(np.sin(x))
    w = diffspec(v)
    err = np.max(np.abs(w - np.cos(x) * v))
    Ns.append(N)
    errs.append(err)

fig, ax = plt.subplots(figsize=(6, 4))
ax.semilogy(Ns, errs, "o-", label=r"$\max_j |Du - u'|$")
ax.axhline(np.finfo(float).eps, ls="--", c="gray", label="machine epsilon")
ax.set_xlabel("N (grid points)")
ax.set_ylabel("maximum error")
ax.set_title(r"Spectral derivative of $e^{\sin x}$")
ax.legend()

out = Path(__file__).resolve().parents[1] / "figures" / "convergence_exp_sin.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print(f"saved {out}")