from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from dmhd.spectral import diffspec
from dmhd.timestep import rk4_step

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


T = 1.0
LAMBDAS = -np.arange(1, 6)
NS = np.arange(10, 160, 10)
HS = T / NS


def integrate(f, y0, n):
    """n RK4 steps from t=0 to t=T."""
    h = T / n
    y = y0
    for i in range(n):
        y = rk4_step(f, i * h, y, h)
    return y


# Test 1: y' = lam*y, y(0)=1  ->  y(T) = exp(lam*T); relative error
err_exp = np.zeros((len(LAMBDAS), len(NS)))
for a, lam in enumerate(LAMBDAS):
    exact = np.exp(lam * T)
    for b, n in enumerate(NS):
        y = integrate(lambda t, y, lam=lam: lam * y, 1.0, n)
        err_exp[a, b] = abs(y - exact) / exact

# Test 2: y' = t^3, y(0)=0  ->  y(T) = T^4/4; RK4 is exact for cubics in t
err_t3 = np.array([abs(integrate(lambda t, y: t**3, 0.0, n) - T**4 / 4) for n in NS])

slopes = [np.polyfit(np.log(HS), np.log(e), 1)[0] for e in err_exp]
for lam, s in zip(LAMBDAS, slopes):
    print(f"lambda = {lam:+d}: fitted order {s:.3f}")
print(f"y'=t^3 error: max {err_t3.max():.2e}")

# ---- figure ----
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
blues = plt.cm.Blues(np.linspace(0.40, 0.95, len(LAMBDAS)))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)

for color, lam, e, s in zip(blues, LAMBDAS, err_exp, slopes):
    ax1.loglog(HS, e, "o-", color=color, lw=1.5, ms=4)
    ax1.annotate(rf"$\lambda={lam}$  (slope {s:.2f})", (HS[0], e[0]), xytext=(6, 0),
                 textcoords="offset points", fontsize=8, color=MUTED, va="center")

ref = err_exp[0, -1] * (HS / HS[-1]) ** 4 * 0.2
ax1.loglog(HS, ref, "--", color="gray", lw=1, label=r"slope 4 ($\propto h^4$)")
ax1.set_xlabel("step size $h = T/n$")
ax1.set_ylabel(r"relative error at $T=1$")
ax1.set_title(r"RK4 on $y'=\lambda y$: fourth-order convergence", fontsize=10, color=INK)
ax1.legend(loc="lower right", fontsize=8, frameon=False)
ax1.set_xlim(right=HS[0] * 3)

nonzero = err_t3 > 0  # exact zeros cannot be drawn on a log axis
ax2.semilogy(NS[nonzero], err_t3[nonzero], "o", color=blues[-2], ms=5,
             label=r"$|y_n - T^4/4|$ (exact zeros not shown)")
ax2.axhline(np.finfo(float).eps, ls="--", c="gray", lw=1, label="machine epsilon")
ax2.set_xlabel("n (steps)")
ax2.set_ylabel("absolute error at $T=1$")
ax2.set_title(r"RK4 on $y'=t^3$: exact up to roundoff", fontsize=10, color=INK)
ax2.set_ylim(1e-18, 1e-14)
ax2.legend(loc="lower right", fontsize=8, frameon=False)

for ax in (ax1, ax2):
    ax.grid(True, which="major", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

out = Path(__file__).resolve().parents[1] / "figures" / "rk4_convergence.png"
fig.savefig(out, dpi=200, bbox_inches="tight")
print(f"saved {out}")