import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from dmhd.solvers.burgers import solve_burgers, u_exact
from matplotlib.ticker import NullFormatter

n = 96
nu = 1
a = 0.5
T = 1
x = 2*np.pi/n*np.arange(n)
steps = np.array([4, 8, 16, 32, 64])
u0 = u_exact(x, 0, nu, a)
errs = np.zeros(len(steps))  # max |u - u_exact(x, T)| per step count
for i, step in enumerate(steps):
    u = solve_burgers(u0, nu, T, step)
    errs[i] = np.max(np.abs(u-u_exact(x, T, nu, a)))


dts = T / steps
orders = np.log2(errs[:-1] / errs[1:])
print("errors:", errs)
print("observed orders:", orders)

fig, ax = plt.subplots()
ax.loglog(dts, errs, "o-", label="IF-RK4")
ax.loglog(dts, errs[-1] * (dts / dts[-1]) ** 4, "k--", label="slope 4")
ax.set_xlabel("dt")
ax.set_xticks(dts)
ax.set_xticklabels([f"T/{s}" for s in steps])
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlabel("dt = T / n_steps")
ax.set_ylabel("max error at T")
ax.set_title(f"Burgers dt-convergence (n={n}, nu={nu}, a={a}, T={T})")
ax.legend()
out = Path(__file__).resolve().parent.parent / "figures" / "burgers_dt_convergence.png"
fig.savefig(out, dpi=150, bbox_inches="tight")

slope = np.polyfit(np.log(dts[-3:]), np.log(errs[-3:]), 1)[0]
print("fitted slope (last 3 points):", slope)