"""Refit the attached-biomass profile Cb(x) using colloid filtration theory
(CFT), calibrated against the real 0-52 cm data, instead of the raw cubic
polynomial currently used in MICP_Simulation_T3B.py.

Where the "real data" comes from: Razbani et al. (2025) fit a third-degree
polynomial to their measured attached-biomass density (their Fig. 2, itself
derived from Martinez et al. 2013's unattached-bacteria profile scaled to
total retained cells). We don't have their raw digitized points, but we DO
have the polynomial that reproduces them at the 10 real cell centers - so
those 10 (x, Cb) pairs, evaluated from the cubic at the actual calibrated
cell positions, are used here as the calibration target. This is the same
real information Razbani et al. fit to, just re-fit with a mechanistic
functional form instead of an unconstrained polynomial.

CRITICAL DIRECTION NOTE: biomass was injected TOP-DOWN (Razbani et al.,
citing Martinez et al. 2013), while cementation solution (what TRANSPORT
models) flows BOTTOM-UP. Cell 1 = bottom = cementation inlet; the top of
the 52 cm column (x=52) is where biomass was actually introduced. So the
correct filtration distance is measured FROM THE TOP, not from cell 1 -
this is why Cb *increases* with cell number even though filtration theory
predicts *decay* away from an injection point: cell 10 (x=49.4) is close
to the biomass injection point (distance ~2.6 cm), and cell 1 (x=2.6) is
far from it (distance ~49.4 cm).

A single-site exponential decay from that top point fits poorly (R^2=0.58)
- the real profile's growth ratio between cells is not constant, which
rules out a clean single-rate exponential. A two-site (bi-exponential)
model - a standard CFT extension for exactly this kind of profile
(Bradford et al.; Tufenkji & Elimelech) - fits well (R^2=0.977).
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

CELL_L = 5.2
N_CELLS = 10
COLUMN_LEN = N_CELLS * CELL_L  # 52 cm

x = np.array([i * CELL_L - CELL_L / 2 for i in range(1, N_CELLS + 1)])  # from bottom (cementation inlet)
x_top = COLUMN_LEN - x                                                  # from top (biomass injection point)


def cb_cubic(x):
    """The original Razbani et al. (2025) polynomial fit, evaluated at the
    real calibrated cell centers - our stand-in for their raw data."""
    return (5.66e-1 * x ** 3 + 1.1e2 * x ** 2 + 4.67e3 * x - 5.62e3) / 1e6


Cb_data = cb_cubic(x)


def r2(y, pred):
    ss_res = np.sum((y - pred) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    return 1 - ss_res / ss_tot


# --- single-site CFT: Cb = Cb0 * exp(-lambda * x_top) ---
A = np.vstack([x_top, np.ones_like(x_top)]).T
slope, intercept = np.linalg.lstsq(A, np.log(Cb_data), rcond=None)[0]
lam_single, Cb0_single = -slope, np.exp(intercept)
pred_single = Cb0_single * np.exp(-lam_single * x_top)
r2_single = r2(Cb_data, pred_single)


# --- two-site CFT: Cb = A1*exp(-l1*x_top) + A2*exp(-l2*x_top) ---
def biexp(xt, A1, l1, A2, l2):
    return A1 * np.exp(-l1 * xt) + A2 * np.exp(-l2 * xt)


popt, _ = curve_fit(biexp, x_top, Cb_data, p0=[0.7, 0.05, 0.3, 15], maxfev=20000)
A1, l1, A2, l2 = popt
pred_biexp = biexp(x_top, *popt)
r2_biexp = r2(Cb_data, pred_biexp)

print("Calibration points (from the real 0-52cm cubic, at the 10 actual cell centers):")
for xi, xti, cbi in zip(x, x_top, Cb_data):
    print(f"  x={xi:5.1f} cm (from bottom)  x_top={xti:5.1f} cm (from biomass injection)  Cb={cbi:.4f}")

print(f"\nSingle-site CFT:  Cb0={Cb0_single:.4f}  lambda={lam_single:.5f} /cm  "
      f"(filtration length {1/lam_single:.1f} cm)  R^2={r2_single:.4f}  <- poor fit")
print(f"Two-site CFT:     A1={A1:.4f} lambda1={l1:.5f} /cm (slow site, length {1/l1:.1f} cm)")
print(f"                  A2={A2:.4f} lambda2={l2:.5f} /cm (fast site, length {1/l2:.2f} cm)")
print(f"                  R^2={r2_biexp:.4f}  <- good fit, adopted")

# ---------------------------------------------------------------------------
# What this means for extending the column - THIS IS NOT A FREE CHOICE.
# x_top = 0 is a physical location: the actual top of Martinez et al.'s real
# 52 cm column, where biomass injection stopped. Two different "extensions"
# are physically different scenarios, not two guesses at the same thing:
#
#   (a) Extending BELOW the original bottom (more column before cell 1):
#       x_top keeps INCREASING past 52 cm - the fitted decay continues
#       naturally and Cb keeps falling toward zero. Well-posed.
#   (b) Extending ABOVE the original top (more column past cell 10, which
#       is what "expand to 1 m" has meant everywhere else in this project):
#       those cells are BEYOND where biomass was ever introduced. Sliding
#       x_top negative and evaluating exp(+lambda*|x_top|) would make Cb
#       grow again - not because more bacteria are there, but because nothing
#       in a purely mathematical extrapolation knows the injection zone
#       ended. Physically the sane assumption for "same biomass inoculation,
#       longer column of otherwise identical sand" is Cb -> 0 beyond x=52,
#       since no bacteria were ever introduced past that point.
# ---------------------------------------------------------------------------

print("\nExtension beyond the calibrated 52 cm domain depends on which physical "
      "scenario is intended (see docstring/comments) - this script does not "
      "assume one; it only refits the validated 0-52cm range.")

fig, ax = plt.subplots(figsize=(6.5, 6))
x_fine = np.linspace(0, COLUMN_LEN, 200)
xtop_fine = COLUMN_LEN - x_fine
ax.plot(cb_cubic(x_fine), x_fine, "k-", lw=1, alpha=0.5, label="Original cubic (Razbani et al. 2025)")
ax.plot(Cb_data, x, "ko", ms=7, label="Calibration points (10 real cell centers)")
ax.plot(biexp(xtop_fine, *popt), x_fine, "b--", lw=2,
        label=f"Two-site CFT refit (R²={r2_biexp:.3f})")
ax.plot(Cb0_single * np.exp(-lam_single * xtop_fine), x_fine, "g:", lw=1.5,
        label=f"Single-site CFT (R²={r2_single:.3f})")
ax.set_xlabel("Cb (attached biomass, OD-equivalent)")
ax.set_ylabel("Distance from bottom / cementation inlet (cm)")
ax.legend(fontsize=8.5)
ax.grid(alpha=0.4)
fig.tight_layout()
fig.savefig("MICP_biomass_CFT_refit.png", dpi=150)
print("saved MICP_biomass_CFT_refit.png")
plt.show()
