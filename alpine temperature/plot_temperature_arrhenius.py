"""Compare temperature-dependence of ureolysis rate from different literature sources
(Wang et al. 2023, Sun et al. 2018), fit Arrhenius, and extrapolate into the alpine
range (4-12 C) relevant for this project. Produces two panels:
  1. Arrhenius plot (ln rate vs 1/T) with both fits, alpine zone shaded.
  2. Predicted R_max vs temperature, scaled from the model's baseline (7e-5 mol/L/s/OD
     at 25 C), using each fitted Ea.
"""
import numpy as np
import matplotlib.pyplot as plt

R = 8.314  # J/mol/K

# --- literature data ---
T_wang_C = np.array([4, 20, 35, 50])
r_wang = np.array([0.018, 0.038, 0.056, 0.135])  # M/h/OD (CaCO3 precipitation rate)

T_sun_C = np.array([15, 20, 25, 30])
r_sun = np.array([3.3e-6, 6.7e-6, 8.3e-6, 10.0e-6])  # mol/L/s (pure urease activity)

datasets = {
    "Wang 2023 (precipitation rate)": (T_wang_C, r_wang, "tab:blue"),
    "Sun 2018 (urease activity)": (T_sun_C, r_sun, "tab:orange"),
}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

fits = {}
for label, (T_C, r, color) in datasets.items():
    T_K = T_C + 273.15
    x = 1 / T_K
    y = np.log(r)
    slope, intercept = np.polyfit(x, y, 1)
    Ea = -slope * R
    fits[label] = (slope, intercept, Ea, color)

    ax1.scatter(x, y, color=color, zorder=3, label=f"{label} (data)")
    x_fit = np.linspace(1 / (55 + 273.15), 1 / (2 + 273.15), 100)
    ax1.plot(x_fit, slope * x_fit + intercept, color=color, linewidth=1.5)

    # highlight the extrapolated (alpine, 4-12C) portion as dashed
    x_alpine = np.linspace(1 / (12 + 273.15), 1 / (4 + 273.15), 50)
    ax1.plot(x_alpine, slope * x_alpine + intercept, color=color, linestyle="--", linewidth=2.5)

# shade the alpine zone
ax1.axvspan(1 / (12 + 273.15), 1 / (4 + 273.15), color="gray", alpha=0.15)
ax1.text(1 / (8 + 273.15), ax1.get_ylim()[1] * 0 + max(y) , "zona alpina\n4-12°C",
         ha="center", va="top", fontsize=9, color="gray")

ax1.set_xlabel("1/T (1/K)")
ax1.set_ylabel("ln(velocità di reazione)")
ax1.set_title("Arrhenius plot: due fonti, due pendenze diverse\n(tratteggio = estrapolazione, non misurato)")
ax1.legend(fontsize=8, loc="lower left")
ax1.invert_xaxis()  # so temperature increases left-to-right, more intuitive
ax1.grid(alpha=0.3)

# --- panel 2: predicted R_max vs T for the actual model parameter ---
Rmax_ref = 7e-5  # mol/L/s/OD, T3B baseline
T_ref_K = 25 + 273.15

T_range_C = np.linspace(2, 32, 200)
T_range_K = T_range_C + 273.15

for label, (slope, intercept, Ea, color) in fits.items():
    ratio = np.exp(-Ea / R * (1 / T_range_K - 1 / T_ref_K))
    Rmax_pred = Rmax_ref * ratio
    ax2.plot(T_range_C, Rmax_pred, color=color, label=f"{label.split(' (')[0]}: Ea={Ea/1000:.1f} kJ/mol")

ax2.axvspan(4, 12, color="gray", alpha=0.15)
ax2.axhline(Rmax_ref, color="black", linestyle=":", linewidth=1, label="R_max baseline (25°C)")
ax2.set_xlabel("Temperatura (°C)")
ax2.set_ylabel("R_max previsto (mol/L/s/OD)")
ax2.set_title("R_max scalato al variare della temperatura")
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3)

fig.tight_layout()
fig.savefig("temperature_arrhenius_comparison.png", dpi=150)
print("saved temperature_arrhenius_comparison.png")
for label, (slope, intercept, Ea, color) in fits.items():
    print(label, "Ea =", Ea / 1000, "kJ/mol")
