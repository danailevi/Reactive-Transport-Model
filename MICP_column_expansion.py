"""Hybrid depth-extension of the Test 3B model (0.52 m) toward 1.5 m, using
Wang et al. (2026) U5 column data ("Large-scale experimental investigation
of biotreated sand column using different grouting pipe configurations").

Why "hybrid" and not a straight extended PHREEQC run: the attached-biomass
profile Cb(x) is only calibrated over 0-52 cm (Razbani et al. 2025, Fig. 2).
Extrapolating that cubic fit further either blows up unboundedly (raw) or
requires an arbitrary "hold it flat" assumption with zero empirical support
(see MICP_Simulation_T3B_depth.py) - either way, PHREEQC output beyond 52 cm
is not a validated prediction of THIS system, since Martinez et al.'s column
was only ever run at 52 cm.

U5 is a *different* experimental system (different sand, bacteria, pipe
injection geometry, treatment recipe) than the Test 3B model, so it cannot
directly validate this model's magnitude at depth. What it DOES provide is a
real, measured example of how ureolytic MICP treatment actually declines
with distance past its own peak in a comparable meter-scale column. That
decline SHAPE (normalized to U5's own peak, not its absolute values) is
borrowed here and anchored to our own validated 52 cm peak, giving an
extension that is at least grounded in one real spatial-decline observation
rather than an arbitrary flat assumption. This is still not a validated
prediction of Test 3B at 1.5 m - it is the most defensible approximation
available without new attached-biomass data spanning that depth.
"""

import phreeqpy.iphreeqc.phreeqc_com as phreeqc_mod
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

phreeqc = phreeqc_mod.IPhreeqc()
phreeqc.load_database(r"C:\Program Files\USGS\IPhreeqcCOM 3.8.6-17100\database\phreeqc.dat")

EXTEND_TO_CM = 150.0  # match U5's column length

# ---------------------------------------------------------------------------
# 1. Validated near-field: the real, calibrated 10-cell / 52 cm Test 3B model.
#    No cap, no extrapolation - Cb(x) stays entirely inside its fitted domain.
# ---------------------------------------------------------------------------
initial_conditions = """
SOLUTION_MASTER_SPECIES
    Urea Urea  0 CO(NH2)2 60.06
     Amm   AmmH+  0.0  N    14.0067

SOLUTION_SPECIES
 Urea = Urea
            log_k  0
            Vm  14.5
    AmmH+ = AmmH+
        log_k 0
    AmmH+ = Amm + H+
        log_k -9.23
       Cl- + H+ = HCl
        log_k 0.7
        delta_h 0 kcal
        -gamma 0 0.06
    Cl- + Ca+2 = CaCl+
        log_k 0.3
        delta_h 0 kcal
        -gamma 0 0.06
    2Cl- + Ca+2 = CaCl2
        log_k .66
        delta_h 0 kcal
        -gamma 0 0.06

KNOBS
    -iterations 500
    -convergence_tolerance 1e-9
    -tolerance 1e-14
    -step_size 100
    -pe_step_size 5
    -diagonal_scale true
END

RATES
Ureolysis
-START
11 cell_l = 5.2
12 t_d    = 288000*1.2
15 initial_por = .37
16 Cb = (5.66e-1*(cell_no*cell_l-cell_l/2)^3+1.1e2*(cell_no*cell_l-cell_l/2)^2+4.67e3*(cell_no*cell_l-cell_l/2)-5.62e3)/1e6
18 decay_ratio = exp(-TOTAL_TIME/t_d)
20 rate = Cb*soln_vol*get_por(cell_no)/initial_por*PARM(1)*MOL("Urea")/(MOL("Urea")+PARM(2))
30 moles = rate * TIME * decay_ratio
40 save moles
-END

Calcite
-START
20 initial_por = .37
21 V_cell = 102
22 d_calcite = 36.9
23 pow_a = 2
25 A0 = 1.
30 si_cc = SI("Calcite")
70 Area = A0 * (initial_por/get_por(cell_no))^pow_a * soln_vol
80 rate = Area * (parm(1) + parm(2)*act("H+")) * (1-10^SI("Calcite"))
100 moles = rate * TIME * get_por(cell_no)/initial_por
110 change_por(initial_por - d_calcite*kin("Calcite")/V_cell, cell_no)
200 SAVE moles
-END

SOLUTION 0
Urea 1.015
Amm .0567
Ca 2
Cl 4.0567
C .0038 as Ca Calcite 0.0
Na .0038
-units mol/kgw
-pH   7  charge_balance
water 38.7e-3
END

SOLUTION 1-10
    -units    mol/kgw
    -pH   7  charge_balance
    water 38.7e-3
END

Kinetics 1-10
    Ureolysis
         -m0 100
        -formula AmmH+ 2 CO3-2 1 Urea -1 H2O -2
        -parms 7.0e-5 0.0173
    Calcite
        -m0   0
        -parms 1.55e-6 .501
END

TRANSPORT
   cells 10
   shifts 178
   time_step 1055 seconds
   flow_direction forward
   boundary_conditions flux constant
   length 0.052
   diffusion_coefficient 2.7e-10
   porosities 10*.37
   -dispersivities 10*.005

SELECTED_OUTPUT
    -step    true
     -h    false
    -sim    false
    -soln   true
     -pe     false
    -totals false
    -time   true
    -pH true
    -tot  Ca  Urea Amm
    -si Calcite
    -k Calcite
    -dist false
    -state false
END
"""

phreeqc.run_string(initial_conditions)
output = phreeqc.get_selected_output_array()
table = pd.DataFrame(output[1:], columns=output[0])

N_CELLS = 10
V_cell = 102e-3  # L
l_cell = 5.2     # cm
shift_m = 178

final = table[table["step"] == shift_m]
final = final[(final["soln"] >= 1) & (final["soln"] <= N_CELLS)].sort_values("soln")

heights = np.arange(1, N_CELLS + 1) * l_cell - l_cell / 2   # cm from inlet, validated range
calcite = final["k_Calcite"].values / V_cell                 # mol/L

final.to_csv(f"T3B_calcite_shift{shift_m}.csv", index=False)

# ---------------------------------------------------------------------------
# 2. U5 real data: find its own peak, and fit an exponential decay to its
#    post-peak decline, normalized to ITS OWN peak (shape only, not magnitude).
# ---------------------------------------------------------------------------
u5 = pd.read_csv("U5_depth_averaged.csv")
u5["calcite_molL"] = u5["CCC_areaweighted_pct"] * 0.156      # CCC% -> mol/L bulk
u5 = u5.sort_values("depth_from_top_cm").reset_index(drop=True)

u5_peak_idx = u5["calcite_molL"].idxmax()
u5_peak_depth = u5.loc[u5_peak_idx, "depth_from_top_cm"]
u5_peak_val = u5.loc[u5_peak_idx, "calcite_molL"]

post_peak = u5[u5["depth_from_top_cm"] > u5_peak_depth]
delta_u5 = post_peak["depth_from_top_cm"].values - u5_peak_depth
frac_u5 = post_peak["calcite_molL"].values / u5_peak_val

# ln(frac) = -delta/L, fit forced through the origin (delta=0 -> frac=1 by definition)
slope = np.sum(delta_u5 * np.log(frac_u5)) / np.sum(delta_u5 ** 2)
decay_length = -1 / slope
pred = np.exp(slope * delta_u5)
ss_res = np.sum((frac_u5 - pred) ** 2)
ss_tot = np.sum((frac_u5 - frac_u5.mean()) ** 2)
r2 = 1 - ss_res / ss_tot

print(f"U5 peak: {u5_peak_val:.3f} mol/L at {u5_peak_depth:.0f} cm")
print(f"Post-peak decay length fitted from U5: {decay_length:.0f} cm (R^2={r2:.3f}, "
      f"n={len(delta_u5)} points beyond the peak)")

# ---------------------------------------------------------------------------
# 3. Hybrid extension: anchor OUR validated peak, apply U5's decay SHAPE.
# ---------------------------------------------------------------------------
our_peak_depth = heights[-1]
our_peak_val = calcite[-1]

ext_depths = np.arange(our_peak_depth, EXTEND_TO_CM + 1e-9, 1.0)
ext_calcite = our_peak_val * np.exp(-(ext_depths - our_peak_depth) / decay_length)

print(f"Our validated peak: {our_peak_val:.3f} mol/L at {our_peak_depth:.1f} cm")
print(f"Hybrid extension at {EXTEND_TO_CM:.0f} cm: {ext_calcite[-1]:.3f} mol/L "
      f"({ext_calcite[-1]/our_peak_val*100:.1f}% of our peak)")

# ---------------------------------------------------------------------------
# 4. Plot: validated model (solid) + hybrid tail (dashed) + U5 reference.
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.5, 7))
ax.plot(calcite, heights, "ko-", label="Simulation (validated, 0-52 cm)", ms=6)
ax.plot(ext_calcite, ext_depths, "b--",
        label=f"Hybrid extension (U5 decay shape, L={decay_length:.0f} cm)")
ax.plot(u5["calcite_molL"], u5["depth_from_top_cm"], "r^-",
        label="U5, Wang et al. 2026 (different system)", ms=6, alpha=0.8)
ax.axvline(0, color="gray", lw=0.5)
ax.set_xlabel("Calcite (mol/L)")
ax.set_ylabel("Distance from inlet (cm)")
ax.legend(fontsize=9)
ax.grid(alpha=0.4)
fig.tight_layout()
fig.savefig("MICP_column_expansion_hybrid.png", dpi=150)
print("saved MICP_column_expansion_hybrid.png")
plt.show()
