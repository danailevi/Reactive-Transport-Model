"""Maximum treatment depth of MICP via substrate (urea) depletion.

Question: in a hypothetical, much deeper soil column, how far down would a
single continuous-flow MICP treatment actually reach before the injected
urea is consumed?

Method: cells 1-10 use the real, calibrated attached-biomass profile Cb(x)
from Razbani et al. (2025), fit to measured data over 0-52 cm. Beyond cell
10, Cb is held FLAT at its highest calibrated value (at x=49.4 cm) rather
than extrapolated - this is a deliberate best-case upper bound (sustained
peak bacterial density at every depth), not a claim about real biomass
beyond 52 cm. The column is run for ~6 biomass-decay time constants
(t_d ~ 4 days) so ureolysis has essentially finished everywhere the front
has reached, giving the true asymptotic depth profile rather than a
still-rising transient.

Depth is assessed via CUMULATIVE urea consumed per cell (k_Ureolysis) and
total calcite precipitated, not instantaneous urea concentration - once
biomass is exhausted, flow-through concentration relaxes back toward the
inlet value everywhere, so a concentration snapshot alone is misleading
about how far treatment actually reached.
"""

import phreeqpy.iphreeqc.phreeqc_com as phreeqc_mod
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DATABASE = r"C:\Program Files\USGS\IPhreeqcCOM 3.8.6-17100\database\phreeqc.dat"

N_CELLS = 40            # 40 * 5.2 cm = 208 cm hypothetical column
CELL_L = 5.2             # cm
X_CAP = 49.4             # cm; cell-10 centre, edge of the calibrated Cb domain
V_CELL_L = 102e-3        # L
TIME_STEP = 1055         # s
T_D = 288000 * 1.2       # biomass decay constant, s (~4.0 days)
N_TD = 6                 # simulate out to N_TD * t_d
SHIFTS = int(np.ceil(N_TD * T_D / TIME_STEP))
UREA_INLET = 1.015       # mol/kgw


def build_input():
    return f"""
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
11 cell_l = {CELL_L}
12 t_d    = {T_D}
13 x_cap  = {X_CAP}
15 initial_por = .37
16 x = cell_no*cell_l-cell_l/2
17 IF x > x_cap THEN x = x_cap
18 Cb = (5.66e-1*x^3+1.1e2*x^2+4.67e3*x-5.62e3)/1e6
20 decay_ratio = exp(-TOTAL_TIME/t_d)
22 rate = Cb*soln_vol*get_por(cell_no)/initial_por*PARM(1)*MOL("Urea")/(MOL("Urea")+PARM(2))
30 moles = rate * TIME * decay_ratio
40 save moles
-END

Calcite
-START
20 initial_por = .37
21 V_cell = {V_CELL_L * 1000}
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
Urea {UREA_INLET}
Amm .0567
Ca 2
Cl 4.0567
C .0038 as Ca Calcite 0.0
Na .0038
-units mol/kgw
-pH   7  charge_balance
water 38.7e-3
END

SOLUTION 1-{N_CELLS}
    -units    mol/kgw
    -pH   7  charge_balance
    water 38.7e-3
END

Kinetics 1-{N_CELLS}
    Ureolysis
         -m0 100
        -formula AmmH+ 2 CO3-2 1 Urea -1 H2O -2
        -parms 7.0e-5 0.0173
    Calcite
        -m0   0
        -parms 1.55e-6 .501
END

TRANSPORT
   cells {N_CELLS}
   shifts {SHIFTS}
   time_step {TIME_STEP} seconds
   flow_direction forward
   boundary_conditions flux constant
   length 0.052
   diffusion_coefficient 2.7e-10
   porosities {N_CELLS}*.37
   -dispersivities {N_CELLS}*.005
   -punch_frequency 20

USER_PUNCH
headings por
 10 punch get_por(cell_no)

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


def run():
    phreeqc = phreeqc_mod.IPhreeqc()
    phreeqc.load_database(DATABASE)

    print(f"{N_CELLS} cells * {CELL_L} cm = {N_CELLS*CELL_L:.0f} cm hypothetical column")
    print(f"Cb(x): real calibrated cubic for x<={X_CAP}cm, held flat at its peak "
          f"value beyond that (best-case, not extrapolated)")
    print(f"Simulating {N_TD} t_d = {SHIFTS} shifts ({SHIFTS*TIME_STEP/86400:.1f} days)")

    phreeqc.run_string(build_input())
    out = phreeqc.get_selected_output_array()
    df = pd.DataFrame(out[1:], columns=out[0])
    return df


def analyze(df):
    """Depth profile from cumulative calcite (a real, monotonic tracer of how
    much reaction has actually happened at each depth). NOTE: PHREEQC's
    k_Ureolysis is NOT usable as "moles of urea consumed" here - it tracks
    the reaction extent of the Ureolysis KINETICS placeholder, whose -m0 100
    is an arbitrary large number chosen only so the reactant never
    numerically runs out; it saturates near its own m0 ceiling at nearly
    every cell regardless of depth and does not represent real urea mass.
    Real urea mass balance would need to be read from -tot Urea
    concentrations with a proper inflow/outflow integral, which this script
    does not attempt.
    """
    final_step = int(df["step"].max())
    final = (df[(df["step"] == final_step) & (df["soln"] >= 1) & (df["soln"] <= N_CELLS)]
             .sort_values("soln"))

    heights = np.arange(1, N_CELLS + 1) * CELL_L - CELL_L / 2
    calcite = final["k_Calcite"].values / V_CELL_L              # mol/L
    urea_frac = final["Urea(mol/kgw)"].values / UREA_INLET      # instantaneous, NOT cumulative

    print(f"\nFinal state at t = {final_step*TIME_STEP/86400:.1f} days:")
    peak_idx = np.argmax(calcite)
    print(f"  calcite peaks at depth {heights[peak_idx]:.1f} cm ({calcite[peak_idx]:.3f} mol/L)")
    print(f"  calcite at the deepest simulated cell ({heights[-1]:.1f} cm): "
          f"{calcite[-1]:.3f} mol/L ({calcite[-1]/calcite[peak_idx]*100:.1f}% of peak)")

    # exponential fit to the post-peak decline, to extrapolate where
    # precipitation would drop to a target fraction of its peak
    post_peak = heights > heights[peak_idx]
    if post_peak.sum() >= 3:
        x = heights[post_peak]
        y = calcite[post_peak]
        slope, _ = np.polyfit(x, np.log(y), 1)  # y = exp(intercept)*exp(slope*x)
        decay_depth = -1 / slope
        print(f"  post-peak decline fits ~exponential with characteristic "
              f"length {decay_depth:.0f} cm (fit over {x.min():.0f}-{x.max():.0f} cm, "
              "extrapolate with caution)")
        for target_frac in (0.5, 0.1):
            depth_at_frac = heights[peak_idx] + decay_depth * np.log(1 / target_frac)
            print(f"    depth where calcite ~ falls to {target_frac*100:.0f}% of peak "
                  f"(extrapolated): {depth_at_frac:.0f} cm")

    return heights, calcite, urea_frac


if __name__ == "__main__":
    df = run()
    heights, calcite, urea_frac = analyze(df)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True)

    ax1.plot(heights, calcite, "ko-", ms=4)
    ax1.axvline(X_CAP, color="red", ls="--", alpha=0.5,
                label=f"edge of calibrated domain ({X_CAP:.0f} cm)")
    ax1.set_ylabel("Calcite precipitation (mol/L)")
    ax1.grid()
    ax1.legend()

    ax2.plot(heights, urea_frac, "bo-", ms=4)
    ax2.axvline(X_CAP, color="red", ls="--", alpha=0.5)
    ax2.set_xlabel("Depth (cm)")
    ax2.set_ylabel("Urea remaining (fraction of inlet, instantaneous)")
    ax2.grid()

    plt.suptitle(f"Substrate-depletion depth check\n(best-case flat biomass beyond "
                 f"{X_CAP:.0f} cm; calcite is the reliable cumulative signal)")
    plt.tight_layout()
    plt.savefig("T3B_depth_substrate_depletion.png", dpi=150)

    out_df = pd.DataFrame({"depth_cm": heights, "calcite_mol_L": calcite,
                            "urea_frac_remaining": urea_frac})
    out_df.to_csv("T3B_depth_substrate_depletion.csv", index=False)

    plt.show()
