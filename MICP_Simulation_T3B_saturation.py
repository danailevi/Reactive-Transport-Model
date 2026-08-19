"""Time evolution of calcite precipitation for the extended (98.8 cm) Test 3B column.

Answers: at what point does calcite precipitation in a given cell stop
increasing, and what happens after that?

Two competing mechanisms could set the plateau:
  1. Biomass decay - ureolysis is scaled by decay_ratio = exp(-TOTAL_TIME/t_d)
     with t_d ~ 4 days, so new carbonate production should become negligible
     after ~3-6 t_d regardless of cell position.
  2. Porosity clogging - the reactive-surface-area term
     Area = A0*(initial_por/get_por(cell_no))^2 has no floor; if a cell's
     porosity is driven to zero (or below) by rapid local precipitation, this
     term blows up or the model becomes non-physical, which would happen
     *before* the biomass-decay plateau.

We track three cells over ~6 biomass-decay time constants: cell 1 (bottom,
injection, low attached-biomass), cell 10 (mid-column, within the original
0-52 cm calibrated range), and cell 19 (top/outlet, deep in the
extrapolated - and, per MICP_Simulation_T3B.py, unbounded - biomass profile).
"""

import phreeqpy.iphreeqc.phreeqc_com as phreeqc_mod
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DATABASE = r"C:\Program Files\USGS\IPhreeqcCOM 3.8.6-17100\database\phreeqc.dat"

N_CELLS = 19          # 19 * 5.2 cm = 98.8 cm, matches MICP_Simulation_T3B.py
CELL_L = 5.2           # cm
V_CELL_L = 102e-3      # L
TIME_STEP = 1055       # s
T_D = 288000 * 1.2     # biomass decay constant, s (~4.0 days)
N_TD = 6               # simulate out to N_TD * t_d
SHIFTS = int(np.ceil(N_TD * T_D / TIME_STEP))
TRACK_CELLS = (1, 10, 19)


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
   -punch_cells {' '.join(str(c) for c in TRACK_CELLS)}
   -punch_frequency 1

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
    -k Ureolysis Calcite
    -dist false
    -state false
END
"""


def run():
    phreeqc = phreeqc_mod.IPhreeqc()
    phreeqc.load_database(DATABASE)

    print(f"t_d = {T_D/86400:.2f} days; simulating {N_TD} t_d = "
          f"{SHIFTS} shifts ({SHIFTS*TIME_STEP/86400:.1f} days)")
    print(f"tracking cells {TRACK_CELLS} (bottom=1/injection, top={N_CELLS}/outlet)")

    phreeqc.run_string(build_input())
    out = phreeqc.get_selected_output_array()
    df = pd.DataFrame(out[1:], columns=out[0])
    return df


def analyze(df):
    """For each tracked cell: calcite(t), porosity(t), and where each
    mechanism (biomass-decay plateau vs porosity clogging) kicks in."""
    results = {}
    for cell in TRACK_CELLS:
        sub = df[df["soln"] == cell].sort_values("step")
        t_days = sub["time"].values / 86400
        calcite = sub["k_Calcite"].values / V_CELL_L
        por = sub["por"].values

        # Plateau: first time step where calcite has stopped changing by
        # more than 0.5% of its final value, and stays that way.
        final_calcite = calcite[-1]
        if final_calcite > 0:
            frac_of_final = calcite / final_calcite
            plateaued = np.where(frac_of_final >= 0.995)[0]
            t_plateau = t_days[plateaued[0]] if len(plateaued) else np.nan
        else:
            t_plateau = np.nan

        # Clogging: first time porosity reaches (or would go below) zero.
        clogged = np.where(por <= 1e-6)[0]
        t_clog = t_days[clogged[0]] if len(clogged) else np.nan

        results[cell] = dict(t_days=t_days, calcite=calcite, por=por,
                              t_plateau=t_plateau, t_clog=t_clog,
                              final_calcite=final_calcite,
                              final_por=por[-1] if len(por) else np.nan)

        print(f"\ncell {cell}:")
        print(f"  final calcite = {final_calcite:.3f} mol/L, "
              f"final porosity = {results[cell]['final_por']:.3f}")
        print(f"  reaches 99.5% of its final value at t = {t_plateau:.2f} days"
              if not np.isnan(t_plateau) else "  never plateaus within simulated window")
        print(f"  porosity hits ~0 at t = {t_clog:.2f} days"
              if not np.isnan(t_clog) else "  porosity never reaches 0 within simulated window")
    return results


if __name__ == "__main__":
    df = run()
    results = analyze(df)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8), sharex=True)

    t_ref = np.linspace(0, N_TD * T_D / 86400, 200)
    decay_ratio = np.exp(-t_ref * 86400 / T_D)

    for cell in TRACK_CELLS:
        r = results[cell]
        ax1.plot(r["t_days"], r["calcite"], marker="o", ms=3, label=f"cell {cell}")
        ax2.plot(r["t_days"], r["por"], marker="o", ms=3, label=f"cell {cell}")
        if not np.isnan(r["t_clog"]):
            ax1.axvline(r["t_clog"], color="C" + str(TRACK_CELLS.index(cell)),
                        ls=":", alpha=0.5)

    ax1b = ax1.twinx()
    ax1b.plot(t_ref, decay_ratio, "k--", alpha=0.4, label="decay_ratio (ureolysis envelope)")
    ax1b.set_ylabel("decay_ratio = exp(-t/t_d)")
    ax1b.set_ylim(0, 1.05)

    ax1.set_ylabel("Calcite precipitation (mol/L)")
    ax1.grid()
    ax1.legend(loc="upper left")
    ax1b.legend(loc="upper right")

    ax2.axhline(0, color="red", ls="--", alpha=0.5, label="porosity = 0")
    ax2.set_xlabel("Time (days)")
    ax2.set_ylabel("Porosity")
    ax2.grid()
    ax2.legend()

    plt.tight_layout()
    plt.savefig("T3B_saturation_timeseries.png", dpi=150)

    final = df[df["step"] == df["step"].max()]
    final.to_csv("T3B_saturation_final.csv", index=False)
    df.to_csv("T3B_saturation_full_timeseries.csv", index=False)

    plt.show()
