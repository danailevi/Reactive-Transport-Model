"""MICP Test 3B with configurable grid resolution.

Reproduces Fig. 7 of Razbani et al. (2025): grid-convergence check by
running the same physical problem at 10 and 20 cells and overlaying the
calcite profiles.

Everything that depends on the grid is derived from `n_cells`:
    cell_l    = H / n                    (cm)
    V_cell    = 102 mL * (10 / n)        (bulk volume per cell)
    water     = 38.7 g * (10 / n)        (pore water per cell)
    time_step = 1055 s * (10 / n)        (residence time per cell)
    shifts    = total_time / time_step   (keeps total treatment time fixed)

The biomass polynomial uses (cell_no*cell_l - cell_l/2), which is the cell
CENTRE in cm, so it evaluates to physically-correct values regardless of n.
"""

import phreeqpy.iphreeqc.phreeqc_com as phreeqc_mod
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


DATABASE = r"C:\Program Files\USGS\IPhreeqcCOM 3.8.6-17100\database\phreeqc.dat"


def build_input(n_cells, column_height_cm=52.0, total_time_s=187_790):
    """Return (phreeqc_input_string, cell_l_cm, V_cell_L, shifts)."""
    cell_l    = column_height_cm / n_cells
    V_cell_mL = 102.0 * (10.0 / n_cells)
    water_kg  = 38.7e-3 * (10.0 / n_cells)
    time_step = 1055.0 * (10.0 / n_cells)
    shifts    = int(round(total_time_s / time_step))

    header = """
SOLUTION_MASTER_SPECIES
    Urea    Urea    0 CO(NH2)2  60.06
    Amm     AmmH+   0.0    N    14.0067

SOLUTION_SPECIES
Urea = Urea
    log_k 0
    Vm 14.5
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
"""

    body = f"""
RATES
Ureolysis
-START
11 cell_l = {cell_l}
12 t_d = 288000
15 initial_por = .38
16 Cb = (5.66e-1*(cell_no*cell_l-cell_l/2)^3+1.1e2*(cell_no*cell_l-cell_l/2)^2+4.67e3*(cell_no*cell_l-cell_l/2)-5.62e3)/1e6
18 decay_ratio = exp(-TOTAL_TIME/t_d)
20 rate = Cb*soln_vol*get_por(cell_no)/initial_por*PARM(1)*MOL("Urea")/(MOL("Urea")+PARM(2))
30 moles = rate * TIME * decay_ratio
40 save moles
-END

Calcite
-START
20 initial_por = .38
21 V_cell = {V_cell_mL}
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
Urea .05
Amm .0567
Ca .05
Cl .1567
C .0038 as Ca Calcite 0.0
Na .0038
-units mol/kgw
-pH 7 charge_balance
water {water_kg}
END

SOLUTION 1-{n_cells}
    -units    mol/kgw
    -pH  7  charge_balance
    water {water_kg}
END

Kinetics 1-{n_cells}
    Ureolysis
        -m0 100
        -formula AmmH+ 2 CO3-2 1 Urea -1 H2O -2
        -parms 7.0e-5 0.0173
    Calcite
        -m0 0
        -parms 1.55e-6 .501
END

TRANSPORT
    cells {n_cells}
    shifts {shifts}
    time_step {time_step} seconds
    flow_direction forward
    boundary_conditions flux constant
    length {cell_l / 100}
    diffusion_coefficient 1.0e-9
    porosities {n_cells}*.38
    -dispersivities {n_cells}*.01

USER_PUNCH
headings por
 10 punch get_por(cell_no)

SELECTED_OUTPUT
    -step true
    -h false
    -sim false
    -soln true
    -pe false
    -totals false
    -time true
    -pH true
    -tot Ca Urea Amm
    -si Calcite
    -k Ureolysis Calcite
    -dist false
    -state false
END
"""
    return header + body, cell_l, V_cell_mL * 1e-3, shifts


def run(n_cells):
    """Run one grid resolution; return (heights_cm, calcite_mol_per_L)."""
    phreeqc = phreeqc_mod.IPhreeqc()
    phreeqc.load_database(DATABASE)

    inp, cell_l, V_cell_L, shifts = build_input(n_cells)
    print(f"  n={n_cells:>3d} | cell_l={cell_l:.3f} cm | "
          f"shifts={shifts} | time_step={1055*10/n_cells:.1f} s")

    phreeqc.run_string(inp)
    out = phreeqc.get_selected_output_array()
    df  = pd.DataFrame(out[1:], columns=out[0])

    final_step = int(df["step"].max())
    final = (df[(df["step"] == final_step)
                & (df["soln"] >= 1)
                & (df["soln"] <= n_cells)]
             .sort_values("soln"))

    heights = np.arange(1, n_cells + 1) * cell_l - cell_l / 2
    calcite = final["k_Calcite"].values / V_cell_L
    return heights, calcite


if __name__ == "__main__":
    print("Grid-convergence test (Razbani et al. 2025, Fig. 7):")
    h10, c10 = run(5)
    h20, c20 = run(30)

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(c10, h10, "ko", ms=8,          label="Simulation, 5 cells")
    ax.plot(c20, h20, "k+", ms=10, mew=2,  label="Simulation, 30 cells")
    ax.set_xlabel("Calcite precipitation (mol/L)")
    ax.set_ylabel("Column height (cm)")
    ax.set_xlim(0, 0.6)
    ax.set_ylim(0, 52)
    ax.grid()
    ax.legend()
    plt.tight_layout()
    plt.savefig("T3B_grid_convergence.png", dpi=150)
    plt.show()
