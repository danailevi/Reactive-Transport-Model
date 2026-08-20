"""Proof-of-concept 2D extension of the MICP model via operator splitting.

Architecture: PHREEQC still does ONLY the local chemistry (ureolysis +
calcite precipitation), exactly as in test_3B.in's RATES blocks. Python
does the transport: at each time step we write one PHREEQC batch-reaction
input covering every cell in a 2D grid, run it once, read back the new
concentrations, then move mass between neighbouring cells ourselves
(advection in x = flow direction, dispersion in x and y) before writing
the next step's input.

Simplifications vs. the validated 1D model (flagged clearly, not hidden):
  - No porosity feedback (get_por/change_por are TRANSPORT-only intrinsics
    in PHREEQC and aren't available in plain batch REACTION mode).
  - Biomass density (Cb) and reactive area are held at fixed uniform
    values per cell, rather than the fitted spatial polynomial used in
    the real 1D column.
  - Transport coefficients (advection/dispersion fractions per step) are
    illustrative, not calibrated to a real Peclet number.
This is meant to demonstrate the COUPLING MECHANISM, not to be a
publication-ready 2D result.

Injection is a POINT SOURCE on the left edge (not the whole edge), so the
plume has to spread sideways as well as advect forward -- the one thing a
1D column can never show.
"""
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PHREEQC = "/home/cacchina_fra/iGEM/phreeqc_install/bin/phreeqc"

NX, NY = 5, 5                 # grid size
INJECT = (0, 2)                # (x, y) of the point source, left edge, middle row
N_STEPS = 20
DT = 1055                      # s, same time step as test_3B.in
ADVECT_FRAC = 0.5              # fraction of a cell's solute pushed to the next cell in x each step
DISPERSE_FRAC = 0.12           # fraction mixed with each y-neighbour each step (lateral spreading)

RMAX, KM = 7.0e-5, 0.0173      # same as test_3B.in
CB = 3.0e-2                    # uniform biomass density (arbitrary fixed value, mol-scale placeholder)
K1, K2, AREA = 1.55e-6, 0.501, 1.0

SPECIES_BLOCK = """SOLUTION_MASTER_SPECIES
    Urea Urea  0 CO(NH2)2 60.06

SOLUTION_SPECIES
    Urea = Urea
            log_k  0
            Vm  14.5
    Cl- + H+ = HCl
            log_k   0.7
            delta_h 0       kcal
            -gamma  0       0.06
    Cl- + Ca+2 = CaCl+
            log_k   0.3
            delta_h 0       kcal
            -gamma  0       0.06
    2Cl- + Ca+2 = CaCl2
            log_k   .66
            delta_h 0       kcal
            -gamma  0       0.06

RATES
Ureolysis
-START
10 rate = PARM(3)*soln_vol*PARM(1)*MOL("Urea")/(MOL("Urea")+PARM(2))
20 moles = rate * TIME
30 SAVE moles
-END

Calcite
-START
10 rate = PARM(3)*(PARM(1)+PARM(2)*act("H+"))*(1-10^SI("Calcite"))
20 moles = rate * TIME
30 SAVE moles
-END
"""

INLET = dict(Urea=0.05, Amm=0.0567, Ca=0.05, Cl=0.1567, C=0.0038, Na=0.0038)
FRESH = dict(Urea=0.0, Amm=0.0, Ca=0.0, Cl=0.0, C=0.0, Na=0.0)

# state[x, y] -> dict of concentrations (mol/kgw)
state = [[dict(FRESH) for _ in range(NY)] for _ in range(NX)]
cum_calcite = np.zeros((NX, NY))

snapshots = {}

def cell_id(x, y):
    return x * NY + y + 1

SELECTED_OUTPUT_BLOCK = """SELECTED_OUTPUT
    -step true
    -soln true
    -pH true
    -tot Ca Urea Amm C
    -k Ureolysis Calcite
END
"""

def write_input(step):
    lines = [SPECIES_BLOCK, SELECTED_OUTPUT_BLOCK]
    for x in range(NX):
        for y in range(NY):
            cid = cell_id(x, y)
            c = INLET if (x, y) == INJECT else state[x][y]
            lines.append(f"""SOLUTION {cid}
Urea {c['Urea']:.6e}
Amm {c['Amm']:.6e}
Ca {c['Ca']:.6e}
Cl {c['Cl']:.6e}
C {c['C']:.6e}
Na {c['Na']:.6e}
-units mol/kgw
-pH 7 charge_balance
water 38.7e-3
END
""")
            lines.append(f"""KINETICS {cid}
    Ureolysis
        -m0 100
        -formula AmmH+ 2 CO3-2 1 Urea -1 H2O -2
        -parms {RMAX} {KM} {CB}
    Calcite
        -m0 {cum_calcite[x, y]:.6e}
        -parms {K1} {K2} {AREA}
-steps {DT} seconds in 1 step
""")
    fname = f"prototype_2d_step{step}.in"
    with open(fname, "w") as f:
        f.write("\n".join(lines))
    return fname


def run_step(step):
    fname = write_input(step)
    outname = fname.replace(".in", ".out")
    subprocess.run([PHREEQC, fname, outname, "phreeqc.dat"], capture_output=True, text=True)
    with open(outname, errors="ignore") as f:
        content = f.read()
    errors = content.count("not converged") + content.count("Model failed")
    df = pd.read_csv("selected_output_1.sel", sep="\t")
    df.columns = [c.strip() for c in df.columns]
    return df, errors


def transport_step():
    """Advect in +x, disperse in y. Explicit, first-order, illustrative."""
    new_state = [[dict(state[x][y]) for y in range(NY)] for x in range(NX)]
    for key in FRESH:
        # advection: push ADVECT_FRAC of each cell's concentration into the next cell in x
        for y in range(NY):
            for x in range(NX - 1, 0, -1):
                incoming = ADVECT_FRAC * state[x - 1][y][key]
                new_state[x][y][key] = (1 - ADVECT_FRAC) * state[x][y][key] + incoming
            new_state[0][y][key] = (1 - ADVECT_FRAC) * state[0][y][key]
        # lateral dispersion: mix a bit with each y-neighbour
        base = [[new_state[x][y][key] for y in range(NY)] for x in range(NX)]
        for x in range(NX):
            for y in range(NY):
                total = base[x][y]
                n = 0
                if y > 0:
                    total += DISPERSE_FRAC * (base[x][y - 1] - base[x][y])
                    n += 1
                if y < NY - 1:
                    total += DISPERSE_FRAC * (base[x][y + 1] - base[x][y])
                    n += 1
                new_state[x][y][key] = total
    return new_state


total_errors = 0
for step in range(N_STEPS):
    df, errors = run_step(step)
    total_errors += errors
    for x in range(NX):
        for y in range(NY):
            cid = cell_id(x, y)
            row = df[df["soln"] == cid].iloc[-1]
            for key, col in [("Urea", "Urea"), ("Ca", "Ca"), ("Amm", "Amm"), ("C", "C")]:
                state[x][y][key] = max(row[col], 0.0)
            # k_Calcite is the CUMULATIVE moles present (we fed the previous
            # cumulative in as -m0), so this REPLACES cum_calcite rather than
            # adding to it -- adding was the bug that caused runaway growth.
            cum_calcite[x, y] = max(row["k_Calcite"], 0.0)
    state = transport_step()
    if step in (2, 9, 19):
        snapshots[step] = cum_calcite.copy()
    print(f"step {step:2d}  total errors so far: {total_errors}  "
          f"max calcite so far: {cum_calcite.max():.5f} mol")

print(f"\nDone. Total convergence errors across all {N_STEPS} steps x {NX*NY} cells: {total_errors}")

fig, axes = plt.subplots(1, len(snapshots), figsize=(5 * len(snapshots), 4.5))
vmax = max(s.max() for s in snapshots.values())
for ax, (step, grid) in zip(axes, snapshots.items()):
    im = ax.imshow(grid.T, origin="lower", cmap="viridis", vmin=0, vmax=vmax,
                    extent=[0, NX, 0, NY], aspect="equal")
    ax.set_title(f"step {step}")
    ax.set_xlabel("x (flow direction)")
    ax.plot(INJECT[0] + 0.5, INJECT[1] + 0.5, "r*", markersize=14, label="injection point")
axes[0].set_ylabel("y (lateral)")
axes[0].legend(loc="upper right", fontsize=8)
fig.colorbar(im, ax=axes, label="Cumulative calcite precipitated (mol)", shrink=0.85)
fig.suptitle("2D prototype: calcite plume spreading from a point injection", fontsize=12)
fig.savefig("prototype_2d_result.png", dpi=150, bbox_inches="tight")
print("saved prototype_2d_result.png")
