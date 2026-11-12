"""2D prototype v2: adds the advisor's temperature and pH guidance to the
operator-splitting architecture from prototype_2d.py (kept unmodified as the
"bare architecture" reference).

New physics added to the RATES -> Ureolysis BASIC code, on top of the same
Michaelis-Menten kinetics used throughout this project:

  1. Temperature (Arrhenius), same functional form as the 1D alpine sweep:
        Tfactor = exp(-Ea/R * (1/(TC+273.15) - 1/Tref))
     Ea = 41200 J/mol, the advisor's recommended MIDPOINT between the two
     literature brackets (Wang 2023: 31.2 kJ/mol, Sun 2019: 51.7 kJ/mol),
     rather than picking one source over the other.

  2. pH inhibition (new -- not present anywhere in the model before this):
        pH_factor = exp(-(PH-pH_opt)^2 / (2*sigma^2))
     A Gaussian "empirical inhibition factor" (per the advisor's explicit
     recommendation over a hard on/off threshold), pH_opt=8.0, calibrated
     so activity = 20% of maximum at pH 6.5 (Omoregie et al., single
     literature source -- see caveat below).

Two scenarios are run and compared, both at a uniform alpine temperature
of 8C (midpoint of the 4-12C target), both with identical transport:
  - "neutral": every cell's pH is charge-balanced (same treatment as
    prototype_2d.py and the whole 1D model) -- a non-acidic control.
  - "acidic": every cell's pH is FIXED at 5.5 (representing acidic alpine
    soil pore water) instead of charge-balanced. Applied grid-wide, not
    just at the injection point, because the simplified Python-side
    transport only carries Urea/Ca/Amm/Cl/C/Na between cells -- an
    acidity signature set only at the injection cell would not propagate
    downstream in this simplified bookkeeping. A grid-wide fixed soil pH
    is also the more physically appropriate choice: acidity here
    represents a property of the alpine soil water the fluid moves
    through, not just the injected fluid itself.

CAVEAT ON THE pH CURVE: pH_opt=8.0 and the 20%-at-6.5 anchor come from a
single literature source (Omoregie et al., cultural-conditions study on
S. pasteurii) found via search-engine summary; a second source was sought
to cross-check it but was blocked by a CAPTCHA and not obtained. Treat the
resulting pH_factor values at 4.5-6 (the actual alpine target range, far
from the pH 6.5-9 the source data covers) as an illustrative extrapolation,
not a validated prediction -- this needs a directly-sourced pH-activity
curve (ideally with data points inside 4.5-6 itself) before being trusted
quantitatively.

Run from inside this folder. The PHREEQC command defaults to phreeqc.bat
on Windows and phreeqc elsewhere; set the PHREEQC environment variable to
override it. The database is the patched ../phreeqc.dat in the repo root.
"""
import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PHREEQC = os.environ.get("PHREEQC", "phreeqc.bat" if os.name == "nt" else "phreeqc")
DATABASE = os.path.join("..", "phreeqc.dat")

NX, NY = 5, 5
INJECT = (0, 2)
N_STEPS = 20
DT = 1055
ADVECT_FRAC = 0.5
DISPERSE_FRAC = 0.12

RMAX, KM = 7.0e-5, 0.0173
CB = 3.0e-2
K1, K2, AREA = 1.55e-6, 0.501, 1.0

ALPINE_TEMP_C = 8          # uniform grid temperature, midpoint of 4-12C target
EA_R = 41200                # J/mol, advisor's midpoint between Wang/Sun
TREF = 298.15
PH_OPT = 8.0
SIGMA_PH = 0.836             # fit to 20% activity at pH 6.5 (Omoregie et al.)

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
5  Ea = %(ea)s
6  Tref = %(tref)s
7  Tfactor = exp(-Ea/8.314*(1/(TC+273.15)-1/Tref))
8  pH_opt = %(ph_opt)s
9  sigma_pH = %(sigma_ph)s
10 pH_factor = exp(-(PH-pH_opt)^2/(2*sigma_pH^2))
15 rate = PARM(3)*soln_vol*PARM(1)*Tfactor*pH_factor*MOL("Urea")/(MOL("Urea")+PARM(2))
20 moles = rate * TIME
30 SAVE moles
-END

Calcite
-START
10 rate = PARM(3)*(PARM(1)+PARM(2)*act("H+"))*(1-10^SI("Calcite"))
20 moles = rate * TIME
30 SAVE moles
-END
""" % dict(ea=EA_R, tref=TREF, ph_opt=PH_OPT, sigma_ph=SIGMA_PH)

INLET = dict(Urea=0.05, Amm=0.0567, Ca=0.05, Cl=0.1567, C=0.0038, Na=0.0038)
FRESH = dict(Urea=0.0, Amm=0.0, Ca=0.0, Cl=0.0, C=0.0, Na=0.0)

SELECTED_OUTPUT_BLOCK = """SELECTED_OUTPUT
    -step true
    -soln true
    -pH true
    -tot Ca Urea Amm C
    -k Ureolysis Calcite
END
"""


def cell_id(x, y):
    return x * NY + y + 1


def write_input(step, scenario, tag):
    ph_line = "-pH 7 charge_balance" if scenario == "neutral" else "-pH 5.5"
    lines = [SPECIES_BLOCK, SELECTED_OUTPUT_BLOCK]
    for x in range(NX):
        for y in range(NY):
            cid = cell_id(x, y)
            c = INLET if (x, y) == INJECT else STATE[x][y]
            lines.append(f"""SOLUTION {cid}
Urea {c['Urea']:.6e}
Amm {c['Amm']:.6e}
Ca {c['Ca']:.6e}
Cl {c['Cl']:.6e}
C {c['C']:.6e}
Na {c['Na']:.6e}
-units mol/kgw
{ph_line}
-temp {ALPINE_TEMP_C}
water 38.7e-3
END
""")
            lines.append(f"""KINETICS {cid}
    Ureolysis
        -m0 100
        -formula AmmH+ 2 CO3-2 1 Urea -1 H2O -2
        -parms {RMAX} {KM} {CB}
    Calcite
        -m0 {CUM_CALCITE[x, y]:.6e}
        -parms {K1} {K2} {AREA}
-steps {DT} seconds in 1 step
""")
    fname = f"prototype_2d_v2_{tag}_step{step}.in"
    with open(fname, "w") as f:
        f.write("\n".join(lines))
    return fname


def run_step(step, scenario, tag):
    fname = write_input(step, scenario, tag)
    outname = fname.replace(".in", ".out")
    subprocess.run([PHREEQC, fname, outname, DATABASE], capture_output=True, text=True)
    with open(outname, errors="ignore") as f:
        content = f.read()
    errors = content.count("not converged") + content.count("Model failed")
    df = pd.read_csv("selected_output_1.sel", sep="\t")
    df.columns = [c.strip() for c in df.columns]
    return df, errors


def transport_step():
    new_state = [[dict(STATE[x][y]) for y in range(NY)] for x in range(NX)]
    for key in FRESH:
        for y in range(NY):
            for x in range(NX - 1, 0, -1):
                incoming = ADVECT_FRAC * STATE[x - 1][y][key]
                new_state[x][y][key] = (1 - ADVECT_FRAC) * STATE[x][y][key] + incoming
            new_state[0][y][key] = (1 - ADVECT_FRAC) * STATE[0][y][key]
        base = [[new_state[x][y][key] for y in range(NY)] for x in range(NX)]
        for x in range(NX):
            for y in range(NY):
                total = base[x][y]
                if y > 0:
                    total += DISPERSE_FRAC * (base[x][y - 1] - base[x][y])
                if y < NY - 1:
                    total += DISPERSE_FRAC * (base[x][y + 1] - base[x][y])
                new_state[x][y][key] = total
    return new_state


results = {}
for scenario, tag in [("neutral", "neu"), ("acidic", "acid")]:
    STATE = [[dict(FRESH) for _ in range(NY)] for _ in range(NX)]
    CUM_CALCITE = np.zeros((NX, NY))
    total_errors = 0
    for step in range(N_STEPS):
        df, errors = run_step(step, scenario, tag)
        total_errors += errors
        for x in range(NX):
            for y in range(NY):
                cid = cell_id(x, y)
                row = df[df["soln"] == cid].iloc[-1]
                for key, col in [("Urea", "Urea"), ("Ca", "Ca"), ("Amm", "Amm"), ("C", "C")]:
                    STATE[x][y][key] = max(row[col], 0.0)
                CUM_CALCITE[x, y] = max(row["k_Calcite"], 0.0)
        STATE = transport_step()
        print(f"[{scenario}] step {step:2d}  errors so far: {total_errors}  "
              f"max calcite so far: {CUM_CALCITE.max():.6f} mol")
    results[scenario] = CUM_CALCITE.copy()
    print(f"[{scenario}] DONE. total errors = {total_errors}, "
          f"total calcite = {CUM_CALCITE.sum():.6f} mol, max cell = {CUM_CALCITE.max():.6f} mol\n")

# --- comparison plot ---
fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
vmax = max(g.max() for g in results.values())
for ax, scenario in zip(axes, ["neutral", "acidic"]):
    im = ax.imshow(results[scenario].T, origin="lower", cmap="viridis", vmin=0, vmax=vmax,
                    extent=[0, NX, 0, NY], aspect="equal")
    ax.set_title(f"{scenario} pH scenario\n(total calcite = {results[scenario].sum():.4f} mol)")
    ax.set_xlabel("x (flow direction)")
    ax.plot(INJECT[0] + 0.5, INJECT[1] + 0.5, "r*", markersize=14)
axes[0].set_ylabel("y (lateral)")
fig.colorbar(im, ax=axes, label="Cumulative calcite precipitated (mol)", shrink=0.85)
fig.suptitle(f"2D prototype v2: temperature (Ea={EA_R/1000:.1f} kJ/mol, {ALPINE_TEMP_C}C) "
             f"+ pH inhibition, neutral vs. acidic alpine soil", fontsize=11)
fig.tight_layout()
fig.savefig("prototype_2d_v2_ph_comparison.png", dpi=150)
print("saved prototype_2d_v2_ph_comparison.png")

ratio = results["acidic"].sum() / results["neutral"].sum() * 100
print(f"\nAcidic-scenario total calcite is {ratio:.1f}% of the neutral-scenario total.")
