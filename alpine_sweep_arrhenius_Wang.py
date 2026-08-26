"""Generate and run test_3B at several alpine temperatures (4-10 C), applying an
Arrhenius correction (Ea=31.2 kJ/mol, from Wang et al. 2023) to R_max in the
Ureolysis RATES block. Does NOT modify the original test_3B.in (25 C baseline).
Collects total calcite precipitated at the final step for each temperature and
plots the trend.
"""
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

PHREEQC = r"C:\Program Files\USGS\phreeqc-3.8.6-17100-x64\bin\phreeqc.bat"

with open("test_3B.in", "r") as f:
    template = f.read()

# --- 1. add the Arrhenius temperature correction to the Ureolysis rate ---
old_rates = """11 cell_l = 5.2
12 t_d = 288000
15 initial_por = .38
16 Cb = (5.66e-1*(cell_no*cell_l-cell_l/2)^3+1.1e2*(cell_no*cell_l-cell_l/2)^2+4.67e3*(cell_no*cell_l-cell_l/2)-5.62e3)/1e6
18 decay_ratio = exp(-TOTAL_TIME/t_d)
20 rate = Cb*soln_vol*get_por(cell_no)/initial_por*PARM(1)*MOL("Urea")/(MOL("Urea")+PARM(2))"""

new_rates = """11 cell_l = 5.2
12 t_d = 288000
13 Ea = 31200
14 Tref = 298.15
15 initial_por = .38
16 Cb = (5.66e-1*(cell_no*cell_l-cell_l/2)^3+1.1e2*(cell_no*cell_l-cell_l/2)^2+4.67e3*(cell_no*cell_l-cell_l/2)-5.62e3)/1e6
17 Tfactor = exp(-Ea/8.314*(1/(TC+273.15)-1/Tref))
18 decay_ratio = exp(-TOTAL_TIME/t_d)
20 rate = Cb*soln_vol*get_por(cell_no)/initial_por*PARM(1)*Tfactor*MOL("Urea")/(MOL("Urea")+PARM(2))"""

assert old_rates in template, "RATES block not found as expected"
template = template.replace(old_rates, new_rates)

results = []
for T in range(4, 11):
    variant = template.replace(
        '-pH     7       charge_balance\nwater 38.7e-3',
        f'-pH     7       charge_balance\n-temp   {T}\nwater 38.7e-3',
    )
    variant = variant.replace(
        '-pH        7       charge_balance\n     water 38.7e-3',
        f'-pH        7       charge_balance\n     -temp      {T}\n     water 38.7e-3',
    )
    fname = f"test_3B_alpine_{T}.in"
    with open(fname, "w") as f:
        f.write(variant)

    outname = f"test_3B_alpine_{T}.out"
    subprocess.run(
        [PHREEQC, fname, outname, "phreeqc.dat"],
        capture_output=True, text=True,
    )

    df = pd.read_csv("selected_output_1.sel", sep="\t")
    df.columns = [c.strip() for c in df.columns]
    cells = df[(df["soln"] >= 1) & (df["soln"] <= 10)]
    last_step = cells[cells["step"] == cells["step"].max()]
    total_calcite = last_step["k_Calcite"].sum()
    n_errors = 0
    with open(outname, errors="ignore") as f:
        content = f.read()
        n_errors = content.count("not converged") + content.count("Model failed")
    results.append({"T": T, "total_calcite_mol": total_calcite, "errors": n_errors})
    print(f"T={T}C  total calcite={total_calcite:.4f} mol  errors={n_errors}")

results_df = pd.DataFrame(results)
results_df.to_csv("alpine_sweep_results.csv", index=False)

BLUE_LINE   = "#2a78d6"   # sequential step 450
BLUE_FILL   = "#cde2fb"   # sequential step 100
INK_PRIMARY = "#0b0b0b"
INK_SECOND  = "#52514e"
INK_MUTED   = "#898781"
GRID        = "#e1e0d9"
SURFACE     = "#fcfcfb"

plt.rcParams["font.family"] = ["Segoe UI", "DejaVu Sans", "sans-serif"]

fig, ax = plt.subplots(figsize=(7.2, 6.4), dpi=200)
fig.patch.set_facecolor(SURFACE)
ax.set_facecolor(SURFACE)

x = results_df["T"]
y = results_df["total_calcite_mol"]

ax.fill_between(x, 0, y, color=BLUE_FILL, zorder=1)
ax.plot(x, y, "-", color=BLUE_LINE, lw=2.4, zorder=2)
ax.plot(x, y, "o", color=BLUE_LINE, ms=9,
        markeredgecolor=SURFACE, markeredgewidth=1.6, zorder=3)

ax.set_xlim(x.min(), x.max())
ax.set_ylim(0, y.max() * 1.15)

ax.set_xlabel("Temperature (°C)", fontsize=12, color=INK_SECOND, labelpad=10)
ax.set_ylabel("Total calcite precipitation (mol)", fontsize=12, color=INK_SECOND, labelpad=10)
ax.tick_params(colors=INK_MUTED, labelsize=10.5, length=0)

ax.grid(True, color=GRID, lw=1, zorder=0)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)
for spine in ("left", "bottom"):
    ax.spines[spine].set_color(GRID)
    ax.spines[spine].set_linewidth(1.2)

fig.suptitle("Effect of Temperature on Total Cementation", fontsize=17, fontweight="bold",
             color=INK_PRIMARY, x=0.065, ha="left", y=0.975)
ax.set_title("Alpine-range sweep, Arrhenius-corrected Rₘₐₓ (Ea=31.2 kJ/mol, Wang 2023)",
             fontsize=11.5, color=INK_MUTED, loc="left", pad=14)

plt.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("alpine_temperature_sweep.png", dpi=300, facecolor=SURFACE)
print("saved alpine_temperature_sweep.png")
