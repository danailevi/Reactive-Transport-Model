"""Run the alpine-adapted Test 3B (R_max Arrhenius + t_d biomass-survival
correction) at a representative cold temperature, and plot the same kind of
calcite precipitation profile chart as plot_results.py (profile along the
column at 5 time points), and the original 25 C baseline, as two separate
figures that share a y-axis range so they stay directly comparable.

Scenario used: T = 4 C, R_max Ea = 31.2 kJ/mol (Wang 2023),
t_d survival Ea_d = 48 kJ/mol (Q10 ~ 2, the more conservative of the two
brackets tested in run_alpine_sweep_td.py).
"""
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

PHREEQC = r"C:\Program Files\USGS\phreeqc-3.8.6-17100-x64\bin\phreeqc.bat"
T_ALPINE = 4
EA_RMAX = 31200
EA_D = 48000

with open("../test_3B.in") as f:
    base_template = f.read()

old_rates = """11 cell_l = 5.2
12 t_d = 288000
15 initial_por = .38
16 Cb = (5.66e-1*(cell_no*cell_l-cell_l/2)^3+1.1e2*(cell_no*cell_l-cell_l/2)^2+4.67e3*(cell_no*cell_l-cell_l/2)-5.62e3)/1e6
18 decay_ratio = exp(-TOTAL_TIME/t_d)
20 rate = Cb*soln_vol*get_por(cell_no)/initial_por*PARM(1)*MOL("Urea")/(MOL("Urea")+PARM(2))"""

new_rates = f"""11 cell_l = 5.2
12 t_d0 = 288000
13 Ea = {EA_RMAX}
14 Ead = {EA_D}
15 Tref = 298.15
16 initial_por = .38
17 Cb = (5.66e-1*(cell_no*cell_l-cell_l/2)^3+1.1e2*(cell_no*cell_l-cell_l/2)^2+4.67e3*(cell_no*cell_l-cell_l/2)-5.62e3)/1e6
18 Tfactor = exp(-Ea/8.314*(1/(TC+273.15)-1/Tref))
19 t_d = t_d0 * exp(Ead/8.314*(1/(TC+273.15)-1/Tref))
20 decay_ratio = exp(-TOTAL_TIME/t_d)
21 rate = Cb*soln_vol*get_por(cell_no)/initial_por*PARM(1)*Tfactor*MOL("Urea")/(MOL("Urea")+PARM(2))"""

assert old_rates in base_template
alpine_template = base_template.replace(old_rates, new_rates)
alpine_template = alpine_template.replace(
    "-pH     7       charge_balance\nwater 38.7e-3",
    f"-pH     7       charge_balance\n-temp   {T_ALPINE}\nwater 38.7e-3",
).replace(
    "-pH        7       charge_balance\n     water 38.7e-3",
    f"-pH        7       charge_balance\n     -temp      {T_ALPINE}\n     water 38.7e-3",
)

with open("test_3B_alpine_profile.in", "w") as f:
    f.write(alpine_template)


def run_and_load(infile, outfile):
    subprocess.run([PHREEQC, infile, outfile, r"..\phreeqc.dat"],
                    capture_output=True, text=True)
    df = pd.read_csv("selected_output_1.sel", sep="\t")
    df.columns = [c.strip() for c in df.columns]
    with open(outfile, errors="ignore") as fh:
        content = fh.read()
    errors = content.count("not converged") + content.count("Model failed")
    return df, errors


df_25, err_25 = run_and_load(r"..\test_3B.in", "test_3B_baseline.out")
df_alpine, err_alpine = run_and_load("test_3B_alpine_profile.in", "test_3B_alpine_profile.out")
print(f"25C baseline: errors={err_25}")
print(f"4C alpine (R_max+t_d): errors={err_alpine}")

# shared y-axis range so the two separate figures stay visually comparable
cal_max = max(
    df[(df["soln"] >= 1) & (df["soln"] <= 10)]["k_Calcite"].max()
    for df in (df_25, df_alpine)
)

for df, title, fname in [
    (df_25, "25°C — lab baseline (original)", "alpine_profile_baseline.png"),
    (df_alpine, "4°C — alpine-adapted (R_max + t_d)", "alpine_profile_adapted.png"),
]:
    fig, ax = plt.subplots(figsize=(6, 5))
    cells = df[(df["soln"] >= 1) & (df["soln"] <= 10)]
    steps = sorted(cells["step"].unique())
    picks = [steps[i] for i in (0, len(steps)//4, len(steps)//2, 3*len(steps)//4, len(steps)-1)]
    for step in picks:
        sub = cells[cells["step"] == step].sort_values("dist_x") if "dist_x" in cells.columns else cells[cells["step"] == step]
        # dist_x may be named differently; fall back to soln order if absent
        x = sub["dist_x"] if "dist_x" in sub.columns else sub["soln"]
        ax.plot(x, sub["k_Calcite"], marker="o", markersize=3, label=f"step {step}")
    ax.set_xlabel("Distance along column (m)")
    ax.set_ylabel("Calcite precipitated (mol)")
    ax.set_ylim(0, cal_max * 1.05)
    ax.set_title(title)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(fname, dpi=150)
    print(f"saved {fname}")
