"""Alpine temperature sweep with a biomass-survival (t_d) correction.

Motivation
----------
The earlier sweep (run_alpine_sweep_compare.py) only slowed the *rate* of
ureolysis in the cold via an Arrhenius factor on R_max. That reproduces the
intuition "cold -> slower -> less cement", but it contradicts the experimental
finding of Peng, Liu & Achal (2019): sand columns treated at 10 C end up
*stronger* than at 30 C, because the bacteria stay active far longer in the cold.

In the model, "how long the bacteria stay active" is the decay time t_d, used as
    decay_ratio = exp(-TOTAL_TIME / t_d)
which multiplies the ureolysis rate. Until now t_d was a fixed 288000 s at every
temperature. Biological decay is itself a temperature-dependent process, so it
should slow down in the cold exactly like every other rate — i.e. t_d should get
*longer* as temperature drops. We apply the same Arrhenius form, but inverted:

    t_d(T) = t_d0 * exp( +Ea_d/R * (1/T - 1/Tref) )

With Ea_d = 0 this reduces to the constant-t_d case, i.e. the previous
R_max-only sweep, which we keep as the reference curve.

Scenarios
---------
R_max Arrhenius is fixed to Wang 2023 (Ea = 31.2 kJ/mol) in all runs, so the only
thing changing between curves is the t_d correction:
    * Ea_d = 0      -> R_max only  (reference, = previous behaviour)
    * Ea_d = 48 kJ/mol -> decay slows with Q10 ~= 2
    * Ea_d = 76 kJ/mol -> decay slows with Q10 ~= 3

Q10 ~= 2-3 is the usual range for microbial maintenance/decay. The two values
bracket a plausible response rather than claiming one exact number.

Does NOT modify the original test_3B.in (25 C baseline).
"""
import subprocess
import math
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PHREEQC = "/home/cacchina_fra/iGEM/phreeqc_install/bin/phreeqc"
EA_RMAX = 31200          # Wang 2023, J/mol  (fixed for all runs)
TOTAL_TIME = 178 * 1055  # shifts * time_step, s  (for reporting decay at end)

with open("test_3B.in") as f:
    base_template = f.read()

old_rates = """11 cell_l = 5.2
12 t_d = 288000
15 initial_por = .38
16 Cb = (5.66e-1*(cell_no*cell_l-cell_l/2)^3+1.1e2*(cell_no*cell_l-cell_l/2)^2+4.67e3*(cell_no*cell_l-cell_l/2)-5.62e3)/1e6
18 decay_ratio = exp(-TOTAL_TIME/t_d)
20 rate = Cb*soln_vol*get_por(cell_no)/initial_por*PARM(1)*MOL("Urea")/(MOL("Urea")+PARM(2))"""

new_rates = """11 cell_l = 5.2
12 t_d0 = 288000
13 Ea = {ea}
14 Ead = {ead}
15 Tref = 298.15
16 initial_por = .38
17 Cb = (5.66e-1*(cell_no*cell_l-cell_l/2)^3+1.1e2*(cell_no*cell_l-cell_l/2)^2+4.67e3*(cell_no*cell_l-cell_l/2)-5.62e3)/1e6
18 Tfactor = exp(-Ea/8.314*(1/(TC+273.15)-1/Tref))
19 t_d = t_d0 * exp(Ead/8.314*(1/(TC+273.15)-1/Tref))
20 decay_ratio = exp(-TOTAL_TIME/t_d)
21 rate = Cb*soln_vol*get_por(cell_no)/initial_por*PARM(1)*Tfactor*MOL("Urea")/(MOL("Urea")+PARM(2))"""

SCENARIOS = {
    "R_max only (t_d const)": {"ead": 0,     "tag": "td0",  "color": "tab:gray",  "ls": "--"},
    "+ t_d survival (Q10~2)": {"ead": 48000, "tag": "td48", "color": "tab:blue",  "ls": "-"},
    "+ t_d survival (Q10~3)": {"ead": 76000, "tag": "td76", "color": "tab:green", "ls": "-"},
}

def run_variant(template, T):
    variant = template.replace(
        "-pH     7       charge_balance\nwater 38.7e-3",
        f"-pH     7       charge_balance\n-temp   {T}\nwater 38.7e-3",
    ).replace(
        "-pH        7       charge_balance\n     water 38.7e-3",
        f"-pH        7       charge_balance\n     -temp      {T}\n     water 38.7e-3",
    )
    fname = f"/tmp/td_{T}.in"
    outname = f"/tmp/td_{T}.out"
    with open(fname, "w") as f:
        f.write(variant)
    subprocess.run([PHREEQC, fname, outname, "phreeqc.dat"],
                   capture_output=True, text=True)
    df = pd.read_csv("selected_output_1.sel", sep="\t")
    df.columns = [c.strip() for c in df.columns]
    cells = df[(df["soln"] >= 1) & (df["soln"] <= 10)]
    last = cells[cells["step"] == cells["step"].max()]
    with open(outname, errors="ignore") as fh:
        content = fh.read()
    errors = content.count("not converged") + content.count("Model failed")
    return last["k_Calcite"].sum(), errors

# baseline: unmodified original at 25 C
subprocess.run([PHREEQC, "test_3B.in", "/tmp/base.out", "phreeqc.dat"],
               capture_output=True, text=True)
dfb = pd.read_csv("selected_output_1.sel", sep="\t")
dfb.columns = [c.strip() for c in dfb.columns]
cb = dfb[(dfb["soln"] >= 1) & (dfb["soln"] <= 10)]
baseline = cb[cb["step"] == cb["step"].max()]["k_Calcite"].sum()
print(f"[baseline 25 C, unmodified] total calcite = {baseline:.4f} mol\n")

all_results = {}
for label, info in SCENARIOS.items():
    template = base_template.replace(
        old_rates, new_rates.format(ea=EA_RMAX, ead=info["ead"]))
    rows = []
    for T in range(4, 11):
        calcite, errors = run_variant(template, T)
        td = 288000 * math.exp(info["ead"] / 8.314 * (1 / (T + 273.15) - 1 / 298.15))
        active_end = math.exp(-TOTAL_TIME / td)
        rows.append({"T": T, "total_calcite_mol": calcite,
                     "pct_of_25C": 100 * calcite / baseline,
                     "t_d_s": td, "active_frac_end": active_end, "errors": errors})
        print(f"[{label}] T={T:2d}C  calcite={calcite:.4f} mol "
              f"({100*calcite/baseline:4.0f}% of 25C)  t_d={td/1e3:6.0f}k s  "
              f"active@end={active_end:.2f}  err={errors}")
    df = pd.DataFrame(rows)
    df.to_csv(f"alpine_td_{info['tag']}.csv", index=False)
    all_results[label] = df
    print()

# ---- plot ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.6))
for label, info in SCENARIOS.items():
    df = all_results[label]
    ax1.plot(df["T"], df["total_calcite_mol"], marker="o",
             color=info["color"], ls=info["ls"], label=label)
    ax2.plot(df["T"], df["pct_of_25C"], marker="o",
             color=info["color"], ls=info["ls"], label=label)
ax1.axhline(baseline, color="black", ls=":", lw=1, label="25 C baseline")
ax1.set_xlabel("Temperature (C)"); ax1.set_ylabel("Total calcite at end (mol)")
ax1.set_title("Cement produced vs temperature"); ax1.legend(fontsize=8); ax1.grid(alpha=.3)
ax2.axhline(100, color="black", ls=":", lw=1)
ax2.set_xlabel("Temperature (C)"); ax2.set_ylabel("% of 25 C cementation")
ax2.set_title("Cold performance relative to lab temperature")
ax2.legend(fontsize=8); ax2.grid(alpha=.3)
fig.suptitle("Effect of a temperature-dependent biomass-survival time (t_d)\n"
             "R_max Arrhenius fixed to Wang 2023 (Ea=31.2 kJ/mol)", fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig("alpine_td_survival_sweep.png", dpi=130)
print("Saved plot -> alpine_td_survival_sweep.png")
