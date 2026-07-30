"""Same alpine temperature sweep (4-10 C) as run_alpine_sweep.py, but run twice:
once with Ea=31.2 kJ/mol (Wang et al. 2023) and once with Ea=51.7 kJ/mol
(Sun et al. 2019), so the two estimates can be compared directly on one plot.
Does NOT modify the original test_3B.in (25 C baseline).
"""
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

PHREEQC = "/mnt/c/home/cacchina_fra/iGEM/phreeqc_install/bin/phreeqc"

with open("test_3B.in", "r") as f:
    base_template = f.read()

old_rates = """11 cell_l = 5.2
12 t_d = 288000
15 initial_por = .38
16 Cb = (5.66e-1*(cell_no*cell_l-cell_l/2)^3+1.1e2*(cell_no*cell_l-cell_l/2)^2+4.67e3*(cell_no*cell_l-cell_l/2)-5.62e3)/1e6
18 decay_ratio = exp(-TOTAL_TIME/t_d)
20 rate = Cb*soln_vol*get_por(cell_no)/initial_por*PARM(1)*MOL("Urea")/(MOL("Urea")+PARM(2))"""

new_rates_template = """11 cell_l = 5.2
12 t_d = 288000
13 Ea = {ea}
14 Tref = 298.15
15 initial_por = .38
16 Cb = (5.66e-1*(cell_no*cell_l-cell_l/2)^3+1.1e2*(cell_no*cell_l-cell_l/2)^2+4.67e3*(cell_no*cell_l-cell_l/2)-5.62e3)/1e6
17 Tfactor = exp(-Ea/8.314*(1/(TC+273.15)-1/Tref))
18 decay_ratio = exp(-TOTAL_TIME/t_d)
20 rate = Cb*soln_vol*get_por(cell_no)/initial_por*PARM(1)*Tfactor*MOL("Urea")/(MOL("Urea")+PARM(2))"""

SOURCES = {
    "Wang 2023": {"ea": 31200, "tag": "wang", "color": "tab:blue"},
    "Sun 2019": {"ea": 51700, "tag": "sun", "color": "tab:green"},
}

all_results = {}

for label, info in SOURCES.items():
    template = base_template.replace(
        old_rates, new_rates_template.format(ea=info["ea"])
    )

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
        fname = f"test_3B_alpine_{info['tag']}_{T}.in"
        with open(fname, "w") as f:
            f.write(variant)

        outname = f"test_3B_alpine_{info['tag']}_{T}.out"
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
        print(f"[{label}] T={T}C  total calcite={total_calcite:.4f} mol  errors={n_errors}")

    all_results[label] = pd.DataFrame(results)
    all_results[label].to_csv(f"alpine_sweep_results_{info['tag']}.csv", index=False)

# --- baseline (25 C, unmodified original) ---
subprocess.run([PHREEQC, "test_3B.in", "test_3B.out", "phreeqc.dat"],
               capture_output=True, text=True)
df = pd.read_csv("selected_output_1.sel", sep="\t")
df.columns = [c.strip() for c in df.columns]
cells = df[(df["soln"] >= 1) & (df["soln"] <= 10)]
last_step = cells[cells["step"] == cells["step"].max()]
baseline_calcite = last_step["k_Calcite"].sum()
print(f"[baseline 25C] total calcite={baseline_calcite:.4f} mol")

# --- combined plot ---
fig, ax = plt.subplots(figsize=(7, 5))
for label, info in SOURCES.items():
    df_r = all_results[label]
    ax.plot(df_r["T"], df_r["total_calcite_mol"], "o-", color=info["color"],
             label=f"{label}: Ea={info['ea']/1000:.1f} kJ/mol")
ax.axhline(baseline_calcite, color="tab:red", linestyle="--",
           label="baseline laboratorio 25°C (originale)")
ax.set_xlabel("Temperatura (°C)")
ax.set_ylabel("Calcite totale precipitata a fine simulazione (mol)")
ax.set_title("Confronto: due stime diverse per il rallentamento al freddo")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("alpine_temperature_sweep_compare.png", dpi=150)
print("saved alpine_temperature_sweep_compare.png")
