"""Generate and run test_3B at several alpine temperatures (4-10 C), applying an
Arrhenius correction (Ea=31.2 kJ/mol, from Wang et al. 2023) to R_max in the
Ureolysis RATES block. Does NOT modify the original test_3B.in (25 C baseline).
Collects total calcite precipitated at the final step for each temperature and
plots the trend.
"""
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

PHREEQC = "/mnt/c/home/cacchina_fra/iGEM/phreeqc_install/bin/phreeqc"

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

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(results_df["T"], results_df["total_calcite_mol"], "o-", color="tab:blue")
ax.set_xlabel("Temperatura (°C)")
ax.set_ylabel("Calcite totale precipitata a fine simulazione (mol)")
ax.set_title("Effetto della temperatura sulla cementazione totale\n(Ea=31.2 kJ/mol, Wang 2023)")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("alpine_temperature_sweep.png", dpi=150)
print("saved alpine_temperature_sweep.png")
