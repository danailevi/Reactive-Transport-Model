"""Plot PHREEQC Test 3A (stop-flow) results, matching Figures 4-6 of Razbani et al. 2024."""
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

PHREEQC = r"C:\Program Files\USGS\phreeqc-3.8.6-17100-x64\bin\phreeqc.bat"

# regenerate selected_output_1.sel from Test 3A itself: it is a generic filename
# that any other .in run (e.g. the alpine 3B profile) overwrites, so never trust
# whatever is currently on disk.
subprocess.run([PHREEQC, r"..\test_3A.in", "test_3A.out", r"..\phreeqc.dat"],
               capture_output=True, text=True)

df = pd.read_csv("selected_output_1.sel", sep="\t")
df.columns = [c.strip() for c in df.columns]

# keep only the stop-flow phase (solution 8, the "D outlet" cell)
stop = df[df["soln"] == 8].copy()

# normalization constants (inlet concentrations), same as Razbani's own script
C_u_0 = 0.05
C_ca_0 = 0.05
C_amm_0 = 0.0567

t_hours = (stop["time"] - stop["time"].iloc[0]) / 3600.0

fig_c, ax_c = plt.subplots(figsize=(6, 4.5))
ax_c.plot(t_hours, stop["Ca"] / C_ca_0, "k:", label="[Calcium]")
ax_c.plot(t_hours, stop["Urea"] / C_u_0, "k-", label="[Urea]")
ax_c.plot(t_hours, stop["Amm"] / C_amm_0, "b-", label="[Ammonium]")
ax_c.set_xlabel("Time of stop-flow (h)")
ax_c.set_ylabel("Normalized concentration")
ax_c.set_title("Test 3A - stop-flow, cell D (outlet)")
ax_c.legend(fontsize=8)
ax_c.grid(alpha=0.3)
fig_c.tight_layout()
fig_c.savefig("results_plot_3A_concentration.png", dpi=150)

fig_ph, ax_ph = plt.subplots(figsize=(6, 4.5))
ax_ph.plot(t_hours, stop["pH"], "r-")
ax_ph.set_xlabel("Time of stop-flow (h)")
ax_ph.set_ylabel("pH")
ax_ph.set_title("pH during stop-flow")
ax_ph.grid(alpha=0.3)
fig_ph.tight_layout()
fig_ph.savefig("results_plot_3A_pH.png", dpi=150)

print("saved results_plot_3A_concentration.png and results_plot_3A_pH.png")
