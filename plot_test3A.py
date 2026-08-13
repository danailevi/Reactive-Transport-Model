"""Plot PHREEQC Test 3A (stop-flow) results, matching Figures 4-6 of Razbani et al. 2024."""
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("selected_output_1.sel", sep="\t")
df.columns = [c.strip() for c in df.columns]

# keep only the stop-flow phase (solution 8, the "D outlet" cell)
stop = df[df["soln"] == 8].copy()

# normalization constants (inlet concentrations), same as Razbani's own script
C_u_0 = 0.05
C_ca_0 = 0.05
C_amm_0 = 0.0567

t_hours = (stop["time"] - stop["time"].iloc[0]) / 3600.0

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

axes[0].plot(t_hours, stop["Ca"] / C_ca_0, "k:", label="[Calcium]")
axes[0].plot(t_hours, stop["Urea"] / C_u_0, "k-", label="[Urea]")
axes[0].plot(t_hours, stop["Amm"] / C_amm_0, "b-", label="[Ammonium]")
axes[0].set_xlabel("Tempo dall'inizio dello stop-flow (h)")
axes[0].set_ylabel("Concentrazione normalizzata")
axes[0].set_title("Test 3A - stop-flow, cella D (outlet)")
axes[0].legend(fontsize=8)
axes[0].grid(alpha=0.3)

axes[1].plot(t_hours, stop["pH"], "r-")
axes[1].set_xlabel("Tempo dall'inizio dello stop-flow (h)")
axes[1].set_ylabel("pH")
axes[1].set_title("pH durante lo stop-flow")
axes[1].grid(alpha=0.3)

fig.tight_layout()
fig.savefig("results_plot_3A.png", dpi=150)
print("saved results_plot_3A.png")
