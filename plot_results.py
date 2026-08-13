"""Plot PHREEQC transport results from selected_output_1.sel (test_3B.in)."""
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("selected_output_1.sel", sep="\t")
df.columns = [c.strip() for c in df.columns]

# keep only the 10 transported cells (drop boundary solution 0 and any extra row)
cells = df[(df["soln"] >= 1) & (df["soln"] <= 10)].copy()

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

# 1) CaCO3 (calcite) precipitation profile along the column, at several times
steps_to_plot = sorted(cells["step"].unique())
picks = [steps_to_plot[i] for i in (0, len(steps_to_plot)//4, len(steps_to_plot)//2,
                                     3*len(steps_to_plot)//4, len(steps_to_plot)-1)]
for step in picks:
    sub = cells[cells["step"] == step].sort_values("dist_x")
    axes[0].plot(sub["dist_x"], sub["k_Calcite"], marker="o", markersize=3,
                 label=f"step {step}")
axes[0].set_xlabel("Distance along column (m)")
axes[0].set_ylabel("Calcite precipitated (mol)")
axes[0].set_title("CaCO3 precipitation profile")
axes[0].legend(fontsize=8)

# 2) pH profile at the final step
last = cells[cells["step"] == cells["step"].max()].sort_values("dist_x")
axes[1].plot(last["dist_x"], last["pH"], marker="o", color="tab:red")
axes[1].set_xlabel("Distance along column (m)")
axes[1].set_ylabel("pH")
axes[1].set_title(f"pH profile at final step ({int(cells['step'].max())})")

fig.tight_layout()
fig.savefig("results_plot.png", dpi=150)
print("saved results_plot.png")
