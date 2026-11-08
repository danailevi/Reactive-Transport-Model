# MICP Reactive Transport Model

A 1D reactive transport model of microbially induced calcium carbonate precipitation (MICP) in a sand column, built in [PHREEQC](https://www.usgs.gov/software/phreeqc-version-3). The baseline model reproduces Razbani et al.'s published Test 3B experiment; the rest of this repository extends it with sensitivity analysis, concentration variation, low-temperature (alpine) adaptation, a mechanistic biomass term, and a 2D prototype.

Full write-up: [iGEM UZH Wiki — Reactive Transport Model](https://2026.igem.wiki/uzurich/reactive-transport-model)

## Repository structure

| Path | Contents |
|---|---|
| `test_3A.in`, `test_3B.in` | Canonical PHREEQC input files (stop-flow and continuous-flow experiments) that everything else is built from or validated against. |
| `phreeqc.dat` | PHREEQC thermodynamic database, patched with the `AmmH+` species (see below). |
| `MICP_Simulation_T3B.py` | Main model: Test 3B reproduction with alpine soil properties (silty sand) incorporated. |
| `Razbani reproduction/` | Reproduction of Razbani et al.'s original Test 3A/3B results and a grid-convergence check — the validation baseline everything else builds on. |
| `Sensitivity/` | ±20% sensitivity sweeps for the three most influential parameters: attached biomass (`Cb`), half-saturation constant (`Km`), and biomass decay time (`t_d`). |
| `concentration variation/` | Urea and CaCl₂ concentration sweeps, and an initial-pH sweep. |
| `alpine temperature/` | Low-temperature (4–10°C) adaptation: Arrhenius correction on `R_max`, temperature-dependent biomass survival, and literature comparison. |
| `biomass/` | Mechanistic (colloid filtration theory) refit of the attached-biomass profile, as an alternative to the original cubic fit. |
| `2D prototype/` | Early 2D extension of the model (PHREEQC input decks only, stepped manually). |

## Requirements

- [PHREEQC 3](https://www.usgs.gov/software/phreeqc-version-3) (tested with 3.8.6-17100), both the IPhreeqcCOM interface and the standalone CLI (`phreeqc.bat`)
- [phreeqpy](https://www.phreeqpy.com/) (Python interface to IPhreeqcCOM)
- `numpy`, `pandas`, `matplotlib`, `scipy`

Scripts that drive PHREEQC through `phreeqpy`/IPhreeqcCOM hardcode the database path as installed by the Windows PHREEQC installer (`C:\Program Files\USGS\IPhreeqcCOM <version>\database\phreeqc.dat`); scripts that shell out to the CLI hardcode `phreeqc.bat`'s install path the same way. Update those paths at the top of a script if your install location differs.

## Running the scripts

Each subfolder is self-contained and expects to be run **from inside that folder**, e.g.:
```
cd "Sensitivity"
python MICP_Simulation_T3B_Cb_sensitivity.py
```
Scripts that reference the root-level `test_3A.in`/`test_3B.in`/`phreeqc.dat` (anything that shells out to `phreeqc.bat`) do so via a relative `..\` path — this only works with backslashes, since `phreeqc.bat`'s internal `copy` step doesn't resolve forward-slash `../` paths correctly on Windows. Each script writes its own output (`.csv`, `.png`) into its own folder.

## The physical picture

- Sand column, injected from the bottom with **urea + Ca²⁺ + chloride** and from the top with a solution rich in **bacteria** (*S. pasteurii*), which settle and remain attached to the grains.
- The bacteria hydrolyze urea (**ureolysis**) → producing ammonium and carbonate.
- The carbonate meets calcium → precipitates **calcite (CaCO₃)**, which cements the sand.
- The calcite that deposits **reduces porosity**, slowing down the subsequent reaction at that point (feedback).

**Overall reaction:**
$$
\mathrm{CO(NH_2)_2} + 2\,\mathrm{H_2O} \longrightarrow 2\,\mathrm{NH_4^+} + \mathrm{CO_3^{2-}}
$$

## Structure of the `test_3B.in` code

| Block | What it does | Lines |
|---|---|---|
| `SOLUTION_MASTER_SPECIES` / `SOLUTION_SPECIES` | "Vocabulary": teaches PHREEQC what urea is and how ions like Cl⁻/Ca²⁺ behave | 1–20 |
| `RATES → Ureolysis` | Rate formula for urea breakdown (depends on biomass and available urea) | 22–33 |
| `RATES → Calcite` | Precipitation rate formula + porosity update | 35–47 |
| `SOLUTION 0` | Recipe of the injected liquid (urea, Ca, ammonium, pH) | 49–59 |
| `SOLUTION 1-10` | What's in the tube at the start (only water) | 61–65 |
| `KINETICS 1-10` | Turns on both reactions in each of the column's 10 slices | 67–75 |
| `TRANSPORT` | Makes the liquid flow: 10 cells, 178 shifts, ~2 simulated days total | 77–86 |
| `SELECTED_OUTPUT` | What to write to file for plotting afterwards | 88–96 |

## The main formulas

**Ureolysis kinetics (Michaelis-Menten, weighted by biomass):**
$$
R_u = C_b \, V_{soln} \, \frac{\phi(x)}{\phi_0} \, R_{max} \, \frac{[\mathrm{Urea}]}{[\mathrm{Urea}] + K_m}
$$

**Biomass decay over time:**
$$
D(t) = \exp\left(-\frac{t}{t_d}\right)
$$

**Calcite precipitation/dissolution:**
$$
R_p = S\left(k_1 + k_2\, a_{\mathrm{H^+}}\right)\left(1 - \Omega\right), \qquad \Omega = \frac{a_{CO_3^{2-}} \cdot a_{Ca^{2+}}}{K}
$$

**Reactive area (shrinks as pores clog up):**
$$
S = S_0 \left(\frac{\phi}{\phi_0}\right)^{p}
$$

## Parameter table (Test 3B, from the paper — Table 3)

| Parameter | Value | Source | Corresponds to |
|---|---|---|---|
| $R_{max}$ | $7\times10^{-5}$ mol/L/s/OD | Wang & Nackenhorst (2020) | `PARM(1)` in `Kinetics → Ureolysis` |
| $K_m$ | 0.0173 mol/kgH₂O (17.3 mM) | Krajewska (2009) | `PARM(2)` |
| $S_0$ | 1 m²/L | Martinez et al. (2014) | `-parms` in `Kinetics → Calcite` |
| $p$ | 2 | assumed | `pow_a` in the `RATES` code |
| $t_d$ | 288,000 s | Fauriel & Laloui (2012) | `t_d` in the `RATES` code |
| Dispersivity | 0.005 m | Wang & Nackenhorst (2020) | ⚠️ our file has `0.01` — to be checked |

> [!warning] Note on K_m
> There's also a more recent value, $K_m = 305$ mM (Lauchnor et al. 2015, whole-cell kinetics), much higher than 17.3 mM. Razbani chose 17.3 mM to stay faithful to Martinez et al.'s original calibration. This isn't an error in our file, but it's a point worth discussing/citing if we decide to update it.

## Sensitivity analysis — where to focus effort

From the paper:
- **Low sensitivity**: initial reactive area $S_0$ (above a ~5 m²/L threshold) and exponent $p$ — precipitation isn't the limiting factor, so refining these parameters doesn't change the result much.
- **Sensitive**: $K_m$ (±40% → significant differences), biomass density (±20%), $t_d$ (±30%).

Confirmed by our own ±20% sweeps in `Sensitivity/` (peak calcite shift from baseline):
- $C_b$ (attached biomass): −8.0% / +5.4% — the most influential parameter
- $K_m$: −4.4% / +5.0%
- $t_d$: −4.5% / +2.9%

**Practical conclusion**: when adapting the model for the Alps, it's worth focusing on how cold changes **ureolysis ($R_{max}$, $K_m$)** and **bacterial survival ($t_d$)**, not precipitation kinetics — and future lab time is best spent re-measuring the attached-biomass profile, since $C_b$ dominates the sensitivity.

## Adaptations beyond the original Razbani et al. model

- **Sensitivity analysis** (`Sensitivity/`) — quantifying how much each kinetic parameter actually matters.
- **Concentration variation** (`concentration variation/`) — urea/CaCl₂ sweeps at wet-lab-relevant concentrations, and an initial-pH sweep.
- **Low (alpine) temperature adaptation** (`alpine temperature/`) — an Arrhenius correction on $R_{max}$ (the `TC` variable, available in PHREEQC's `RATES` BASIC) plus temperature-dependent biomass survival; not considered in the original paper at all.
- **Alpine soil properties** (silty sand) — incorporated directly into `MICP_Simulation_T3B.py`.
- **Mechanistic attached-biomass term** (`biomass/`) — a colloid-filtration-theory refit of $C_b(x)$, as an alternative to the original unconstrained cubic fit.
- **2D model prototype** (`2D prototype/`) — early extension beyond the original paper's 1D column.

Still outside scope: unsaturated flow conditions (would need the Faeli et al. approach) and the acidic pH typical of alpine soils (4.5–6), both left for future work.

## Bug fixed in our `phreeqc.dat`

The `AmmH+ = Amm + H+` reaction didn't have its own `-log_k` (it was missing, so it was read as 0 instead of the real pKa ~9.25). Fixed by adding:
```
AmmH+ = Amm + H+
    -log_k -9.252
    -delta_h 12.48 kcal
    -analytic 0.6322 -0.001225 -2835.76
```
This resolved the "Model failed to converge" error we were getting every time.
## Useful literature for the temperature adaptation

- No $E_a$ published specifically for *S. pasteurii*. Reasonable proxy: jack bean urease, $E_a \approx 24$ kJ/mol (range 18–35 kJ/mol across various studies). Our alpine sweeps additionally compare two MICP-specific estimates: Wang et al. (2023), $E_a \approx 31.2$ kJ/mol, and Sun et al. (2019), $E_a \approx 51.7$ kJ/mol.
- [Peng, Liu & Achal (2019), PLoS ONE](https://pmc.ncbi.nlm.nih.gov/articles/PMC6581288/) — sand columns at 10/15/20/25/30°C: at 10°C the final mechanical strength is **2.35× higher** than at 30°C, because the crystals grow more slowly and adhere better (>100 μm at 10°C vs <10μm at 30°C).
- [Lauchnor et al. (2015)](https://academic.oup.com/jambio/article/118/6/1321/6717359) — base kinetics ($V_{max}$=200 mmol/L/h, $K_m$=305 mM, $k_1$=0.35 h⁻¹ for the first-order model).

## Main sources

Razbani, M.A., Røyne, A., Jettestuen, E. (2025). *A One-Dimensional Reactive Transport Model for Microbially Induced Calcium Carbonate Precipitation (MICP) in Phreeqc*. Geomicrobiology Journal, 42(1), 19–29.
DOI: [10.1080/01490451.2024.2425295](https://doi.org/10.1080/01490451.2024.2425295)
Code: [Razbani 2024, Zenodo](https://doi.org/10.5281/zenodo.13886156)

Reference experiment that the model reproduces: Martinez et al. (2013), *J Geotech Geoenviron Eng* 139(4):587–598 — a 52 cm long sand column, tests "3A" (stop-flow) and "3B" (continuous flow). Our file is called `test_3B.in` because it exactly reproduces **Test 3B** of that experiment.

[Phreeqpy](https://www.phreeqpy.com/) is the Python interface to PHREEQC used throughout. The full description of each PHREEQC keyword block is in the [PHREEQC manual](https://pubs.usgs.gov/tm/06/a43/).
