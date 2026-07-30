---
tags: [iGEM, dry-lab, MICP, PHREEQC]
created: 2026-07-09
---

# Modello MICP — note di riferimento

> [!note] Cos'è
> Note per capire `test_3B.in`, il modello PHREEQC che simula la cementazione batterica (MICP) in una colonna di sabbia. Basato sul paper originale di Razbani, l'autore del codice che stiamo usando come punto di partenza.

## Fonte principale

Razbani, M.A., Røyne, A., Jettestuen, E. (2025). *A One-Dimensional Reactive Transport Model for Microbially Induced Calcium Carbonate Precipitation (MICP) in Phreeqc*. Geomicrobiology Journal, 42(1), 19–29.
DOI: [10.1080/01490451.2024.2425295](https://doi.org/10.1080/01490451.2024.2425295)
Codice: [Razbani 2024, Zenodo](https://doi.org/10.5281/zenodo.13886156)

Esperimento di riferimento che il modello riproduce: Martinez et al. (2013), *J Geotech Geoenviron Eng* 139(4):587–598 — colonna di sabbia lunga 52 cm, test "3A" (stop-flow) e "3B" (flusso continuo). Il nostro file si chiama `test_3B.in` perché riproduce esattamente il **Test 3B** di quell'esperimento.

## Il quadro fisico

- Colonna di sabbia, iniettata dal basso con **urea + Ca²⁺ + cloruro** e dall'alto con una soluzione ricca di **batteri** (*S. pasteurii*), che si depositano e restano attaccati ai granelli.
- I batteri idrolizzano l'urea (**ureolisi**) → producono ammonio e carbonato.
- Il carbonato incontra il calcio → precipita **calcite (CaCO₃)**, che cementa la sabbia.
- La calcite che si deposita **riduce la porosità**, rallentando la reazione successiva in quel punto (feedback).

**Reazione complessiva:**
$$
\mathrm{CO(NH_2)_2} + 2\,\mathrm{H_2O} \longrightarrow 2\,\mathrm{NH_4^+} + \mathrm{CO_3^{2-}}
$$

## Struttura del codice `test_3B.in`

| Blocco | Cosa fa | Righe |
|---|---|---|
| `SOLUTION_MASTER_SPECIES` / `SOLUTION_SPECIES` | "Vocabolario": insegna a PHREEQC cos'è l'urea e come si comportano ioni come Cl⁻/Ca²⁺ | 1–20 |
| `RATES → Ureolysis` | Formula della velocità di rottura dell'urea (dipende da biomassa e urea disponibile) | 22–33 |
| `RATES → Calcite` | Formula della velocità di precipitazione + aggiornamento della porosità | 35–47 |
| `SOLUTION 0` | Ricetta del liquido iniettato (urea, Ca, ammonio, pH) | 49–59 |
| `SOLUTION 1-10` | Cosa c'è nel tubo all'inizio (solo acqua) | 61–65 |
| `KINETICS 1-10` | Accende le due reazioni in ognuna delle 10 fettine della colonna | 67–75 |
| `TRANSPORT` | Fa scorrere il liquido: 10 celle, 178 shift, ~2 giorni simulati totali | 77–86 |
| `SELECTED_OUTPUT` | Cosa scrivere su file per fare i grafici dopo | 88–96 |

## Le formule principali

**Cinetica di ureolisi (Michaelis-Menten, pesata su biomassa):**
$$
R_u = C_b \, V_{soln} \, \frac{\phi(x)}{\phi_0} \, R_{max} \, \frac{[\mathrm{Urea}]}{[\mathrm{Urea}] + K_m}
$$

**Decadimento della biomassa nel tempo:**
$$
D(t) = \exp\left(-\frac{t}{t_d}\right)
$$

**Precipitazione/dissoluzione della calcite:**
$$
R_p = S\left(k_1 + k_2\, a_{\mathrm{H^+}}\right)\left(1 - \Omega\right), \qquad \Omega = \frac{a_{CO_3^{2-}} \cdot a_{Ca^{2+}}}{K}
$$

**Area reattiva (si riduce quando i pori si intasano):**
$$
S = S_0 \left(\frac{\phi}{\phi_0}\right)^{p}
$$

## Tabella parametri (Test 3B, dal paper — Tabella 3)

| Parametro | Valore | Fonte | Corrisponde a |
|---|---|---|---|
| $R_{max}$ | $7\times10^{-5}$ mol/L/s/OD | Wang & Nackenhorst (2020) | `PARM(1)` in `Kinetics → Ureolysis` |
| $K_m$ | 0.0173 mol/kgH₂O (17.3 mM) | Krajewska (2009) | `PARM(2)` |
| $S_0$ | 1 m²/L | Martinez et al. (2014) | `-parms` in `Kinetics → Calcite` |
| $p$ | 2 | assunto | `pow_a` nel codice `RATES` |
| $t_d$ | 288 000 s | Fauriel & Laloui (2012) | `t_d` nel codice `RATES` |
| Dispersività | 0.005 m | Wang & Nackenhorst (2020) | ⚠️ il nostro file ha `0.01` — da controllare |

> [!warning] Nota sul K_m
> Esiste anche un valore più recente, $K_m = 305$ mM (Lauchnor et al. 2015, cinetica su cellula intera), molto più alto di 17.3 mM. Razbani ha scelto 17.3 mM per restare fedele alla calibrazione originale di Martinez et al. Non è un errore nel nostro file, ma è un punto da discutere/citare se decidiamo di aggiornarlo.

## Analisi di sensibilità (dal paper) — dove concentrare gli sforzi

- **Poco sensibile**: area reattiva iniziale $S_0$ (sopra una soglia ~5 m²/L) ed esponente $p$ — la precipitazione non è il fattore limitante, quindi affinare questi parametri non cambia molto il risultato.
- **Sensibile**: $K_m$ (±40% → differenze rilevanti), densità di biomassa (±20%), $t_d$ (±30%).
- **Conclusione pratica**: quando adatteremo il modello alle Alpi, conviene concentrarsi su come il freddo cambia **ureolisi ($R_{max}$, $K_m$)** e **sopravvivenza dei batteri ($t_d$)**, non sulla cinetica di precipitazione.

## Cosa NON è nel modello originale (nostro lavoro da aggiungere)

- ❌ **Temperatura** — il paper non la considera mai. Serve aggiungere una correzione tipo Arrhenius su $R_{max}$ nel BASIC delle `RATES` (variabile `TC` disponibile in PHREEQC).
- ❌ Condizioni non sature (serve l'approccio di Faeli et al., come da roadmap di luglio/agosto)
- ❌ pH acido tipico dei suoli alpini (4.5–6) — da testare separatamente

## Bug risolto nel nostro `phreeqc.dat`

La reazione `AmmH+ = Amm + H+` non aveva un `-log_k` proprio (mancava, veniva letto come 0 invece del vero pKa ~9.25). Corretto aggiungendo:
```
AmmH+ = Amm + H+
    -log_k -9.252
    -delta_h 12.48 kcal
    -analytic 0.6322 -0.001225 -2835.76
```
Questo risolveva l'errore "Model failed to converge" che uscivamo ogni volta.

## Letteratura utile per l'adattamento alla temperatura

- Nessun $E_a$ pubblicato specifico per *S. pasteurii*. Proxy ragionevole: ureasi jack bean, $E_a \approx 24$ kJ/mol (range 18–35 kJ/mol tra vari studi).
- [Peng, Liu & Achal (2019), PLoS ONE](https://pmc.ncbi.nlm.nih.gov/articles/PMC6581288/) — colonne di sabbia a 10/15/20/25/30°C: a 10°C la resistenza meccanica finale è **2.35× più alta** che a 30°C, perché i cristalli crescono più lentamente e aderiscono meglio (>100 μm a 10°C vs <10μm a 30°C).
- [Lauchnor et al. (2015)](https://academic.oup.com/jambio/article/118/6/1321/6717359) — cinetica base ($V_{max}$=200 mmol/L/h, $K_m$=305 mM, $k_1$=0.35 h⁻¹ per il modello del prim'ordine).

## Prossimi passi (dal nostro TODO di progetto)

- [ ] Leggere Razbani et al. 2024 per intero, annotando le equazioni chiave *(fatto con queste note)*
- [ ] Riprodurre Test 3A (stop-flow) con lo stesso codice
- [ ] Test di convergenza griglia (raddoppiare celle 10→20)
- [ ] Trovare $R_{max}$ a bassa temperatura (4–12°C) per *S. pasteurii* o *Psychrobacter*
- [ ] Implementare correzione Arrhenius su $R_{max}$ nelle `RATES`
