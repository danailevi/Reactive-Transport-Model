"""Simulation of the Test 3B."""

import phreeqpy.iphreeqc.phreeqc_com as phreeqc_mod
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

phreeqc = phreeqc_mod.IPhreeqc()
phreeqc.load_database(r"C:\Program Files\USGS\IPhreeqcCOM 3.8.6-17100\database\phreeqc.dat") 

initial_conditions = """
#-------------------------------------------
SOLUTION_MASTER_SPECIES
    Urea Urea  0 CO(NH2)2 60.06 
     Amm   AmmH+  0.0  N    14.0067
    
    
SOLUTION_SPECIES 
 Urea = Urea
            log_k  0 
            Vm  14.5
    AmmH+ = AmmH+
        log_k 0
    AmmH+ = Amm + H+
        log_k -9.23        
       Cl- + H+ = HCl
        log_k 0.7
        delta_h 0 kcal
        -gamma 0 0.06
    Cl- + Ca+2 = CaCl+
        log_k 0.3
        delta_h 0 kcal
        -gamma 0 0.06
    2Cl- + Ca+2 = CaCl2
        log_k .66
        delta_h 0 kcal
        -gamma 0 0.06
#---------------------------------------------

KNOBS
    -iterations 500
    -convergence_tolerance 1e-9
    -tolerance 1e-14
    -step_size 100
    -pe_step_size 5
    -diagonal_scale true
END

RATES
Ureolysis
-START
11 cell_l = 5.2         # cell length, cm
12 t_d    = 288000*1.2      # biomass decay constant, s
15 initial_por = .38    # initial porosity ε₀
16 Cb = (5.66e-1*(cell_no*cell_l-cell_l/2)^3+1.1e2*(cell_no*cell_l-cell_l/2)^2+4.67e3*(cell_no*cell_l-cell_l/2)-5.62e3)/1e6 # OD, attached biomass profile
18 decay_ratio = exp(-TOTAL_TIME/t_d)
20 rate = Cb*soln_vol*get_por(cell_no)/initial_por*PARM(1)*MOL("Urea")/(MOL("Urea")+PARM(2))
30 moles = rate * TIME * decay_ratio # * get_por(cell_no)/initial_por
40 save moles
50 REM PARM(1) = maximum urease activity (mol/l)      #Rmax
60 REM PARM(2) = half saturation constant             #Km
-END

Calcite
-START
20 initial_por = .38
21 V_cell = 102 # ml
22 d_calcite = 36.9 # cm³/mol
23 pow_a = 2
25 A0 = 1. # m2/l
30 si_cc = SI("Calcite")
70 Area = A0 * (initial_por/get_por(cell_no))^pow_a * soln_vol
80 rate = Area * (parm(1) + parm(2)*act("H+")) * (1-10^SI("Calcite"))
100 moles = rate * TIME * get_por(cell_no)/initial_por 
110 change_por(initial_por - d_calcite*kin("Calcite")/V_cell,cell_no) 
200 SAVE moles
-END
#-------------------------------------------------------------------------------

SOLUTION 0
Urea .05
Amm .0567
Ca 0.05
Cl .1567
C .0038 as Ca Calcite 0.0
Na .0038
-units 	mol/kgw
-pH 	7	charge_balance
water 38.7e-3 # kg
END

SOLUTION 1-10
     -units 	mol/kgw
     -pH 	7	charge_balance
     water 38.7e-3 # kg
END

Kinetics 1-10 
    Ureolysis
         -m0 100  
        -formula AmmH+ 2 CO3-2 1 Urea -1 H2O -2
        -parms 7.0e-5 0.0173 #[Rmax] [Km]
    Calcite
        -m0   0# mol/L
        -parms 1.55e-6 .501 
END

TRANSPORT 
   cells 10
   shifts 178 # 50 h
   time_step 1055 seconds # 16 min
   flow_direction forward
   boundary_conditions flux constant 
   length 0.052
   diffusion_coefficient 1.0e-9
   porosities 10*.38
   -dispersivities 10*.01

USER_PUNCH
headings por 
 10 punch get_por(cell_no)

SELECTED_OUTPUT
    -step    true
     -h    false
    -sim    false
    -soln   true
     -pe     false
    -totals false
    -time   true
    -pH true
    -tot  Ca  Urea Amm
    -si Calcite
    -k Ureolysis Calcite
    -dist false
    -state false
END
"""
np.set_printoptions(precision=2)
pd.set_option("display.precision", 3)


phreeqc.run_string(initial_conditions)
output = phreeqc.get_selected_output_array()
table = pd.DataFrame(output[1:], columns=output[0])
print("\n")
print(table)

n_cells = 10  # total number of cells
V_cell = 102e-3  # l
t_step = 1055  # time steps
shift_m = 178
row_0 = 12  # cell number
t = table.loc[4:]["time"]
C_u_0 = 0.05  # initial urea concentration molar
C_ca_0 = 0.05  # initial calcium concentration molar
por_0 = 0.38  # initial porosity
l_cell = 5.2  # cm

final = table[table["step"] == 178]
final = final[(final["soln"] >= 1) & (final["soln"] <= 10)].sort_values("soln")

heights = np.arange(1, n_cells + 1) * l_cell - l_cell / 2   # cell centres (cm)
calcite = final["k_Calcite"].values / V_cell                # mol/L

final.to_csv(
    "T3B_calcite_shift178.csv", 
    index=False
)

fig, ax = plt.subplots()
ax.plot(calcite, heights, "ko", label="Simulation", ms=8)
ax.set_xlabel("Calcite precipitation (mol/L)")
ax.set_ylabel("Column height (cm)")
ax.set_xlim(0, 0.6)
ax.set_ylim(0, 52)
ax.grid()
ax.legend()


plt.show()
