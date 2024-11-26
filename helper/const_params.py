# -*- coding: utf-8 -*-
"""
Created on Tue Nov  7 09:17:43 2023

@author: smtobroc
"""
import numpy as np
import pandas as pd


# =============================================================================
# KINETICS
# =============================================================================
# #X: haldane kinetic from do-mpc model  
# mu_m  = 0.02
# K_m	  = 0.05
# K_i	  = 5.0
# v_par = 0.004
# Y_p	  = 1.2
# Y_x=0.4
# Y_x=0.8

#params marischen
k1=3.19*10 
k2=0.05
k3=1.83
k4=18.2

k_xyl=3.02
k_lac=8.93
k_xat=34.9

#dummy
k6=0.1
K2=1 #KEINE ANGABE MARISCHEN?
kX=18.2
k_aKG=34
Y_s=0.4


  
# =============================================================================
# MOLECULAR WEIGHTS NCBI Pubchem
# =============================================================================
#METABOLISM    
M_Xyl=150.13 #g/mol NCBI
M_Xlac=148.11 #g/mol NCBI
M_Xat=165.12 #g/mol  NCBI
M_aKG=146.10 #g/mol  NCBI
M_Lys=146.19  #g/mol  NCBI
M_HLys=162.19 #g/mol  NCBI (3,4,5? weight anyway the same)
M_O2=32 #g/mol  NCBI
M_CO2=44.01 #g/mol  NCBI

#pH
M_H3PO4=98.00 #g/mol  NCBI
M_NH4OH=35.05 #g/mol  NCBI
#M_Na2HPO4=141.96 #g/mol  NCBI
M_Na2HPO4_2H2O=177.99#g/mol SIGMA ALDIRICH https://www.sigmaaldrich.com/DE/de/product/sial/71645
M_KH2PO4=136.09 #g/mol  NCBI 
M_NaCl=58.44 #g/mol  NCBI
M_NH4Cl=53.49 #g/mol  NCBI
M_MgSO4=120.37 #g/mol  NCBI

M_PO4=92 #g/mol
M_PO4=31 #MARISCHEN 
# =============================================================================
# PURER PHOSPHOR!
# =============================================================================

#M_FeSO4=151.91 #g/mol  NCBI
#TRACE ELEMENTS

#M_NaOH=40.00 #g/mol  NCBI
M_NH3=17.03 #g/mol  NCBI


# =============================================================================
# CONCENTRATION
# =============================================================================
cm_Na2HPO4_2H2O=8.5 #g/L  J_Hand
cm_KH2PO4=10 #g/L  J_Hand
cm_NaCl=0.5 #g/L  J_Hand
cm_NH4Cl=2 #g/L  J_Hand

#c_Trace=1e-3 #L/L

c_MgSO4=1*2e-3 #mol/L  J_Hand
c_FeSO4=0.032e-3 #mol/L  J_Hand
c_Na2HP=cm_Na2HPO4_2H2O/M_Na2HPO4_2H2O #mol/L
c_KH2P= cm_KH2PO4/M_KH2PO4
c_NH4Cl= cm_NH4Cl/M_NH4Cl
c_NaCl= cm_NaCl/M_NaCl
   
#MOL/L
c_H3PO4=1.74 # mol/L, 10 #PROZENT J_Hand, https://www.sedgeochem.uni-bremen.de/acidconc.html
c_NH4OH=2.87 # mol/L, WIKI, 5 PROZENT  J_Hand
#  
c_F_P = 1.74*1e3 #mM #H3PO4
c_F_N = 2.87*1e3 #mM #NH4OH

c_Fe2=0.032 #mM Cofactor=const? 


# =============================================================================
# pH ACID / BASE DISSOCIATION CONSTANT
# =============================================================================
#equilibrium constants: https://chem.libretexts.org/Ancillary_Materials/Reference/Reference_Tables/Equilibrium_Constants/E2._Base_Dissociation_Constants_at_25C
# find something for 30 °C ?
#UNITS: K_W mol^2*L^(-2)
#all the oter mol/L=M, concentrations in equations mM CONSIDER!

#Waterproduct
K_W=1e-14           #water
#Base
K_Am=1.75e-5         #Ammonia Base -> Kb-value  #1.82e-5 at 30^C

#MARISCHEN, NH3, H2CO3 f(T), but T=const in incubation!?, AND SALINITY, MAYBE LATER

#Acid
K_C_1=4.5e-7         #carbon dioxide 1: H2CO3-> H+ + HCO3-
K_C_2=4.7e-11        #carbon dioxide 2: HCO3-> H+ + CO3--
K_A_1=6.9e-3         #phosphoric acid 1: H3PO4-> H2PO4- + H+
K_A_2=6.2e-8         #phosphoric acid 2:  H2PO4 -> HPO4--  + H+ 
K_A_3=4.8e-13        #phosphoric acid

K_Xat_1=1.99e-4

#K_S_1= STRONG
#CHECK! SO4-- input =A_tot... formula should remain the same
K_S_2=1e-2 #Sulfuric acid, HSO4- -> SO4-- + H+, BUT MgSO4 / FeSO4 input!!

#if ph set, not modelled
pH_set=6.6
x_set=10**(-pH_set)
# =============================================================================
# AERIATION
# =============================================================================
F_Air=3 #L/h
x_O_Air_in=0.2094
x_C_Air_in=0.0003
x_C_Air_in=420*1e-6


p_=101325 #Pa
R=8.3145 #J/molK
T=273.15+30 #K

    
V_M=R*T/p_

V_vessel=0.25 #L
V_L_max=0.23 #due to foaming....
# =============================================================================
# ### EHER 0.23 reality!!
# =============================================================================

c_O2_Air_in=x_O_Air_in*p_/(R*T) #mM =mol/m^3 
c_CO2_Air_in=x_C_Air_in*p_/(R*T) #mM =mol/m^3 

#book, 25°C, H_PC 
H_O=790.6 #UNIT!!
H_C=29.7 #UNIT!!!
    
#30°C
K_H_O=1700
K_H_C=2400

import math
H_O_b=1/H_O
H_O_b30=H_O_b*math.exp(K_H_O*(1/T-1/298.15))
H_O_30=1/H_O_b30

H_C_b=1/H_C
H_C_b30=H_C_b*math.exp(K_H_C*(1/T-1/298.15))
H_C_30=1/H_C_b30


O2_L_max=1e3*1*x_O_Air_in/H_O_30

# =============================================================================
# INITIALS
# =============================================================================
# States
V_L_0=0.2 #L

X_0 = 0.041*3  # g/l OD=0.2*3 #? DACHTE DASbox 0.6 OD, nicht 0.2 beimpft? Thesis 0.6, data 0.2?
V_0 = V_L_0  # L
Xyl_0 = 133.22  # mM=mmol/l, 20g/L=133.22mM
Xlac_0 = 0  # mM
Xat_0 = 0  # mM
aKG_0 = 0.004 # mM #if Haldane model aKG_0 >0.004/1.2, otherwise first step negative..
Lys_0 = 100  # mM # Handke 100 mM, aber nicht fedbatch, recommendation <100
HLys_0 = 0  # mM
O2_0_Air = c_O2_Air_in  # mM
CO2_0_Air = c_CO2_Air_in  # mM

O2_L_0=1e3*1*x_O_Air_in/H_O_30
CO2_L_0=1e3*1*x_C_Air_in/H_C_30 #total CO2, including HCO3-, CO3--

V_Air_0=V_vessel-V_L_0

P_0 = 1e3*(c_Na2HP+c_KH2P)  # mM
N_0 = 1e3*c_NH4Cl  # mM


k_O2_max=F_Air*V_M/(1/H_O_30*V_L_0) #1/h
k_CO2_max=F_Air*V_M/(1/H_C_30*V_L_0)
# =============================================================================
# MOLAR DENSITY!! instead of V_M
# =============================================================================



# =============================================================================
# PANDAS TRY LATER
# =============================================================================
# data = {
#   "M": [150.13, 148.11, 165.12, 146.10, 146.19, 162.19, 32, 44.01, 98.00, 35.05, 177.99, 136.09, 58.44, 53.49, 120.37], #g/mol
#   #"concentration [g/L]": [50, 40, 45]
# }

# df = pd.DataFrame(data, index = ["Xyl", "Xlac", "Xat", 'aKG', 'Lys', 'HLys', 'O2', 'CO2', 'H3PO4', 'NH4OH','Na2HPO4_2H2O','KH2PO4', 'NaCl', 'NH4Cl', 'MgSO4' ])

# print(df) 

