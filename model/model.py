
import numpy as np
from casadi import *
from casadi.tools import *
import pdb
import sys
import os
rel_do_mpc_path = os.path.join('..','..')
sys.path.append(rel_do_mpc_path)
import do_mpc
from helper.const_params import *

def template_model(symvar_type='SX'):
    """
    --------------------------------------------------------------------------
    template_model: Variables / RHS / AUX
    --------------------------------------------------------------------------
    """
    model_type = 'continuous' # either 'discrete' or 'continuous'
    model = do_mpc.model.Model(model_type, symvar_type)
    
    # =============================================================================
    # STATES
    # =============================================================================
    X = model.set_variable('_x',  'X')  # Biomass in g/L
    Xyl = model.set_variable('_x',  'Xyl')  # Xylose concentration in mM, fed
    
    Xlac = model.set_variable('_x',  'Xlac') # Xylonolacton concentration in mM 
    Xat = model.set_variable('_x',  'Xat') # Xylonat concentration in mM

    O2_L = model.set_variable('_x',  'O2_L') # O2_Liquid concentration in mM
    O2_A= model.set_variable('_x',  'O2_A') # O2_Air concentration in mM
    #CO2_L_tot=model.set_variable('_x',  'CO2_L_tot') # total CO2_Liquid concentration (CO2_aq, HCO3-, CO3--) in mM
    #CO2_A= model.set_variable('_x',  'CO2_A') # CO2_Air concentration in mM

    V_L = model.set_variable('_x',  'V_L')  # Reactor volume in L
    
    Lys = model.set_variable('_x',  'Lys') # Lysine concentration in mM, fed
    HLys = model.set_variable('_x',  'HLys') # Hydroxylysine concentration in mM

# =============================================================================
#     #INPUTS:
# =============================================================================
    #pH_inp = model.set_variable('_u', 'pH_inp')   #Substrate feed Xylonat
    
    #INPUT FEEDING XYLOSE
    F_XylN= model.set_variable('_u',  'F_XylN')     #ph control base NH4OH
    


# =============================================================================
#     PARAMS

    #BIOMASS PRODUCTION
    kX= model.set_variable('_p',  'kX')
    Y_s= model.set_variable('_p',  'Y_s')
    Ks_Xyl_X= model.set_variable('_p',  'Ks_Xyl_X')
   
    #PRODUCT SYNTHESIS
    kP= model.set_variable('_p',  'kP')
    Km_P_Xyl= model.set_variable('_p',  'Km_P_Xyl')
    Km_P_Lys= model.set_variable('_p',  'Km_P_Lys')

    
    #INTERMEDIATES CONVERSION
    k1= model.set_variable('_p',  'k1')
    k2= model.set_variable('_p',  'k2')
    k3= model.set_variable('_p',  'k3')
    Km_Xyl_Xat= model.set_variable('_p',  'Km_Xyl_Xat')
    Km_Xlac= model.set_variable('_p',  'Km_Xlac')
    K2= model.set_variable('_p',  'K2')   
    
    
    #GAS PHASE EXCHANGE
    kla_O2= model.set_variable('_p',  'kla_O2')
    #kla_CO2= model.set_variable('_p',  'kla_CO2')  
    
    # pH_min= model.set_variable('_p',  'pH_min')
    # pH_max= model.set_variable('_p',  'pH_max')
    # pH_opt= model.set_variable('_p',  'pH_opt')  
   
    f_XP= model.set_variable('_p',  'f_XP')
    

    c_F_Xyl=500/M_Xyl*1e3 #mM

    # =============================================================================
    # VOLUME
    # =============================================================================    
      
    V_Air=V_vessel-V_L
    model.set_expression('V_Air', V_Air)

    F_in=F_XylN
    
    D=F_in/V_L #Dilution
    model.set_expression('D', D)
    
    F_out=DM(0)
# =============================================================================
#     KINETICS
# =============================================================================
    mu=kX*Xyl/(Ks_Xyl_X+Xyl)
    rS=X*mu/Y_s
    rX=mu*X

    r1=k1*X*Xyl/(Km_Xyl_Xat+Xyl)
    r2=k2*(Xlac-Xat/K2)
    r3=k3*X*(Xlac-Xat/K2)/(Xlac-Xat/K2+Km_Xlac)
    
    rP=X*kP*Lys*Xyl/(Km_P_Xyl*Xyl+Km_P_Lys*Lys+Lys*Xyl)

# =============================================================================
# pH
# =============================================================================
 
    
    # x_=10**(-pH_inp)
    # model.set_expression('x_', x_)
    

    # k_pH0=((pH_inp-pH_min)*(pH_inp-pH_max))/((pH_inp-pH_min)*(pH_inp-pH_max)-(pH_inp-pH_opt)**2)
    
    # k_pH = if_else(k_pH0 > 0, k_pH0, 0)
    # model.set_expression('k_pH', k_pH)
    

    # # =============================================================================
    # #     pH-dependency
    # # =============================================================================

    #CO2 dissociation water
    # CO2_L_aq=CO2_L_tot/(1+K_C_1/xsim_+K_C_1*K_C_2/(xsim_*xsim_))
    # HCO3_L_aq=CO2_L_aq*K_C_1/xsim_
    # CO3_L_aq=HCO3_L_aq*K_C_2/xsim_
    #CO2_L_aq=CO2_L_tot*0.98
    #CO2_L_aq=CO2_L_tot/(1+K_C_1/x_+K_C_1*K_C_2/(x_*x_))
    # HCO3_L_aq=CO2_L_aq*K_C_1/x_
    # CO3_L_aq=HCO3_L_aq*K_C_2/x_
    
    #model.set_expression('CO2_L_aq', CO2_L_aq)
    # model.set_expression('HCO3_L_aq', HCO3_L_aq)
    # model.set_expression('CO3_L_aq', CO3_L_aq)
    
# =============================================================================
# Gas phase
# =============================================================================

    
    V_M=R*T/p_ #molar volume,  V_M=RT/p_, T=303.15 K, p_ in Pa
    
    c_O2_in=x_O_Air_in*p_/(R*T) #mol/m3 = mmol/L
    #c_CO2_in=x_C_Air_in*p_/(R*T) #mol/m3 = mmol/L
    
    x_O2_A=O2_A*R*T/p_
    #x_CO2_A=CO2_A*R*T/p_
    
    model.set_expression('x_O2_A', x_O2_A)
    #model.set_expression('x_CO2_A', x_CO2_A)
    
    #DO
    O2_L_max=1e3*1*x_O_Air_in/H_O_30 #mM
    DO=O2_L/O2_L_max
    model.set_expression('DO', DO)
    
    #OTR, CTR
    OTR=kla_O2*(1e3*1*x_O2_A/H_O_30-O2_L)
    #CTR=kla_CO2*(CO2_L_aq-1e3*1*x_CO2_A/H_C_30)
    
    
    model.set_expression('OTR', OTR)
   # model.set_expression('CTR', CTR)
# =============================================================================
#     CONVERSION X_LYS
# =============================================================================
    c_Lys_0=100
    #n_Lys_0=c_Lys_0*V_L_0 #20mmol BATCH
    
    #V_L_0_FB=V_L_0
    V_L_0_FB=0.2
    n_Lys_0=c_Lys_0*V_L_0_FB #FEDBATCH

    X_Lys=1-(Lys*V_L/n_Lys_0)
    model.set_expression('X_Lys', X_Lys)

    
    
# =============================================================================
# ODE#s STATES
# =============================================================================
    model.set_rhs('X', -D*X+(rX+f_XP*rP))
    model.set_rhs('Xyl',-D*Xyl +c_F_Xyl*F_XylN/V_L-(rS+r1+rP))

    model.set_rhs('Xlac', -D*Xlac+ (r1-r2+r3))
    model.set_rhs('Xat', -D*Xat+(r2-r3))
    
    #model.set_rhs('O2_L', -(5*rS+rP)+OTR-O2_L*D+O2_L_0*F_in/V_L)
   # model.set_rhs('CO2_L_tot', (5*rS+rP)-CTR-CO2_L_tot*D+CO2_L_0*F_in/V_L)
    
    model.set_rhs('O2_L', -(5*rS+rP)+OTR+F_in/V_L*(O2_L_0-O2_L))
    #model.set_rhs('CO2_L_tot', (5*rS+rP)-CTR+F_in/V_L*(CO2_L_0-CO2_L_tot))
    
    # model.set_rhs('O2_L', -(5*rS+rP)+OTR-O2_L*D+O2_L_0*D)
    # model.set_rhs('CO2_L_tot', (5*rS+rP)-CTR-CO2_L_tot*D+CO2_L_0*D)
    
    
    model.set_rhs('O2_A', 1/V_Air*(F_Air*(c_O2_in-O2_A)-OTR*V_L)) 
   # model.set_rhs('CO2_A', 1/V_Air*(F_Air*(c_CO2_in-CO2_A)+CTR*V_L)) 
    
    model.set_rhs('V_L', F_in-F_out) 

    model.set_rhs('Lys', -D*Lys -rP)
    model.set_rhs('HLys', -D*HLys+rP) 


    # Build the model
    model.setup()
    
    return model  
    
    
