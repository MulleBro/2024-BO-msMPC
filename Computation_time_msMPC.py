import time
import matplotlib.gridspec as gridspec
import do_mpc
import numpy as np
import matplotlib.pyplot as plt
from casadi import *
from casadi.tools import *
import pdb
import sys
import os
rel_do_mpc_path = os.path.join('..', '..')
sys.path.append(rel_do_mpc_path)

from helper.const_params import *
from helper.aux_functions import *

from model.model import template_model
from mpc.MPC_config_soft_multistage import template_mpc
from scipy import optimize


""" User settings: """
n_sim=72
dt=1


"""
Get configured do-mpc modules:
"""
model = template_model()
mpc = template_mpc(model)
estimator = do_mpc.estimator.StateFeedback(model)
model_list(model)
#%%
# =============================================================================
# INITIALS
# =============================================================================
x0=np.zeros([model.n_x])
x0[get_ind_x(model,'X')]=0.05
x0[get_ind_x(model,'Xyl')]=0
x0[get_ind_x(model,'Xlac')]=0
x0[get_ind_x(model,'Xat')]=0
x0[get_ind_x(model,'O2_L')]=O2_L_max
x0[get_ind_x(model,'O2_A')]=x_O_Air_in/V_M
x0[get_ind_x(model,'V_L')]=0.2
x0[get_ind_x(model,'Lys')]=100
x0[get_ind_x(model,'HLys')]=0

#params
P_df=pd.read_csv('./p_est/res/p_opt.csv')  
P=np.array(P_df['val'])

#Inputs
#F_Xyl
u0=np.array([0]).reshape(-1,1) #ml/h
#%%
# =============================================================================
# SIMULATOR
# =============================================================================

simulator = do_mpc.simulator.Simulator(model)

params_simulator = {
    'integration_tool': 'idas',
    'abstol': 1e-10,
    'reltol': 1e-10,
    't_step': dt
}

simulator.set_param(**params_simulator)


# =============================================================================
# INITIALIZATION
# =============================================================================

p_num = simulator.get_p_template()
p_names=model_names(model)[1]
def p_fun(t_now):
    for i in range(model.n_p):
        p_num[p_names[i]]=P[get_ind_p(model,p_names[i])]
    return p_num

simulator.set_p_fun(p_fun)
simulator.setup()
simulator.reset_history()
simulator.x0 = x0
simulator.u0 = u0
simulator.set_initial_guess() 
estimator.x0 = x0
mpc.x0=x0
mpc.u0=u0
mpc.set_initial_guess() #initial guess using x0, u0


# =============================================================================
# MAIN LOOP MPC
# =============================================================================
tt=np.zeros(n_sim)
for k in range(n_sim):
    tic=time.time()

    u0 = mpc.make_step(x0)
    tt[k]=time.time()-tic
    y_next = simulator.make_step(u0)
    x0 = estimator.make_step(y_next)
    
# =============================================================================
#    print mean computation time 
# =============================================================================
print('mean computation time of multistge MPC: '+str(np.round(tt.mean(),4))+' s')    
    


