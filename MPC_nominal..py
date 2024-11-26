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
from mpc.MPC_config_soft import template_mpc


""" User settings: """
n_sim=72
dt=1



store_results=False
filename='./results/mpc_nom.pkl' #NEED TO UPDATE FILENAME IF NEW RESULTS STORED

"""
Get configured do-mpc modules:
"""
model = template_model()
mpc = template_mpc(model)
estimator = do_mpc.estimator.StateFeedback(model)

model_list(model) #name all states, inputs, parameters and auxillary expressions
#%%
# =============================================================================
# INITIALS
# =============================================================================
#states
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

#parameters
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


# INITIALIZATION

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

#%%

# =============================================================================
# MAIN LOOP MPC
# =============================================================================


for k in range(n_sim):

    u0 = mpc.make_step(x0)
    y_next = simulator.make_step(u0)
    x0 = estimator.make_step(y_next)

if store_results:
    do_mpc.data.save_results([mpc, simulator], 'mpc_nom')




#%%
# =============================================================================
# Loading saved results
# =============================================================================

#data impot preparation
import pickle

def load_pickle(path_to_file):
    with open(path_to_file, 'rb') as f:
        data = pickle.load(f)
    return data


# Data preparation
alldata=load_pickle(filename)
mpc_d=alldata.get("mpc")
simulator_d=alldata.get("simulator")


#%%

# =============================================================================
# linear interpolation of t_batch (99% conversion of lysine) and obtained product titer
# =============================================================================

X_Lys=mpc.data._aux[:,get_ind_aux(model, 'X_Lys')] #conversion
H_Lys=mpc.data._x[:,get_ind_x(model, 'HLys')] #product titer

#t_batch
t_99=np.argmax(X_Lys>=0.99)
X99=X_Lys[t_99]
X99_1=X_Lys[t_99-1]
t_99_int=t_99-1+(0.99-X99_1)/(X99-X99_1) 

#product titer at t_batch
P99=H_Lys[t_99]
P99_1=H_Lys[t_99-1]
P99_int=P99_1+(t_99_int-(t_99-1))/1* (P99-P99_1)

print('production time t_batch of nomnial mpc: '+str(np.round(t_99_int,2))+' h')
print('obtained product titer at t_batch: '+str(np.round(P99_int,2))+' mmol/l')