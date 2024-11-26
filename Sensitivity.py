

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
from mpc.MPC_config_soft_sen import template_mpc



model = template_model()
mpc = template_mpc(model)

model_list(model)
dt=1

timepoint=42

x_names=model_names(model)[0]
p_names=model_names(model)[1]
u_names=model_names(model)[3]
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


# INITIALIZATION, FUNCTIONS FOR DIFFERENT MODES AND SETUP

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
mpc.x0=x0
mpc.u0=u0
mpc.set_initial_guess() #initial guess using x0, u0

#%%
# =============================================================================
# MPC OPEN LOOP PREDICTION
# =============================================================================
nlp_diff =do_mpc.differentiator.DoMPCDifferentiator(mpc) #INITIALIZATION
n_sim=1

for i in range(n_sim):  
    ui = mpc.make_step(x0)
    x0 = simulator.make_step(ui)
    
    dx_dp_num, dlam_dp_num = nlp_diff.differentiate() #calculating parametric sensitivity at current optimal solution


store_results=False
if store_results:
    do_mpc.data.save_results([mpc], 'mpc_sensitivity')



#%%
# =============================================================================
# COMPUTATION sensitivity
# =============================================================================

optx_mpc=mpc._opt_x
#mpc._opt_x['_x']
J_mpc=mpc._nlp_obj

#SYMBOLIC
JAC_Jmpc=jacobian(J_mpc,optx_mpc) #symbolic but in J datapoints not symbolic
JAC_Jmpc_fun=Function('JAC_J_fun', [optx_mpc],[JAC_Jmpc])
opt_x_Jmpc = mpc._opt_x(JAC_Jmpc.T) #['_x', time_step, scenario, collocation_point, _x_name]


# #NOT SYMBOLIC
optx_mpc_dat_all=mpc.data._opt_x_num
optx_mpc_dat=optx_mpc_dat_all[-1].reshape(-1,1)
Jmpc_res=JAC_Jmpc_fun(optx_mpc_dat)

dJ_dp=Jmpc_res@dx_dp_num
dJ_dp_full=dJ_dp.full()

dJ_dp_est=dJ_dp_full[0,model.n_x:dJ_dp_full.shape[1]-model.n_u]
dJ_dp_est_scal=dJ_dp_est*P



dict_dJdp_mpc = {'p_name': p_names, 'val': dJ_dp_est_scal} 
df_dJdp_mpc = pd.DataFrame(dict_dJdp_mpc)

#%%
# =============================================================================
# PLOTTING
# =============================================================================

x=np.linspace(0,model.n_p-1,model.n_p)
p_name2=[r'$\mu_{\text{max}}$','$Y_{S}$',' $K_{S}$',r'$v_{\text{max},P}$','$K_{M,3}$','$K_{M,4}$',r'$v_{\text{max},1}$',r'$v_{\text{max},2}$',r'$v_{\text{max},3}$','$K_{M,1}$','$K_{M,2}$','$K_2$','$kla$','$k_{PX}$']


df_dJdp_abs=df_dJdp_mpc.copy()
df_dJdp_abs['val']=df_dJdp_abs['val'].abs()

df_dJdp_abs_pp=df_dJdp_abs.copy()
df_dJdp_abs['p_name']=p_name2

df_dJdp_abs_sort=df_dJdp_abs.sort_values('val', ascending=False)
df_dJdp_abs_pp_sort=df_dJdp_abs_pp.sort_values('val', ascending=False)



import helper.config_mpl
import importlib
importlib.reload(helper.config_mpl)



fig, ax=plt.subplots(1,1)
ax.bar(x,df_dJdp_abs_sort['val'])
plt.xticks(x, df_dJdp_abs_sort['p_name'], rotation=90)  
plt.ticklabel_format(axis='y', style='sci', scilimits=(4,4))

ax.set_xlabel(r'$\bm{\tilde{\Theta}}$ ')
ax.set_ylabel(r'$ \text{abs}(\partial J_{\text{MPC}} / \partial \bm{\tilde{\Theta}} \cdot \bm{\tilde{\Theta}})$')

plt.setp(ax.spines.values(), lw=0.5);

plt.subplots_adjust(left=0.15, bottom=0.32, right=0.99, top=0.92)


