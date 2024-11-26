
# =============================================================================
# #import
# =============================================================================
import pickle
import do_mpc
import matplotlib.pyplot as plt
import numpy as np

from helper.const_params import *
from helper.aux_functions import *

from model.model import template_model

model=template_model()


def load_pickle(path_to_file):
    with open(path_to_file, 'rb') as f:
        data = pickle.load(f)
    return data


# =============================================================================
# DEFINE PICKLE FILES IN WHICH TRAJECTORIES ARE SAVED
# =============================================================================
n_runs=100

folder='results_100_manu0'


alll0=['./results/results_100_manu0/100_nom.pkl']

for i in range(1,10):
    alll0.append('./results/results_100_manu0/00'+str(i)+'_100_nom.pkl')

for i in range(10,n_runs):
    alll0.append('./results/results_100_manu0/0'+str(i)+'_100_nom.pkl')
  

# =============================================================================
# PLOTTING NOMINAL MPC UNDER UNCERTAINTS, 100 PARAMETER SETS, LOOPING OVER PICKLE FILES
# =============================================================================


import helper.config_mpl
import importlib
importlib.reload(helper.config_mpl)
    
for j in range(len(alll0)):
    filename0=alll0[j]

    alldata0=load_pickle(filename0)
    mpc0=alldata0.get("mpc")
    sim_time=mpc0._x.shape[0]
    
    
    if j==0:
        figDO, axDO = plt.subplots(2,1, sharex=True)
 
    axDO[1].plot(mpc0._aux[:,get_ind_aux(model, 'DO')]*100,alpha=0.15, color='b') 


    
axDO[1].hlines(30, 0, int(np.max(mpc0._time)), linestyle='--', color='r')
axDO[1].set_ylabel(r'$\bm{\Theta}_{\mathcal{U}}$''\n'r'$DO [\%]$', rotation=90, ha="center")
axDO[-1].set_xlabel(r'$t\ [\text{h}]$')

        
#%%
# =============================================================================
# PLOTTING NOMINAL MPC UNDER NOMINL CONDITIONS
# =============================================================================
filename_nom='./results/mpc_nom.pkl'
alldata_nom=load_pickle(filename_nom)
mpc_nom=alldata_nom.get("mpc")

axDO[0].set_ylabel(r'$\bm{\Theta}_{\text{nom}}$''\n'r'$DO [\%]$', rotation=90, ha="center")
axDO[0].plot(mpc_nom._aux[:43,get_ind_aux(model, 'DO')]*100,alpha=1, color='b') 
axDO[0].hlines(30, 0, int(np.max(mpc0._time)), linestyle='--', color='r')
plt.setp(axDO[0].spines.values(), lw=0.5);
plt.setp(axDO[1].spines.values(), lw=0.5);
plt.subplots_adjust(left=0.18, bottom=0.19, right=0.99, top=0.98,wspace=0.22, hspace=0.2)
plt.show()

#%%
# =============================================================================
# LINEAR INTERPOLATION FOR NOMINAL MPC UNDER NOMINAL CONDIIONS: t_batch (99% CONVERSION) AND PRODUCT TITER
# =============================================================================
X_Lys=mpc_nom._aux[:,get_ind_aux(model, 'X_Lys')]
t_99_2a=np.argmax(X_Lys>=0.99)
t_99_2=int(mpc_nom._time[t_99_2a][0])

X99=X_Lys[t_99_2]
X99_1=X_Lys[t_99_2-1]

t_99=t_99_2-1+(0.99-X99_1)/(X99-X99_1) 

HLys_99=mpc_nom._x[t_99_2-1,get_ind_x(model, 'HLys')]+(t_99-(t_99_2-1))*(mpc_nom._x[t_99_2,get_ind_x(model, 'HLys')]-mpc_nom._x[t_99_2-1,get_ind_x(model, 'HLys')])

 #control 99% 
XLys_99=mpc_nom._aux[t_99_2-1,get_ind_aux(model, 'X_Lys')]+(t_99-(t_99_2-1))*(mpc_nom._aux[t_99_2,get_ind_aux(model, 'X_Lys')]-mpc_nom._aux[t_99_2-1,get_ind_aux(model, 'X_Lys')])


print('production time t_batch of nomnial mpc: '+str(np.round(t_99,2))+' h')
print('obtained product titer at t_batch: '+str(np.round(HLys_99,2))+' mmol/l')