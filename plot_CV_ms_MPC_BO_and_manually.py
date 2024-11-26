
# =============================================================================
# #import
# =============================================================================
import pickle
import do_mpc
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
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


alll0=['./results/results_100_manu8/100_nom8.pkl']

for i in range(1,10):
    alll0.append('./results/results_100_manu8/00'+str(i)+'_100_nom8.pkl')

for i in range(10,n_runs):
    alll0.append('./results/results_100_manu8/0'+str(i)+'_100_nom8.pkl')
    
    
    
alllBO=['./results/results_100_BO/100_nomBO.pkl']

for i in range(1,10):
    alllBO.append('./results/results_100_BO/00'+str(i)+'_100_nomBO.pkl')

for i in range(10,n_runs):
    alllBO.append('./results/results_100_BO/0'+str(i)+'_100_nomBO.pkl')
    

# =============================================================================
# PLOTTING MULTI-STAGE MPC, TUNED VIA BO AND MANUALLY, 100 PARAMETER SETS, LOOPING OVER PICKLE FILES
# =============================================================================


import helper.config_mpl
import importlib
importlib.reload(helper.config_mpl)
    
    
for j in range(len(alll0)):
    filename0=alll0[j]
    filenameBO=alllBO[j]

    alldata0=load_pickle(filename0)
    mpc0=alldata0.get("mpc")
    
    alldataBO=load_pickle(filenameBO)
    mpcBO=alldataBO.get("mpc")
 
    sim_time=mpc0._x.shape[0]
    
    
    if j==0:
        figDO, axDO = plt.subplots(2,2, sharex='col')
 
    axDO[0,0].plot(mpc0._aux[:,get_ind_aux(model, 'DO')]*100,alpha=0.15, color='b') 
    axDO[1,0].plot(mpcBO._aux[:,get_ind_aux(model, 'DO')]*100,alpha=0.15, color='b') 

    axDO[0,1].plot(mpc0._aux[:,get_ind_aux(model, 'DO')]*100,alpha=0.15, color='b') 
    axDO[1,1].plot(mpcBO._aux[:,get_ind_aux(model, 'DO')]*100,alpha=0.15, color='b') 
    
axDO[0,0].hlines(30, 0, int(np.max(mpc0._time)), linestyle='--', color='r')
axDO[1,0].hlines(30, 0, int(np.max(mpc0._time)), linestyle='--', color='r')
axDO[0,1].hlines(30, 0, int(np.max(mpc0._time)), linestyle='--', color='r')
axDO[1,1].hlines(30, 0, int(np.max(mpc0._time)), linestyle='--', color='r')




axDO[0,0].set_ylabel(r'$\bm{\Delta_{\tilde{\Theta},\text{man}}}$''\n'r'$DO [\%]$', rotation=90, ha="center")
axDO[1,0].set_ylabel(r'$\bm{\Delta_{\tilde{\Theta},\text{BO}}}$''\n'r'$DO [\%]$', rotation=90, ha="center")


axDO[-1,0].set_xlabel(r'$t\ [\text{h}]$')
axDO[-1,1].set_xlabel(r'$t\ [\text{h}]$')

y0=29.5
y1=40
x0=22
x1=40
axDO[0,1].set_ylim(y0,y1)
axDO[1,1].set_ylim(y0,y1)
axDO[0,1].set_xlim(x0,x1)
plt.setp(axDO[0,0].spines.values(), lw=0.5);
plt.setp(axDO[0,1].spines.values(), lw=0.5);
plt.setp(axDO[1,0].spines.values(), lw=0.5);
plt.setp(axDO[1,1].spines.values(), lw=0.5);


#create rectangle for zooming in
axDO[0,0].hlines(y0, x0, x1, linestyle='--', color='k')
axDO[0,0].hlines(y1, x0, x1, linestyle='--', color='k')
axDO[0,0].vlines(x0, y0, y1, linestyle='--', color='k')
axDO[0,0].vlines(x1, y0, y1, linestyle='--', color='k')

axDO[1,0].hlines(y0, x0, x1, linestyle='--', color='k')
axDO[1,0].hlines(y1, x0, x1, linestyle='--', color='k')
axDO[1,0].vlines(x0, y0, y1, linestyle='--', color='k')
axDO[1,0].vlines(x1, y0, y1, linestyle='--', color='k')

 
dx=16
dy=20
dy=0


axDO[1,0].annotate('', (x1,(y0+y1)/2), xytext=(x1+dx,(y0+y1)/2+dy), xycoords='data', textcoords=None, arrowprops=dict(arrowstyle= '<|-', color='k',  lw=1, ls='-'), annotation_clip=None)
axDO[0,0].annotate('', (x1,(y0+y1)/2), xytext=(x1+dx,(y0+y1)/2+dy), xycoords='data', textcoords=None, arrowprops=dict(arrowstyle= '<|-', color='k', lw=1, ls='-'), annotation_clip=None)


plt.subplots_adjust(left=0.19, bottom=0.18, right=0.98, top=0.96, wspace=0.22, hspace=0.2 )