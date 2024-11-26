
import numpy as np
from casadi import *
from casadi.tools import *
import pdb
import sys
import os
rel_do_mpc_path = os.path.join('..','..')
sys.path.append(rel_do_mpc_path)
import do_mpc


from helper.aux_functions import *
from helper.const_params import *

import pandas as pd



def template_mpc(model):
    """
    --------------------------------------------------------------------------
    template_mpc: tuning parameters
    --------------------------------------------------------------------------
    """
    mpc = do_mpc.controller.MPC(model)
    
    max_it=3000
    
            
    #max_it=1000
    #max_it=2000
    
    del_mult_proz=np.array(([7.68],[0.60],[7.58]))
    
        
     
    
    setup_mpc = {
        'n_horizon': 12,
        'n_robust': 1,
        'open_loop': 0,
        't_step': 1,
        'state_discretization': 'collocation',
        'collocation_type': 'radau',
        'collocation_deg': 2,
        'collocation_ni': 2,
        'store_full_solution': True,
        #Use MA27 linear solver in ipopt for faster calculations:
        'nlpsol_opts': {'ipopt.linear_solver': 'ma27','ipopt.max_iter':max_it,'ipopt.print_level':3}#,'ipopt.acceptable_tol':5*1e-2, 'ipopt.tol':5*1e-2}
    }
        
        
     
    mpc.set_param(**setup_mpc)
    
     
    w_P=1
    w_inp_pen=1
      
    mterm = -1e2*model.x['HLys']
      
    # #stage cost
      
    lterm = -1e2*model.x['HLys']+w_inp_pen*(model.u['F_XylN'])
      
      
    mpc.set_objective(mterm=mterm, lterm=lterm)
    mpc.set_rterm( F_XylN=1)
    # =============================================================================
    #   Boundaries
    # =============================================================================
    P_df=pd.read_csv('./p_est/res/p_opt.csv')   
    P=np.array(P_df['val'])
      
    #States
    x_names=model_names(model)[0]
    for i in range(len(x_names)):
        mpc.bounds['lower', '_x', x_names[i]] = 0.0
      
    mpc.bounds['lower', '_x', 'V_L'] = 90*1e-3 #L
    #mpc.bounds['lower', '_x', 'O2_L'] = O2_L_max*0.3 #mM #DO MINIMUM 0.3
    #mpc.bounds['lower','_aux','DO'] = 0.3
      
      
    mpc.bounds['upper', '_x', 'V_L'] = 230*1e-3 # L
    #mpc.bounds['upper', '_x', 'Xyl'] = 120.0 #mM
      
      
      
    mpc.bounds['lower','_u','F_XylN'] = 0.0
    mpc.bounds['upper','_u','F_XylN'] = 5*1e-3 # L/h
      
      
    mpc.set_nl_cons('g_DO', -(model.x['O2_L']/O2_L_max), ub=-0.3, soft_constraint=True, penalty_term_cons=1e8)
      
    mpc.set_nl_cons('g_Xyl', (model.x['Xyl']), ub=120, soft_constraint=True, penalty_term_cons=1e8)
      
    
    P_df=pd.read_csv('./p_est/res/p_opt.csv')   
    P=np.array(P_df['val'])
    
    
      
    #%%
      
      
            
    P_1=np.copy(P) #a+ b+ c+
    P_1[get_ind_p(model,'kP')]=P[get_ind_p(model,'kP')]*(1+del_mult_proz[1]/100)
    P_1[get_ind_p(model,'kX')]=P[get_ind_p(model,'kX')]*(1+del_mult_proz[0]/100)
    P_1[get_ind_p(model,'Ks_Xyl_X')]=P[get_ind_p(model,'Ks_Xyl_X')]*(1+del_mult_proz[2]/100)
    
    
    P_2=np.copy(P) #a- b- c+
    P_2[get_ind_p(model,'kP')]=P[get_ind_p(model,'kP')]*(1-del_mult_proz[1]/100)
    P_2[get_ind_p(model,'kX')]=P[get_ind_p(model,'kX')]*(1-del_mult_proz[0]/100)
    P_2[get_ind_p(model,'Ks_Xyl_X')]=P[get_ind_p(model,'Ks_Xyl_X')]*(1+del_mult_proz[2]/100)
    
    P_3=np.copy(P) #a- b+ c+
    P_3[get_ind_p(model,'kP')]=P[get_ind_p(model,'kP')]*(1-del_mult_proz[1]/100)
    P_3[get_ind_p(model,'kX')]=P[get_ind_p(model,'kX')]*(1+del_mult_proz[0]/100)
    P_3[get_ind_p(model,'Ks_Xyl_X')]=P[get_ind_p(model,'Ks_Xyl_X')]*(1+del_mult_proz[2]/100)
    
    P_4=np.copy(P) #a+ b- c+
    P_4[get_ind_p(model,'kP')]=P[get_ind_p(model,'kP')]*(1+del_mult_proz[1]/100)
    P_4[get_ind_p(model,'kX')]=P[get_ind_p(model,'kX')]*(1-del_mult_proz[0]/100)
    P_4[get_ind_p(model,'Ks_Xyl_X')]=P[get_ind_p(model,'Ks_Xyl_X')]*(1+del_mult_proz[2]/100)
    
    
    P_5=np.copy(P) #a+ b+ c-
    P_5[get_ind_p(model,'kP')]=P[get_ind_p(model,'kP')]*(1+del_mult_proz[1]/100)
    P_5[get_ind_p(model,'kX')]=P[get_ind_p(model,'kX')]*(1+del_mult_proz[0]/100)
    P_5[get_ind_p(model,'Ks_Xyl_X')]=P[get_ind_p(model,'Ks_Xyl_X')]*(1-del_mult_proz[2]/100)
    
    P_6=np.copy(P) #a- b- c-
    P_6[get_ind_p(model,'kP')]=P[get_ind_p(model,'kP')]*(1-del_mult_proz[1]/100)
    P_6[get_ind_p(model,'kX')]=P[get_ind_p(model,'kX')]*(1-del_mult_proz[0]/100)
    P_6[get_ind_p(model,'Ks_Xyl_X')]=P[get_ind_p(model,'Ks_Xyl_X')]*(1-del_mult_proz[2]/100)
    
    P_7=np.copy(P) #a- b+ c-
    P_7[get_ind_p(model,'kP')]=P[get_ind_p(model,'kP')]*(1-del_mult_proz[1]/100)
    P_7[get_ind_p(model,'kX')]=P[get_ind_p(model,'kX')]*(1+del_mult_proz[0]/100)
    P_7[get_ind_p(model,'Ks_Xyl_X')]=P[get_ind_p(model,'Ks_Xyl_X')]*(1-del_mult_proz[2]/100)
    
    P_8=np.copy(P) #a+ b- c-
    P_8[get_ind_p(model,'kP')]=P[get_ind_p(model,'kP')]*(1+del_mult_proz[1]/100)
    P_8[get_ind_p(model,'kX')]=P[get_ind_p(model,'kX')]*(1-del_mult_proz[0]/100)
    P_8[get_ind_p(model,'Ks_Xyl_X')]=P[get_ind_p(model,'Ks_Xyl_X')]*(1-del_mult_proz[2]/100)
    
       
       
       
     
    n_combinations = 9
    p_template = mpc.get_p_template(n_combinations)
    p_template['_p',0] = P
    p_template['_p',1] = P_1
    p_template['_p',2] = P_2
    p_template['_p',3] = P_3
    p_template['_p',4] = P_4
    p_template['_p',5] = P_5
    p_template['_p',6] = P_6
    p_template['_p',7] = P_7
    p_template['_p',8] = P_8
      
    
    
    def p_fun(t_now):
        return p_template
    
    mpc.set_p_fun(p_fun)
    
    mpc.setup()
    
    return mpc