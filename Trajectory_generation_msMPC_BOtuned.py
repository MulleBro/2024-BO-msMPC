from bayes_opt import BayesianOptimization
from bayes_opt import UtilityFunction
from bayes_opt.logger import JSONLogger
from bayes_opt.event import Events
from bayes_opt.util import load_logs

import multiprocessing as mp
from functools import partial
from do_mpc.tools import load_pickle
import sys
sys.path.append('../../../../')
import do_mpc
import time

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

from model.model import template_model

sys.path.append('../../../../')
import multiprocessing as mp
from functools import partial
from do_mpc.tools import load_pickle
from sklearn.gaussian_process.kernels import RBF, Matern


# =============================================================================
# INPUTS
# =============================================================================
n_inits = 1 #initialization points per dimension
n_iter=0 #BO iterations after initialization
kappa_value=1.96

n_runs=100 # number of investigated sets of parameters (drawn from uniform distribution)

# =============================================================================
# INITIALIZATION POINTS manually set to BO result
# =============================================================================
#kX
p_inits=np.array([8.61300587])

#kP
p_inits1=np.array([0.03263058])

#Ks_Xyl_X
p_inits2=np.array([6.00862537])

# =============================================================================
# LOAD SEEDS FOR SAME SETS OF PARAMETERS FROM UNIFORM DISTRIBUTION
# =============================================================================
plan_seeds = load_pickle('./sample_plan/BO_plan_seeds_100.pkl')

seeds=np.zeros([int(plan_seeds[-1]['id'])+1])

for i in range(seeds.shape[0]):
    seeds[i]=plan_seeds[i]['this_is_my_seed']



# =============================================================================
#     INTIALIZATION OPTIMIZER AND UTILITY_FUNCTION
# =============================================================================

p_branch=['kX', 'kP', 'Ks_Xyl_X']

optimizer = BayesianOptimization(
    f=None,
    pbounds={'kX': (lb_x, ub_x), 'kP': (lb_x1, ub_x1), 'Ks_Xyl_X': (lb_x2, ub_x2)},
    verbose=2,
    random_state=1,
)

optimizer.set_gp_params(kernel=RBF())

utility = UtilityFunction(kind="ucb", kappa=kappa_value)

logger = JSONLogger(path="./BO/logs/BO_log_manualBO")
optimizer.subscribe(Events.OPTIMIZATION_STEP, logger)

#%%
# =============================================================================
# MPC FUNCTION FOR BAYESIAN OPTIMIZATION
# =============================================================================
def mpc_(del_mult_proz,this_is_my_seed):
        
        #n_horizon=10.3
        """ User settings: """
        #n_sim=72 # steps, nicht t_total
        n_sim=72
        dt=1
        np.random.seed(this_is_my_seed)
        
        model = template_model()
        def template_mpc(model):
            mpc = do_mpc.controller.MPC(model)
            

            max_it=2000
            
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
            mpc.bounds['upper', '_x', 'V_L'] = 230*1e-3 # L
            mpc.bounds['lower','_u','F_XylN'] = 0.0
            mpc.bounds['upper','_u','F_XylN'] = 5*1e-3 # L/h


            mpc.set_nl_cons('g_DO', -(model.x['O2_L']/O2_L_max), ub=-0.3, soft_constraint=True, penalty_term_cons=1e8)
            mpc.set_nl_cons('g_Xyl', (model.x['Xyl']), ub=120, soft_constraint=True, penalty_term_cons=1e8)
       
            
            P_df=pd.read_csv('./p_est/res/p_opt.csv')   
            P=np.array(P_df['val'])
            
            

            #%%

# =============================================================================
#         multi-stage branching
# =============================================================================
                    
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
    
           
        """
        Get configured do-mpc modules:
        """
        
        mpc = template_mpc(model)
        #simulator = template_simulator(model)
        estimator = do_mpc.estimator.StateFeedback(model)
        model_list(model)

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
                

        P_df=pd.read_csv('./p_est/res/p_opt.csv')  
        P=np.array(P_df['val'])
        
        #Inputs
        #F_Xyl
        u0=np.array([ 0]).reshape(-1,1) #ml/h

        
        # =============================================================================
        # SIMULATOR
        # =============================================================================
        
        simulator = do_mpc.simulator.Simulator(model)
        
        params_simulator = {
            'integration_tool': 'idas',
            'abstol': 1e-14,
            'reltol': 1e-14,
            't_step': dt
        }
        
        simulator.set_param(**params_simulator)
        
        
        p_num = simulator.get_p_template()
        p_names=model_names(model)[1]
        

        #const random value
        for i in range(model.n_p):
            p_num[p_names[i]]=np.random.uniform(0.95,1.05)*P[get_ind_p(model,p_names[i])]
        

        def p_fun(t_now):
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
        mpc.set_initial_guess() #initial guess using x0, z0, u0
        
        
        # =============================================================================
        # MAIN LOOP MPC
        # =============================================================================
        
        
        t = time.time()
        failure=0
        # do stuff
        
        u0 = mpc.make_step(x0)
        y_next = simulator.make_step(u0)
        x0 = estimator.make_step(y_next)
        
            
        while (np.max(mpc.data._aux[:,get_ind_aux(model, 'X_Lys')])<0.99 and mpc.data._aux[:,get_ind_aux(model, 'X_Lys')].shape[0]<n_sim):
        
            try:
                u0 = mpc.make_step(x0)
                y_next = simulator.make_step(u0)
                x0 = estimator.make_step(y_next)
                
            except RuntimeError:
                failure=1
                break
            
        elapsed = time.time() - t  
         

    # =============================================================================
    #     CONSTRAINT VIOLATION X
    # =============================================================================
        if failure!=1:
            n_mpc=mpc.data._aux[:,get_ind_aux(model, 'X_Lys')].shape[0]
            lbx=np.zeros([n_mpc,model.n_x]) 
            #lbx[:,get_ind_x(model, 'V_L')]=90*1e-3 #LEAVE OUT _> CALC EASIER
            lbx[:,get_ind_x(model, 'O2_L')]=O2_L_max*0.3
            
            ubx=np.ones([n_mpc,model.n_x])*1e10
            ubx[:,get_ind_x(model, 'Xyl')] = 120.0 #mM
            
            
            lbx_viol=mpc.data._x<lbx
            n_viol_lb=np.sum(lbx_viol)
            
            ubx_viol=mpc.data._x>ubx
            n_viol_ub=np.sum(ubx_viol)
            
            n_viol=n_viol_lb+n_viol_ub
            
            #value penalized
            vlbx_viol=mpc.data._x-lbx
            v_viol_lb=vlbx_viol[lbx_viol]/0.3*100
            viol_lb=np.sum(np.abs(v_viol_lb))

           
            vubx_viol=ubx-mpc.data._x
            v_viol_ub=vubx_viol[ubx_viol]/120*100
            viol_ub=np.sum(np.abs(v_viol_ub))
            
            v_total=viol_ub+viol_lb
            
        # =============================================================================
        #     TIME UNTIL HLYS CONCENTRATION SUFFICIENT
        # =============================================================================
            X_Lys=mpc.data._aux[:,get_ind_aux(model, 'X_Lys')]
            t_99=np.argmax(X_Lys>=0.99)
            if t_99!=0:
                t_99=int(mpc.data._time[t_99][0])
                #ADDED INTEROPLATION
                X99=X_Lys[t_99]
                X99_1=X_Lys[t_99-1]

                t_99=t_99-1+(0.99-X99_1)/(X99-X99_1) 
            else:
                t_99=100 #penalty if goal not reached
                
            
        # =============================================================================
        #             SAVING TRAJECTORIES
        # =============================================================================
            do_mpc.data.save_results([mpc, simulator], '100_nomBO')

    
        # =============================================================================
        # COST
        # =============================================================================

            
            JJ_=-(t_99+(1e3*v_total)+np.sum((np.array(del_mult_proz))))

        else: # filtering out runs with runtime errors
            JJ_=np.nan
        return JJ_,t_99,n_viol,viol_ub,viol_lb,del_mult_proz
    
    
#%%

# =============================================================================
# BAYESIAN OPTIMIZATION
# =============================================================================
delli=np.zeros([len(p_branch)])


# =============================================================================
# INITIALIZATION PHASE
# =============================================================================
jjj=-1
for i in range(p_inits.shape[0]):
    for j in range(p_inits1.shape[0]):
        for k in range(p_inits2.shape[0]):

                jjj+=1

            
                delli[0]=p_inits[i]
                delli[1]=p_inits1[j]
                delli[2]=p_inits2[k]

                
                dict_p = {}
                for dd in range(len(p_branch)):
                    dict_p[p_branch[dd]] = delli[dd]
                    
                print(dict_p)
                
                
                sp = do_mpc.sampling.SamplingPlanner()
                sp.set_param(overwrite = True)
                sp.data_dir = './BO/sample_plan_3D_c/'
                
                

                
                
                sp.set_sampling_var('del_mult_proz', lambda: [delli[0],delli[1],delli[2]])
                sp.set_sampling_var('this_is_my_seed', lambda: int(seeds[-1]))
                
                for FF in range(n_runs-1):
                    sp.add_sampling_case(del_mult_proz=np.array([delli[0],delli[1],delli[2]]), this_is_my_seed=int(seeds[FF]))
                
                plan = sp.gen_sampling_plan(n_samples=1)
                
                sp.export('BO_plan_3D_c_'+str(jjj))

                sampler = do_mpc.sampling.Sampler(plan)
                sampler.set_param(overwrite=True)
                sampler.data_dir ='./BO/sample_results_3D_c_'+str(jjj)+'/'
                    
                    
                
                  
                
                  
                sampler.set_sample_function(mpc_)
                 
                with mp.Pool(processes=55) as pool:
                    p = pool.map(sampler.sample_idx, list(range(sampler.n_samples)))
                  
 
                  
                # DataHandling
                dh = do_mpc.sampling.DataHandler(plan)
                  
                dh.data_dir = './BO/sample_results_3D_c_'+str(jjj)+'/'
                
                dh_df=pd.DataFrame(dh[:])
                dh_df2=dh_df.dropna()
                
                
                resi=np.zeros(dh_df2.shape[0])
                
                for j_ind in range(resi.shape[0]):
                    resi[j_ind]=dh_df2.loc[:, 'res'][j_ind][0]
                
                if dh_df2.empty==False:
                    res_mean=resi.mean()
                    
                    print(dh_df2)
                    print(res_mean)
                    
                    optimizer.register(
                        params=dict_p,
                        target=res_mean,
                    )





# =============================================================================
# OPTIMIZATION PHASE: NOT PRESENT, only generating trajectories
# =============================================================================

