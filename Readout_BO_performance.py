from do_mpc.tools import load_pickle
import numpy as np


folder='BO_full'

# =============================================================================
# DATA PREPARATION
# =============================================================================

plan_B0 = load_pickle('./BO/'+str(folder)+'/sample_plan_3D_c/BO_plan_3D_c_0.pkl')
seeds=np.zeros([int(plan_B0[-1]['id'])+1])

for i in range(seeds.shape[0]):
    seeds[i]=plan_B0[i]['this_is_my_seed']
   

res_BB=np.zeros([6,seeds.shape[0]])
res_BB_branch=np.zeros([3,seeds.shape[0]])


# =============================================================================
# READING DATA FROM BEST BO RESULT
# =============================================================================
num_BO=88
for i in range(10):
    res_BB[0:5,i] = load_pickle('./BO/'+str(folder)+'/sample_res_3D_c_opt'+str(num_BO)+'/sample_00'+str(i)+'.pkl')[0:5]
    res_BB_branch[:,i] = load_pickle('./BO/'+str(folder)+'/sample_res_3D_c_opt'+str(num_BO)+'/sample_00'+str(i)+'.pkl')[-1]
    res_BB[5,i]=np.sum(res_BB_branch[:,i]**2)
    
for i in range(10,seeds.shape[0]):
    res_BB[0:5,i] = load_pickle('./BO/'+str(folder)+'/sample_res_3D_c_opt'+str(num_BO)+'/sample_0'+str(i)+'.pkl')[0:5]
    res_BB_branch[:,i] = load_pickle('./BO/'+str(folder)+'/sample_res_3D_c_opt'+str(num_BO)+'/sample_0'+str(i)+'.pkl')[-1]
    res_BB[5,i]=np.sum(res_BB_branch[:,i]**2)


runs_viol_BB=np.where(res_BB[2,:]>0)[0].shape[0]
viol_mean_BB=res_BB[2,:].mean()
viol_sd_BB=res_BB[2,:].std()
t99_mean_BB=res_BB[1,:].mean()
t99_sd_BB=res_BB[1,:].std()

res_BB_branches=np.array([res_BB_branch[0][0],res_BB_branch[1][0],res_BB_branch[2][0]])



v_total_BB=res_BB[3,:]+res_BB[4,:] #ub + lb
v_viol_mean_BB=v_total_BB.mean()
v_viol_sd_BB=v_total_BB.std()


print()
print('branching:'+str(res_BB_branches))
#print('mean cost:'+str(res_BB[0].mean()))
print('runs with violations:'+str(runs_viol_BB)+'  %:'+str(np.round((runs_viol_BB/seeds.shape[0]*100),2)))
print('mean number violations per run:'+str(viol_mean_BB.round(2))+'+-'+str(viol_sd_BB.round(2)))
print('mean value violation value per run:'+str(v_viol_mean_BB))
print('mean t99:'+str(t99_mean_BB.round(2))+'+-'+str(t99_sd_BB.round(2)))

