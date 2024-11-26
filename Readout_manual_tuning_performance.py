from do_mpc.tools import load_pickle
import numpy as np

folder='BO_manual_multi'

    
# =============================================================================
# ALL 0
# =============================================================================

plan_B0 = load_pickle('./BO/'+str(folder)+'/sample_plan/BO_plan_0.pkl')
seeds=np.zeros([int(plan_B0[-1]['id'])+1])

for i in range(seeds.shape[0]):
    seeds[i]=plan_B0[i]['this_is_my_seed']
   

res_B0=np.zeros([6,seeds.shape[0]])
res_B0_branch=np.zeros([3,seeds.shape[0]])


for i in range(10):
    res_B0[0:5,i] = load_pickle('./BO/'+str(folder)+'/sample_results_0/sample_00'+str(i)+'.pkl')[0:5]
    res_B0_branch[:,i] = load_pickle('./BO/'+str(folder)+'/sample_results_0/sample_00'+str(i)+'.pkl')[-1]
    res_B0[5,i]=np.sum(res_B0_branch[:,i]**2)
    
for i in range(10,seeds.shape[0]):
    res_B0[0:5,i] = load_pickle('./BO/'+str(folder)+'/sample_results_0/sample_0'+str(i)+'.pkl')[0:5]
    res_B0_branch[:,i] = load_pickle('./BO/'+str(folder)+'/sample_results_0/sample_0'+str(i)+'.pkl')[-1]
    res_B0[5,i]=np.sum(res_B0_branch[:,i]**2)


runs_viol_B0=np.where(res_B0[2,:]>0)[0].shape[0]
viol_mean_B0=res_B0[2,:].mean()
viol_sd_B0=res_B0[2,:].std()
t99_mean_B0=res_B0[1,:].mean()
t99_sd_B0=res_B0[1,:].std()

res_B0_branches=np.array([res_B0_branch[0][0],res_B0_branch[1][0],res_B0_branch[2][0]])



v_total_B0=res_B0[3,:]+res_B0[4,:] #ub + lb
v_viol_mean_B0=v_total_B0.mean()
v_viol_sd_B0=v_total_B0.std()


print()
print('branching:'+str(res_B0_branches))
print('runs with violations:'+str(runs_viol_B0)+'  %:'+str(np.round((runs_viol_B0/seeds.shape[0]*100),2)))
print('mean number violations per run:'+str(viol_mean_B0.round(2))+'+-'+str(viol_sd_B0.round(2)))
print('mean value violation value per run:'+str(v_viol_mean_B0))
print('mean t99:'+str(t99_mean_B0.round(2))+'+-'+str(t99_sd_B0.round(2)))



# =============================================================================
# ALL 3
# =============================================================================

res_B3=np.zeros([6,seeds.shape[0]])
res_B3_branch=np.zeros([3,seeds.shape[0]])


for i in range(10):
    res_B3[0:5,i] = load_pickle('./BO/'+str(folder)+'/sample_results_1/sample_00'+str(i)+'.pkl')[0:5]
    res_B3_branch[:,i] = load_pickle('./BO/'+str(folder)+'/sample_results_1/sample_00'+str(i)+'.pkl')[-1]
    res_B3[5,i]=np.sum(res_B3_branch[:,i]**2)
    
for i in range(10,seeds.shape[0]):
    res_B3[0:5,i] = load_pickle('./BO/'+str(folder)+'/sample_results_1/sample_0'+str(i)+'.pkl')[0:5]
    res_B3_branch[:,i] = load_pickle('./BO/'+str(folder)+'/sample_results_1/sample_0'+str(i)+'.pkl')[-1]
    res_B3[5,i]=np.sum(res_B3_branch[:,i]**2)


runs_viol_B3=np.where(res_B3[2,:]>0)[0].shape[0]
viol_mean_B3=res_B3[2,:].mean()
viol_sd_B3=res_B3[2,:].std()
t99_mean_B3=res_B3[1,:].mean()
t99_sd_B3=res_B3[1,:].std()

res_B3_branches=np.array([res_B3_branch[0][0],res_B3_branch[1][0],res_B3_branch[2][0]])

v_total_B3=res_B3[3,:]+res_B3[4,:] #ub + lb
v_viol_mean_B3=v_total_B3.mean()
v_viol_sd_B3=v_total_B3.std()


print()
print('branching:'+str(res_B3_branches))
print('runs with violations:'+str(runs_viol_B3)+'  %:'+str(np.round((runs_viol_B3/seeds.shape[0]*100),2)))
print('mean number violations per run:'+str(viol_mean_B3.round(2))+'+-'+str(viol_sd_B3.round(2)))
print('mean value violation value per run:'+str(v_viol_mean_B3))
print('mean t99:'+str(t99_mean_B3.round(2))+'+-'+str(t99_sd_B3.round(2)))



# =============================================================================
# ALL 6
# =============================================================================


res_B6=np.zeros([6,seeds.shape[0]])
res_B6_branch=np.zeros([3,seeds.shape[0]])


for i in range(10):
    res_B6[0:5,i] = load_pickle('./BO/'+str(folder)+'/sample_results_2/sample_00'+str(i)+'.pkl')[0:5]
    res_B6_branch[:,i] = load_pickle('./BO/'+str(folder)+'/sample_results_2/sample_00'+str(i)+'.pkl')[-1]
    res_B6[5,i]=np.sum(res_B6_branch[:,i]**2)
    
for i in range(10,seeds.shape[0]):
    res_B6[0:5,i] = load_pickle('./BO/'+str(folder)+'/sample_results_2/sample_0'+str(i)+'.pkl')[0:5]
    res_B6_branch[:,i] = load_pickle('./BO/'+str(folder)+'/sample_results_2/sample_0'+str(i)+'.pkl')[-1]
    res_B6[5,i]=np.sum(res_B6_branch[:,i]**2)


runs_viol_B6=np.where(res_B6[2,:]>0)[0].shape[0]
viol_mean_B6=res_B6[2,:].mean()
viol_sd_B6=res_B6[2,:].std()
t99_mean_B6=res_B6[1,:].mean()
t99_sd_B6=res_B6[1,:].std()

res_B6_branches=np.array([res_B6_branch[0][0],res_B6_branch[1][0],res_B6_branch[2][0]])

v_total_B6=res_B6[3,:]+res_B6[4,:] #ub + lb
v_viol_mean_B6=v_total_B6.mean()
v_viol_sd_B6=v_total_B6.std()

print()
print('branching:'+str(res_B6_branches))
print('runs with violations:'+str(runs_viol_B6)+'  %:'+str(np.round((runs_viol_B6/seeds.shape[0]*100),2)))
print('mean number violations per run:'+str(viol_mean_B6.round(2))+'+-'+str(viol_sd_B6.round(2)))
print('mean value violation value per run:'+str(v_viol_mean_B6))
print('mean t99:'+str(t99_mean_B6.round(2))+'+-'+str(t99_sd_B6.round(2)))


# =============================================================================
# ALL 7
# =============================================================================
folder='BO_manual_7'

res_B7=np.zeros([6,seeds.shape[0]])
res_B7_branch=np.zeros([3,seeds.shape[0]])


for i in range(10):
    res_B7[0:5,i] = load_pickle('./BO/'+str(folder)+'/sample_results/sample_00'+str(i)+'.pkl')[0:5]
    res_B7_branch[:,i] = load_pickle('./BO/'+str(folder)+'/sample_results/sample_00'+str(i)+'.pkl')[-1]
    res_B7[5,i]=np.sum(res_B7_branch[:,i]**2)
    
for i in range(10,seeds.shape[0]):
    res_B7[0:5,i] = load_pickle('./BO/'+str(folder)+'/sample_results/sample_0'+str(i)+'.pkl')[0:5]
    res_B7_branch[:,i] = load_pickle('./BO/'+str(folder)+'/sample_results/sample_0'+str(i)+'.pkl')[-1]
    res_B7[5,i]=np.sum(res_B7_branch[:,i]**2)


runs_viol_B7=np.where(res_B7[2,:]>0)[0].shape[0]
viol_mean_B7=res_B7[2,:].mean()
viol_sd_B7=res_B7[2,:].std()
t99_mean_B7=res_B7[1,:].mean()
t99_sd_B7=res_B7[1,:].std()

res_B7_branches=np.array([res_B7_branch[0][0],res_B7_branch[1][0],res_B7_branch[2][0]])

v_total_B7=res_B7[3,:]+res_B7[4,:] #ub + lb
v_viol_mean_B7=v_total_B7.mean()
v_viol_sd_B7=v_total_B7.std()


print()
print('branching:'+str(res_B7_branches))
print('runs with violations:'+str(runs_viol_B7)+'  %:'+str(np.round((runs_viol_B7/seeds.shape[0]*100),2)))
print('mean number violations per run:'+str(viol_mean_B7.round(2))+'+-'+str(viol_sd_B7.round(2)))
print('mean value violation value per run:'+str(v_viol_mean_B7))
print('mean t99:'+str(t99_mean_B7.round(2))+'+-'+str(t99_sd_B7.round(2)))

# =============================================================================
# ALL 8
# =============================================================================
folder='BO_manual_8'

res_B8=np.zeros([6,seeds.shape[0]])
res_B8_branch=np.zeros([3,seeds.shape[0]])


for i in range(10):
    res_B8[0:5,i] = load_pickle('./BO/'+str(folder)+'/sample_results/sample_00'+str(i)+'.pkl')[0:5]
    res_B8_branch[:,i] = load_pickle('./BO/'+str(folder)+'/sample_results/sample_00'+str(i)+'.pkl')[-1]
    res_B8[5,i]=np.sum(res_B8_branch[:,i]**2)
    
for i in range(10,seeds.shape[0]):
    res_B8[0:5,i] = load_pickle('./BO/'+str(folder)+'/sample_results/sample_0'+str(i)+'.pkl')[0:5]
    res_B8_branch[:,i] = load_pickle('./BO/'+str(folder)+'/sample_results/sample_0'+str(i)+'.pkl')[-1]
    res_B8[5,i]=np.sum(res_B8_branch[:,i]**2)


runs_viol_B8=np.where(res_B8[2,:]>0)[0].shape[0]
viol_mean_B8=res_B8[2,:].mean()
viol_sd_B8=res_B8[2,:].std()
t99_mean_B8=res_B8[1,:].mean()
t99_sd_B8=res_B8[1,:].std()

res_B8_branches=np.array([res_B8_branch[0][0],res_B8_branch[1][0],res_B8_branch[2][0]])

v_total_B8=res_B8[3,:]+res_B8[4,:] #ub + lb
v_viol_mean_B8=v_total_B8.mean()
v_viol_sd_B8=v_total_B8.std()


print()
print('branching:'+str(res_B8_branches))
print('runs with violations:'+str(runs_viol_B8)+'  %:'+str(np.round((runs_viol_B8/seeds.shape[0]*100),2)))
print('mean number violations per run:'+str(viol_mean_B8.round(2))+'+-'+str(viol_sd_B8.round(2)))
print('mean value violation value per run:'+str(v_viol_mean_B8))
print('mean t99:'+str(t99_mean_B8.round(2))+'+-'+str(t99_sd_B8.round(2)))



