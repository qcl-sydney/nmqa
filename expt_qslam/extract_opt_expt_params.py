import numpy as np

datakey=1
IDX1_SHP = 30
IDX2_SHP = 30 
path = '/scratch/QCL_RG/qslamdatapaper_v3/'

meta_ssim_pairs_all = []
meta_empr_pairs_all = []

for idx_job_array in range(1, IDX1_SHP*IDX2_SHP + 1):
    
    idx_1, idx_2 = np.unravel_index(idx_job_array - 1 , (IDX1_SHP, IDX2_SHP) )
    param_index = str(idx_1)+'_'+str(idx_2)
    filename = '2019_Jun_qslam_exptdata_'+str(datakey)+'_param_'+str(param_index)+'.npz'
    dataobj = np.load(path+filename)
    
    meta_ssim_pairs_all.append(dataobj["meta_ssim_pairs"])
    meta_empr_pairs_all.append(dataobj["meta_empr_pairs"])
    

np.savez(path+'2019_Jun_qslam_exptdata_collatedresults',
         meta_ssim_pairs_all=meta_ssim_pairs_all,
         meta_empr_pairs_all=meta_empr_pairs_all)

for idx_algo in range(2):
    
    opt_ssim_idx = np.argsort(np.mean(meta_ssim_pairs_all[:, idx_algo, -3:], axis=1))
    opt_empr_idx = np.argsort(np.mean(meta_empr_pairs_all[:, idx_algo, -3:], axis=1))
    
    
    print('Numerical optimisation in the high msmt regime')
    
    print('SSIM for algorithm: ', idx_algo)
    print('Top ten parameters:', opt_ssim_idx[0:10])
    print('Best parameter:', opt_ssim_idx[0])
    print('Min Error for last 3 max-iter datapoints: ', np.mean(meta_ssim_pairs_all[opt_ssim_idx[0], idx_algo, -3:], axis=1))
    print()
    print('SSIM for algorithm: ', idx_algo)
    print('Top ten parameters:', opt_empr_idx[0:10])
    print('Best parameter:',opt_empr_idx[0])
    print('with min Error for last 3  max-iter datapoints: ', np.mean(meta_empr_pairs_all[opt_empr_idx[0], idx_algo, -3:], axis=1))
    print()
    print()
    
