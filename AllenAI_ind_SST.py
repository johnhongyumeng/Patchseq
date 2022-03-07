# -*- coding: utf-8 -*-
"""
Created on Sun Jan 23 16:51:32 2022
Analyze the Allen's data to generate correct AI, and rABC. Just measure for every 
qualified cells. Generate a CSV file as output
Please Notice this code needs to run under py3.6. Some functions are not capatible with 3.8, like spike_df
@author: John
"""

#%% import  

# use py36, not py3
from ipfx.dataset.create import create_ephys_data_set
from ipfx.utilities import drop_failed_sweeps
from ipfx.epochs import get_stim_epoch
from ipfx.stim_features import get_stim_characteristics
from ipfx.feature_extractor import SpikeFeatureExtractor
from ipfx import spike_train_features
import statistics

import numpy as np
import pandas as pd
import scipy as sp

import matplotlib.pyplot as plt
import seaborn as sns


#%%
# read in dataframe containing metadata and nwb-related information
pv_metadata_nwbs = pd.read_csv("L23_sst_metadata_nwbs.csv", index_col=0)
# generate a dataframe containing the unique index and nwb file name
pv_nwb_filelist = pv_metadata_nwbs[['index', 'file_name']]    #+ '%0D'
# add '%0D' at the end of each nwb files' names
#pv_nwb_filelist.loc[:, 'file_name'] = pv_metadata_nwbs[['file_name']] + '%0D' # for later sweep selecting

pv_metadata_nwbs['file_name']=pv_metadata_nwbs['file_name'].str.replace('nwb','nwb%0D')
index_col = pv_metadata_nwbs[['index']]
nwb_filelist = pv_metadata_nwbs[['file_name']]
n_cells = len(nwb_filelist)

#%% Now plotting for one cell as an example. Later generate a loop to do the job.
i_cell=14;
cell_nwb_filename = 'SSTL23/'+nwb_filelist['file_name'][i_cell]
data_set = create_ephys_data_set(nwb_file = cell_nwb_filename)
drop_failed_sweeps(data_set)     # The function of these following lines are not making sense to me yet. 
long_square_table = data_set.filtered_sweep_table(stimuli=data_set.ontology.long_square_names) # get sweep table for Long Square sweeps
long_square_sweeps = data_set.sweep_set(long_square_table.sweep_number)
long_square_sweeps.select_epoch("recording")
long_square_sweeps.align_to_start_of_epoch("experiment")
sweep=long_square_sweeps.sweeps[0]
n_sweep=len(long_square_sweeps.sweeps)


#%%  Test on the single spike
i_sweep=13
sweep=long_square_sweeps.sweeps[i_sweep]
t_start, duration, i_amp, _, _ = get_stim_characteristics(sweep.i, sweep.t)
spfx = SpikeFeatureExtractor(start = t_start, end = (t_start+duration), min_peak = 20,dv_cutoff=15)
spike_df = spfx.process(sweep.t, sweep.v, sweep.i)

# testing
plt.figure(figsize=(12,6))
plt.plot(sweep.t,sweep.v)

vmax_sweep=max(sweep.v)
Vthr_sweep=max(vmax_sweep-20,0)
flag_spike=(sweep.v[0:-1]<Vthr_sweep) & (sweep.v[1:]>=Vthr_sweep) & (sweep.t[:-1]>t_start) & (sweep.t[:-1]<(t_start+duration))# get the spikes
t_spike= sweep.t[np.append(flag_spike,False)]- t_start
if len(t_spike)<10:
    print('Not enough spike for the highest injection')
else:
    t_isi=np.diff(t_spike)
    FI_vec=1/t_isi       
    # for testing the figure.
    plt.figure(figsize=(6,6))
    plt.plot(t_spike[:-1],FI_vec)


#spike_df = spfx.process(sweep.t, sweep.v, sweep.i)
#%% Define functions     
def model_func(t, A, K, C):
    return A * np.exp(K * t) + C

def fit_exp_nonlinear(t, y,p0=[-5, 0.01, 25]):
    opt_parms, parm_cov = sp.optimize.curve_fit(model_func, t, y, p0, maxfev=10000)
    A, K, C = opt_parms
    return A, K, C

#%%

A, K, C = fit_exp_nonlinear(t_spike[:-1], FI_vec, [2*FI_vec[-1], -1, FI_vec[-1]]) 
t_axis=np.linspace(0.0,1.0,num=10001,endpoint=True)
fit_y = model_func(t_axis, A, K, C)
plt.figure(figsize=(6,6))
plt.plot(t_spike[:-1],FI_vec)
plt.plot(t_axis, fit_y, '-',linewidth=1) 

#%%
AI_sweep=1-fit_y[-1]/fit_y[0]
AI_sweep
ABC_vec=-fit_y+ (fit_y[-1]+  (fit_y[0]-fit_y[-1])*(1-t_axis))
rABC=sum(ABC_vec)/len(ABC_vec)/(fit_y[0]-fit_y[-1])*2
rABC

