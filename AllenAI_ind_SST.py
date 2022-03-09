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
i_cell=6;
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
i_sweep=17
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

#%% Try AllenSDK code to extract half-width

from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
sweep_ext = EphysSweepFeatureExtractor(t=sweep.t, v=sweep.v, i=sweep.i, start=t_start+0.0, end=t_start+duration+0.0)
sweep_ext.process_spikes()
#sweep_ext.process(keys="width")

#print("Avg spike threshold: %.01f mV" % sweep_ext.spike_feature("threshold_v").mean())
print("Avg spike width: %.02f ms" %  (1e3 * np.nanmean(sweep_ext.spike_feature("width"))))

plt.figure(figsize=(12,6))
plt.plot(sweep.t,sweep.v)
plt.plot(sweep.t,sweep.i)

#plt.xlim([1.17,1.175])
#plt.xlim([1.5,2.0])


sweep_ext.spike_feature("width")
np.nanmean(sweep_ext.spike_feature("width"))
#%% Now let's try to extract input resistance. Seems there isn't a way that just calculate the input resistance.
# Let me do it. Also for the fitting. Ignoring the sweep selection for now.

flag_dep=(sweep.t>=t_start+duration/2)  & (sweep.t<t_start+duration)
flag_rest=(sweep.t>=t_start+duration*1.5/2)  & (sweep.t<t_start+2*duration)
v_dep= sweep.v[flag_dep]
V_rest= sweep.v[flag_rest]
R_sweep=(v_dep.mean()-V_rest.mean())/i_amp*1000


#%% Fitting the timescale
flag_relax=(sweep.t>t_start+duration)  & (sweep.t<=t_start+duration+0.1)  # fitting by using the 100 ms.
ind_all=np.arange(len(sweep.t))
ind_vec=ind_all[flag_relax]
t_relax=sweep.t[flag_relax]-(t_start+duration)
v_relax=sweep.v[flag_relax]-sweep.v[ind_vec[-1]]

At, Kt, Ct = fit_exp_nonlinear(t_relax, v_relax, [-10, -10, 0]) 
fit_relax= model_func(t_relax,At,Kt,Ct)
plt.figure(figsize=(6,6))
plt.plot(t_relax,v_relax)
plt.plot(t_relax, fit_relax, '-',linewidth=1) 
plt.ylabel('injection relax' )

Err_fit= np.sqrt(np.mean(np.abs(fit_relax-v_relax)**2) )       

tau_sweep= -1000/Kt