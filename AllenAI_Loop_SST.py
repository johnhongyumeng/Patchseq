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
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor

import numpy as np
import pandas as pd
import scipy as sp

import matplotlib.pyplot as plt
import seaborn as sns

import os

#%% Define functions     
def model_func(t, A, K, C):
    return A * np.exp(K * t) + C

def fit_exp_nonlinear(t, y,p0=[-5, 0.01, 25]):
    opt_parms, parm_cov = sp.optimize.curve_fit(model_func, t, y, p0, maxfev=10000)
    A, K, C = opt_parms
    return A, K, C

def ReLU_func(x,A,B):
    return (x-A)*((x-A)>0).astype(int) *B

def fit_relu(x,y,p0=[100,1]):
    opt_parms, parm_cov = sp.optimize.curve_fit(ReLU_func, x, y, p0, maxfev=10000)
    A, B = opt_parms
    return A, B

#%% ONLY FOR THE FIRST TIME!
#Output_df=pv_nwb_filelist
#Output_df['rABC']=newdata.tolist()



#%% Get the file names
#Output_df=pd.read_csv("Output_df.csv",index_col=0)
Output_df=pd.read_csv("Output_df_sst.csv")

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
#for i_cell in [5,6]:
for i_cell in range(n_cells):
    print(i_cell)
    cell_flag=True
    
    TempPath='SSTL23/Cell'+str(i_cell)+'/'
    if os.path.isdir(TempPath)==False:  # To save figures
        os.makedirs(TempPath)
    
    t_spike_max=[]
    FI_vec_max=[]
    
    cell_nwb_filename = 'SSTL23/'+nwb_filelist['file_name'][i_cell]
    data_set = create_ephys_data_set(nwb_file = cell_nwb_filename)
    drop_failed_sweeps(data_set)     # The function of these following lines are not making sense to me yet. 
    long_square_table = data_set.filtered_sweep_table(stimuli=data_set.ontology.long_square_names) # get sweep table for Long Square sweeps
    long_square_sweeps = data_set.sweep_set(long_square_table.sweep_number)
    long_square_sweeps.select_epoch("recording")
    long_square_sweeps.align_to_start_of_epoch("experiment")
    n_sweep=len(long_square_sweeps.sweeps)
    
    
    
    
    #%% Initialize
    sweep=long_square_sweeps.sweeps[n_sweep-1]
    
    
    t_start, duration, i_amp, _, _ = get_stim_characteristics(sweep.i, sweep.t)
    vmax_sweep=max(sweep.v)
    Vthr_sweep=max(vmax_sweep-20,-20)
    
    # In case the cell don't have enough spikes.
    Output_df.at[i_cell,'maxInj']=i_amp
    flag_spike=(sweep.v[0:-1]<Vthr_sweep) & (sweep.v[1:]>=Vthr_sweep) & (sweep.t[:-1]>t_start) & (sweep.t[:-1]<(t_start+duration))# get the spikes
    Output_df.at[i_cell,'maxRate']=len(flag_spike)
    
    
    
    I_vec=np.full(n_sweep,np.nan)
    nspike_vec=np.full(n_sweep,np.nan)
    rABC_vec=np.full(n_sweep,np.nan)
    AI_vec= np.full(n_sweep,np.nan)
    Err_FI_vec= np.full(n_sweep,np.nan)
    
    R_vec=np.full(n_sweep,np.nan)
    Tau_vec=np.full(n_sweep,np.nan)
    hw_all=np.empty(0)
    
    sweep2=sweep;
    sweep0=long_square_sweeps.sweeps[0]
    sweep1=sweep; 
    
    # Initialize the fitting of the passive parameters.
    At0=-10
    Kt0=-100
    Ct0=0
    #%%  Test on the single spike
    for i_sweep in range(n_sweep-1,-1,-1):
        #print(i_sweep)
        sweep=long_square_sweeps.sweeps[i_sweep]
        _, _, i_amp, _, _ = get_stim_characteristics(sweep.i, sweep.t)
    
        flag_spike=(sweep.v[0:-1]<Vthr_sweep) & (sweep.v[1:]>=Vthr_sweep) & (sweep.t[:-1]>t_start) & (sweep.t[:-1]<(t_start+duration))# get the spikes
        t_spike= sweep.t[np.append(flag_spike,False)]- t_start
        I_vec[i_sweep]=i_amp
        nspike_vec[i_sweep]=len(t_spike)
        
        if (i_sweep==n_sweep-1) & (len(t_spike)<10) :
            print('Not enough spike for the highest injection') 
            cell_flag=False
    
            break
        elif len(t_spike)>5: # fitting for at least 5 spikes.
    
            t_isi=np.diff(t_spike)
            FI_vec=1/t_isi       
        # for testing the figure.
        # plt.figure(figsize=(6,6))
        # plt.plot(t_spike[:-1],FI_vec)
    
    
        # Sanity test on fitting
    
            A, K, C = fit_exp_nonlinear(t_spike[:-1], FI_vec, [2*FI_vec[-1], -0.1, FI_vec[-1]]) 
            t_axis=np.linspace(0.0,1.0,num=10001,endpoint=True)
            fit_y = model_func(t_axis, A, K, C)
            #plt.figure(figsize=(6,6))
            #plt.plot(t_spike[:-1],FI_vec)
            #plt.plot(t_axis, fit_y, '-',linewidth=1) 
        
            # get the readout, 
            AI_sweep=1-fit_y[-1]/fit_y[0]
            ABC_vec=-fit_y+ (fit_y[-1]+  (fit_y[0]-fit_y[-1])*(1-t_axis))
            rABC=sum(ABC_vec)/len(ABC_vec)/(fit_y[0]-fit_y[-1])*2
            AI_vec[i_sweep]=AI_sweep
            rABC_vec[i_sweep]=rABC
            
            #including the error of the fitting now
            fit_FI=model_func(t_spike[:-1], A, K, C)
            Err_FI=np.sqrt(np.mean((np.abs(fit_FI-FI_vec)/fit_FI)**2) )
            Err_FI_vec[i_sweep]=Err_FI
            
            sweep1=sweep
            t_spike_rheo=t_spike[:-1]
            FI_vec_rheo=FI_vec
            
            if (i_sweep==n_sweep-1):
                t_spike_max=t_spike[:-1]
                FI_vec_max=FI_vec
                fig=plt.figure(figsize=(6,6))
                plt.plot(t_spike[:-1],FI_vec)
                plt.plot(t_axis, fit_y, '-',linewidth=1) 
                fig.savefig(TempPath+'SampleFit'+'.jpg',bbox_inches='tight')       
        elif i_amp<0:                
            flag_dep=(sweep.t>=t_start+duration/2)  & (sweep.t<t_start+duration)
            flag_rest=(sweep.t>=t_start+duration*1.5/2)  & (sweep.t<t_start+2*duration)
            v_dep= sweep.v[flag_dep]
            V_rest= sweep.v[flag_rest]
            R_sweep=(v_dep.mean()-V_rest.mean())/i_amp*1000
            R_vec[i_sweep]=R_sweep
            
            # Fitting the timescale
            flag_relax=(sweep.t>t_start+duration)  & (sweep.t<=t_start+duration+0.1)  # fitting by using the 100 ms.
            ind_all=np.arange(len(sweep.t))
            ind_vec=ind_all[flag_relax]
            t_relax=sweep.t[flag_relax]-(t_start+duration)
            v_relax=sweep.v[flag_relax]-sweep.v[ind_vec[-1]]
            
            At, Kt, Ct = fit_exp_nonlinear(t_relax, v_relax, [At0, Kt0, Ct0]) 
            
            
            fit_relax= model_func(t_relax,At,Kt,Ct)
            Err_fit= np.sqrt(np.mean(np.abs(fit_relax-v_relax)**2) )       
            if Err_fit< 0.15:
                At0=At
                Kt0=Kt
                Ct0=Ct
                tau_sweep= -1000/Kt
                Tau_vec[i_sweep]=tau_sweep    
            else:
                print('Warning! Passive fitting error= %.02f too large, skip.' % (Err_fit))    
    #        plt.figure(figsize=(6,6))
    #        plt.plot(t_relax,v_relax)
    #        plt.plot(t_relax, fit_relax, '-',linewidth=1) 
            
            
        if len(t_spike)>0 :
            sweep_ext = EphysSweepFeatureExtractor(t=sweep.t, v=sweep.v, i=sweep.i, start=t_start+0.0, end=t_start+duration+0.0)
            sweep_ext.process_spikes()
            sweep_hw=sweep_ext.spike_feature("width")
            hw_all=np.append(hw_all,sweep_hw)
                 
    #%% generate the figure output for each cell
    fig=plt.figure(figsize=(16, 12))  # Checking for adaptation index.         
    plt.subplot(411)
    plt.plot(I_vec, nspike_vec,'.',markersize=3)
    if cell_flag:
        Rheo, k=  fit_relu(I_vec,nspike_vec,[100,1])   
        fit_rate=ReLU_func(I_vec, Rheo,k)
        plt.plot(I_vec,fit_rate,'-',linewidth=1)
    
    
    plt.xlabel('Injection (nA)')
    plt.ylabel('Firing rate (sp/s)')
    
    plt.subplot(412)
    plt.plot(I_vec, AI_vec)
    plt.ylabel('Adaptation Index')
    
    plt.subplot(413)
    plt.plot(I_vec, rABC_vec)
    plt.ylabel('ratio of Area_Between_curve')
    
    plt.subplot(414)
    plt.plot(I_vec, R_vec)
    plt.plot(I_vec, Tau_vec,'.-',markersize=5)
    
    plt.ylabel('R and tau')

    
    fig.savefig(TempPath+'Summary'+'.jpg',bbox_inches='tight')
    
    #%% Relu fitting to get an understanding of the cell
    
    fig, ax=plt.subplots(figsize=(12,6))
    plt.plot(sweep0.t-t_start,sweep0.v)
    plt.plot(sweep1.t-t_start,sweep1.v)
    plt.plot(sweep2.t-t_start,sweep2.v)
    ax.set_xlim([-0.1,1.1])
    fig.savefig(TempPath+'Sweeps'+'.jpg',bbox_inches='tight')
    
    if cell_flag:
        fig, ax=plt.subplots(figsize=(6,6))
        plt.plot(t_spike_max,FI_vec_max,'.',markersize=5)
        plt.plot(t_spike_rheo,FI_vec_rheo,'.',markersize=5)
        ax.set_xlim([-0.1,1.1])
        fig.savefig(TempPath+'FIs'+'.jpg',bbox_inches='tight')
    
    
    #%% Output the data into Output_df. Make sure do not overwrite what's saved before.
    if cell_flag:
        indout=np.argmax(nspike_vec)
        Output_df.at[i_cell,'AI']=AI_vec[indout]
        Output_df.at[i_cell,'rABC']=rABC_vec[indout]
        Output_df.at[i_cell,'rheo']=Rheo
        Output_df.at[i_cell,'Slope']=k
        Output_df.at[i_cell,'maxInj']=I_vec[indout]
        Output_df.at[i_cell,'maxRate']=nspike_vec[indout]
        Output_df.at[i_cell,'hw']=hw_all.mean()*1e3
        Output_df.at[i_cell,'R']=np.nanmean(R_vec)
        Output_df.at[i_cell,'Tau']=np.nanmean(Tau_vec)
        Output_df.at[i_cell,'AI_err']=Err_FI_vec[indout]
        if (I_vec[indout]>Rheo+100) |(I_vec[indout]>Rheo*1.5) | (nspike_vec[indout]>40):
            Output_df.at[i_cell,'Flag_select']=1
        else:
            Output_df.at[i_cell,'Flag_select']=0    
    else:
        Output_df.at[i_cell,'Flag_select']=0    

# Contains the output AI, rABC rheo Slope maxInj maxRate Flag_select
Output_df.to_csv('Output_df_sst.csv', index=False)

