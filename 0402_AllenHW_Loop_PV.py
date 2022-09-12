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
import statistics

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

#%% ONLY FOR THE FIRST TIME! Adding new features to the existing dataframe
"""
pv_metadata_nwbs = pd.read_csv("0309_pv_metadata_nwbs.csv", index_col=0)
# pv_metadata_nwbs['file_name']=pv_metadata_nwbs['file_name'].str.replace('nwb','nwb%0D')

pv_nwb_filelist = pv_metadata_nwbs[['index', 'file_name']]    #+ '%0D'
df_allPV=pv_nwb_filelist


df_allPV= df_allPV.reindex(columns = df_allPV.columns.tolist() + [
                    'HW',
                    'HWstd',
                    'Slope_rise',
                    'Slope_drop',
                    'V_reset',
                    'V_thr',
                    'V_max',
                    'R',
                    'Tau',
                           ])

df_allPV= df_allPV.astype( {
                    'HW': float,
                    'HWstd': float,
                    'Slope_rise':float,
                    'Slope_drop':float,
                    'V_reset': float,
                    'V_thr': float,
                    'V_max': float,
                    'R': float,
                    'Tau': float,
                   })
df_allPV.to_csv('0402_df_allPV_0_218.csv', index=True)


##Output_df=pv_nwb_filelist
#Output_df['rABC']=newdata.tolist()
"""

#%% Get the file names
#Output_df=pd.read_csv("Output_df.csv",index_col=0)
df_allPV=pd.read_csv("0403_df_allPV_0_217_219_498_500_551_553_602_604_619_621_622_624_833_835_864.csv")



#index_col = pv_metadata_nwbs[['index']]
nwb_filelist = df_allPV[['file_name']]
n_cells = len(nwb_filelist)

#%% Now plotting for one cell as an example. Later generate a loop to do the job.
#for i_cell in range(n_cells):
for i_cell in range(866, n_cells): # an error on cell 218, t_start = None # same error on cell 499 # same for 552 # same for 620, 623, 834, 865
    
    print(i_cell)
    cell_nwb_filename = nwb_filelist['file_name'][i_cell]
    data_set = create_ephys_data_set(nwb_file = cell_nwb_filename)
    drop_failed_sweeps(data_set)     # The function of these following lines are not making sense to me yet. 
    long_square_table = data_set.filtered_sweep_table(stimuli=data_set.ontology.long_square_names) # get sweep table for Long Square sweeps
    long_square_sweeps = data_set.sweep_set(long_square_table.sweep_number)
    long_square_sweeps.select_epoch("recording")
    long_square_sweeps.align_to_start_of_epoch("experiment")
    sweep=long_square_sweeps.sweeps[0]
    n_sweep=len(long_square_sweeps.sweeps)
    
    
    TempPath='0402_PV_HW/Cell'+str(i_cell)+'/'
    if os.path.isdir(TempPath)==False:  # To save figures
        os.makedirs(TempPath)
        
    #%% Initialize
    sweep=long_square_sweeps.sweeps[n_sweep-1]
    
    t_start, duration, i_amp, _, _ = get_stim_characteristics(sweep.i, sweep.t)
    vmax_sweep=max(sweep.v)
    Vthr_sweep=max(vmax_sweep-20,-20)
    
    df_cell =  pd.DataFrame({'Inj': pd.Series(dtype='float'),
                       'Ind_AP': pd.Series(dtype='int'),
                        'HW': pd.Series(dtype='float'),
                        'Slope_rise': pd.Series(dtype='float'),
                        'Slope_drop': pd.Series(dtype='float'),
                        'V_reset': pd.Series(dtype='float'),
                        'V_thr': pd.Series(dtype='float'),
                        'V_max': pd.Series(dtype='float'),
                       })
    
    
    R_vec=np.full(n_sweep,np.nan)
    Tau_vec=np.full(n_sweep,np.nan)
    # Initialize the fitting of the passive parameters.
    At0=-10
    Kt0=-120
    Ct0=0
    #%%  Test on the single spike
    for i_sweep in range(n_sweep-1,-1,-1):
            
        sweep=long_square_sweeps.sweeps[i_sweep]
        t_start, duration, i_amp, _, _ = get_stim_characteristics(sweep.i, sweep.t)
        vmax_sweep=max(sweep.v)
        Vthr_sweep=max(vmax_sweep-20,-20)
        
        # testing
    #    plt.figure(figsize=(12,6))
    #    plt.plot(sweep.t,sweep.v)
        
        vmax_sweep=max(sweep.v)
        Vthr_sweep=max(vmax_sweep-20,-20)
        
        # for debug
        # print(t_start) # for i_cell = 218 the 3rd sweep in the loop, its t_start is None. skip this cell by hand (20022/04/02)
        #
        
        flag_spike=(sweep.v[0:-1]<Vthr_sweep) & (sweep.v[1:]>=Vthr_sweep) & (sweep.t[:-1]>t_start) & (sweep.t[:-1]<(t_start+duration))# get the spikes
        t_spike= sweep.t[np.append(flag_spike,False)]- t_start
        
        
        # Analyzing the HW based on my own code from Matlab. Analyze_AP and Get_APs funcitons 
        xaxis_dummy=np.arange(len(sweep.t))
        spike_ind=xaxis_dummy[np.append(flag_spike,False)]
        
        dtstep=sweep.t[1]-sweep.t[0]   # calculate the timestep. 
        t_window= 0.004   # 4ms as the time window to calculate the AP. 
        Ind_window= int(round(t_window/dtstep))
        shift_t= 0.001
        Ind_shift= int(round(shift_t/dtstep))
        
        flag_slope_thr=20000  # mV/s  The threshold for detecting the AP.
        
        for i_AP in range(len(t_spike)):
    #        print(i_AP)
            TempVolt=sweep.v[spike_ind[i_AP]:spike_ind[i_AP]+Ind_window];
            TempInd_max=np.argmax(TempVolt)
            Tuned_ind=spike_ind[i_AP]+TempInd_max
            APwaveform=sweep.v[Tuned_ind-Ind_shift:Tuned_ind+Ind_window-Ind_shift];
        #    xtest=xaxis_dummy[Tuned_ind-Ind_shift:Tuned_ind+Ind_window-Ind_shift]
            # Now start writing the Analyze code for AP. Maybe write in a function later
            dAPdt=(APwaveform[1:]-APwaveform[:-1])/dtstep   # Notice the length is Ind_window-1
            Slope_rise= max(dAPdt[:Ind_shift])
            Slope_drop=min(dAPdt[Ind_shift:])
            V_max=max(APwaveform)
            V_reset=min(APwaveform[Ind_shift:])
            
            # Using interpolation to calculate the V_threshold
            flag_vthr=(dAPdt[:-1]<=flag_slope_thr)  & (dAPdt[1:]>flag_slope_thr)
            ind_vthr=[itemp for itemp, x in enumerate(flag_vthr) if x]
            if len(ind_vthr)==0:
                 print('Warning! The firing threshold is not found!!')
            else:
               a1= APwaveform[ind_vthr[0]]; a2=APwaveform[ind_vthr[0]+1];
               b1= dAPdt[ind_vthr[0]]; b2=dAPdt[ind_vthr[0]+1];
               V_thr=a1+ (flag_slope_thr-b1)*(a2-a1)/(b2-b1);
               
               V_hw=1/2*(V_thr+V_max)
               
               flag_rise=(APwaveform[:Ind_shift]<=V_hw)  & (APwaveform[1:Ind_shift+1]>V_hw)
               ind_rise=[itemp for itemp, x in enumerate(flag_rise) if x]
               dt_rise=dtstep*(V_hw-APwaveform[ind_rise[0]]) /(APwaveform[ind_rise[0]+1]-APwaveform[ind_rise[0]])    
               
               flag_drop=(APwaveform[Ind_shift:-1]>=V_hw)  & (APwaveform[Ind_shift+1:]<V_hw)
               ind_drop=[itemp+Ind_shift for itemp, x in enumerate(flag_drop) if x]   #+Ind_shift
               
               # for cell 603 index error bug
               # print(len(APwaveform))
               # print(ind_drop)
               # print('itemp, Ind_shift, flag_drop')
               # print(enumerate(flag_drop))
               # print(Ind_shift)
               # print(flag_drop)
               # for cell 603, one loop's ind_drop is empty because all x in itemp, x in enumerate(flag_drop) are False
               
               dt_drop=dtstep*(V_hw-APwaveform[ind_drop[0]]) /(APwaveform[ind_drop[0]+1]-APwaveform[ind_drop[0]])  # index error for cell 603 here
               
               if len(ind_rise)==0 | len(ind_drop)==0:
                   print('Warning! The rise voltage or drop voltage is not found!!')
               else:          
                   V_HW=(ind_drop[0]-ind_rise[0])*dtstep+ dt_drop-dt_rise
                   df_temp={'Inj':i_amp,
                       'Ind_AP': i_AP,
                        'HW': V_HW*1000,
                        'Slope_rise': Slope_rise/1000,
                        'Slope_drop': Slope_drop/1000,
                        'V_reset': V_reset,
                        'V_thr': V_thr,
                        'V_max': V_max,
                       }
                   df_cell=df_cell.append(df_temp,ignore_index= True)
        
        if i_amp<0:                
    
            
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
        
        
    # Now generate some output figure
    
    fig, axes=plt.subplots(figsize=(16, 12),nrows=4, ncols=1)  # Checking for adaptation index.
    #ind_df=np.arange(len(df_cell.index) )
             
    df_cell.plot(y="HW",use_index=True,ax=axes[0])
    #plt.ylabel('Half width (ms)')
    
    df_cell.plot(y={"Ind_AP",'Inj'},use_index=True,ax=axes[1])
    #plt.ylabel('Index of AP')
    
    #plt.subplot(413)
    df_cell.plot(y={'V_reset','V_thr','V_max'},use_index=True,ax=axes[2])
    #plt.ylabel('Reset, Thr, Max of Voltage (mV)')
    
    #plt.subplot(414)
    df_cell.plot(y={'Slope_rise','Slope_drop'},use_index=True,ax=axes[3])
    
    #plt.ylabel('R and tau')
    
    
    fig.savefig(TempPath+'HWSummary'+'.jpg',bbox_inches='tight')
     
    df_cell.to_csv(TempPath+'HWandAP.csv', index=False)  # Saving the cell.
    
    
    df_cell_1st=df_cell[df_cell['Ind_AP']==0]
    df_cell_latter=df_cell[df_cell['Ind_AP']!=0]
    
    ave_col=df_cell_latter.mean(axis=0)
    std_col=df_cell_latter.std(axis=0)
    
    df_allPV.at[i_cell,'HW']=ave_col['HW']
    df_allPV.at[i_cell,'HWstd']=std_col['HW']
    df_allPV.at[i_cell,'Slope_rise']=ave_col['Slope_rise']
    df_allPV.at[i_cell,'Slope_drop']=ave_col['Slope_drop']
    df_allPV.at[i_cell,'V_reset']=ave_col['V_reset']
    df_allPV.at[i_cell,'V_thr']=ave_col['V_thr']
    df_allPV.at[i_cell,'V_max']=ave_col['V_max']
    df_allPV.at[i_cell,'R']=np.nanmean(R_vec)
    df_allPV.at[i_cell,'Tau']=np.nanmean(Tau_vec)




# Contains the output AI, rABC rheo Slope maxInj maxRate Flag_select
df_allPV.to_csv('0403_df_allPV_0_217_219_498_500_551_553_602_604_619_621_622_624_833_835_864_866_end.csv', index=True)
