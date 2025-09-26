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
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib
matplotlib.use('agg')  # shut down the inline figure genreatation, but not in the "plot" pane. Adding a "matplotlibrc" at %USERPROFILE%/.matplotlib to shutdown all the figures.
#matplotlib.use('Qt5Agg') # generate the inline figures
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.io as io
import scipy.signal as signal
from scipy.optimize import curve_fit
from scipy import integrate
import logging
import importlib
import os
import warnings

from pynwb import NWBHDF5IO
from allensdk.ephys import ephys_extractor as efex
from allensdk.ephys import ephys_features as ft

from ipfx.dataset.create import create_ephys_data_set
from ipfx.utilities import drop_failed_sweeps
from ipfx.stim_features import get_stim_characteristics


sns.set_style()

import matplotlib as mpl


#from sklearn import linear_model
#ransac = linear_model.RANSACRegressor()

import scipy as sp
import re


#os.chdir('./ephy_temperature/0308_ephy_temp_nwbs')

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

def get_time_voltage_current_currindex0(nwb):
    df = nwb.sweep_table.to_dataframe()
    voltage = np.zeros((len(df['series'][0][0].data[:]), int((df.shape[0]+1)/2)))
    time = np.arange(len(df['series'][0][0].data[:]))/df['series'][0][0].rate
    voltage[:, 0] = df['series'][0][0].data[:]
    current_initial = df['series'][1][0].data[12000]*df['series'][1][0].conversion
    curr_index_0 = int(-current_initial/20) # index of zero current stimulation
    current = np.linspace(current_initial, (int((df.shape[0]+1)/2)-1)*20+current_initial, \
                         int((df.shape[0]+1)/2))
    for i in range(curr_index_0):   # Find all voltage traces from minimum to 0 current stimulation
        voltage[:, i+1] = df['series'][0::2][(i+1)*2][0].data[:]
    for i in range(curr_index_0, int((df.shape[0]+1)/2)-1):   # Find all voltage traces from 0 to highest current stimulation
        voltage[:, i+1] = df['series'][1::2][i*2+1][0].data[:]
    voltage[:, curr_index_0] = df.loc[curr_index_0*2][0][0].data[:]    # Find voltage trace for 0 current stimulation
    return time, voltage, current, curr_index_0


#%% ONLY FOR THE FIRST TIME! Adding new features to the existing dataframe
"""
ephy_temp_meta = pd.read_csv("0411Names.csv", names=['file_name'])
ephy_temp_meta['Cell'] = [s.replace("-", "_") for s in [
    "".join(item) for item in [
        re.findall('_cell-(.*)_icephys.nwb', i_str) for i_str in ephy_temp_meta['file_name']
    ]
]]


pv_metadata_nwbs = ephy_temp_meta

pv_nwb_filelist = pv_metadata_nwbs[['Cell', 'file_name']]    #+ '%0D'
df_allPV=pv_nwb_filelist
nwb_filelist = df_allPV[['file_name']]
n_cells = len(nwb_filelist)
"""
# ## Collect filenames and directories
# root_dir="C:\AllenOtherCells"
# cell_data = []
# for root, dirs, files in os.walk(root_dir):
#     for file in files:
#         cell_data.append({
#             "directory": root,
#             "filename": file
#         })

# # Create DataFrame
# df = pd.DataFrame(cell_data)

# # Show first few rows
# print(df.head())

# # (Optional) Save to CSV
# df.to_csv("cells.csv", index=False)



PV_names=pd.read_csv('./cells.csv')
df_allPV=PV_names[['filename']]

nwb_filelist = PV_names[['filename']]
n_cells = len(nwb_filelist)

Path_lists=PV_names[['directory']]

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
#df_allPV.to_csv('0309_room_temp_fitting.csv', index=True)


##Output_df=pv_nwb_filelist
#Output_df['rABC']=newdata.tolist()
#"""

#%% Get the file names
#Output_df=pd.read_csv("Output_df.csv",index_col=0)
#df_allPV=pd.read_csv("0309_room_temp_fitting.csv")



#index_col = pv_metadata_nwbs[['index']]
#i_cell_vec=np.array([1])

#%% Now plotting for one cell as an example. Later generate a loop to do the job.
#for i_cell in range(115,n_cells): # i_cell = 0,2, 4, 33, 48, 50, 99, 104, 136, 143, 159, 160, 166, 189, 190, 196, 206, 207, 243, 245, 251, 256, 258, 260, 261, 262, 270, 273 have bug "list index out of range" for line around 254: dt_drop=dtstep*(V_hw-APwaveform[ind_drop[0]]) /(APwaveform[ind_drop[0]+1]-APwaveform[ind_drop[0]])  
# 0: low quqlity; 2:L2/3 IT_3;50: IT; 99: Sst; 104: low quality                 284, 287, 314, 321, 324, 328, 340, 379
#for i_cell in i_cell_vec:
# for i_cell in range(0,n_cells):
for i_cell in range(n_cells):    
    print(i_cell)
    cell_flag=True
#    print(df_allPV['Cell'][i_cell])
    fpath = Path_lists['directory'][i_cell]+'\\'+nwb_filelist['filename'][i_cell]
    
    '''
    io_ = NWBHDF5IO(fpath, 'r', load_namespaces=True)
    nwb = io_.read()
    time, voltage, current, curr_index_0 = get_time_voltage_current_currindex0(nwb)
    
    n_sweep = voltage.shape[1]
    t_start = 0.1 # defaulted by Scala et. al.
    t_end = 0.7 # from Scala default
    duration = t_end - t_start
    '''
    try:

        data_set = create_ephys_data_set(nwb_file = fpath)
        drop_failed_sweeps(data_set)     # The function of these following lines are not making sense to me yet. 
        long_square_table = data_set.filtered_sweep_table(stimuli=data_set.ontology.long_square_names) # get sweep table for Long Square sweeps
        long_square_sweeps = data_set.sweep_set(long_square_table.sweep_number)
        long_square_sweeps.select_epoch("recording")
        long_square_sweeps.align_to_start_of_epoch("experiment")
        sweep=long_square_sweeps.sweeps[0]
        n_sweep=len(long_square_sweeps.sweeps)    
        sweep=long_square_sweeps.sweeps[n_sweep-1]
        t_start, duration, i_amp, _, _ = get_stim_characteristics(sweep.i, sweep.t)
    
        TempPath='All_figures/Cell'+str(i_cell)+'/'
        if os.path.isdir(TempPath)==False:  # To save figures
            os.makedirs(TempPath)
            
        #%% Initialize
    
    
        df_cell =  pd.DataFrame({'Inj': pd.Series(dtype='float'),
                           'Ind_AP': pd.Series(dtype='int'),
                            'HW': pd.Series(dtype='float'),
                            'Slope_rise': pd.Series(dtype='float'),
                            'Slope_drop': pd.Series(dtype='float'),
                            'V_reset': pd.Series(dtype='float'),
                            'V_thr': pd.Series(dtype='float'),
                            'V_max': pd.Series(dtype='float'),
                            't_spike': pd.Series(dtype='float'),
                           })
        
        df_cell_output =  pd.DataFrame({'Inj': pd.Series(dtype='float'),
                           'Ind_AP': pd.Series(dtype='int'),
                            'HW': pd.Series(dtype='float'),
                            'Slope_rise': pd.Series(dtype='float'),
                            'Slope_drop': pd.Series(dtype='float'),
                            'V_reset': pd.Series(dtype='float'),
                            'V_thr': pd.Series(dtype='float'),
                            'V_max': pd.Series(dtype='float'),
                            't_spike': pd.Series(dtype='float'),
                           })    
        R_vec=np.full(n_sweep,np.nan)
        Tau_vec=np.full(n_sweep,np.nan)
        Inj_vec=np.full(n_sweep,np.nan)
        nspike_vec=np.full(n_sweep,np.nan)
        rABC_vec=np.full(n_sweep,np.nan)
        AI_vec= np.full(n_sweep,np.nan)
        Err_FI_vec= np.full(n_sweep,np.nan)
        
        # Things added to calculate the adaptation index
        sweep0=long_square_sweeps.sweeps[0]  # Hyperpolarized step. I never expect this can be different.
        sweep1_flag=False                     # Rheobase step  
        sweep1_flag2=False
        sweep2_flag=False                     # Max firing rate step`
        
        
        # Initialize the fitting of the passive parameters.
        At0=-10
        Kt0=-120
        Ct0=0
        indout=n_sweep-1
        #%%  Test on the single spike
        maxSpikenum=0;
        for i_sweep in range(n_sweep):
    #        print(i_sweep)
    #    for i_sweep in range(28,29):
    #    for i_sweep in range(n_sweep-1,-1,-1):
            '''
            i_amp = current[i_sweep]
            sweep_v = voltage[:, i_sweep]
            sweep_t = time
            '''
            sweep=long_square_sweeps.sweeps[i_sweep]
            t_start, duration, i_amp, _, _ = get_stim_characteristics(sweep.i, sweep.t)
            sweep_v = sweep.v
            sweep_t = sweep.t
    #        print(sweep_t.shape)  Can not believe the resolution are different between pre-rheobase and post-rheobase
            vmax_sweep=max(sweep_v)
            vup_sweep=min(vmax_sweep,20)
            vlow_sweep=max(-40,np.median(sweep_v))
            Vthr_sweep=max(vup_sweep*0.8+vlow_sweep*0.2,0)
            # V_step=sweep.v()
            df_sweep =  pd.DataFrame({'Inj': pd.Series(dtype='float'),
                       'Ind_AP': pd.Series(dtype='int'),
                        'HW': pd.Series(dtype='float'),
                        'Slope_rise': pd.Series(dtype='float'),
                        'Slope_drop': pd.Series(dtype='float'),
                        'V_reset': pd.Series(dtype='float'),
                        'V_thr': pd.Series(dtype='float'),
                        'V_max': pd.Series(dtype='float'),
                        't_spike': pd.Series(dtype='float'),
                       })        
        #   sweep=long_square_sweeps.sweeps[i_sweep]
        #   t_start, duration, i_amp, _, _ = get_stim_characteristics(sweep.i, sweep.t)
        #    vmax_sweep=max(sweep.v)
        #    Vthr_sweep=max(vmax_sweep-20,-20)
            
            # testing
        #    plt.figure(figsize=(12,6))
        #    plt.plot(sweep_t,sweep_v)
        #    plt.title('Ind sweep'+str(i_sweep))
        #    vmax_sweep=max(sweep.v)
        #    Vthr_sweep=max(vmax_sweep-20,-20)
            
            # Previous check for spike. Now adding the voltage check check.
            flag_spike=(sweep_v[0:-1]<Vthr_sweep) & (sweep_v[1:]>=Vthr_sweep) & (sweep_t[:-1]>t_start) & (sweep_t[:-1]<(t_start+duration))# get the spikes
            t_spike= sweep_t[np.append(flag_spike,False)]- t_start
            t_spike_flag= np.full(len(t_spike),True)
    
    
            if len(t_spike)==0 and np.abs(i_amp)<200 and np.abs(i_amp)>15:                
        
                
                flag_dep=(sweep_t>=t_start+duration/2)  & (sweep_t<t_start+duration)
                flag_rest=(sweep_t>=t_start+duration*1.5/2)  & (sweep_t<t_start+2*duration)
                v_dep= sweep_v[flag_dep]
                V_rest= sweep_v[flag_rest]
                R_sweep=(v_dep.mean()-V_rest.mean())/i_amp*1000
                R_vec[i_sweep]=R_sweep
                
                # Fitting the timescale
                flag_relax=(sweep_t>t_start+duration)  & (sweep_t<=t_start+duration+0.1)  # fitting by using the 100 ms.
                ind_all=np.arange(len(sweep_t))
                ind_vec=ind_all[flag_relax]
                t_relax=sweep_t[flag_relax]-(t_start+duration)
                v_relax=sweep_v[flag_relax]-sweep_v[ind_vec[-1]]
                
                try:
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
                except:
                    print('Warning! Passive fitting fails. Skip')
                        
    
    
            # Analyzing the HW based on my own code from Matlab. Analyze_AP and Get_APs funcitons 
            xaxis_dummy=np.arange(len(sweep_t))
            spike_ind=xaxis_dummy[np.append(flag_spike,False)]
            
            dtstep=sweep_t[1]-sweep_t[0]   # calculate the timestep. 
            t_window= 0.005  # 4ms as the time window to calculate the AP. 
            Ind_window= int(round(t_window/dtstep))
            shift_t= 0.002
            Ind_shift= int(round(shift_t/dtstep))
    
    
            flag_slope_thr=20000  # mV/s  The threshold for detecting the AP.
            
            for i_AP in range(0,len(t_spike)):
        #        print(i_AP)
                if i_AP< len(t_spike)-1:
                    Ind_extend= min(Ind_window+spike_ind[i_AP],spike_ind[i_AP+1])
                else:
                    Ind_extend=Ind_window+spike_ind[i_AP]
    
    
                TempVolt=sweep_v[spike_ind[i_AP]:Ind_extend]
                TempInd_max=np.argmax(TempVolt)
                Tuned_ind=spike_ind[i_AP]+TempInd_max
                APwaveform=sweep_v[Tuned_ind-Ind_shift:Tuned_ind+Ind_window-Ind_shift];
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
                    print('Warning! This AP is a false-positive')
                    t_spike_flag[i_AP]=False
                    # For debugging only. Shut down later. 
                    #fig, axes=plt.subplots(3,1,figsize=(15,8))
                    #axes[0].plot(APwaveform)
                    #axes[0].set_title('TimeAt'+str(sweep_t[spike_ind[i_AP]]))
                    #axes[1].plot(dAPdt)
                    #axes[1].set_title('TimeAt'+str(sweep_t[spike_ind[i_AP]]))                
                    #axes[2].plot(sweep_t,sweep_v)
                    #axes[2].set_title('Vthr_sweep'+str(Vthr_sweep))
                    #plt.show()        
                    
                #    input('Paused for debugging. Go check the figure. Type anything to continue..') 
                     
                else:
                   a1= APwaveform[ind_vthr[0]]; a2=APwaveform[ind_vthr[0]+1];
                   b1= dAPdt[ind_vthr[0]]; b2=dAPdt[ind_vthr[0]+1];
                   V_thr=a1+ (flag_slope_thr-b1)*(a2-a1)/(b2-b1);
                   
                   V_hw=1/2*(V_thr+V_max)
                   
                   flag_rise=(APwaveform[:Ind_shift]<=V_hw)  & (APwaveform[1:Ind_shift+1]>V_hw)
                   ind_rise=[itemp for itemp, x in enumerate(flag_rise) if x]
                   
                   try:
                       dt_rise=dtstep*(V_hw-APwaveform[ind_rise[0]]) /(APwaveform[ind_rise[0]+1]-APwaveform[ind_rise[0]])
                       flag_drop=(APwaveform[Ind_shift:-1]>=V_hw)  & (APwaveform[Ind_shift+1:]<V_hw)
                       ind_drop=[itemp+Ind_shift for itemp, x in enumerate(flag_drop) if x]   #+Ind_shift
                       dt_drop=dtstep*(V_hw-APwaveform[ind_drop[0]]) /(APwaveform[ind_drop[0]+1]-APwaveform[ind_drop[0]])
                       
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
                                't_spike': sweep_t[Tuned_ind] -t_start,
                               }
                           df_sweep=df_sweep.append(df_temp,ignore_index= True)
                   except:
                       print("Error! Plot the sweep!!!")
                       fig, axes=plt.subplots(figsize=(12,6))
                       axes.plot(sweep_t,sweep_v) 
                       fig.savefig(TempPath+'ErrorSweepInj='+str(i_amp)+'Vthr_sweep'+str(Vthr_sweep)+'.jpg',bbox_inches='tight')
                       plt.close()
    
            t_spike_proofed=t_spike[t_spike_flag]
            Inj_vec[i_sweep]=i_amp
            nspike_vec[i_sweep]=len(t_spike_proofed)
            if len(t_spike_proofed)<=50:   # 50(AP/s)*0.6=30 
                df_cell=df_cell.append(df_sweep,ignore_index= True)
            df_cell_output=df_cell_output.append(df_sweep,ignore_index= True)
            # Get the rheobase sweep
            if len(t_spike_proofed)>0 and sweep1_flag==False:
                sweep1=sweep
                sweep1_flag=True
            # Get the sweep with the maximum number of APs. Stop the scanning.
    
            if len(t_spike_proofed)>maxSpikenum:
                maxSpikenum=len(t_spike)
                last_spike_t=t_spike[-1]
                #sweep2=voltage[:, i_sweep]
                sweep2=sweep
                
                indout=i_sweep
            '''        
            elif len(t_spike_proofed)<(maxSpikenum-2):
                sweep2_flag=True
                break       
            '''
            
            # Now introducing the AI calculation
            
            if len(t_spike_proofed)>5: # fitting for at least 5 spikes.
        
                t_isi=np.diff(t_spike_proofed)
                FI_vec=1/t_isi
                # for debugging
                #plt.figure(figsize=(6,6))
                #plt.subplot(211)
                #plt.plot(sweep_t,sweep_v)
                #plt.subplot(212)
                #plt.plot(t_spike_proofed[:-1],FI_vec)
                try:
                  
                    K0_vec=np.array([-10,-10,-100,-100,-1,-0.1,-0.1])
                    AI0_vec=np.array([0.5,1,0.1,-0.1,0.1,-0.1,0.1])
                    for K0,AIamp in zip(K0_vec,AI0_vec):     # This version is faster. For a more detailed version, see *_getrheo
                        try:
                            A, K, C = fit_exp_nonlinear(t_spike_proofed[:-1], FI_vec, [AIamp*FI_vec[-1], K0, FI_vec[-1]])
                            fit_FI=model_func(t_spike_proofed[:-1], A, K, C)
                            Err_FI=np.sqrt(np.mean((np.abs(fit_FI-FI_vec)/fit_FI)**2) ) 
                            if Err_FI<0.5:
                                break
                        except:
                            continue
                        
                    t_axis=np.linspace(0.0,1.0,num=10001,endpoint=True)
                    fit_y = model_func(t_axis, A, K, C)
                    #plt.figure(figsize=(6,6))
                    #plt.plot(t_spike[:-1],FI_vec)
                    #plt.plot(t_axis, fit_y, '-',linewidth=1) 
                
                    # get the readout, 
                    AI_sweep=1-fit_y[-1]/fit_y[0]
    #                ABC_vec=-fit_y+ (fit_y[-1]+  (fit_y[0]-fit_y[-1])*(0.6-t_axis)/0.6)
                    ABC_vec=-fit_y+ (fit_y[-1]+  (fit_y[0]-fit_y[-1])*(1-t_axis)/1.0)
    
                    rABC=sum(ABC_vec)/len(ABC_vec)/(fit_y[0]-fit_y[-1])*2
                    AI_vec[i_sweep]=AI_sweep
                    rABC_vec[i_sweep]=rABC
                    
                    #including the error of the fitting now
    
                    Err_FI_vec[i_sweep]=Err_FI
                except:
                    print("Error! The adaptation curve doesn't fit correctly!!!")
                    AI_vec[i_sweep]=np.nan
                    fig=plt.figure(figsize=(6,6))
                    plt.subplot(211)
                    plt.plot(sweep_t,sweep_v)
                    plt.subplot(212)
                    plt.plot(t_spike_proofed[:-1],FI_vec)
                    fig.savefig(TempPath+'FitErrorSweep'+str(i_sweep)+'.jpg',bbox_inches='tight')
                    plt.close()
    
                if sweep1_flag2==False:
                    t_spike_rheo=t_spike_proofed[:-1]
                    FI_vec_rheo=FI_vec
                    sweep1_flag2=True
                    
                t_spike_max=t_spike_proofed[:-1]
                FI_vec_max=FI_vec
    #            fig=plt.figure(figsize=(6,6))
    #                plt.plot(t_spike[:-1],FI_vec)
    #                plt.plot(t_axis, fit_y, '-',linewidth=1) 
    #                fig.savefig(TempPath+'SampleFit'+'.jpg',bbox_inches='tight')       
    
        
        if maxSpikenum<10 :
                print('Not enough spike for the highest injection') 
                cell_flag=False
            
                
           
    
        # Now generate some output figure
        if df_cell.shape[0]>0:
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
            plt.close()
    
        if df_cell_output.shape[0]>0:
        
            df_cell_output.to_csv(TempPath+'HWandAP.csv', index=False)  # Saving the cell.
        
        
        df_cell_1st=df_cell[df_cell['Ind_AP']==0]
        df_cell_latter=df_cell[df_cell['Ind_AP']!=0]
        
        ave_col=df_cell_latter.mean(axis=0)
        std_col=df_cell_latter.std(axis=0)
    
        #%% generate the figure output for each cell
        fig=plt.figure(figsize=(16, 12))  # Checking for adaptation index.         
        plt.subplot(411)
        plt.plot(Inj_vec[0:indout], nspike_vec[0:indout],'.',markersize=3)
        if cell_flag:
            Rheo, k=  fit_relu(Inj_vec[0:indout],nspike_vec[0:indout],[100,1])   
            fit_rate=ReLU_func(Inj_vec[0:indout], Rheo,k)
            plt.plot(Inj_vec[0:indout],fit_rate,'-',linewidth=1)
        
        
        plt.xlabel('Injection (nA)')
        plt.ylabel('Firing rate (sp/s)')
        
        plt.subplot(412)
        plt.plot(Inj_vec, AI_vec)
        plt.ylabel('Adaptation Index')
        
        plt.subplot(413)
        plt.plot(Inj_vec, rABC_vec)
        plt.ylabel('ratio of Area_Between_curve')
        
        plt.subplot(414)
        plt.plot(Inj_vec, R_vec)
        plt.plot(Inj_vec, Tau_vec,'.-',markersize=5)
        
        plt.ylabel('R and tau')
        fig.savefig(TempPath+'Summary'+'.jpg',bbox_inches='tight')
        plt.close()
    
        #%% Relu fitting to get an understanding of the cell
        
        fig, ax=plt.subplots(figsize=(12,6))
        plt.plot(sweep0.t-t_start,sweep0.v)
        plt.plot(sweep1.t-t_start,sweep1.v)
        plt.plot(sweep2.t-t_start,sweep2.v)
        ax.set_xlim([-0.1,1.1])
        fig.savefig(TempPath+'Sweeps'+'.jpg',bbox_inches='tight')
        plt.close()
        if cell_flag:
            fig, ax=plt.subplots(figsize=(6,6))
            plt.plot(t_spike_max,FI_vec_max,'.',markersize=5)
            plt.plot(t_spike_rheo,FI_vec_rheo,'.',markersize=5)
            ax.set_xlim([-0.1,0.65])
            fig.savefig(TempPath+'FIs'+'.jpg',bbox_inches='tight')
            plt.close()
    
    
    
        
        df_allPV.at[i_cell,'HW']=ave_col['HW']
        df_allPV.at[i_cell,'HWstd']=std_col['HW']
        df_allPV.at[i_cell,'Slope_rise']=ave_col['Slope_rise']
        df_allPV.at[i_cell,'Slope_drop']=ave_col['Slope_drop']
        df_allPV.at[i_cell,'V_reset']=ave_col['V_reset']
        df_allPV.at[i_cell,'V_thr']=ave_col['V_thr']
        df_allPV.at[i_cell,'V_max']=ave_col['V_max']
        df_allPV.at[i_cell,'LastSpike']=last_spike_t
        df_allPV.at[i_cell,'MaxSpikeNum']=maxSpikenum
        df_allPV.at[i_cell,'R']=np.nanmean(R_vec)
        df_allPV.at[i_cell,'Tau']=np.nanmean(Tau_vec)
    
        if cell_flag:
            df_allPV.at[i_cell,'AI']=AI_vec[indout]
            df_allPV.at[i_cell,'rABC']=rABC_vec[indout]
            df_allPV.at[i_cell,'rheo']=Rheo
            df_allPV.at[i_cell,'Slope']=k
            df_allPV.at[i_cell,'maxInj']=Inj_vec[indout]
            df_allPV.at[i_cell,'maxRate']=nspike_vec[indout]
            df_allPV.at[i_cell,'AI_err']=Err_FI_vec[indout]
            df_allPV.at[i_cell,'Flag_select']=1
        else:
            df_allPV.at[i_cell,'Flag_select']=0    
    except:
        print("Error! The loading process failed!!!")
        
from datetime import datetime
date_str = datetime.now().strftime('%m%d%y')
filename = f'{date_str}_AllenV1PV_fitting_runall.csv'
df_allPV.to_csv(filename, index=True)

# Contains the output AI, rABC rheo Slope maxInj maxRate Flag_select
#df_allPV.to_csv('041823_AllenV1PV_fitting_runall.csv', index=True)
#df_allPV.to_csv('041823_AllenV1_fitting_runto1664.csv', index=True)
