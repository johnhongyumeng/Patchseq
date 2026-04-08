# -*- coding: utf-8 -*-
"""
Created on Wed May 28 09:29:02 2025
The script is designed to test the individual cell. Meant to merge into the 
major loop later
@author: yawnt
"""

# use py36, not py3
import pandas as pd
import numpy as np
import seaborn as sns
#import matplotlib
#matplotlib.use('agg')  # shut down the inline figure genreatation, but not in the "plot" pane. Adding a "matplotlibrc" at %USERPROFILE%/.matplotlib to shutdown all the figures.
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.io as io
import scipy.signal as signal
from scipy.optimize import curve_fit
from scipy import integrate
import scipy.stats as stats

import logging
import importlib
import os
import warnings
from lmfit import Parameters, Minimizer, report_fit
from lmfit.models import ExponentialModel, LinearModel, ConstantModel



import logging
# Set up the logger to log to a file
logging.basicConfig(filename='warningM1all.log', level=logging.WARNING)


#from pynwb import NWBHDF5IO
#from allensdk.ephys import ephys_extractor as efex
#from allensdk.ephys import ephys_features as ft

sns.set_style()

import matplotlib as mpl


#from sklearn import linear_model
#ransac = linear_model.RANSACRegressor()

import scipy as sp
import re
from tqdm import tqdm
import os

def ReLU_func(x,A,B):
    return (x-A)*((x-A)>0).astype(int) *B

def fit_relu(x,y,p0=[100,1]):
    opt_parms, parm_cov = sp.optimize.curve_fit(ReLU_func, x, y, p0, maxfev=10000)
    A, B = opt_parms
    return A, B

def eval_fit(model_func, params, t):
    vals = [params[name].value for name in model_func.__code__.co_varnames[1:model_func.__code__.co_argcount]]
    return model_func(t, *vals)


def double_exp_func(t,A1,A2,k1,k2,C):
    # A1>0 ,A2 can be both, K1, K2<0 C>0 # k1 is faster. 
    # suggest value for k1 [-20, -100] equals timescale [10 50]ms
    #  k2 [-1, -20] equals timescale [50 to 1000]ms
    return A1*np.exp(-k1*t)+A2*np.exp(-k2*t)+C

def residual_de(params, t, data):  # de for double_exp_func
    A1 = params['A1']
    A2 = params['A2']
    k1 = params['k1']
    k2 = params['k2']
    C  = params['C']
    model = double_exp_func(t, A1, A2, k1, k2, C)
    return model - data

def single_exp_func(t, A2,  k2, C):
    return A2*np.exp(-k2*t)+C

def residual_se(params, t, data):  # se for single_exp_func
    #A1 = params['A1']
    A2 = params['A2']
    #k1 = params['k1']
    k2 = params['k2']
    C  = params['C']
    model = single_exp_func(t, A2, k2, C)
    return model - data

def classic_adap_func(t, A2, k2, C):
    return A2*np.exp(-k2*t)+C

def residual_ca(params, t, data):  # ca for classic_adap_func
    #A1 = params['A1']
    A2 = params['A2']
    #m  = params['m']
    #k1 = params['k1']
    k2 = params['k2']
    C  = params['C']
    model = classic_adap_func(t, A2, k2, C)
    return model - data

def multi_fitting(t, y, initials=None):
    # Rewrite all the fitting functions. Fit to three.
    # double_exp_func, single_exp_func,

    if initials is None:
        initials = np.array([500, 50,  100, 10, 2])  
        # Fast_amplitude, slow_amplitude, const, fast_decay, slow_decay, constant

    params_de= Parameters()
    params_de.add('A1', value=initials[0], min=0, max=2000)
    params_de.add('A2', value=initials[1], min=0, max=100)
    params_de.add('k1', value=initials[2], min=50, max=500)   # fast decay
    params_de.add('k2', value=initials[3], min=1, max=50)     # slow decay
    params_de.add('C', value=initials[4], min=0)

    params_se= Parameters()
    params_se.add('A2', value=initials[1], min=-100, max=500)
    params_se.add('k2', value=initials[3], min=1, max=50)     # slow decay
    params_se.add('C', value=initials[4], min=0)




    minimizer_de = Minimizer(residual_de, params_de, fcn_args=(t, y))
    result_de  = minimizer_de.minimize()

    minimizer_se = Minimizer(residual_se, params_se, fcn_args=(t, y))
    result_se  = minimizer_se.minimize()

    n_data = len(t)
    n_params = len(result_de.params)
    df = n_data - n_params
    
    # Extract fitted value and standard error
    A1_val = result_de.params['A1'].value
    A1_err = result_de.params['A1'].stderr
    t_stat = A1_val / A1_err
    pvalue_ob = 2 * (1 - stats.t.cdf(abs(t_stat), df))
    if np.isnan(pvalue_ob):
        pvalue_ob = 1.0
    return result_de, result_se, pvalue_ob
    # return result_de,  result_ca

def plot_multi_fits(t, y, result_de, result_se, title="Multi-model fitting"):
    # Generate fitted curves
    t_axis=np.linspace(0.0,1.0,num=10001,endpoint=True)

    fit_de = double_exp_func(
        t_axis,
        result_de.params['A1'].value,
        result_de.params['A2'].value,
        result_de.params['k1'].value,
        result_de.params['k2'].value,
        result_de.params['C'].value,
    )
    fit_se = single_exp_func(
        t_axis,
        result_se.params['A2'].value,
        result_se.params['k2'].value,
        result_se.params['C'].value,
    )

    fit_de_ob = single_exp_func(
        t_axis,
        result_de.params['A1'].value,
        result_de.params['k1'].value,
        0,
    )
    fit_de_ad = single_exp_func(
        t_axis,
        result_de.params['A2'].value,
        result_de.params['k2'].value,
        result_de.params['C'].value,
    )


    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(t, y, 'o', markersize=3, label='Data', alpha=0.5)
    plt.plot(t_axis, fit_de, '-', linewidth=2, label='Double Exp Fit')
    plt.plot(t_axis, fit_de_ob, '-', linewidth=2, label='OB_part')
    plt.plot(t_axis, fit_de_ad, '-', linewidth=2, label='Adap_part')
    plt.plot(t_axis, fit_se, '--', linewidth=2, label='Single Exp Fit')
    plt.xlabel('Time')
    plt.ylabel('Response')
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    
def model_func(t, A, K, C):
    return A * np.exp(K * t) + C

def fit_exp_nonlinear(t, y,p0=[-5, 0.01, 25]):
    bounds = ([-np.inf, -1, 0], [0, 1, 50])  # lower and upper bounds for A, K, C
    opt_parms, parm_cov = sp.optimize.curve_fit(model_func, t, y, p0,bounds=bounds, maxfev=10000)
    A, K, C = opt_parms
    return A, K, C



def fit_2exps_drop(t, y,p0=[500, 50, 0, -50, -10, 2]):
    bounds = ([0, -500, -50, -100,-20,0], [1000, 500, 50, -20,-1,500])  # lower and upper bounds for A, K, C
    opt_parms, parm_cov = sp.optimize.curve_fit(double_exp_func, t, y, p0,bounds=bounds, maxfev=10000)
    A1,A2,A3,k1,k2,C = opt_parms
    return A1,A2,A3,k1,k2,C


# Residuals with L1-like penalty
def residual_with_L1(params, t, y, l1_weight=1e-3):
    A1 = params['A1']
    A2 = params['A2']
    A3 = params['A3']
    k1 = params['k1']
    k2 = params['k2']
    C  = params['C']

    model = double_exp_func(t, A1, A2, A3, k1, k2, C)
    residual = y - model

    # L1 penalty term encourages sparsity
    l1_penalty = l1_weight * (abs(A1) + abs(A2) + abs(A3))
    return residual + l1_penalty

# Improved fitting function
def fit_2exps_with_l1(t, y, initials=None, l1_weight=1e-3):
    # Default initial values
    if initials is None:
        initials = np.array([500, 50, 0, -50, -10, 2])

    # Unpack initials
    A1_0, A2_0, A3_0, k1_0, k2_0, C_0 = initials

    # Define parameter set
    params = Parameters()
    params.add('A1', value=A1_0, min=0, max=1000)
    params.add('A2', value=A2_0, min=-500, max=500)
    params.add('A3', value=A3_0, min=-50, max=50)
    params.add('k1', value=k1_0, min=-500, max=-20)  # 2 ms to 50 ms
    params.add('k2', value=k2_0, min=-20, max=-1)    # 50 ms to 1000 ms 
    params.add('C',  value=C_0,  min=0, max=500)

    # Create Minimizer with the penalized residuals
    minner = Minimizer(residual_with_L1, params, fcn_args=(t, y, l1_weight))
    result = minner.minimize()

    # Print fitting report
    report_fit(result)

    return result

ephy_temp_meta = pd.read_csv("0308_room_temp_nwb_names.csv", names=['file_name'])
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
duration=0.6

from datetime import datetime
date_str = datetime.now().strftime('%m%d%y')
filename = f'./testoutput/{date_str}_M1_addon_no1AP.csv'


#SST
# n_cells = 1823    # PV: 888 SST: 1823
# halfcook_path='./SST_all_figures/Cell'
# filename = f'./testoutput/{date_str}_V1_SST_addon.csv'

#PV
# n_cells = 888    # PV: 888 SST: 1823
# halfcook_path='./PV_all_figures/Cell'
# filename = f'./testoutput/{date_str}_V1_PV_addon.csv'

df_cells =  pd.DataFrame({'Measured Inj': pd.Series(dtype='float'),
                          'LastSpike': pd.Series(dtype='float'),
                          'InitRate': pd.Series(dtype='float'),
                          'rheo1st': pd.Series(dtype='float'),
                         'rheo_fit': pd.Series(dtype='float'),
                         'rheo_diff':  pd.Series(dtype='int'),
                         'slope': pd.Series(dtype='float'),
                         'error_slope':  pd.Series(dtype='int'),
                         'n_fitsweep':  pd.Series(dtype='int'),
                         'Inj_max': pd.Series(dtype='float'),
                         'AI' : pd.Series(dtype='float'),
                         'rABC' : pd.Series(dtype='float'),
                         'AIerr': pd.Series(dtype='float'),
                         'K': pd.Series(dtype='float'),
                         'A':  pd.Series(dtype='float'),
                         'C':  pd.Series(dtype='float'),                            
                         'pvalue_ob': pd.Series(dtype='float'),
                         'Kob_de': pd.Series(dtype='float'),
                         'Aob_de':  pd.Series(dtype='float'),
                         'Kad_de': pd.Series(dtype='float'),
                         'Aad_de':  pd.Series(dtype='float'),
                         'C_de':  pd.Series(dtype='float'),               
                         'AIerr_de': pd.Series(dtype='float'),
                         'OBnum':  pd.Series(dtype='int'),   # Onset bursting number
                         'OBdrop':  pd.Series(dtype='float'),
                         'Detect_OBt':  pd.Series(dtype='float'),
                    })




#i_cell_vec= np.array([408,336,192,477,440])
# i_cell_vec= np.array([0])

# for i_cell in tqdm(range(10)):
for i_cell in tqdm(range(n_cells)):
# i_cell_vec= np.array([59])
# i_cell_vec= np.array([1569,247,59,617,97]) # Cells show strong OB

# for i_cell in i_cell_vec:
    
    #print(i_cell)
    TempPath='0310_room_figure/Cell'+str(i_cell)+'/'

    # TempPath=halfcook_path+str(i_cell)+'/'

    df_sweeps =  pd.DataFrame({'Inj': pd.Series(dtype='float'),
                             'nAP': pd.Series(dtype='int'),
                             'start_i': pd.Series(dtype='int'),
                             'end_i': pd.Series(dtype='int'),
                             'LastSpike':  pd.Series(dtype='int'),
                             'InitRate':  pd.Series(dtype='int'),
                             'AI' : pd.Series(dtype='float'),
                             'rABC' : pd.Series(dtype='float'),
                             'AIerr': pd.Series(dtype='float'),
                             'K': pd.Series(dtype='float'),
                             'A':  pd.Series(dtype='float'),
                             'C':  pd.Series(dtype='float'),                            
                             'pvalue_ob': pd.Series(dtype='float'),
                             'Kob_de': pd.Series(dtype='float'),
                             'Aob_de':  pd.Series(dtype='float'),
                             'Kad_de': pd.Series(dtype='float'),
                             'Aad_de':  pd.Series(dtype='float'),
                             'C_de':  pd.Series(dtype='float'),               
                             'AIerr_de': pd.Series(dtype='float'),
                             'OBnum':  pd.Series(dtype='int'),   # Onset bursting number
                             'OBdrop':  pd.Series(dtype='float'),
                             'Detect_OBt':  pd.Series(dtype='float'),
                             'OBrate_rec':  pd.Series(dtype='float'),
                    })



    

    try:
        # Extract essentials from the HWandAP.csv
        df_APs= pd.read_csv(TempPath+'HWandAP.csv')
        n_APs= df_APs.shape[0]
        ind_sweep=  0
        ind_all_APs= 0
        while ind_all_APs< n_APs:
            if df_APs.at[ind_all_APs,'Ind_AP']==0:
                df_sweeps.at[ind_sweep,'Inj']=df_APs.at[ind_all_APs,'Inj']
                df_sweeps.at[ind_sweep,'start_i']=ind_all_APs
                if ind_all_APs>0:
                    df_sweeps.at[ind_sweep-1,'end_i']=ind_all_APs-1;
                    df_sweeps.at[ind_sweep-1,'nAP']=ind_all_APs-df_sweeps.at[ind_sweep-1,'start_i']
                ind_sweep+=1
            ind_all_APs +=1
        if ind_all_APs == n_APs:
            df_sweeps.at[ind_sweep-1,'end_i']=ind_all_APs-1;
            df_sweeps.at[ind_sweep-1,'nAP']=ind_all_APs-df_sweeps.at[ind_sweep-1,'start_i']
        
        # For some reason, I thought the soring is bad for M1 data but needed for V1 data. 
        # It is reflected in the results if run wiht sorted. Anyway, I am not using fitted rheobase. 
        # ending of regist the start and end of the APs. I need to change here
        # From the logic, I think I need the sort.  But somehow I remember there are two type of data,
        # such that I cannot do that
        #df_sweeps=df_sweeps.sort_values('Inj')   # Removing this line is the key comparing to getrheo.py
        df_sweeps = df_sweeps.reset_index(drop=True)
        
        Inj_vec=df_sweeps['Inj'].to_numpy()
        rate_vec=df_sweeps['nAP'].to_numpy()

        maxRate_ind= rate_vec.argmax(axis=0)    
        Ind= np.where( (Inj_vec>0)& (rate_vec>0))[0]
        Ind=Ind[Ind<=maxRate_ind]
        rheo1st= min(Inj_vec[Ind])   #  if sourted:  Inj_vec[Ind[0]] 
        df_cells.at[i_cell,'rheo1st']=rheo1st
        
        Ind2AP=np.where( (Inj_vec>0) & (rate_vec>1) & (rate_vec<=30))[0]
        Ind2AP=Ind2AP[Ind2AP<=maxRate_ind]
        df_cells.at[i_cell,'n_fitsweep']= len(Ind2AP)

        if len(Ind2AP)==0:
            print('Error: ' +str(i_cell)+ 'No valid sweeps.')
            df_cells.at[i_cell,'rheo_fit']= np.nan
            df_cells.at[i_cell,'slope']= np.nan
            df_cells.at[i_cell,'rheo_diff']= np.nan
            df_cells.at[i_cell,'OBnum']= np.nan
            df_cells.at[i_cell,'OBdrop']= np.nan

        else:
            

            n_sweep=len(Inj_vec)
            #Ind_last_temp= max(Ind[-1],Ind2AP[0]+4)
            #Ind_last= min(maxRate_ind,Ind2AP[-1])+1   # Using maxRate_ind instead.
            Ind_last=Ind2AP[-1]+1
            
            Inj_vec_fit=Inj_vec[Ind[0]:Ind_last]
            rate_vec_fit=rate_vec[Ind[0]:Ind_last]
            
            Inj_vec_more= np.array([-40,-30,-20,-10])+rheo1st
            rate_vec_more= np.array([0,0,0,0])
            
            Inj_vec_fit=np.concatenate((Inj_vec_more,Inj_vec_fit))
            rate_vec_fit=np.concatenate((rate_vec_more,rate_vec_fit))
            
            fig=plt.figure(figsize=(8, 6))
            plt.plot(Inj_vec,rate_vec,'.',markersize=15)
        
            try:
                rheo0_vec=np.array([-20,-10,0,0,10,10,20])+rheo1st
                k0_vec=np.array([0,0.2,0.8,0.1,0.8,0.1,0.8,1])

                err_min=100
                k=np.nan
                Rheo=np.nan
                fit_rate=np.nan
                for rheo0,k0 in zip(rheo0_vec,k0_vec):
                    try:
                        Rheo_temp, k_temp=  fit_relu(Inj_vec_fit,rate_vec_fit,[rheo0,k0])  
                        fit_rate_temp=ReLU_func(Inj_vec_fit, Rheo_temp,k_temp)
                        err_temp=np.sqrt(np.mean((np.abs(fit_rate_temp-rate_vec_fit))**2)) 
                        if (err_temp<err_min) & (Rheo_temp>0):
                            Rheo=Rheo_temp
                            k=k_temp
                            fit_rate=fit_rate_temp
                            err_min=err_temp
                    except:
                        continue
                plt.plot(Inj_vec_fit,fit_rate,'-.',linewidth=1,markersize=5)

            except:
                Rheo=np.nan
                k=np.nan

            df_cells.at[i_cell,'rheo_fit']= Rheo-rheo1st
            df_cells.at[i_cell,'rheo_diff']= Rheo
            df_cells.at[i_cell,'slope']= k/duration  
            df_cells.at[i_cell,'error_slope']= err_min 
            
            fig.savefig(TempPath+'fIfitRedo'+'.jpg',bbox_inches='tight')
        max_nAP=rate_vec[maxRate_ind]
        df_cells.at[i_cell,'Inj_max']=Inj_vec[maxRate_ind]
        df_cells.at[i_cell,'nAPmax']=max_nAP
        df_cells.at[i_cell,'nSweep']=maxRate_ind
        
         # Now let's do some fitting of AI
        for i_sweep in range(maxRate_ind+1):
            if df_sweeps.at[i_sweep,'nAP']> 5: #Used to be at least 5AP, now change to 10AP for fitting 6 parameters.
                                                # Change to 5 because the if is shorted. Let's check this
                t_spikes= df_APs.loc[df_sweeps.at[i_sweep,'start_i']:df_sweeps.at[i_sweep,'end_i'],'t_spike'].values
                df_sweeps.at[i_sweep,'LastSpike']=t_spikes[-1]
                t_spikes=t_spikes[1:]   
                t_isi= np.diff(t_spikes)
                df_sweeps.at[i_sweep,'InitRate']=1/t_isi[0]

                # drop the spikes if t_isi<0.002 (2 ms 500Hz)
                flag_longenough= t_isi>= 0.001
                t_spikes=t_spikes[np.append(flag_longenough,True)]
                t_isi= np.diff(t_spikes)
                if not all(flag_longenough):
                    print('Warning, some spikes are too close to each other.Icell'+str(i_cell)+ 'i_sweep'+str(i_sweep))    
                    logging.warning('Warning, some spikes are too close to each other.Icell'+str(i_cell)+ 'i_sweep'+str(i_sweep))
                FI_vec=1/t_isi
                try:    
                    K0_vec=np.array([10,10,10,10])
                    AI0_vec=np.array([10,20,0.1,1])
        
                    err_min= 0.5
                    err_de_min=0.5
                    Best_Fitted_Result = None
                    Best_Fitted_de_Result = None
                    p_ob_best=1
                    for K0,AIamp in zip(K0_vec,AI0_vec):
                        try:
                            test_inits= np.array([500, AIamp, 100, K0, 2])
                            Fitted_de, Fitted_se,p_ob= multi_fitting(t_spikes[:-1], FI_vec,initials=test_inits)   
                            '''
                            report_fit(Fitted_de)
                            report_fit(Fitted_se)
                            plot_multi_fits(t_spikes[:-1], FI_vec, Fitted_de, Fitted_se)
                            '''
                            y_fit=eval_fit(single_exp_func, Fitted_se.params, t_spikes[:-1])
                            Err_FI = np.sqrt(np.mean( (1-FI_vec/y_fit)**2   ))
                            fit_sanity =  eval_fit(single_exp_func, Fitted_se.params, 0)
                            
        
                            #fit_sanity= model_func(0, Atemp, Ktemp, Ctemp)                   
                            if (Err_FI<err_min) & (fit_sanity>=0):
                                Best_Fitted_Result=Fitted_se
                                err_min=Err_FI
        
                            y_fit_de=eval_fit(double_exp_func, Fitted_de.params, t_spikes[:-1])
                            Err_FI_de = np.sqrt(np.mean( (1-FI_vec/y_fit_de)**2   ))
                            
        
                            if (Err_FI_de<err_de_min):
                                Best_Fitted_de_Result=Fitted_de
                                err_de_min=Err_FI_de
                                p_ob_best=p_ob
                                AIamp_best=AIamp
                            '''
                            report_fit(Best_Fitted_de_Result)
                            report_fit(Best_Fitted_Result)
                            plot_multi_fits(t_spikes[:-1], FI_vec, Fitted_de, Fitted_se)
                            '''
                        except:
                            continue
                    '''
                    plot_multi_fits(t_spikes[:-1], FI_vec, Best_Fitted_de_Result, Best_Fitted_Result)
                    '''
                    t_axis=np.linspace(0.0,duration,int(duration/0.0001+1),endpoint=True)
                    fit_y =  eval_fit(single_exp_func, Best_Fitted_Result.params, t_axis)
                    AI_sweep=1-fit_y[-1]/fit_y[0]
                    ABC_vec=-fit_y+ (fit_y[-1]+  (fit_y[0]-fit_y[-1])*(duration-t_axis)/duration)
                    rABC=sum(ABC_vec)/len(ABC_vec)/(fit_y[0]-fit_y[-1])*2
                    
                    if i_sweep==maxRate_ind:
                        fig=plt.figure(figsize=(4,3))
                        plt.plot(t_spikes[:-1],FI_vec,markersize=10)
                        plt.plot(t_axis, fit_y, '.-',linewidth=3) 
                        fig.savefig(TempPath+'MaxSweepAIfit'+'.jpg',bbox_inches='tight')
        
                    # For the ISI:
                    num_ob=0
                    for isi in t_isi:
                        if isi<0.01:
                            num_ob+=1
                        else:
                            break
                    in_window = (t_spikes[1:] >= 0.03) & (t_spikes[1:] <=0.07)
                    if np.any(in_window):
                        FI_avg_window = np.mean(FI_vec[in_window])
                        OBdrop=1- FI_avg_window/FI_vec[0]
                    else:
                        OBdrop=np.nan    
                        
                    # get the readout, 
        
                     
                    df_sweeps.at[i_sweep,'AI']=AI_sweep
                    df_sweeps.at[i_sweep,'rABC']=rABC
                    df_sweeps.at[i_sweep,'AIerr']=err_min
        
                    #including the error of the fitting now
                    df_sweeps.at[i_sweep,'A']=Best_Fitted_Result.params['A2'].value
                    df_sweeps.at[i_sweep,'K']=Best_Fitted_Result.params['k2'].value
                    df_sweeps.at[i_sweep,'C']=Best_Fitted_Result.params['C'].value
                    df_sweeps.at[i_sweep,'pvalue_ob']=p_ob_best
                    df_sweeps.at[i_sweep,'Aob_de']=Best_Fitted_de_Result.params['A1'].value
                    df_sweeps.at[i_sweep,'Kob_de']=Best_Fitted_de_Result.params['k1'].value
                    df_sweeps.at[i_sweep,'Aad_de']=Best_Fitted_de_Result.params['A2'].value            
                    df_sweeps.at[i_sweep,'Kad_de']=Best_Fitted_de_Result.params['k2'].value
                    df_sweeps.at[i_sweep,'C_de']=Best_Fitted_de_Result.params['C'].value     
                    df_sweeps.at[i_sweep,'AIerr_de']=err_de_min   
                    # Old staff
                    Aad_de=Best_Fitted_de_Result.params['A2'].value           
                    Kad_de=Best_Fitted_de_Result.params['k2'].value
                    C_de=Best_Fitted_de_Result.params['C'].value  
                    eval_init_rate=single_exp_func(t_spikes[0],Aad_de,Kad_de,C_de)
                    
                    OBrate_rec= df_sweeps.at[i_sweep,'InitRate']-eval_init_rate
                    df_sweeps.at[i_sweep,'OBrate_rec']=OBrate_rec
                    
                    df_sweeps.at[i_sweep,'OBnum']=num_ob
                    df_sweeps.at[i_sweep,'OBdrop']=OBdrop
        
                except:
                    df_sweeps.at[i_sweep,'AI']=np.nan
                    df_sweeps.at[i_sweep,'rABC']=np.nan
                    df_sweeps.at[i_sweep,'AIerr']=np.nan
        
                    #including the error of the fitting now
                    df_sweeps.at[i_sweep,'A']=np.nan
                    df_sweeps.at[i_sweep,'K']=np.nan
                    df_sweeps.at[i_sweep,'C']=np.nan
                    df_sweeps.at[i_sweep,'pvalue_ob']=np.nan
                    df_sweeps.at[i_sweep,'Aob_de']=np.nan
                    df_sweeps.at[i_sweep,'Kob_de']=np.nan
                    df_sweeps.at[i_sweep,'Aad_de']=np.nan           
                    df_sweeps.at[i_sweep,'Kad_de']=np.nan
                    df_sweeps.at[i_sweep,'C_de']=np.nan   
                    df_sweeps.at[i_sweep,'AIerr_de']=np.nan  
                    # Old staff
                    df_sweeps.at[i_sweep,'OBnum']=np.nan
                    df_sweeps.at[i_sweep,'OBdrop']=np.nan
                            
        # Just for plotting the dFI_vec curve     
        fig=plt.figure(figsize=(6, 6))
               
        for i_sweep in range(maxRate_ind+1):
            if df_sweeps.at[i_sweep,'nAP']> 5: #at least have 5 ISIs.
                t_spikes= df_APs.loc[df_sweeps.at[i_sweep,'start_i']:df_sweeps.at[i_sweep,'end_i'],'t_spike'].values
                df_sweeps.at[i_sweep,'LastSpike']=t_spikes[-1]
                t_isi= np.diff(t_spikes)
                # drop the spikes if t_isi<0.002 (2 ms 500Hz)
                flag_longenough= t_isi>= 0.001
                t_spikes=t_spikes[np.append(flag_longenough,True)]
                t_isi= np.diff(t_spikes)
                if not all(flag_longenough):
                    print('Warning, some spikes are too close to each other.Icell'+str(i_cell)+ 'i_sweep'+str(i_sweep))    
                    logging.warning('Warning, some spikes are too close to each other.Icell'+str(i_cell)+ 'i_sweep'+str(i_sweep))
                FI_vec=1/t_isi
                dFI_vec= np.diff(FI_vec)/(t_spikes[2:]-t_spikes[:-2])*2    # FI(2.5)-FI(1.5)/( (t3-t1)/2)
                plt.plot(t_spikes[1:-1],dFI_vec)
                plt.plot(t_spikes[1:-1],dFI_vec,'.k')
                small_dFI=dFI_vec>-500
                if np.any(small_dFI):
                    OBt_ind= np.argmax(small_dFI)
                    if OBt_ind==0:
                        df_sweeps.at[i_sweep,'Detect_OBt']= 0
                    else:
                        Detect_OBt= t_spikes[OBt_ind]
                        df_sweeps.at[i_sweep,'Detect_OBt']= Detect_OBt
                else:
                    Detect_OBt= t_spikes[-1]
                    df_sweeps.at[i_sweep,'Detect_OBt']= Detect_OBt
                    
                
        plt.xlabel('Spike Time')
        plt.ylabel('dFI/dt (Hz/s)')
        plt.xlim([0, duration])
        fig.savefig(TempPath+'dFIcurve1s'+'.jpg',bbox_inches='tight')
        plt.xlim([0, 0.1])
        fig.savefig(TempPath+'dFIcurve'+'.jpg',bbox_inches='tight')
                            
                    
        df_sweeps.to_csv(TempPath+'Sweeps.csv',index=True)      
        
        if max_nAP>=6:
            dAdap_vec=df_sweeps['AI'].to_numpy()
            LastSpike_vec=df_sweeps['LastSpike'].to_numpy()
            # legit_sweeps= (Inj_vec>rheo1st*1.5) & (Inj_vec<rheo1st*4) & ~np.isnan(dAdap_vec) &  (LastSpike_vec>0.5)  
            legit_sweeps= (Inj_vec>rheo1st*1.5) & (Inj_vec<rheo1st*4) & (dAdap_vec>-0.2) &  (LastSpike_vec>0.5)  
            if np.any(legit_sweeps):
                hero_ind= np.argmax(legit_sweeps)
                df_cells.at[i_cell,'hero or Max']= 'hero'
            else:
                hero_ind= maxRate_ind
                df_cells.at[i_cell,'hero or Max']= 'Max'
                
            df_cells.at[i_cell,'Measured Inj']=df_sweeps.at[hero_ind,'Inj']
            df_cells.at[i_cell,'LastSpike']=df_sweeps.at[hero_ind,'LastSpike']
            df_cells.at[i_cell,'AI']=df_sweeps.at[hero_ind,'AI']
            df_cells.at[i_cell,'rABC']=df_sweeps.at[hero_ind,'rABC']
            df_cells.at[i_cell,'AIerr']=df_sweeps.at[hero_ind,'AIerr']
            df_cells.at[i_cell,'K']=df_sweeps.at[hero_ind,'K']
            df_cells.at[i_cell,'A']=df_sweeps.at[hero_ind,'A']
            df_cells.at[i_cell,'C']=df_sweeps.at[hero_ind,'C']
            
            pvalue_vec=df_sweeps['pvalue_ob'].to_numpy()

            
            if np.any(legit_sweeps):
                legit_indices = np.where(legit_sweeps)[0]
                pvalue_legit_vec=pvalue_vec[legit_sweeps]
                min_pvalue_ind= np.argmin(pvalue_legit_vec)
                legit_length=len(legit_indices) 
                orig_pvalue_ind= legit_indices[min_pvalue_ind]
            else:
                orig_pvalue_ind= maxRate_ind
                legit_length=1

            df_cells.at[i_cell,'InitRate']=df_sweeps.at[orig_pvalue_ind,'InitRate']
            df_cells.at[i_cell,'pvalue_ob']=df_sweeps.at[orig_pvalue_ind,'pvalue_ob']*legit_length
            df_cells.at[i_cell,'Kob_de']=df_sweeps.at[orig_pvalue_ind,'Kob_de']
            df_cells.at[i_cell,'Aob_de']=df_sweeps.at[orig_pvalue_ind,'Aob_de']
            df_cells.at[i_cell,'Kad_de']=df_sweeps.at[orig_pvalue_ind,'Kad_de']
            df_cells.at[i_cell,'Aad_de']=df_sweeps.at[orig_pvalue_ind,'A']
            df_cells.at[i_cell,'C_de']=df_sweeps.at[orig_pvalue_ind,'C_de']
            df_cells.at[i_cell,'AIerr_de']=df_sweeps.at[orig_pvalue_ind,'AIerr_de']
            df_cells.at[i_cell,'OBrate_rec']=df_sweeps.at[orig_pvalue_ind,'OBrate_rec']

            
            
            df_cells.at[i_cell,'OBnum']=df_sweeps.at[orig_pvalue_ind,'OBnum']
            df_cells.at[i_cell,'OBdrop']=df_sweeps.at[orig_pvalue_ind,'OBdrop']
            df_cells.at[i_cell,'Detect_OBt']=df_sweeps.at[orig_pvalue_ind,'Detect_OBt']
        
            fig=plt.figure(figsize=(6, 6))
            plt.subplot(211)
            t_spikes= df_APs.loc[df_sweeps.at[hero_ind,'start_i']:df_sweeps.at[hero_ind,'end_i'],'t_spike'].values
            t_isi= np.diff(t_spikes)
            flag_longenough= t_isi>= 0.002
            t_spikes=t_spikes[np.append(flag_longenough,True)]
            t_isi= np.diff(t_spikes)
            FI_vec=1/t_isi
        
            # Still st needed here to calculate the FI_vec.
        
            plt.plot(t_spikes[:-1],FI_vec)
            t_axis=np.linspace(0.0,1.0,num=int(duration/0.0001+1),endpoint=True)
            if df_cells.at[i_cell,'AIerr_de']<0.5:
                A1 =  df_sweeps.at[hero_ind,'Aob_de']
                A2 =  df_sweeps.at[hero_ind,'Aad_de']
                k1 =  df_sweeps.at[hero_ind,'Kob_de']
                k2 =  df_sweeps.at[hero_ind,'Kad_de']
                C  =  df_sweeps.at[hero_ind,'C_de']
                fit_y_de = double_exp_func(t_axis, A1,A2,k1,k2,C)
                plt.plot(t_axis, fit_y_de, '-',linewidth=1,label='fitted_double_exp') 
            if df_cells.at[i_cell,'AIerr']<0.5:
                A=    df_sweeps.at[hero_ind,'A']
                K=    df_sweeps.at[hero_ind,'K']
                C=    df_sweeps.at[hero_ind,'C']
                fit_y = single_exp_func(t_axis, A, K, C)
                plt.plot(t_axis, fit_y, '-',linewidth=1,label='fitted_single_exp') 
            plt.legend()    
            plt.subplot(212)
            plt.plot(df_sweeps.loc[0:maxRate_ind,'Inj'].values,df_sweeps.loc[0:maxRate_ind,'AI'].values, label='AI')
            plt.plot(df_sweeps.loc[0:maxRate_ind,'Inj'].values,df_sweeps.loc[0:maxRate_ind,'rABC'].values, label='rABC')
            plt.plot(df_sweeps.loc[0:maxRate_ind,'Inj'].values,df_sweeps.loc[0:maxRate_ind,'AIerr'].values, label='AIerr')
            plt.plot(df_sweeps.loc[0:maxRate_ind,'Inj'].values,df_sweeps.loc[0:maxRate_ind,'AI'].values,'.k')
            plt.plot(df_sweeps.loc[0:maxRate_ind,'Inj'].values,df_sweeps.loc[0:maxRate_ind,'rABC'].values,'.k')
            plt.plot(df_sweeps.loc[0:maxRate_ind,'Inj'].values,df_sweeps.loc[0:maxRate_ind,'AIerr'].values,'.k')
            plt.legend()    
            fig.savefig(TempPath+'AdapfitatHero'+'.jpg',bbox_inches='tight')
                    
    #find the rheo1st.
    #Inj_vec= df_sweeps[['file_name']]
            
    
    except FileNotFoundError:
        print('Error: ' +str(i_cell)+ 'Cell File not found.')
        df_cells.at[i_cell,'rheo1st']=np.nan
        df_cells.at[i_cell,'rheo_fit']= np.nan
        df_cells.at[i_cell,'Inj_max']= np.nan
        df_cells.at[i_cell,'slope']= np.nan
        df_cells.at[i_cell,'AI_cell']= np.nan
        df_cells.at[i_cell,'nAPmax']=np.nan
        df_cells.at[i_cell,'nSweep']= np.nan
        df_cells.at[i_cell,'K']=np.nan
        df_cells.at[i_cell,'A']=np.nan
        df_cells.at[i_cell,'C']=np.nan

    plt.close('all')
    
'''
    except:
        print('Error: may come from the min firing rate>30 for every sweep. Marked with maxnAP=-100')
        df_cells.at[i_cell,'rheo1st']=np.nan
        df_cells.at[i_cell,'rheo_fit']= np.nan
        df_cells.at[i_cell,'Inj_max']= np.nan
        df_cells.at[i_cell,'slope']= np.nan
        df_cells.at[i_cell,'AI_cell']= np.nan
        df_cells.at[i_cell,'nAPmax']=-100
        df_cells.at[i_cell,'nSweep']= np.nan
'''
logging.shutdown()


'''
filename = f'./testoutput/{date_str}_till{i_cell}V1_SST_addon.csv'
'''

df_cells.to_csv(filename,index=True)

    
    
    
    
    
    
    
    
    
