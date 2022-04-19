# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 11:35:24 2021

@author: Yael
"""
import matplotlib.pyplot as plt
import numpy as np
#import time
#import seaborn as sns
#from scipy.interpolate import make_interp_spline, BSpline
#from scipy.interpolate import interp1d



#PARAMETERS
#define parameters
# N = 10000
#N = 10607700 # Wuhan
# N = 308000000 # USA
# Z = 3.5 #incubation period, days from exposed to infected 3.5

# Dp = 3.5 #days for preclinical to become clinical, 3.5 (6 days until symptom onset+3 days until covid test results)
# Dc = 10.4 #days for clinical to recover 10.4
# Ds = 7 #days 

N = 8804190 # USA
#N=8000000
Z = 5 #incubation period, days from exposed to infected 5
#β = 0.85 #transmission rate
α1 = 0.65 #wildtype virus, fraction become clinical, (0.6)
Dp = 3 #days for preclinical to become clinical, 3 
Dc = 14 #days for clinical to recover 13
Ds = 6 #days for subclinical to recover, 6

α2 = 0.05 #attenuated virus, fraction become clinical, 0.05

μ = 0.6 #factor by which subclinical is less contagious than preclinical (number is correct for symptomatic?)
#μ2 = 0.3 #factor by which subclinical is less contagious than preclinical (number is correct for symptomatic?)
#β = 0.85 #transmission rate


a = 1
#Δ = 6 #lag current *reported* number of cases and actual infectious individuals
t = np.arange(1000)



def simulate_with_smart_NPI(Z, Dp, Dc, Ds, βs_range, α1, α2, μ, a, Δ,start,on_off=[],thresh=-1):
 
    E0, Is0 = 500, 1500
    #E0, Is0 = 2000, 0
    
    S0 = N - E0 - Is0
    βs = np.zeros(t.size) 
    y = np.zeros((t.size, 11))
    #.     S   E   Ip  Ic   Is   R  / E    Ip  Ic   Is   R
    y[0] = S0, E0/2, 0, 0, Is0/2, 0, E0/2, 0, 0, Is0/2, 0
    
    
    Ic_curr_week_sum = 0
    Ic_last_week_sum = 0 
    β_ = βs_range.max()
    βs[0] = β_
    R_0 = 0
    
    for t_ in t[1:]:
        S, E, Ip, Ic, Is, R, EE, IIp, IIc, IIs, RR = y[t_ - 1]
        #start the count only from start
        
        if len(on_off)>0:
                if (t_>on_off[0] and t_<on_off[1]) or (t_>on_off[2] and t_<on_off[3]):
                        β_ = βs_range.min()
                else:
                        β_ = βs_range.max()
        
        
        
        if t_==start:
                Ic_last_week_sum = sum(y[max(t_-Δ,0):t_,3])+sum(y[max(t_-Δ,0):t_,8])
                β_ = βs_range.max()
        if t_>start:
                if (t_-start)%Δ==0: 
                    #IcΔ = y[t_ - Δ - 1, 3] # 3 is the index of Ic and 8 of IIc
                    Ic_curr_week_sum = sum(y[t_-Δ:t_,3])+sum(y[t_-Δ:t_,8])
                    if Ic_last_week_sum>0:
                        
                        if thresh==-1:
                                R_0=(Ic_curr_week_sum/Ic_last_week_sum)**(4/7)
                                β_ = βs_range.min() if R_0>1 else βs_range.max()

                        else:
                                #R_0 = (Ic_curr_week_sum/Δ) / (thresh*N)
                                mean_Ic = (Ic_curr_week_sum/Δ)
                                if β_ == βs_range.min():
                                        if mean_Ic < 0.01*N:
                                                β_ = βs_range.max()
                                else: #if beta = max beta
                                        if mean_Ic>thresh*N:
                                                β_ = βs_range.min()
                        
                    #elif Ic_curr_week_sum>0:    
                    #    β_ = βs_range.min()
                        
                    Ic_last_week_sum = Ic_curr_week_sum
    
            
        #β_ = β - a * IcΔ / N
        βs[t_] = β_
        dy = [
            -β_ * S * (Ip + μ * Is + IIp + μ * IIs) / N, # S
            # wt
            β_ * S * (Ip + μ * Is) / N - E / Z,   # E
            α1 * E / Z - Ip / Dp,                 # Ip
            Ip / Dp - Ic / Dc,                    # Ic
            (1 - α1) * E / Z - Is / Ds,           # Is
            Is / Ds + Ic / Dc,                    # R
            # mut
            β_ * S * (IIp + μ * IIs) / N - EE / Z, # E
            α2 * EE / Z - IIp / Dp,                # Ip
            IIp / Dp - IIc / Dc,                   # Ic
            (1 - α2) * EE / Z - IIs / Ds,          # Is
            IIs / Ds + IIc / Dc,                   # R      
        ]
        y[t_] = y[t_ - 1] + dy
    return y, βs




def plot_timeline_with_npi_subplot(axs,axs_col, βs_range, on_off_=[300,300,300,300],labels=[]):
        
        #plt.close('all')
        plt.rcParams.update({'font.size': 16})

        delta = 1
        start_ = 1000
        
        
        #LEFT PANEL
        #on_off_=[300,300,300,300]
        
        #RIGHT PANEL
        #on_off_=[31,110,300,300]
        
        y, βs = simulate_with_smart_NPI(Z, Dp, Dc, Ds, 
                                        βs_range, α1, α2,
                                        μ=0.65, a=1, Δ=delta,start=start_,on_off=on_off_)
       
      # 0, 1, 2,  3,  4,  5, 6,   7,   8,   9,  10
        S, E, Ip, Ic, Is, R, EE, IIp, IIc, IIs, RR = y.T/N
                
        #plt.plot(βs,'y')
        
        
        #R,RR=5,10
        #first,second=5,10
        #first,second = (2,3,4),8
        
        #sum of I: 2:5, 7:10
        
        #sum of Ic, Is: 3:5, 8:10

        #R,RR: 5:6, 10:11
        #E,EE: 1:2, 6:7
        #fig, axs = plt.subplots(3,figsize=(6,6),sharex=True,tight_layout=True)
        #fig.suptitle(f'$\Delta={delta}$, start = {start_}, NPI strength = {1-βs_range.min()/βs_range.max()}')
        
         
        #PATCH TO MAKE THE Y-AXIS PERCENTILE INSTEAD OF FRACTION
        N_div = N/100 
        
        
        #SHADE NPI REGION
        if axs_col!=0:
                days = np.array(range(0,len(y.T[3])))
                where_fill = np.bitwise_and(days>on_off_[0],days<=on_off_[1])   
        
                axs[0][axs_col].fill_between( x=days, y1=np.maximum(y.T[5]/N_div,y.T[10]/N_div), y2=0, where=where_fill,
                                             facecolor='gray', color='dimgray',alpha=0.2)
                                                                                    
                axs[1][axs_col].fill_between( x=days, y1=np.maximum(y.T[3]/N_div,y.T[8]/N_div), y2=0, where=where_fill,
                                             facecolor='gray', color='dimgray',alpha=0.2)
        
                axs[2][axs_col].fill_between( x=days, y1=np.maximum(y.T[1]/N_div,y.T[6]/N_div), y2=0, where=where_fill,
                                             facecolor='gray', color='dimgray',alpha=0.2)   
                
                
                axs[1][axs_col].text(x=78,
                            y=0.78,
                            s="NPIs",
                            fontdict=dict(fontsize=11,fontstyle='oblique',fontweight='semibold'))
                
                
                # axs[1][axs_col].text(x=112,
                #             y=0.66,
                #             s="NPIs",
                #             fontdict=dict(fontsize=12,fontstyle='oblique',fontweight='semibold'))
                
                
                #ARROW TO SHOW SWITCHING EVOLUTION
                # axs[0][axs_col].annotate("",xy=(140,0.41*100),xytext=(140,0.32*100),
                #                   arrowprops=dict(arrowstyle="fancy",facecolor='darkgray',
                #                                   alpha=0.5,edgecolor='dimgray'))
        
        
        #####
        
                
        axs[0][axs_col].plot((y.T[5])/N_div,'orange',markevery=on_off_,label='Virulent')
        axs[0][axs_col].plot((y.T[10])/N_div,'purple',label='Attenuated',ls='--')
        
        axs[0][axs_col].text(0.1, 0.9, labels[0], transform=axs[0][axs_col].transAxes,
                fontsize=16, fontweight='bold', va='top', ha='right')
        
        axs[1][axs_col].plot(np.sum(y.T[3:4]/N_div,axis=0),'orange',markevery=on_off_,label='Virulent')
        axs[1][axs_col].plot(np.sum(y.T[8:9]/N_div,axis=0),'purple',label='Attenuated',ls='--')
        
        axs[1][axs_col].text(0.1, 0.9, labels[1], transform=axs[1][axs_col].transAxes,
                fontsize=16, fontweight='bold', va='top', ha='right')
        
        
        axs[2][axs_col].plot(np.sum(y.T[1:2]/N_div,axis=0),'orange',markevery=on_off_,label='Virulent')
        axs[2][axs_col].plot(np.sum(y.T[6:7]/N_div,axis=0),'purple',label='Attenuated',ls='--')
        
        axs[2][axs_col].text(0.1, 0.9, labels[2], transform=axs[2][axs_col].transAxes,
                fontsize=16, fontweight='bold', va='top', ha='right')
        
        if axs_col==0:
                
                axs[0][axs_col].set_ylabel('% Recovered')
                axs[1][axs_col].set_ylabel('% Symptomatic')
                axs[2][axs_col].set_ylabel('% Exposed')
                
                axs[2][axs_col].legend(borderpad=0.1, loc="upper right")
        
               
        axs[2][axs_col].set_xlabel('Days')
        
        # axs[0][axs_col].set_xlim(0,250)
        # axs[1][axs_col].set_xlim(0,250)
        # axs[2][axs_col].set_xlim(0,250)

        plt.setp(axs, xlim=(0,250))
        
        #Y-TICKS
        axs[0][axs_col].set_ylim([0,70])
        y_ticks = np.arange(0,70,20)
        axs[0][axs_col].set_yticks(y_ticks)
        
        #[1][axs_col].set_ylim([0,70])
        y_ticks = np.arange(0,16,5) if axs_col==0 else np.arange(0,8,2)
        axs[1][axs_col].set_yticks(y_ticks)
        
        y_ticks = np.arange(0,16,5) if axs_col==0 else np.arange(0,8,2)
        axs[2][axs_col].set_yticks(y_ticks)
        #####



#fig = plt.figure()
fig, axs = plt.subplots(3,2,figsize=(8,6),sharex=True,tight_layout=True)
        
plot_timeline_with_npi_subplot(axs,0,np.linspace(0.42, 1.2,20),labels=['a','c','e']) 
   
#plot_timeline_with_npi_subplot(axs,1,np.linspace(0.42, 1.2,20),
#                                on_off_=[31,150,300,300],labels=['b','d','f'])        

plot_timeline_with_npi_subplot(axs,1,np.linspace(0.42, 1.2,20),
                                on_off_=[15,110,300,300],labels=['b','d','f'])    

# plot_timeline_with_npi_subplot(axs,1,np.linspace(0.42, 1.2,20),
#                                 on_off_=[10,150,300,300],labels=['b','d','f'])        

plt.savefig('Figure_4.jpeg',dpi=300)
