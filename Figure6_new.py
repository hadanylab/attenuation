# -*- coding: utf-8 -*-
"""
Created on Fri Jan 15 12:58:11 2021

@author: Yael
"""

import matplotlib.pyplot as plt
import numpy as np
#import time
#from matplotlib.lines import Line2D
#import seaborn as sns

#PARAMETERS
#define parameters
# N = 10000
#N =7700 # Wuhan
N = 8804190 # USA
Z = 5 #incubation period, days from exposed to infected 3.5
#β = 1.2 #transmission rate
α1 = 0.65 #wildtype virus, fraction become clinical, 0.2 (actually likely higher, like 0.6)
α2 = 0.65 #attenuated virus, fraction become clinical, 0.05
Dp = 3 #days for preclinical to become clinical, 3.5 (6 days until symptom onset+3 days until covid test results)
Dc = 14 #days for clinical to recover 10.4
Ds = 6 #days for subclinical to recover



μ = 0.6 #factor by which subclinical is less contagious than preclinical (number is correct for symptomatic?)
μ2 = 0.6 #factor by which subclinical is less contagious than preclinical (number is correct for symptomatic?)

d1 = 0.9 #detectability of wt, this is true positive rate. 

a = 1
# = 6 #lag current *reported* number of cases and actual infectious individuals
t = np.arange(1000) 
Δ=0


def simulate_with_detectability(Z, Dp, Dc, Ds, β, α1, α2, μ, a, Δ, p, d1, d2, c, STRAIN=0):
    E0, Is0 = 500, 1500
    S0 = N - E0 - Is0
    βs = np.zeros(t.size) 
    y = np.zeros((t.size, 13))
    dy_arr = np.zeros((t.size, 13))
    perc_pos = np.zeros(t.size)
    detected = np.zeros(t.size)
    #.     S   E   Ip  Ic   Is   R, P  / E    Ip  Ic   Is   R ,PP
    if STRAIN==1:
            y[0] = S0, E0, 0, 0, Is0, 0, 0, 0, 0, 0, 0, 0 ,0
    elif STRAIN==2:
            y[0] = S0, 0, 0, 0, 0, 0, 0,E0, 0, 0, Is0, 0,0
    elif STRAIN==0:    
            y[0] = S0, E0/2, 0, 0, Is0/2, 0, 0, E0/2, 0, 0, Is0/2, 0,0
    
    
    
    βs[0] = β
    
    a=1-c

        
    for t_ in t[1:]:
        S, E, Ip, Ic, Is,  P, R, EE, IIp, IIc, IIs, PP, RR = y[t_ - 1]
        

        p_eff = p

        β_ = β
        βs[t_] = β_
        dy = [
            -β_ * S * (Ip + μ * Is + a*(IIp + μ2 * IIs)) / N, # S
            # wt
            β_ * S * (Ip + μ * Is) / N - E / Z,   # E
            α1 * E / Z - Ip / Dp - p_eff*d1*Ip,       # Ip
            Ip / Dp - Ic / Dc - p_eff*d1*Ic,                    # Ic
            (1 - α1) * E / Z - Is / Ds - p_eff*d1*Is ,          # Is
            p_eff*d1*(Ip+Is+Ic) - P/Dc,          # P
            Is / Ds + Ic / Dc + P/Dc,                    # R
            
            
            # mut
            β_ * S * a*(IIp + μ2 * IIs) / N - EE / Z, # EE
            α2 * EE / Z - IIp / Dp - p_eff*d2*IIp ,                # IIp
            IIp / Dp - IIc / Dc - p_eff*d2*IIc,                   # IIc
            (1 - α2) * EE / Z - IIs / Ds - p_eff*d2*IIs ,          # IIs
            p_eff*d2*(IIp+IIs+IIc) - PP/Dc,          # PP
            IIs / Ds + IIc / Dc + PP/Dc,                   # RR 
            
        ]
        y[t_] = y[t_ - 1] + dy
        
        
        dy_arr[t_] = dy
        perc_pos[t_] = (p_eff*d1*(Ip+Is+Ic)+p_eff*d2*(IIp+IIs+IIc))/(p_eff*(N-P-PP)) 
        detected[t_] = p_eff*d1*(Ip+Is+Ic)+p_eff*d2*(IIp+IIs+IIc) 
                
        
    return y, βs, perc_pos, detected, dy_arr


#get the results for Fig 6, per specific testing rate
#x_range == d2_range, p_ == testing rate, c == cost, β == transmission rate
def get_result_for_testing_rate(x_range,p_,c, β):
        
        res_actual = np.zeros(x_range.size)
        res_confirmed = np.zeros(x_range.size)
        res_perc_pos = np.zeros(x_range.size)
        
               
        #baseline == no test-evasive strain, no cost
        #STRAIN parameter is important
        y1, βs, perc_pos1, detected1, dy_arr = simulate_with_detectability(Z, Dp, Dc, Ds, β, α1, α2,
                                                 μ, a, 0, p_, d1, d1, c=0, STRAIN=1)
        
        
        S1, E1, Ip1, Ic1, Is1, P1, R1, EE1, IIp1, IIc1, IIs1, PP1, RR1 = y1.T
        daily_actual_cases_detectable = E1/Z + EE1/Z
        daily_confirmed_cases_detectable = detected1
        daily_perc_pos_detectable = perc_pos1
        
        peak_t = np.argmax(daily_actual_cases_detectable)
        peak_actual_cases_detectable = [daily_actual_cases_detectable[peak_t]]*len(x_range)
                     
        peak_t = np.argmax(daily_confirmed_cases_detectable)
        peak_confirmed_cases_detectable = [daily_confirmed_cases_detectable[peak_t]]*len(x_range)
         
        peak_t = np.argmax(daily_perc_pos_detectable)
        peak_perc_pos_detectable = [daily_perc_pos_detectable[peak_t]]*len(x_range)
        
        
        for i,d2_ in enumerate(x_range):
                
                
                #REMOVE FOR NOW
                #c_= 0 if d2_==1 else c
                c_ = c
                           
                #no detectable strain, only test-evasive
                #STRAIN parameter is important
                y2, βs, perc_pos2, detected2, dy_arr = simulate_with_detectability(Z, Dp, Dc, Ds, β, α1, α2,
                                                 μ, a, 0, p_, d2_, d2_, c_, STRAIN=2)
                        
                S2, E2, Ip2, Ic2, Is2, P2, R2, EE2, IIp2, IIc2, IIs2, PP2, RR2  = y2.T
                
                daily_actual_cases_evasive = E2/Z+EE2/Z
                daily_confirmed_cases_evasive = detected2
                daily_perc_pos_evasive = perc_pos2


                peak_t = np.argmax(daily_actual_cases_evasive)
                res_actual[i] = daily_actual_cases_evasive[peak_t]
                
                peak_t = np.argmax(daily_confirmed_cases_evasive)
                res_confirmed[i] = daily_confirmed_cases_evasive[peak_t]
                
                peak_t = np.argmax(daily_perc_pos_evasive)
                res_perc_pos[i] = daily_perc_pos_evasive[peak_t]
                
        return peak_actual_cases_detectable, peak_confirmed_cases_detectable, peak_perc_pos_detectable,res_actual, res_confirmed,res_perc_pos



#helper function to plot all sorts of scenarios for testing
def plot_test():
        
          p = 10.9/10**3
          β = 0.42
          
          y1, βs, perc_pos1, detected1, dy_arr = simulate_with_detectability(Z, Dp,
                                                                                   Dc, Ds, β,
                                                                                   α1, α2, μ,
                                                                                   a, Δ, p, 0, 0, c=0.01, STRAIN=2)
          S1, E1, Ip1, Ic1, Is1, R1, EE1, IIp1, IIc1, IIs1, RR1 = y1.T
          daily_actual_cases_detectable = E1/Z + EE1/Z
          #d1 = 1
          #R0 = (β/(Z*d1*p+1))*(α1*Dp/(Dp*d1*p+1))*((1-α1)*Ds*μ/(Ds*d1*p+1))

          plt.plot(daily_actual_cases_detectable,color='blue',
                   marker='*',markevery=[np.argmax(daily_actual_cases_detectable)])
          
          y1, βs, perc_pos1, detected1, dy_arr = simulate_with_detectability(Z, Dp,
                                                                                   Dc, Ds, β,
                                                                                   α1, α2, μ,
                                                                                   a, Δ, p, 0.25, 0.25, c=0.01, STRAIN=2)
          S1, E1, Ip1, Ic1, Is1, R1, EE1, IIp1, IIc1, IIs1, RR1 = y1.T
          daily_actual_cases_detectable = E1/Z + EE1/Z
          plt.plot(daily_actual_cases_detectable,color='green',
                   marker='*',markevery=[np.argmax(daily_actual_cases_detectable)])
          
          y1, βs, perc_pos1, detected1, dy_arr = simulate_with_detectability(Z, Dp,
                                                                                   Dc, Ds, β,
                                                                                   α1, α2, μ,
                                                                                   a, Δ, p, 0.5, 0.5, c=0.01, STRAIN=2)
          S1, E1, Ip1, Ic1, Is1, R1, EE1, IIp1, IIc1, IIs1, RR1 = y1.T
          daily_actual_cases_detectable = E1/Z + EE1/Z
          plt.plot(daily_actual_cases_detectable,color='green',
                   marker='*',markevery=[np.argmax(daily_actual_cases_detectable)])
          
          y1, βs, perc_pos1, detected1, dy_arr = simulate_with_detectability(Z, Dp,
                                                                                   Dc, Ds, β,
                                                                                   α1, α2, μ,
                                                                                   a, Δ, p, 0.9, 0.9, c=0.01, STRAIN=2)
          S1, E1, Ip1, Ic1, Is1, R1, EE1, IIp1, IIc1, IIs1, RR1 = y1.T
          daily_actual_cases_detectable = E1/Z + EE1/Z
          plt.plot(daily_actual_cases_detectable,color='green',
                   marker='*',markevery=[np.argmax(daily_actual_cases_detectable)])
          
          
          
          y1, βs, perc_pos1, detected1, dy_arr = simulate_with_detectability(Z, Dp,
                                                                                   Dc, Ds, β,
                                                                                   α1, α2, μ,
                                                                                   a, Δ, p, 1, 1, c=0, STRAIN=1)
          S1, E1, Ip1, Ic1, Is1, R1, EE1, IIp1, IIc1, IIs1, RR1 = y1.T
          daily_actual_cases_detectable = E1/Z + EE1/Z
          plt.plot(daily_actual_cases_detectable,color='red',marker='*',markevery=[np.argmax(daily_actual_cases_detectable)])
          
          
#plot_test()        
          


def plot_Fig6_timeline(d2_,p_arr,c_, β):
        
        color1 = "black"
        color2 = "forestgreen"
        
        fig, axs = plt.subplots(3,len(p_arr),figsize=(8,6),sharex=True,sharey='row', tight_layout=True)

        for i,p_ in enumerate(p_arr): 
                y1, βs, perc_pos1, detected1, dy_arr = simulate_with_detectability(Z, Dp,
                                                                                   Dc, Ds, β,
                                                                                   α1, α2, μ,
                                                                                   a, Δ, p_, 1, 1, c=0, STRAIN=1)
                
                S1, E1, Ip1, Ic1, Is1, R1, EE1, IIp1, IIc1, IIs1, RR1 = y1.T
                daily_actual_cases_detectable = E1/Z + EE1/Z
                daily_confirmed_cases_detectable = detected1
                
                y2, βs, perc_pos2, detected2, dy_arr = simulate_with_detectability(Z, Dp, 
                                                                                   Dc, Ds, β,
                                                                                   α1, α2,μ,
                                                                                   a, Δ, p_, d2_, d2_, c_, STRAIN=2)
                S2, E2, Ip2, Ic2, Is2, R2, EE2, IIp2, IIc2, IIs2, RR2 = y2.T
                daily_actual_cases_evasive = E2/Z+EE2/Z
                daily_confirmed_cases_evasive = detected2
                
                #fig, axs = plt.subplots(3,1,figsize=(6,6),sharex=True,tight_layout=True)
                plt.setp(axs, xlim=(0,600))
                axs[2][i].set_xlabel('Days')
                
                axs[0][i].plot(t,perc_pos1,label=rf'$d2={1}$',color=color1)
                axs[0][i].plot(t,perc_pos2,label=rf'$d2={d2_}$',color=color2)
                y_ticks = np.arange(0,0.2,0.05)
                axs[0][i].set_yticks(y_ticks)
                axs[0][i].set_yticklabels([f'{"{:.0f}".format(y*100)}%' for y in y_ticks])
        
                
                axs[1][i].set_title(f'Daily tests: {p_*10**3} per 1K people',fontsize=12)
        
                
                axs[1][i].plot(t,daily_confirmed_cases_detectable,label=rf'$d_2={1}$',color=color1)
                axs[1][i].plot(t,daily_confirmed_cases_evasive,label=rf'$d_2={d2_}$',color=color2)
                
                axs[1][i].set_yticks(np.arange(0,4*10**5,10**5))
                #axs[1][i].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
                
                if i==0:
                        axs[0][i].set_ylabel('Daily %\npositive tests')
                        axs[1][i].set_ylabel('Daily confirmed\ncases (per 1K)')
                        axs[2][i].set_ylabel('Daily actual\ncases (per 1K)')
                        
                        axs[1][i].legend(['Detectable','Test-evasive'],loc="upper left",borderpad=0.1)
                
                axs[2][i].plot(daily_actual_cases_detectable,'g',ls='-',color=color1)
                axs[2][i].plot(daily_actual_cases_evasive,'g',ls='-',color=color2)
                axs[2][i].set_yticks(np.arange(0,4*10**6,10**6))
                #axs[2][i].set_ylim(3*10**6,3.8*10**6)
                
                
                #axs[1][i].set_ylim([0,3.5*10**5])
                #y_ticks = axs[1][i].get_yticks()
                y_ticks = np.array([0,3,6,10])*(N/1000)
                axs[1][i].set_yticks(y_ticks)
                axs[1][i].set_yticklabels(["{:.1f}".format((x/N)*1000) for x in y_ticks])
                
                #axs[2][i].set_ylim([0,3.5*10**6])
                y_ticks = axs[2][i].get_yticks()
                axs[2][i].set_yticks(y_ticks)
                axs[2][i].set_yticklabels(["{:.1f}".format((x/N)*1000) for x in y_ticks])
        
                
def plot_Fig6_timeline_one_column(d2_,p_,c_, β):
        
        color1 = "black"
        color2 = "forestgreen"
        
        fig, axs = plt.subplots(2,1,figsize=(6,6),sharex=True,sharey='row', tight_layout=True)

        y1, βs, perc_pos1, detected1, dy_arr = simulate_with_detectability(Z, Dp,
                                                                           Dc, Ds, β,
                                                                           α1, α2, μ,
                                                                           a, 0, p_, 0.9, 0.9, c=0, STRAIN=1)
        
        S1, E1, Ip1, Ic1, Is1, P1, R1, EE1, IIp1, IIc1, IIs1, PP1, RR1 = y1.T
        daily_actual_cases_detectable = E1/Z + EE1/Z
        daily_confirmed_cases_detectable = detected1
        
        y2, βs, perc_pos2, detected2, dy_arr = simulate_with_detectability(Z, Dp, 
                                                                           Dc, Ds, β,
                                                                           α1, α2,μ,
                                                                           a, 0, p_, d2_, d2_, c_, STRAIN=2)
        S2, E2, Ip2, Ic2, Is2, P2, R2, EE2, IIp2, IIc2, IIs2, PP2, RR2= y2.T
        daily_actual_cases_evasive = E2/Z+EE2/Z
        daily_confirmed_cases_evasive = detected2
        
        #fig, axs = plt.subplots(3,1,figsize=(6,6),sharex=True,tight_layout=True)
        plt.setp(axs, xlim=(0,600))
        
  
        axs[0].set_title(f'Daily tests: {p_*10**3} per 1K people',fontsize=12)  
        axs[1].set_xlabel('Days')
        
        axs[1].plot(t,daily_confirmed_cases_detectable,ls='None',label=rf'$d_2={1}$',color=color1,marker='.',markevery=10)
        axs[1].plot(t,daily_confirmed_cases_evasive,ls='None',label=rf'$d_2={d2_}$',color=color2,marker='.',markevery=10)        

        axs[1].set_yticks(np.arange(0,8/6*10**4,10**5))
        #axs[1][i].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        
        axs[1].set_ylabel('Daily confirmed\n new cases (per 1K)')
        axs[0].set_ylabel('Daily actual\nnew cases (per 1K)')
        
        axs[1].legend(['Detectable','Test-evasive'],loc="upper right",borderpad=0.1)
        
        axs[0].plot(t,daily_actual_cases_evasive,ls='None',color=color2,marker='.',markevery=10)
        axs[0].plot(t,daily_actual_cases_detectable,ls='None',color=color1,marker='.',markevery=10)
        #axs[0].bar(t,daily_actual_cases_evasive,width=0.5,color=color2)
        #axs[0].bar(t,daily_actual_cases_detectable,color=color1)
        
        axs[0].set_yticks(np.arange(0,8/6*10**4,10**6))
        #axs[2][i].set_ylim(3*10**6,3.8*10**6)
                
        #axs[1][i].set_ylim([0,3.5*10**5])
        #y_ticks = axs[1].get_yticks()
        y_ticks = np.array([0,0.2,0.4])*(N/1000)
        axs[1].set_yticks(y_ticks)
        axs[1].set_yticklabels(["{:.1f}".format((x/N)*1000) for x in y_ticks])

        
        #axs[2][i].set_ylim([0,3.5*10**6])
        #y_ticks = axs[0].get_yticks()
        y_ticks = np.array([0,2,4])*(N/1000)
        axs[0].set_yticks(y_ticks)
        axs[0].set_yticklabels(["{:.0f}".format((x/N)*1000) for x in y_ticks])
        
        
        #PANEL LETTERS
        axs[0].text(0.05, 1.15, 'b', transform=axs[0].transAxes,
                     fontsize=16, fontweight='bold', va='top', ha='right')
        
        
        axs[1].text(0.05, 1.15, 'd', transform=axs[1].transAxes,
             fontsize=16, fontweight='bold', va='top', ha='right')
        
        #FIGURE FOR % POSITIVE TESTS
        fig,axs2 = plt.subplots(1,1,tight_layout=True)
        axs2.plot(t,perc_pos1,label=rf'$d2={1}$',color=color1)
        axs2.plot(t,perc_pos2,label=rf'$d2={d2_}$',color=color2)
        y_ticks = np.arange(0,0.05,0.01)
        axs2.set_title(f'Daily tests: {p_*10**3} per 1K people',fontsize=12)  
        axs2.set_yticks(y_ticks)
        axs2.set_xlim(0,600)
        axs2.set_yticklabels([f'{"{:.0f}".format(y*100)}%' for y in y_ticks])
        axs2.set_ylabel('Daily % positive tests')
        axs2.legend(['Detectable','Test-evasive'],loc="upper left",borderpad=0.1)


                
        
        

def plot_Fig6(p_arr,c_,β_):
        #plt.figure(figsize=(8, 6),tight_layout=True)
        plt.rcParams.update({'font.size': 14})
        
        #p_range = np.linspace(1/10**3, 0.014,20)
        d2_range = np.append(np.linspace(0, 1-10**-6,30),1)
        
        base_actual, base_confirmed, base_perc_pos, res_actual, res_confirmed, res_perc_pos = get_result_for_testing_rate(d2_range,p_=p_arr[0],c=c_, β=β_)
        base_actual2, base_confirmed2, base_perc_pos2, res_actual2, res_confirmed2,res_perc_pos2 = get_result_for_testing_rate(d2_range,p_=p_arr[1],c=c_, β=β_)
        base_actual3, base_confirmed3, base_perc_pos3, res_actual3, res_confirmed3,res_perc_pos3 = get_result_for_testing_rate(d2_range,p_=p_arr[2],c=c_, β=β_)
        
        
        fig, ax = plt.subplots(2,1,sharex=True,tight_layout=True,figsize=(6,6))
           
        marker_type = 'x'
        
        #FOR SUPP FIG, NO COST
        #marker_type = 'None'

        ax[0].plot(d2_range, base_actual,marker=marker_type, markevery=[6],markersize=8, color='k',ls='-')
        ax[0].plot(d2_range, base_actual2,marker=marker_type, markevery=[14],markersize=8,color='k',ls='--')
        ax[0].plot(d2_range, base_actual3,marker=marker_type, markevery=[21],markersize=8,color='k',ls=':')
        
        marker_type = 'None'

        ax[0].plot(d2_range, res_actual,marker=marker_type,markevery=[6],color='g',ls='-')
        ax[0].plot(d2_range, res_actual2,marker=marker_type,markevery=[14],color='g',ls='--')
        ax[0].plot(d2_range, res_actual3,marker=marker_type,markevery=[21],color='g',ls=':')
        #ax[0].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        ax[0].set_ylabel('Max daily\nactual cases (per 1K)')
        
        
        ax[0].text(0.4,3.4*10**4,f'Daily tests: {"{:.1f}".format(p_arr[0]*10**3)} per 1K people',
                   fontstyle='oblique',fontweight='book',rotation=0,fontsize=11)
        ax[0].text(0.6,3.2*10**4,f'{"{:.1f}".format(p_arr[1]*10**3)} per 1K people',fontstyle='oblique',
                   fontweight='book',rotation=0,fontsize=11)
        ax[0].text(0.1,2.65*10**4,f'{"{:.1f}".format(p_arr[2]*10**3)} per 1K people',fontstyle='oblique',
                   fontweight='book',rotation=0,fontsize=11)
        
        ax[0].set_xlim([0,0.9])
        ax[0].set_ylim(top=3.6*10**4)

        #FOR SUPP FIG, NO COST
        #ax[0].set_ylim(top=1.4*10**6)

        
        
        ax[1].set_xticks(np.arange(0,1,0.1))
        
        

        y_ticks = np.array([3,3.5,4])*(N/1000)
        ax[0].set_yticks(y_ticks)
        ax[0].set_yticklabels(["{:.1f}".format((x/N)*1000) for x in y_ticks])
        

        y_ticks = np.array([0,0.2,0.4])*(N/1000)
        ax[1].set_yticks(y_ticks)
        ax[1].set_yticklabels(["{:.1f}".format((x/N)*1000) for x in y_ticks])
        
        
        marker_type = 'None'
        ax[1].plot(d2_range, base_confirmed,marker=marker_type,color='k',ls='-')
        ax[1].plot(d2_range, base_confirmed2,marker=marker_type,color='k',ls='--')
        ax[1].plot(d2_range, base_confirmed3,marker=marker_type,color='k',ls=':')

        ax[1].plot(d2_range, res_confirmed,marker=marker_type,color='g',ls='-')
        ax[1].plot(d2_range, res_confirmed2,marker=marker_type,color='g',ls='--')
        ax[1].plot(d2_range, res_confirmed3,marker=marker_type,color='g',ls=':')
        
        #ax[1].ticklabel_format(axis="y", style="sci", scilimits=(0,0))
        ax[1].set_ylabel('Max daily\nconfirmed cases (per 1K)')
        ax[1].set_xlabel('Detectability of test-evasive strain (TPR)')

        
        #y-value is from the end of the text when rotation is negative
        ax[1].text(0.01,3.42*10**3,f'Daily tests: {"{:.1f}".format(p_arr[2]*10**3)} per 1K people',
                   fontstyle='oblique',fontweight='book',rotation=0,fontsize=11)
        ax[1].text(0.01,2.2*10**3,f'{"{:.1f}".format(p_arr[1]*10**3)} per 1K people',fontstyle='oblique',
                   fontweight='book',rotation=0,fontsize=11)
        ax[1].text(0.01,1.7*10**3,f'{"{:.1f}".format(p_arr[0]*10**3)} per 1K people',fontstyle='oblique',
                   fontweight='book',rotation=0,fontsize=11)    
        ax[1].set_ylim(bottom=0, top=0.4*10**4)
      
        
        #PANEL LETTERS
        ax[0].text(0.05, 1.15, 'a', transform=ax[0].transAxes,
                     fontsize=16, fontweight='bold', va='top', ha='right')
        
        
        ax[1].text(0.05, 1.15, 'c', transform=ax[1].transAxes,
             fontsize=16, fontweight='bold', va='top', ha='right')
        #ax[1].spines['top'].set_visible(False)

        

        #SUPPLEMENTARY FIGURE FOR % POSTIIVE TESTS
        plt.figure(figsize=(6,5),tight_layout=True)
        
        plt.plot(d2_range,base_perc_pos,color='k',ls='-',label=f'$p$ = {p_arr[0]*10**3}')
        plt.plot(d2_range,base_perc_pos2,color='k',ls='--',label=f'$p$ = {p_arr[1]*10**3}')
        plt.plot(d2_range,base_perc_pos3,color='k',ls=':',label=f'$p$ = {p_arr[2]*10**3}')
        
        plt.plot(d2_range,res_perc_pos,color='g',ls='-',label=f'$p$ = {p_arr[0]*10**3}')
        plt.plot(d2_range,res_perc_pos2,color='g',ls='--',label=f'$p$ = {p_arr[1]*10**3}')
        plt.plot(d2_range,res_perc_pos3,color='g',ls=':',label=f'$p$ = {p_arr[2]*10**3}')
        
        
        plt.text(0.01,0.0425,f'Daily tests: {"{:.1f}".format(p_arr[2]*10**3)} per 1K people',
                   fontstyle='oblique',fontweight='book',rotation=0,fontsize=11)
        plt.text(0.01,0.0385,f'{"{:.1f}".format(p_arr[1]*10**3)} per 1K people',fontstyle='oblique',
                   fontweight='book',rotation=0,fontsize=11)
        plt.text(0.01,0.0335,f'{"{:.1f}".format(p_arr[0]*10**3)} per 1K people',fontstyle='oblique',
                   fontweight='book',rotation=0,fontsize=11)    
        
        
        #plt.legend(borderpad=0.1)
        plt.ylabel('Daily % positive tests')
        #plt.xlabel(r'$\dfrac{d_2}{d_1}$ = relative detectability of test-evasive strain')
        plt.xlabel('Detectability of test-evasive strain')
        x_ticks = np.arange(0,1,0.1)
        plt.xlim([0,0.9])
        plt.ylim(bottom=0)
        plt.xticks(ticks = x_ticks)
        #locs, labels = plt.yticks()
        locs = np.arange(0,0.045,0.01)
        plt.yticks(ticks=locs,labels=[f'{"{:.1f}".format(y*100)}%' for y in locs])
        
        
        




plot_Fig6(p_arr=[4/10**3,6.3/10**3,15/10**3],c_=0.01, β_=0.42)


#plot_Fig6_timeline_one_column(d2_=0.2,p_=6.3/10**3,c_=0.01, β=0.42)


#PAPER HAS THIS:
#plot_Fig6_timeline_one_column(d2_=0.6,p_=15/10**3,c_=0.01, β=0.42)

plt.savefig('Figure_6a.jpeg',dpi=300)

