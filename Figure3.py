# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 11:33:04 2021

@author: Yael
"""


import matplotlib.pyplot as plt
import numpy as np
import time
#import seaborn as sns
#from scipy.interpolate import make_interp_spline, BSpline
from scipy.interpolate import interp1d


N = 8804190 # USA
#N=8000000
Z = 5 #incubation period, days from exposed to infected 5
#β = 0.85 #transmission rate
α1 = 0.65 #wildtype virus, fraction become clinical, (0.6)
Dp = 3 #days for preclinical to become clinical, 3 
Dc = 14 #days for clinical to recover 13
Ds = 6 #days for subclinical to recover, 6

α2 = 0.05 #attenuated virus, fraction become clinical, 0.05
μ = 0 #doesn't work here since we are examining different value of miu

#μ = 0.75 #factor by which subclinical is less contagious than preclinical (number is correct for symptomatic?)
#μ2 = 0.75 #factor by which subclinical is less contagious than preclinical (number is correct for symptomatic?)

a = 1
t = np.arange(500)

def simulate(Z, Dp, Dc, Ds, βs_range, α1, α2, μ, a, Δ, start=1):
        
    μ2 = μ  
    
    E0, Is0 = 500, 1500
    
    S0 = N - E0 - Is0
    βs = np.zeros(t.size) 
    y = np.zeros((t.size, 11))
    #.     S   E   Ip  Ic   Is   R  / E    Ip  Ic   Is   R
    y[0] = S0, E0/2, 0, 0, Is0/2, 0, E0/2, 0, 0, Is0/2, 0
      
    #y[0] =  S0, 499,  0,  0,   1500 , 0,   1, 0,  0,   0,   0
      
    
    βs[0] = βs_range.min()
    
    for t_ in t[1:]:
                     
        β_ = βs_range.min() if t_>=start else βs_range.max()
                
        βs[t_] = β_
        S, E, Ip, Ic, Is, R, EE, IIp, IIc, IIs, RR = y[t_ - 1]
       
        dy = [
            -β_ * S * (Ip + μ * Is + a*(IIp + μ2 * IIs)) / N, # S
            # wt
            β_ * S * (Ip + μ * Is) / N - E / Z,   # E
            α1 * E / Z - Ip / Dp,                 # Ip
            Ip / Dp - Ic / Dc,                    # Ic
            (1 - α1) * E / Z - Is / Ds,           # Is
            Is / Ds + Ic / Dc,                    # R
            # mut
            β_ * S * a*(IIp + μ2 * IIs) / N - EE / Z, # EE
            α2 * EE / Z - IIp / Dp,                # Ip
            IIp / Dp - IIc / Dc,                   # Ic
            (1 - α2) * EE / Z - IIs / Ds,          # Is
            IIs / Ds + IIc / Dc,                   # R      
        ]
        y[t_] = y[t_ - 1] + dy
    return y, βs


def plot_start_point_vs_beta(βs_range,delta=0,ls_='-',color_='green',miu=μ,SMOOTH=False):
        
        plt.rcParams.update({'font.size': 14})

        #starts = np.ones(βs_range.size)*t.max()
        starts = np.zeros(βs_range.size)
       
        
        #no such thing as start=0, min day is 1
        start_range = range(1,t.max(),1)
        
          
        for i,β_ in enumerate(βs_range):

                counter = 0
                for start_ in start_range:
                      
      
                        y,beta = simulate(Z, Dp, Dc, 
                                          Ds, np.linspace(β_,βs_range.max(),20),
                                          α1, α2, miu, a, 0,start=start_)
                  
                        S, E, Ip, Ic, Is, R, EE, IIp, IIc, IIs, RR = y.T/N
                        
                        I0 = (EE[0]+IIc[0]+IIs[0]+IIp[0])/(E[0]+Ic[0]+Is[0]+Ip[0])
                        
                        I_final = RR[-1] / R[-1]
                        
                        s = I_final/I0
                        
                        if s>=1:
                                counter = 0
                                starts[i]=start_
                        elif counter>20: #after 20 negatives break 
                                break
                        else:
                                counter+=1

 
        x_vals = 1-βs_range/βs_range.max()
           
        y_vals = starts
        
        plt.plot(x_vals,y_vals,ls='None', color=color_, marker='.',markersize=5,lw=2)
        
        ############SMOOTH LINES
        
        line=[]
        if SMOOTH:
        #only interpolating on the values larger than start=1 because this 
        #is the part that is quadratic
        
        #add the last element that is y==1 so the interpolation will be nice
        #the values are in reverse (sorted descending)
                relevant_ind = np.nonzero(y_vals>0)
                
                if relevant_ind[0].max()<(len(y_vals)-1):
                        relevant_ind = np.append(relevant_ind[0],relevant_ind[0].max()+1)
                
                pos_y_vals = y_vals[relevant_ind]
                pos_x_vals = x_vals[relevant_ind]
                
                
                #step_ = max(len(pos_y_vals)//5,1)
                step_ = 1
                
                xnew = np.linspace(pos_x_vals.min(), pos_x_vals.max(), 50) 
                
                #PATCH: interpolating on low sampling rate gives a smoother curve
                #and another patch to include the last point
                if len(pos_x_vals)-1 in list(range(0,len(pos_x_vals),step_)):
                        f2 = interp1d(pos_x_vals[0:len(pos_x_vals):step_], 
                                      pos_y_vals[0:len(pos_y_vals):step_],
                                      kind='cubic',fill_value='extrapolate')    
                        
                else:
                        f2 = interp1d(np.append(pos_x_vals[0:len(pos_x_vals):step_],pos_x_vals[-1]), 
                                      np.append(pos_y_vals[0:len(pos_y_vals):step_],pos_y_vals[-1]),
                                      kind='cubic',fill_value='extrapolate')
                
                line, = plt.plot(xnew,f2(xnew),color=color_,lw=2,ls=ls_)
        
        ########## END OF SMOOTH
        

        
        plt.xticks(np.arange(x_vals.min(), x_vals.max(), step=0.1))
        plt.xlabel(r'Impact of NPIs')
        plt.ylabel('Days between outbreak and NPIs')
        plt.ylim(bottom=1.5,top=70)
        #plt.xlim(left=0, right=1-βs_range.min()/βs_range.max())
        if SMOOTH==True:
                plt.xlim(left=0, right=line._x[-1])
        
        
        
        return line


#FIXED NPI

plt.figure(figsize=(8,6))
timestamp = time.time()

smooth_ = True

##NOTE: if weird errors pop out from the smooth, increase the number of data points sampled

beta_min = 0.35


miu1_ = 0.6
l1 = plot_start_point_vs_beta(np.linspace(beta_min,1.2,15),
                          delta=0,ls_='--',color_='indigo', miu=miu1_,SMOOTH=smooth_)

if smooth_:
        plt.fill_between(x=l1._x, y1=l1._y, y2=0, facecolor='indigo', 
                  edgecolor='darkred',alpha=0.2)


miu2_ = 0.65
l2 = plot_start_point_vs_beta(np.linspace(beta_min,1.2,15),
                          delta=0,ls_=':',color_='indigo', miu=miu2_,SMOOTH=smooth_)

if smooth_:
        plt.fill_between(x=l2._x, y1=l2._y, y2=0, facecolor='indigo', 
                         edgecolor='forestgreen',alpha=0.2)


miu3_ = 0.7
l3 = plot_start_point_vs_beta(np.linspace(beta_min,1.2,10),
                          delta=0,ls_='-.',color_='indigo', miu=miu3_,SMOOTH=smooth_)

if smooth_:
        plt.fill_between(x=l3._x, y1=l3._y, y2=0, facecolor='indigo', 
                         edgecolor='midnightblue',alpha=0.2)
        
        
        
        


plt.xlim(0,1-beta_min/1.2)
plt.legend((l1,l2,l3),(rf'$\mu={miu1_}$',rf'$\mu={miu2_}$',rf'$\mu={miu3_}$'),borderpad=0.1,loc='upper left')

elapsed = time.time() - timestamp


#Figure 3 MAIN TEXT
plt.text(0.3,62,'Virulent strain evolves',rotation=0,fontsize=12,fontstyle='oblique',fontweight='semibold')
plt.text(0.12,20,'Attenuated strain\nevolves',rotation=0,fontsize=12,fontstyle='oblique',fontweight='semibold')
###


# #Figure S2 SUPPLEMENTARY
# plt.text(0.2,50,'Virulent strain evolves',rotation=0,fontsize=12,fontstyle='oblique',fontweight='semibold')
# plt.text(0.495,15,'Attenuated strain\nevolves',rotation=0,fontsize=12,fontstyle='oblique',fontweight='semibold')
# ###




print(f"time elapsed: {elapsed}")

plt.savefig('Figure_3.jpeg',dpi=300)
