# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 11:29:50 2021

@author: Yael
"""
import matplotlib.pyplot as plt
import numpy as np
import time




#PARAMETERS
#define parameters
N = 8804190 # USA
#N=8000000
Z = 5 #xin 2021,alene 2021,rai 2021 #latent period, days from exposed to infected (=contagious),
#incubation period is around 5.5-6.5, so latent is estimated by the incubation period - 2 day (zhao)
α1 = 0.65 #sah 2021 #wildtype virus, fraction become isolated (35% truly asymptomatic)
Dp = 3 #days for pre-isolated to become isolated (incubation-latent)
Dc = 14 #days for isolated to leave isolation (14) 
Ds = 6 #days for non-isolated to recover 

α2 = 0.05 #attenuated virus, fraction become isolated, 0.05
#μ = 0.75 #factor by which subclinical is less contagious than preclinical (number is correct for symptomatic?)
#μ2 = 0.75 #factor by which subclinical is less contagious than preclinical (number is correct for symptomatic?)
#β = 0.85 #transmission rate

a = 1
t = np.arange(750)



def simulate(Z, Dp, Dc, Ds, βs_range, α1, α2, μ, a, Δ, start=1):
        
    μ2 = μ  
    
    E0, Is0 = 500, 1500
    
    
    S0 = N - E0 - Is0
    βs = np.zeros(t.size) 
    y = np.zeros((t.size, 11))
    #.     S   E   Ip  Ic   Is   R  / E    Ip  Ic   Is   R
    #y[0] = S0, E0/2, 0, 0, Is0/2, 0, E0/2, 0, 0, Is0/2, 0
    
    #from rarity
    y[0] =  S0, 499,  0,  0,   1500 , 0,   1, 0,  0,   0,   0

    
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


def plot_heatmap(βs,μ_range,s):
#heatmap plot

             
        plt.figure(figsize=(8, 6),tight_layout=True)
        plt.rcParams.update({'font.size': 16})
        
        im = plt.pcolormesh(1-βs/βs.max(), μ_range, s, vmin=0, vmax=2, cmap='PuOr', rasterized=False,shading='auto')   
        cs = plt.contour(1-βs/βs.max(), μ_range, s, vmin=0, vmax=2, colors='k', linewidths=0.6, levels=[1])
        
        plt.xlabel(r'Impact of NPIs')
        
        plt.xticks(np.arange(0,0.8,0.1))
        #plt.ylabel(r'$\mu = \beta_{subclinical}\; /\; \beta_{preclinical}$')
        plt.yticks(np.arange(μ_range.min(),μ_range.max()+10**-6,step=0.1))
        
        #SUPP FIGURE
        #im = plt.pcolormesh(1-βs/βs.max(), μ_range, s, vmin=0.5, vmax=1.5, cmap='PuOr', rasterized=False,shading='auto')   
        #cs = plt.contour(1-βs/βs.max(), μ_range, s, vmin=0.5, vmax=1.5, colors='k', linewidths=0.6, levels=[1])
        # im = plt.pcolormesh(βs, μ_range, s, vmin=0, vmax=2, cmap='PuOr', rasterized=False,shading='auto')   
        # cs = plt.contour(βs, μ_range, s, vmin=0, vmax=2, colors='k', linewidths=0.6, levels=[1])
        
        # plt.xticks(np.arange(0.4,0.9,0.1))
        
        
        ####
        
        
        plt.colorbar(im, label='Selection coefficient for attenuated strain',ticks=np.array([0,0.5,1,1.5,2])*1)

       

        #FOR SUPP FIGURE:
                                
        # plt.ylabel(r'$\alpha_2$'+' = fraction of symptomatic\nin attenuated strain')
        # plt.text(0.52, 0.3, 'Attenuated\nstrain evolves', 
        #           fontdict=dict(color='w',fontstyle='oblique',fontweight='semibold',fontsize=11))
        # plt.text(0.05, 0.3, 'Virulent\nstrain evolves', 
        #           fontdict=dict(color='k',fontstyle='oblique',fontweight='semibold',fontsize=11))                
        # #####


        plt.ylabel('$\mu$ = Relative infectiousness of non-isolated')
        plt.text(0.43, 0.8, 'Attenuated\nstrain evolves', 
                  fontdict=dict(color='w',fontstyle='oblique',fontweight='semibold'))
        plt.text(0.15, 0.5, 'Virulent\nstrain evolves', 
                  fontdict=dict(color='k',fontstyle='oblique',fontweight='semibold'))
        
        return cs

def get_result_for_heatmap_supp(x_range,y_range,start_):
#tested ranges of parameters

        βs_range = np.linspace(0.35, 1.2,20)

        results = [
                    [
                    simulate(Z, Dp, Dc, Ds,
                             βs_range,
                             α1, α2=y, μ=x, a=1, Δ=0, start=start_)[0]
                       
                    for x in x_range
        #             for α2_ in α2s
        #            for a_ in aa #a in detectability is p (testing rate)
                    ]
                #for a_ in aa
                #for α2_ in α2s
                for y in y_range
        ]
        
        # for each _outer_loop_parameter and _inner_loop_parameter, get the last day value of RR/R-1 == (RR-R)/R
        # which is a proxy for which type of virus dominates the population
         
        
        
        y0 = results[0][0][0]
        S, E, Ip, Ic, Is, R, EE, IIp, IIc, IIs, RR = y0
                        
        I0 = (EE+IIc+IIs+IIp)/(E+Ic+Is+Ip)
        
        #s = np.array([[ y[-1, 10] / y[-1, 5] - 1 for y in ys] for ys in results])
        s = np.array([[ (y[-1, 10] / y[-1, 5]) / I0 for y in ys] for ys in results])
        
        
        
        print(s.min(), s.max())
        s[s>2] = 2
        #s[s<-1] = -1
        
        
        return s


def get_result_for_heatmap(βs_range,y_range,start_):
#tested ranges of parameters
     
        results = [
                    [
                    simulate(Z, Dp, Dc, Ds,
                             np.linspace(β_,βs_range.max(),20),
                             α1, α2, μ=x, a=1, Δ=0, start=start_)[0]
                      
                    for β_ in βs_range
        #             for α2_ in α2s
        #            for a_ in aa #a in detectability is p (testing rate)
                    ]
                #for a_ in aa
                #for α2_ in α2s
                for x in y_range
        ]
        
        # for each _outer_loop_parameter and _inner_loop_parameter, get the last day value of RR/R-1 == (RR-R)/R
        # which is a proxy for which type of virus dominates the population
        
        
        y0 = results[0][0][0]
        S, E, Ip, Ic, Is, R, EE, IIp, IIc, IIs, RR = y0
        
        
        #initial ratio of infected between attenuated strain and virulent strain
        I0 = (EE+IIc+IIs+IIp)/(E+Ic+Is+Ip)
        
        s = np.array([[ (y[-1, 10] / y[-1, 5]) / I0 for y in ys] for ys in results])
        print(s.min(), s.max())
        s[s>2] = 2
        
        return s




timestamp = time.time()

βs_range = np.linspace(0.37, 1.2,50)
y_range = np.linspace(0.4,0.9,50)
res = get_result_for_heatmap(βs_range,y_range,1)

#SUPPLEMENTARY FIG

#actually mu range
# βs_range = np.linspace(0.4,0.9,20)
# y_range = np.linspace(0.05,0.59,20)
# res = get_result_for_heatmap_supp(βs_range,y_range,1)
##


plot_heatmap(βs_range,y_range,res)
elapsed = time.time() - timestamp

print(f"time elapsed: {elapsed}")

#plt.savefig('Figure_2.jpeg',dpi=300)
