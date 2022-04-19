# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 11:37:46 2021

@author: Yael
"""
import matplotlib.pyplot as plt
import numpy as np
import time
#from matplotlib.lines import Line2D

#PARAMETERS
#define parameters
# N = 10000
#N = 10607700 # Wuhan
N = 8804190 # USA
#N=8000000
Z = 5 #incubation period, days from exposed to infected 5
#β = 0.85 #transmission rate
α1 = 0.65 #wildtype virus, fraction become clinical, (0.6)
Dp = 3 #days for preclinical to become clinical, 3 
Dc = 14 #days for clinical to recover 13
Ds = 6 #days for subclinical to recover, 6

α2 = 0.65 #attenuated virus, fraction become clinical, 0.05


μ = 0.6 #factor by which subclinical is less contagious than preclinical (number is correct for symptomatic?)
μ2 = 0.6 #factor by which subclinical is less contagious than preclinical (number is correct for symptomatic?)

d1 = 0.9 #detectability of wt, this is true positive rate. 

a = 1
#Δ = 6 #lag current *reported* number of cases and actual infectious individuals
t = np.arange(500) 


def simulate_with_detectability(Z, Dp, Dc, Ds, β, α1, α2, μ, a, Δ, p, d2, c):
    E0, Is0 = 500, 1500
    S0 = N - E0 - Is0
    βs = np.zeros(t.size) 
    y = np.zeros((t.size, 13))
    #perc_pos = np.zeros(t.size)
    #.     S   E   Ip  Ic   Is   P R  / E    Ip  Ic   Is  P R
    y[0] = S0, E0/2, 0, 0, Is0/2, 0,0, E0/2, 0, 0, Is0/2, 0,0
    βs[0] = β
    
    #d1  = 0.8 #detectability of wt, this is true positive rate. 
    #d2 = d1*α2/α1 #detectability of mt; d2<d1
    
    #d2 = d1*c #detectability of mt; d2<d1
    
    a=1-c
    #μ2 = μ*(1-c)
    #μ2 = μ
    #α2 = α1*(1-c)
    
    #p = 1.12*10**6/N #fraction of population that gets tested every time unit (testing rate per day, USA)
    
    for t_ in t[1:]:
        S, E, Ip, Ic, Is, P, R, EE, IIp, IIc, IIs, PP, RR = y[t_ - 1]
        
        
        #people_to_test = p*(N-R-RR)-(Ic+IIc)

        #otherwise undefined behavior
        #print(f'R+RR: {R+RR}, Ic+IIc: {Ic+IIc}, p*(N-R-RR): {p*(N-R-RR)}')
        #assert (people_to_test>=0)
        
        #p_eff = people_to_test/(S+E+EE+Is+IIs+Ip+IIp+R+RR)
        #p_eff = people_to_test/(N-R-RR-Ic-IIc)

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
       
        
    return y, βs

def plot_heatmap(x_range,y_range,s):
#heatmap plot

             
        plt.figure(figsize=(8, 6),tight_layout=True)
        plt.rcParams.update({'font.size': 14})
        im = plt.pcolormesh(1-x_range/x_range.max(), y_range/d1, s, 
                            vmin=0, vmax=2, cmap='PuOr', rasterized=False, shading='auto')
        
        #im = plt.pcolormesh(1-βs/βs.max(), aa , s, vmin=-0.3, vmax=0.3, cmap='RdBu', rasterized=False)
        
        #im = plt.pcolormesh(1-βs/βs.max(), α2s/α1, s, vmin=-0.4, vmax=0.4, cmap='RdBu')
        #im = plt.pcolormesh(aa, μ_range, s, vmin=s.min(), vmax=s.max(), cmap='RdBu')
        # plt.contour(α2s/α1, μ_range, s, vmin=-0.4, vmax=0.4, colors='k', linewidths=0.4)
        #plt.contour(aa, μ_range, s, vmin=-0.4, vmax=0.4, colors='k', linewidths=0.4)
        
        #QUESTION: are there other ways to represent NPI? perhaps the "a" parameter?
        contour = plt.contour(1-x_range/x_range.max(), y_range/d1, s, vmin=0, vmax=2, colors='k', linewidths=0.6, levels=[1])
        #plt.contour(aa, μ_range, s, vmin=-0.3, vmax=0.3, colors='k', linewidths=0.6, levels=[ 0])
        
        # plt.xlabel('Clinical rate of mutant / wildtype')
        # plt.xlabel('Effect of percieved risk on transmission rate')
        #plt.ylabel(r'testing rate (fraction of population tested per day)')
        #plt.xlabel(r'testing rate')
        
        #plt.xticks([0,0.2,0.4,0.6])
        #plt.ylabel(r'$\beta_{subclinical}\; /\; \beta_{preclinical}$')
        plt.xlabel(r'Impact of NPI')
        plt.ylim(0,y_range.max())

        plt.ylabel(r'Detectability of test-evasive strain (TPR)')
        #plt.ylabel(r'Dp/Ds')
        #plt.ylabel(r'a2/a1')
       
        #plt.title(rf'cost is in attenuation, not infectivity, cost of detectability = 0.9')
        plt.colorbar(im, label='Relative fitness of test-evasive strain',ticks=np.array([0,0.5,1,1.5,2])*1)
#plt.text(0.13, 0.6, 'attenuated\nmutant evolves', fontdict={'color':'w'})
        #plt.text(0.02, 0.9, 'attenuated mutant\ndoesn\'t evolve', fontdict={'color':'k'})
        #plt.savefig('../figures/mutator_fitness_α2={}_Ds={}.pdf'.format(α2, Ds))

        return contour



def get_result_for_heatmap(x_range,y_range,c, p_):
        
        #c=0
#tested ranges of parameters
        results = [
                    [
                                  
                     simulate_with_detectability(Z, Dp, Dc, Ds, β_, α1, α2,
                                                 μ, a, 0, p_, d2_, c)[0]
                                    
                     for β_ in x_range

                    ]
         
                for d2_ in y_range
        ]
        
        # for each _outer_loop_parameter and _inner_loop_parameter, get the last day value of RR/R-1 == (RR-R)/R
        # which is a proxy for which type of virus dominates the population
        s = np.array([[ y[-1, 12] / y[-1, 6] for y in ys] for ys in results])
        
        
        print(s.min(),s.max())
        #s[s>1] = 1
        #s[s<-1] = -1
        
        s[s>2] = 2

        
        return s


#helper function to get the contour line coordinates directly
def getContourLine(cs,plot=False):      
        p = cs.collections[0].get_paths()[0]
        v = p.vertices
        x = v[:,0]
        y = v[:,1]     
        if plot:
                plt.plot(x,y,marker='x',linestyle='None')
        return x,y



def plotContourPlot(contours,vals):
        
        plt.figure(figsize=(8, 6),tight_layout=True)
        plt.rcParams.update({'font.size': 14})

        linestyles = ['-','--',':']
        for i,line in enumerate(contours):
                x,y = getContourLine(line)
                      
                plt.plot(x,y,ls = linestyles[i],color = 'black',lw = 3,
                         label = r'p'+f'={"{:.2f}".format(vals[i]*10**3)}')
                
                #aesthetic patch to make sure the colored fill runs to the end 
                #of the plot on the right side
                x_extend = np.arange(round(x[-1],2),0.8,0.05)
                x_fill = np.append(x,x_extend)
                y_fill = np.append(y,np.ones((1,len(x_extend))))
                plt.fill_between(x_fill, y_fill, y2=0, alpha=0.2, color='forestgreen')
        

        plt.text(0.5,0.13,'Low testing rate',rotation=30,fontstyle='oblique')
        plt.text(0.35,0.38,'Medium testing rate',rotation=15,fontstyle='oblique')
        plt.text(0.2,0.6,'High testing rate',rotation=9,fontstyle='oblique')    
        plt.text(0.05,0.1,'Test-evasive strain evolves',rotation=0,fontstyle='oblique',fontweight='semibold')
        plt.text(0.2,0.8,'Detectable strain evolves',rotation=0,fontstyle='oblique',fontweight='semibold')
            
        
        #FOR SUPP FIGURE, COST=0.02
        # plt.text(0.62,0.04,'UK: 6.3',rotation=38,fontstyle='oblique')
        # plt.text(0.2,0.16,'Israel: 10.9 daily tests per 1K people',rotation=20,fontstyle='oblique')    
        # plt.text(0.3,0.1,'Test-evasive strain evolves',rotation=0,fontstyle='oblique',fontweight='semibold')
        # plt.text(0.2,0.5,'Detectable strain evolves',rotation=0,fontstyle='oblique',fontweight='semibold')
        
        
        
        plt.xlabel(r'Impact of NPIs')
        #plt.ylabel(r'$\dfrac{d_2}{d_1}$ = relative detectability of test-evasive strain')
        plt.ylabel('Detectability of test-evasive strain (TPR)')
        
        #data from here Dec 30: https://ourworldindata.org/coronavirus-testing#how-many-tests-are-performed-each-day
        
        plt.xlim(0,1-βs_range.min()/βs_range.max())
        plt.ylim(0,d2_range.max())
        
        
        #MARK FOR THE EXAMPLE IN THE FIGURE CAPTION
        plt.plot(0.5,0.7, marker='X', markersize=10,color='k')
        
        

def runHeatmaps(p_arr,c_):
        contours = []
               
        cost = c_
        res = get_result_for_heatmap(βs_range,d2_range,c=cost, p_=p_arr[0])
        cs = plot_heatmap(βs_range,d2_range,res)
        contours.append(cs)
        
        res = get_result_for_heatmap(βs_range,d2_range,c=cost, p_=p_arr[1])
        cs = plot_heatmap(βs_range,d2_range,res)
        contours.append(cs)
        
        res = get_result_for_heatmap(βs_range,d2_range,c=cost, p_=p_arr[2])
        cs = plot_heatmap(βs_range,d2_range,res)
        contours.append(cs)
        
        return contours


plt.close('all')
timestamp = time.time()

βs_range = np.linspace(0.35, 1.2,20)
d2_range = np.linspace(0, d1,20)


#FIGURE 5: PLOT CONTOUR PL5OT
p_arr = np.array([4/10**3,6.3/10**3,15/10**3]) 


#FOR SUPPLEMENTARY
#p_arr = np.array([15/10**3,20/10**3,30/10**3]) 


c_=0.01
contours = runHeatmaps(p_arr, c_)


plotContourPlot(contours, p_arr)

plt.savefig('Figure_5.jpeg',dpi=300)

