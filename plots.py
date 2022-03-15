#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 01:08:54 2022

@author: ali and autumn
"""
import Graph as ga
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import sem # standard error of the mean: s/sqrt(n)

# function for curve_fit
def quartic(x, a0, a1, a2, a3, a4): # this is an inverse quartic
    return a0 + a1/x + a2/x**2 + a3/x**3 + a4/x**4

class Plots:
    def __init__(self):
        pass
        
    def delta_vs_p(self,p_start = 0.01, p_end = 2, p_step = 0.1, nodes = 50, trials = 10, delta_type = "avg"\
                   ,curve = "quartic"):
# =============================================================================
#         still need to implement ability to specify curve for curve_fit
# =============================================================================
        function_mapper ={
            "quartic": quartic
            }
        
        p_vals = np.arange(p_start, p_end, p_step) # p_end is exclusive
        mean_deltas = [] #stores the mean value of delta for each p_i
        errors = [] # stores the error of the mean delta for each p_i
        
# =============================================================================
#         lists for delta_type = "all"
# =============================================================================
        mean_avg_deltas = []
        errors_avg_deltas = []
        
        mean_min_deltas = []
        errors_min_deltas = []
        
        mean_max_deltas = []
        errors_max_deltas = []
        # IMPERATIVE TO HAVE THE LISTS IN ORDER AVG,MIN,MAX
        mean_deltas_list = [mean_avg_deltas, mean_min_deltas, mean_max_deltas]
        errors_list = [errors_avg_deltas, errors_min_deltas, errors_max_deltas]
        
        for i in range(len(p_vals)):
            p_i = p_vals[i]
            deltas_i = [] # to be used for "avg", "max", "min"
            
            # to be used for "all"
            avg_deltas_i = []
            min_deltas_i = []
            max_deltas_i = []
            deltas_i_list = [avg_deltas_i, min_deltas_i, max_deltas_i]
            
            for j in range(0, trials):
                np.random.seed(j)
                g_i = ga.Graph(nodes= nodes, p=p_i, geodesic=1)
                
                if delta_type == "avg":
                    avg_delta_j = g_i.triangle()[0]
                    deltas_i.append(avg_delta_j)
                elif delta_type == "min":
                    min_delta_j = g_i.triangle()[1]
                    deltas_i.append(min_delta_j)
                elif delta_type == "max":
                    max_delta_j = g_i.triangle()[2]
                    deltas_i.append(max_delta_j)
                
                elif delta_type == "all":
                    
                    for k in range(len(deltas_i_list)):
                        delta_j = g_i.triangle()[k]
                        deltas_i_list[k].append(delta_j)
            
            # mean of the deltas for a given p_i
            if delta_type != "all":
                mean_delta_i = np.mean(deltas_i)
                error_i = sem(deltas_i)
                mean_deltas.append(mean_delta_i)
                errors.append(error_i)
            else:
                for l in range(len(mean_deltas_list)):
                    mean_all_delta_i = np.mean(deltas_i_list[l])
                    error_all_i = sem(deltas_i_list[l])
                    mean_deltas_list[l].append(mean_all_delta_i)
                    errors_list[l].append(error_all_i)
# =============================================================================
#         plots
# =============================================================================
        if delta_type != "all":
            fig, ax = plt.subplots()
            fig.set_size_inches(5, 5)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.set_xlabel('$p$', fontsize=14)
            ax.set_ylabel('{} $\delta$'.format(delta_type.capitalize()), fontsize=14)
            ax.set_xticks(np.arange(p_start,p_end,0.5)) 
            ax.set_xlim(p_start - 0.1, p_end + 0.1)
            ax.set_ylim(-0.5,0.2)
            
            # curve fitting
            popt, pcov = curve_fit(quartic, p_vals, mean_deltas, sigma = errors)
			
			# Find roots of polynomial
            poly = [popt[0], popt[1], popt[2], popt[3], popt[4]]
            roots = np.roots(poly)
            real_roots = [i for i in roots if i.imag == 0]
            x_int = min([i for i in real_roots if i > 0])
            print('x-intercept:', np.real(x_int))
            
            # axis values for the curve fit
            x_vals = np.linspace(p_vals[0], p_vals[-1], 50)
            y_vals = quartic(x_vals, popt[0], popt[1], popt[2], popt[3], popt[4])
            
            #plot
            plt.plot(x_vals, y_vals, linestyle='-', color='C30')

            plt.hlines(0,0.3,2.1, linestyle=':', color='k', alpha=0.5)
            plt.plot(p_vals, mean_deltas, marker='x', linestyle='', color='C30', label='Average $\delta$')
            plt.errorbar(p_vals, mean_deltas, yerr=errors, fmt='o', color='C30', capsize=5)
            
            plt.savefig('ResultPlot1_{}_delta.png'.format(delta_type), dpi=500, bbox_inches='tight')
            plt.show()
            
            print('Mean {} Gromov-Delta for each p:'.format(delta_type), mean_deltas)
            print('Uncertainties on mean {} Gromov-Delta for each p:'.format(delta_type), errors)
			
        else:
            delta_types_list = ["Avg", "Min", "Max"] # for the y axis labels
            fig, axs = plt.subplots(1,3) # axs = (ax_avg, ax_min, ax_max)
            fig.set_size_inches(5, 5)
            for m in range(len(mean_deltas_list)):
                axs[m].spines['right'].set_visible(False)
                axs[m].spines['top'].set_visible(False)
                axs[m].set_xlabel('$p$', fontsize=14)
                axs[m].set_ylabel('{} $\delta$'.format(delta_types_list[m]), fontsize=14)
                axs[m].set_xticks(np.arange(p_start,p_end,0.5)) 
                axs[m].set_xlim(p_start - 0.1, p_end + 0.1)
                axs[m].set_ylim(-0.2,0.2)
                
                # curve fitting
                popt, pcov = curve_fit(quartic, p_vals, mean_deltas_list[m], sigma = errors_list[m])
				
				# Find roots of polynomial
           		poly = [popt[0], popt[1], popt[2], popt[3], popt[4]]
            	roots = np.roots(poly)
            	real_roots = [i for i in roots if i.imag == 0]
            	x_int = min([i for i in real_roots if i > 0])
            	print('x-intercept:', np.real(x_int))
                
                # axis values for the curve fit
                x_vals = np.linspace(p_vals[0], p_vals[-1], 50)
                y_vals = quartic(x_vals, popt[0], popt[1], popt[2], popt[3], popt[4])
                
                #plot
                plt.plot(x_vals, y_vals, linestyle='-', color='C30')

                plt.hlines(0,0.3,2.1, linestyle=':', color='k', alpha=0.5)
                plt.plot(p_vals, mean_deltas_list[m], marker='x', linestyle='', color='C30',\
                         label='{} $\delta$'.format(delta_types_list[m]))
                plt.errorbar(p_vals, mean_deltas_list[m], yerr=errors_list[m], fmt='o', color='C30', capsize=5)
                
                
                print('Mean {} Gromov-Delta for each p:'.format(delta_types_list[m]), mean_deltas_list[m])
                print('Uncertainties on mean {} Gromov-Delta for each p:'.format(delta_types_list[m]), errors_list[m])
            plt.savefig('ResultPlot1_{}_delta.png'.format(delta_type), dpi=500, bbox_inches='tight')
            plt.show()
# =============================================================================
#     need to add unc in x intercept analysis
# =============================================================================
        
        
