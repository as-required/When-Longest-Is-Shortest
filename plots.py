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
from numpy.polynomial.polynomial import polyval

# function for curve_fit
def quartic(x, a0, a1, a2, a3, a4): # this is an inverse quartic
    return a0 + a1/x + a2/x**2 + a3/x**3 + a4/x**4
def arbitrary_poly(x, *params):
    return sum([p*(x**i) for i, p in enumerate(params)])

class Plots:
    def __init__(self):
        pass
        
    def delta_vs_p(self,p_start = 0.01, p_end = 2, p_step = 0.1, nodes = 50, trials = 10, delta_type = "avg"\
                   ,curve = "arbitrary_poly", poly_order = 10, radius = None):
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
# =============================================================================
#                 KEEP GEODESIC OFF TO BE SAFE
# =============================================================================
                if radius is not None:
                    g_i = ga.Graph(nodes= nodes, p=p_i, geodesic=0, auto_radius= 0, radius = radius)
                else:
                    g_i = ga.Graph(nodes= nodes, p=p_i, geodesic=0, auto_radius= 1)
                
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
            popt, pcov = curve_fit(arbitrary_poly, p_vals, mean_deltas, sigma = errors,\
                                   p0 = np.ones(poly_order))
            
            # Find roots of polynomial
            #poly = [popt[0], popt[1], popt[2], popt[3], popt[4]]
            roots = np.roots(popt)
            real_roots = [i for i in roots if i.imag == 0]
            #x_ints = [i for i in real_roots if i > 0]
            print('x-intercept for {}:'.format(delta_type), np.real(real_roots))
            
            # axis values for the curve fit
            x_vals = np.linspace(p_vals[0], p_vals[-1], 50)
            y_vals = polyval(x_vals, popt) #doesn't matter what you pass in as x
            
            
            #plot
            plt.plot(x_vals, y_vals, linestyle='-', color='C30')

            plt.hlines(0,0,2.1, linestyle=':', color='k', alpha=0.5)
            plt.plot(p_vals, mean_deltas, marker='x', linestyle='', color='C30', label='Average $\delta$')
            plt.errorbar(p_vals, mean_deltas, yerr=errors, fmt='o', color='C30', capsize=5)
            
            plt.savefig('ResultPlot1_{}_delta.png'.format(delta_type), dpi=500, bbox_inches='tight')
            plt.show()
            
            #print('Mean {} Gromov-Delta for each p:'.format(delta_type), mean_deltas)
            #print('Uncertainties on mean {} Gromov-Delta for each p:'.format(delta_type), errors)
            
        else:
            delta_types_list = ["Avg", "Min", "Max"] # for the y axis labels
            fig, axs = plt.subplots(3) # axs = (ax_avg, ax_min, ax_max)
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
                popt, pcov = curve_fit(arbitrary_poly, p_vals, mean_deltas_list[m], sigma = errors_list[m]\
                                       , p0 = np.ones(poly_order))
                
                # Find roots of polynomial
                #poly = [popt[0], popt[1], popt[2], popt[3], popt[4]]
                roots = np.roots(popt)
                real_roots = [i for i in roots if i.imag == 0]
                #x_ints = [i for i in real_roots if i > 0]
                print('x-intercept for {}:'.format(delta_types_list[m]), np.real(real_roots))
                
                # axis values for the curve fit
                x_vals = np.linspace(p_vals[0], p_vals[-1], 50)
                y_vals = polyval(x_vals, popt) #doesn't matter what you pass in as x
                x_vals = np.linspace(p_vals[0], p_vals[-1], 50)
                
                
                #plot
                axs[m].plot(x_vals, y_vals, linestyle='-', color='C30')

                axs[m].hlines(0,0,2.1, linestyle=':', color='k', alpha=0.5)
                axs[m].plot(p_vals, mean_deltas_list[m], marker='x', linestyle='', color='C30',\
                         label='{} $\delta$'.format(delta_types_list[m]))
                axs[m].errorbar(p_vals, mean_deltas_list[m], yerr=errors_list[m], fmt='o', color='C30', capsize=5)
                
                
                #print('Mean {} Gromov-Delta for each p:'.format(delta_types_list[m]), mean_deltas_list[m])
                #print('Uncertainties on mean {} Gromov-Delta for each p:'.format(delta_types_list[m]), errors_list[m])
            plt.savefig('ResultPlot1_{}_delta.png'.format(delta_type), dpi=500, bbox_inches='tight')
            plt.show()
# =============================================================================
#     need to add unc in x intercept analysis
# =============================================================================
    def perp_vs_p(self,p_start = 0.01, p_end = 2, p_step = 0.1, nodes = 50, trials = 10, perp_type = "avg"\
                   ,curve = "quartic", poly_order = 10, radius = None):
# =============================================================================
#         still need to implement ability to specify curve for curve_fit
# =============================================================================
        function_mapper ={
            "quartic": quartic
            }
        
        p_vals = np.arange(p_start, p_end, p_step) # p_end is exclusive
        mean_perps = [] #stores the mean value of delta for each p_i
        errors = [] # stores the error of the mean delta for each p_i
        
# =============================================================================
#         lists for perp_type = "all"
# =============================================================================
        mean_avg_perps = []
        errors_avg_perps = []
        
        mean_min_perps = []
        errors_min_perps = []
        
        mean_max_perps = []
        errors_max_perps = []
        # IMPERATIVE TO HAVE THE LISTS IN ORDER AVG,MIN,MAX
        mean_perps_list = [mean_avg_perps, mean_min_perps, mean_max_perps]
        errors_list = [errors_avg_perps, errors_min_perps, errors_max_perps]
        
# =============================================================================
#         Collecting the data and making the graphs
# =============================================================================
        for i in range(len(p_vals)):
            p_i = p_vals[i]
            perps_i = [] # to be used for "avg", "max", "min"
            
            # to be used for "all"
            avg_perps_i = []
            min_perps_i = []
            max_perps_i = []
            perps_i_list = [avg_perps_i, min_perps_i, max_perps_i]
            
            for j in range(0, trials):
                np.random.seed(j)
# =============================================================================
#                 TURN GEODESIC OFF AS DON'T NEED FOR PERP DIST
# =============================================================================
                if radius is not None:
                    g_i = ga.Graph(nodes= nodes, p=p_i, geodesic=0, auto_radius= 0, radius = radius)
                else:
                    g_i = ga.Graph(nodes= nodes, p=p_i, geodesic=0, auto_radius= 1)
                if perp_type == "avg":
                    avg_perp_j = g_i.perp_dist()[0]
                    perps_i.append(avg_perp_j)
                elif perp_type == "min":
                    min_perp_j = g_i.perp_dist()[1]
                    perps_i.append(min_perp_j)
                elif perp_type == "max":
                    max_perp_j = g_i.perp_dist()[2]
                    perps_i.append(max_perp_j)
                
                elif perp_type == "all":
                    
                    for k in range(len(perps_i_list)):
                        perp_j = g_i.perp_dist()[k]
                        perps_i_list[k].append(perp_j)
            
# =============================================================================
#             Collecting statistics
# =============================================================================
            # mean of the perps for a given p_i
            if perp_type != "all":
                mean_perp_i = np.mean(perps_i)
                error_i = sem(perps_i)
                mean_perps.append(mean_perp_i)
                errors.append(error_i)
            else:
                for l in range(len(mean_perps_list)):
                    mean_all_perp_i = np.mean(perps_i_list[l])
                    error_all_i = sem(perps_i_list[l])
                    mean_perps_list[l].append(mean_all_perp_i)
                    errors_list[l].append(error_all_i)
# =============================================================================
#         plots
# =============================================================================
        if perp_type != "all":
            fig, ax = plt.subplots()
            fig.set_size_inches(5, 5)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.set_xlabel('$p$', fontsize=14)
            ax.set_ylabel(r'{} $\Delta_L - \Delta_S$'.format(perp_type.capitalize()), fontsize=14)
            ax.set_xticks(np.arange(p_start,p_end,0.5)) 
            ax.set_xlim(p_start - 0.1, p_end + 0.1)
            ax.set_ylim(-0.5,0.2)
            
            # curve fitting

            popt, pcov = curve_fit(arbitrary_poly, p_vals, mean_perps, sigma = errors,\
                                   p0 = np.ones(poly_order))
            
            # Find roots of polynomial
            #poly = [popt[0], popt[1], popt[2], popt[3], popt[4]]
            roots = np.roots(popt)
            real_roots = [i for i in roots if i.imag == 0]
            #x_ints = [i for i in real_roots if i > 0]
            print('x-intercept for {}:'.format(perp_type), np.real(real_roots))
            
            # axis values for the curve fit
            x_vals = np.linspace(p_vals[0], p_vals[-1], 50)
            y_vals = polyval(x_vals, popt) #doesn't matter what you pass in as x
            
            #plot
            plt.plot(x_vals, y_vals, linestyle='-', color='C30')

            plt.hlines(0,0,2.1, linestyle=':', color='k', alpha=0.5)
            plt.plot(p_vals, mean_perps, marker='x', linestyle='', color='C30', label='Average perps')
            plt.errorbar(p_vals, mean_perps, yerr=errors, fmt='o', color='C30', capsize=5)
            
            plt.savefig('ResultPlot2_{}_perp.png'.format(perp_type), dpi=500, bbox_inches='tight')
            plt.show()
            
            #print('Mean {} Perp Dist for each p:'.format(perp_type), mean_perps)
            #print('Uncertainties on mean {} Perp Dist for each p:'.format(perp_type), errors)
            
        else:
            perp_types_list = ["Avg", "Min", "Max"] # for the y axis labels
            fig, axs = plt.subplots(3) # axs = (ax_avg, ax_min, ax_max)
            fig.set_size_inches(5, 5)
            for m in range(len(mean_perps_list)):
                axs[m].spines['right'].set_visible(False)
                axs[m].spines['top'].set_visible(False)
                axs[m].set_xlabel('$p$', fontsize=14)
                axs[m].set_ylabel(r'{} $\Delta_L - \Delta_S$'.format(perp_types_list[m]), fontsize=14)
                axs[m].set_xticks(np.arange(p_start,p_end,0.5)) 
                axs[m].set_xlim(p_start - 0.1, p_end + 0.1)
                axs[m].set_ylim(-0.2,0.2)
                
                # curve fitting
                popt, pcov = curve_fit(quartic, p_vals, mean_perps_list[m], sigma = errors_list[m])
                
                # Find roots of polynomial
                #poly = [popt[0], popt[1], popt[2], popt[3], popt[4]]
                roots = np.roots(popt)
                real_roots = [i for i in roots if i.imag == 0]
                #x_ints = [i for i in real_roots if i > 0]
                print('x-intercept for {}:'.format(perp_types_list[m]), np.real(real_roots))
                
                # axis values for the curve fit
                x_vals = np.linspace(p_vals[0], p_vals[-1], 50)
                y_vals = polyval(x_vals, popt) #doesn't matter what you pass in as x
                # NEED TO PROPAGATE TO ALL
                
                
                #plot
                axs[m].plot(x_vals, y_vals, linestyle='-', color='C30')

                axs[m].hlines(0,0,2.1, linestyle=':', color='k', alpha=0.5)
                axs[m].plot(p_vals, mean_perps_list[m], marker='x', linestyle='', color='C30',\
                         label='{} $\delta$'.format(perp_types_list[m]))
                axs[m].errorbar(p_vals, mean_perps_list[m], yerr=errors_list[m], fmt='o', color='C30', capsize=5)
                
                
                #print('Mean {} perps for each p:'.format(perp_types_list[m]), mean_perps_list[m])
                #print('Uncertainties on mean {} perps for each p:'.format(perp_types_list[m]), errors_list[m])
            plt.savefig('ResultPlot2_{}_perp.png'.format(perp_type), dpi=500, bbox_inches='tight')
            plt.show()
# =============================================================================
#     need to add unc in x intercept analysis
# =============================================================================

    def length_vs_p(self,p_start = 0.01, p_end = 2, p_step = 0.1, nodes = 50, trials = 10, \
                     poly_order = 10, radius = None, path_type = "metric"):

        
        p_vals = np.arange(p_start, p_end, p_step) # p_end is exclusive
        sh_diffs = [] # Shortest difference list
        lo_diffs = [] # Longest difference list
        sh_errors = []
        lo_errors = []
        

        
        
# =============================================================================
#         Collecting the data and making the graphs
# =============================================================================
        for i in range(len(p_vals)):
            p_i = p_vals[i]
            # to be used for "avg", "max", "min"
            
            sh_diffs_i = []
            lo_diffs_i = []
            
            
            for j in range(0, trials):
                np.random.seed(j)
                if radius is not None:
                    g_i = ga.Graph(nodes= nodes, p=p_i, geodesic=1, auto_radius= 0, radius = radius)
                else:
                    g_i = ga.Graph(nodes= nodes, p=p_i, geodesic=1, auto_radius= 1)
                sh_diff_j = g_i.path_geo_diff()[0]
                sh_diffs_i.append(sh_diff_j)
                lo_diff_j = g_i.path_geo_diff()[1]
                lo_diffs_i.append(lo_diff_j)
            
# =============================================================================
#             Collecting statistics
# =============================================================================
            # mean of the lengths for a given p_i
            sh_diff_i = np.mean(sh_diffs_i)
            sh_error_i = sem(sh_diffs_i)
            sh_diffs.append(sh_diff_i)
            sh_errors.append(sh_error_i)
            
            lo_diff_i = np.mean(lo_diffs_i)
            lo_error_i = sem(lo_diffs_i)
            lo_diffs.append(lo_diff_i)
            lo_errors.append(lo_error_i)
# =============================================================================
#         plots
# =============================================================================
       
        fig, ax = plt.subplots()
        fig.set_size_inches(5, 5)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xlabel('$p$', fontsize=14)
        ax.set_ylabel(r'Difference in $L_p$ length with geodesic', fontsize=14)
        ax.set_xticks(np.arange(p_start,p_end,0.3)) 
        ax.set_xlim(p_start - 0.1, p_end + 0.1)
        ax.set_ylim(0.0,0.001)
        
        # curve fitting for longest and shortest
        # list to store the data for the shortest and longest path respectively
        path_diffs_list = [sh_diffs, lo_diffs]
        path_errors_list = [sh_errors, lo_errors]
        
        popt_list = [] #shortest, then longest
        pcov_list = []
        path_names = ["Shortest", "Longest"]
        
        x_vals_list = []
        y_vals_list = []
        
        for m in range(len(path_diffs_list)):
            popt_m, pcov_m = curve_fit(arbitrary_poly, p_vals, path_diffs_list[m], sigma = path_errors_list[m],\
                                   p0 = np.ones(poly_order))
            
            popt_list.append(popt_m)
            pcov_list.append(pcov_m)
        
            # Find roots of polynomial
        
            roots_m = np.roots(popt_m)
            real_roots_m = [i for i in roots_m if i.imag == 0]
            print('x-intercept for {}:'.format(path_names[m]), np.real(real_roots_m))
        
            # axis values for the curve fit
            x_vals_m = np.linspace(p_vals[0], p_vals[-1], 50)
            y_vals_m = polyval(x_vals_m, popt_m) #doesn't matter what you pass in as x
            
            x_vals_list.append(x_vals_m)
            y_vals_list.append(y_vals_m)
        
        #plot
            plt.plot(x_vals_m, y_vals_m, linestyle='-',\
                     label = "{} Metric Path Fit".format(path_names[m]))
            

        plt.hlines(0,0,2.1, linestyle=':', color='k', alpha=0.5)
        
        plt.errorbar(p_vals, sh_diffs, fmt='x', capsize = 5, yerr = sh_errors,\
                     label='Shortest Metric Path')
        plt.errorbar(p_vals, lo_diffs, fmt='x', capsize = 5, yerr = lo_errors,\
                     label='Longest Metric Path')
        plt.legend()
        
        plt.savefig('ResultPlot3_{}_length.png'.format(path_type), dpi=500, bbox_inches='tight')
        plt.show()
        
        """print('Mean shortest length diff for each p:', sh_diffs)
        print('Uncertainties on mean shortest length Dist for each p:', sh_errors)
        print('Mean longest length diff for each p:', lo_diffs)
        print('Uncertainties on mean longest length Dist for each p:', lo_errors)"""
            
# =============================================================================
#     need to add unc in x intercept analysis
# =============================================================================
        



