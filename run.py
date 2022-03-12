#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 01:56:36 2022

@author: ali and autumn
"""

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import Graph as ga

#%% iniialising graphs
g1 = ga.Graph(nodes = 5, p = 1, geodesic = 1, nx_plot = 1, show_weights= 1)
ghalf = ga.Graph(nodes = 5, p = 0.5, geodesic = 1,nx_plot = 1, show_weights= 1)
g2 = ga.Graph(nodes = 5, p =2, geodesic = 1,  nx_plot = 1, show_weights= 1)

graphs = [[ghalf, 0.5], [g1, 1.0], [g2, 2.0]]

#%% longest\shortest paths
for i in range(len(graphs)):
    
    print("The longest network path for p = {} is".format(graphs[i][1]), graphs[i][0].longest_path(), \
      "with length", graphs[i][0].longest_path_length())

    print("The shortest network path for p = {}  is".format(graphs[i][1]), graphs[i][0].shortest_path(), \
      "with length", graphs[i][0].shortest_path_length())
    print("The longest metric path for p = {} is".format(graphs[i][1]), graphs[i][0].longest_weighted_path(), \
      "with length", graphs[i][0].longest_weighted_path_length())
    print("The shortest metric path for p = {}  is".format(graphs[i][1]), graphs[i][0].shortest_weighted_path(), \
      "with length", graphs[i][0].shortest_weighted_path_length())

# =============================================================================
# - the longest metric path for p = 1/2 is the geodesic from (0,0) to (1,1) 
#   w/ length 4
# - the geodesic is both the shortest and longest metric path for p = 1
# =============================================================================

#%%
ghalf.triangle()
g1.triangle()
g2.triangle()

# We found that:
# 1. delta = 0 for p = 1
# 2. delta > 0 for p = 2
# 3. delta < 0 for p = 1/2
# This agrees with what we expected.

#%%
ghalf.perp_dist()
g1.perp_dist()
g2.perp_dist()

#%%
ghalf.path_geo_diff()
g1.path_geo_diff()
g2.path_geo_diff()


#%% Result Plot 1: Average Gromov-delta vs p

p_vals = []
avg_deltas = []
errors = []

# Loop to iterate over p values
for i in range(1,24):
    
    p_i = (i+1)/12
    p_vals.append(p_i)
    avg_deltas_i = []
    
    # Loop to iterate over random seeds
    for j in range(0,10):
        np.random.seed(j)
        g_i = ga.Graph(nodes=50, p=p_i, geodesic=1)
        avg_delta_j = g_i.triangle()[0]
        avg_deltas_i.append(avg_delta_j)
    
    avg_delta_i = np.mean(avg_deltas_i)
    error_i = np.std(avg_deltas_i)
    avg_deltas.append(avg_delta_i)
    errors.append(error_i)

print('Mean Gromov-Delta for each p:', avg_deltas)
print('Uncertainties on Gromov-Delta for each p:', errors)
    
fig, ax = plt.subplots()
fig.set_size_inches(5, 5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('$p$', fontsize=14)
ax.set_ylabel('Average $\delta$', fontsize=14)
ax.set_xticks(np.arange(0.5,2.01,0.5)) 
ax.set_xlim(0.45,2.1)
ax.set_ylim(-0.15,0.1)

# Curve Fitting

def test(x, a0, a1, a2, a3, a4):
    return a0 + a1/x + a2/x**2 + a3/x**3 + a4/x**4

popt, pcov = curve_fit(test, p_vals, avg_deltas)

x_vals = np.linspace(p_vals[0], p_vals[-1], 50)
y_vals = popt[0] + popt[1]/x_vals + popt[2]/x_vals**2 + popt[3]/x_vals**3 + popt[4]/x_vals**4
plt.plot(x_vals, y_vals, linestyle='-', color='C30')

plt.hlines(0,0.3,2.1, linestyle=':', color='k', alpha=0.5)
plt.plot(p_vals, avg_deltas, marker='o', linestyle='', color='C30', label='Average $\delta$')
plt.errorbar(p_vals, avg_deltas, yerr=errors, fmt='o', color='C30', capsize=5)
#plt.savefig('ResultPlot1.png', dpi=500, bbox_inches='tight')
plt.show()

# Print coefficients
print('a0=', popt[0], 'a1=', popt[1], 'a2=', popt[2], 'a3=', popt[3], 'a4=', popt[4])
# a0= 0.08994396818712974, a1= -0.07027496127410193, a2= -0.028936151547050006
# a3= 0.013985374271716092, a4= -0.003556886743392925

#%%
# Solutions to quartic equation a0*x**4 + a1*x**3 + a2*x**2 + a3*x + a4 = 0 are:
# x1 = -0.55101952015491
# x2 = 0.17188041295944 -0.20749543234584i
# x3 = 0.17188041295944 +0.20749543234584i
# x4 = 0.98857802771145
# These were found using an online quartic equation solver (https://keisan.casio.com/exec/system/1181809416)
# Hence, the x-intercept we are interested in is x4 = 0.98857802771145

x_intercept = 0.98857802771145
print('x-intercept:', x_intercept)

#%% Result Plot 1b: Minimum Gromov-Delta vs p

p_vals = []
min_deltas = []
errors = []

# Loop to iterate over p values
for i in range(1,24):
    
    p_i = (i+3)/12
    p_vals.append(p_i)
    min_deltas_i = []
    
    # Loop to iterate over random seeds
    for j in range(0,4):
        np.random.seed(j)
        g_i = ga.Graph(nodes=50, p=p_i, geodesic=1)
        min_delta_j = g_i.triangle()[1]
        min_deltas_i.append(min_delta_j)
    
    min_delta_i = np.mean(min_deltas_i)
    error_i = np.std(min_deltas_i)
    min_deltas.append(min_delta_i)
    errors.append(error_i)

print('Mean Gromov-Delta for each p:', min_deltas)
print('Uncertainties on Gromov-Delta for each p:', errors)
    
fig, ax = plt.subplots()
fig.set_size_inches(5, 5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('$p$', fontsize=14)
ax.set_ylabel('Minimum $\delta$', fontsize=14)
ax.set_xticks(np.arange(0.5,2.01,0.5)) 
ax.set_xlim(0.2,2.2)

# Curve Fitting

def test(x, a0, a1, a2, a3, a4):
    return a0 + a1/x + a2/x**2 + a3/x**3 + a4/x**4

popt, pcov = curve_fit(test, p_vals, min_deltas)

x_vals = np.linspace(p_vals[0], p_vals[-1], 50)
y_vals = popt[0] + popt[1]/x_vals + popt[2]/x_vals**2 + popt[3]/x_vals**3 + popt[4]/x_vals**4
plt.plot(x_vals, y_vals, linestyle='-', color='C2')

plt.hlines(0,0.3,2.1, linestyle=':', color='k', alpha=0.5)
plt.plot(p_vals, min_deltas, marker='o', linestyle='', color='C2', label='Minimum $\delta$')
plt.errorbar(p_vals, min_deltas, yerr=errors, fmt='o', color='C2', capsize=5)
plt.savefig('ResultPlot1b.png', dpi=500, bbox_inches='tight')
plt.show()

#%% Result Plot 1c: Maximum Gromov-Delta vs p

p_vals = []
max_deltas = []
errors = []

# Loop to iterate over p values
for i in range(1,24):
    
    p_i = (i+5)/12
    p_vals.append(p_i)
    max_deltas_i = []
    
    # Loop to iterate over random seeds
    for j in range(0,10):
        np.random.seed(j)
        g_i = ga.Graph(nodes=50, p=p_i, geodesic=1)
        max_delta_j = g_i.triangle()[2]
        max_deltas_i.append(max_delta_j)
    
    max_delta_i = np.mean(max_deltas_i)
    error_i = np.std(max_deltas_i)
    max_deltas.append(max_delta_i)
    errors.append(error_i)

print('Mean Gromov-Delta for each p:', avg_deltas)
print('Uncertainties on Gromov-Delta for each p:', errors)
    
fig, ax = plt.subplots()
fig.set_size_inches(5, 5)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('$p$', fontsize=14)
ax.set_ylabel('Maximum $\delta$', fontsize=14)
ax.set_xticks(np.arange(0.5,2.01,0.5)) 
ax.set_xlim(0.45,2.2)

# Curve Fitting

def test(x, a0, a1, a2, a3, a4):
    return a0 + a1/x + a2/x**2 + a3/x**3 + a4/x**4

popt, pcov = curve_fit(test, p_vals, max_deltas)

x_vals = np.linspace(p_vals[0], p_vals[-1], 50)
y_vals = popt[0] + popt[1]/x_vals + popt[2]/x_vals**2 + popt[3]/x_vals**3 + popt[4]/x_vals**4
plt.plot(x_vals, y_vals, linestyle='-', color='C3')

plt.hlines(0,0.3,2.1, linestyle=':', color='k', alpha=0.5)
plt.plot(p_vals, max_deltas, marker='o', linestyle='', color='C3', label='Maximum $\delta$')
plt.errorbar(p_vals, max_deltas, yerr=errors, fmt='o', color='C3', capsize=5)
plt.savefig('ResultPlot1c.png', dpi=500, bbox_inches='tight')
plt.show()

#%% Result Plot 2: Average perpendicular distance from geodesic vs p 

p_vals = []
avg_dists = []

# Loop to iterate over p values
for i in range(1,24):
    
    p_i = (i+1)/12
    p_vals.append(p_i)
    avg_dists_i = []
    
    # Loop to iterate over random seeds
    for j in range(0,5):
        np.random.seed(j)
        g_i = ga.Graph(nodes=50, p=p_i, geodesic=1)
        avg_dist_j = g_i.perp_dist()[0]
        avg_dists_i.append(avg_dist_j)
    
    avg_dist_i = np.mean(avg_dists_i)
    error_i = np.std(avg_dists_i)
    avg_dists.append(avg_dist_i)
    errors.append(error_i)
   
fig, ax = plt.subplots()
fig.set_size_inches(5, 5)
ax.set_xlabel('$p$')
ax.set_ylabel('Average Deviation')
ax.set_ylim(0.225,0.230)

plt.plot(p_vals, avg_dists, marker='o', linestyle='')
#plt.savefig('ResultPlot2.png', dpi=500, bbox_inches='tight')
plt.show()

#%% Result Plot 3: Difference in length between geodesic and path vs p

p_vals = []
sh_diffs = [] # Shortest difference list
lo_diffs = [] # Longest difference list
sh_errors = []
lo_errors = []

# Loop to iterate over p values
for i in range(1,24):
    
    p_i = (i+1)/12
    p_vals.append(p_i)
    sh_diffs_i = []
    lo_diffs_i = []
    
    # Loop to iterate over random seeds
    for j in range(0,10):
        np.random.seed(j)
        g_i = ga.Graph(nodes=50, p=p_i, geodesic=1)
        sh_diff_j = g_i.path_geo_diff()[0]
        sh_diffs_i.append(sh_diff_j)
        lo_diff_j = g_i.path_geo_diff()[1]
        lo_diffs_i.append(lo_diff_j)
        
    sh_diff_i = np.mean(sh_diffs_i)
    sh_error_i = np.std(sh_diffs_i)
    sh_diffs.append(sh_diff_i)
    sh_errors.append(sh_error_i)
    lo_diff_i = np.mean(lo_diffs_i)
    lo_error_i = np.std(lo_diffs_i)
    lo_diffs.append(lo_diff_i)
    lo_errors.append(lo_error_i)
    
fig, ax = plt.subplots()
fig.set_size_inches(5, 5)
ax.set_xlabel('$p$')
ax.set_ylabel('$|L_{p, path} - L_{p, geo}|$')

plt.plot(p_vals, sh_diffs, marker='o', linestyle='', color='C30', label='Shortest Metric Path')
plt.plot(p_vals, lo_diffs, marker='o', linestyle='', color='C2', label='Longest Metric Path')
plt.legend()
#plt.savefig('ResultPlot3.png', dpi=500, bbox_inches='tight')
plt.show()
