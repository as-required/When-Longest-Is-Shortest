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
from scipy.optimize import curve_fit

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

#%% test triangle
ghalf.triangle()
g1.triangle()
g2.triangle()

# We found that:
# 1. delta = 0 for p = 1
# 2. delta > 0 for p = 2
# 3. delta < 0 for p = 1/2
# This agrees with what we expected.

#%% test perp_dist to geo
ghalf.perp_dist()
g1.perp_dist()
g2.perp_dist()

# Found that:
# 1. average and maximum difference = 0 for p = 1 
# 2. average and maximum difference > 0 for p = 2
# 3. average and maximum difference < 0 for p = 1/2
# 4. minimum difference = 0 for all p values (as source and sink lie on geodesic)
# This agrees with what we expected:
# The longest path is a better approximation to the geodesic if diff < 0 (we expect this for p < 1)
# The shortest path is a better approximation to the geodesic if diff > 0 (we expect this for p > 1)

#%% test diff of path compared to geo
ghalf.path_geo_diff()
g1.path_geo_diff()
g2.path_geo_diff()
