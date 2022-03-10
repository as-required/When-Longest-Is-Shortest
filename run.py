#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 01:56:36 2022

@author: ali
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

#%% trinagle inequality
g2.triangle()
g1.triangle()
ghalf.triangle()
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

