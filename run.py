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

#%% Testing 
g1 = ga.Graph(nodes = 5, radius = 10.0, p = 1, nx_plot = 1)
ghalf = ga.Graph(nodes = 5, radius = 10.0, p = 1/2, nx_plot = 1)
g2 = ga.Graph(nodes = 5, radius= 10.0, nx_plot = 1) #Euclidean graph

graphs = [[ghalf, 0.5], [g1, 1.0], [g2, 2.0]]

for i in range(len(graphs)):
    
    print("The longest network path for p = {} is".format(graphs[i][1]), graphs[i][0].longest_path(), \
      "with length", graphs[i][0].longest_path_length())

    print("The shortest network path for p = {}  is".format(graphs[i][1]), graphs[i][0].shortest_path(), \
      "with length", graphs[i][0].shortest_path_length())
#%% Triangle inequality
g1.triangle()
ghalf.triangle()
g2.triangle()