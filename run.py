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

g = ga.Graph(nodes = 5, radius= 0.8, nx_plot=1, plt_plot =1)
print("The longest network path is", g.longest_path(), \
      "with length", g.longest_path_length())

print("The shortest network path is", g.shortest_path(), \
      "with length", g.shortest_path_length())
