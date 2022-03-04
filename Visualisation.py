#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 15:00:50 2022

@author: autumn
"""

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

#%% Triangle Inequality Visualisation

np.random.seed(1)

edges = [['x', 'y'], ['y', 'z'], ['x', 'z']]
pos = {'x': (0.0,0.0), 'y': (0.7,0.4), 'z': (1.0,1.0)}
labels = {('x','y'): 'd(x,y)', ('y','z'): 'd(y,z)', ('x','z'): 'd(x,z)'} # Dictionary to store edge labels

G = nx.DiGraph()
G.add_edges_from(edges)

fig, ax = plt.subplots()
fig.set_size_inches(5, 5)
plt.axis('off')

nx.draw_networkx(G, pos, arrows=1, edge_color='black', width=1, linewidths=1, node_size=500, node_color='C30')
nx.draw_networkx_edge_labels(G, pos, edge_labels=labels) # Add edge labels

plt.savefig('TriangleInequality.png', dpi=500)
plt.show()

#%% Geodesic, Longest/Shortest Metric Path, and Longest/Shortest Network Path Visualisation

