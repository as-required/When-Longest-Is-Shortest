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

nx.draw_networkx(G, pos, arrows=1, edge_color='black', width=1, linewidths=1, node_size=500, node_color='C30', \
                 arrowsize=18)
nx.draw_networkx_edge_labels(G, pos, edge_labels=labels) # Add edge labels

plt.savefig('TriangleInequality.png', dpi=500)
plt.show()

#%% Geodesic, Longest Network Path and Shortest Network Path Visualisation
# Blue represents longest network path; red represents shortest network path; grey represents geodesic.

edges1 = [[0,1], [1,2], [2,4], [4,3], [0,5], [5,3]]
pos1 = {0: (0.0,0.0), 1: (0.1,0.5), 2: (0.4,0.6), 3: (1.0,1.0), 4: (0.6,0.9), 5: (0.7,0.2)}
edge_colours1 = ['C30', 'C3', 'C30', 'C30', 'C30', 'C3']

G1 = nx.DiGraph()
G1.add_edges_from(edges1)

fig, ax = plt.subplots()
fig.set_size_inches(5, 5)
plt.axis('off')

nx.draw_networkx(G1, pos1, arrows=1, edge_color=edge_colours1, width=1, linewidths=1, node_size=300, \
                 node_color='black', font_color='black', arrowsize=20)
    
x = np.linspace(0,1)
y = x

plt.plot(x, y, linestyle='--', color='black', alpha=0.6) # Add geodesic

plt.savefig('GeoVsNetworkPaths.png', dpi=500)
plt.show()

#%% Manhattan Space Visualisation

edges2 = [[0,1], [0,2], [2,1], [0,3], [3,4], [4,5], [5,1]]
pos2 = {0: (0.0,0.0), 1: (1.0,1.0), 2: (0.0,1.0), 3: (0.7,0.0), 4: (0.7,0.3), 5: (1.0,0.3)}
edge_colours2 = ['C30', 'C2', 'C3', 'C2', 'C3', 'C3', 'C3']

G2 = nx.DiGraph()
G2.add_edges_from(edges2)

fig, ax = plt.subplots()
fig.set_size_inches(5, 5)
plt.axis('off')

nx.draw_networkx(G2, pos2, arrows=1, edge_color=edge_colours2, width=1, linewidths=1, node_size=500, \
                 node_color='black', font_color='white', arrowsize=20)
    
plt.savefig('ManhattanSpace.png', dpi=500)
plt.show()

