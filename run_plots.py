#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 02:27:04 2022

@author: ali
"""

import Graph as ga
import plots as pls
import numpy as np
import matplotlib.pyplot as plt

plotter = pls.Plots()

#%% delta vs p "all"

plotter.delta_vs_p(p_start = 0.01, p_end = 2, p_step = 0.1, nodes = 50, trials = 10, delta_type = "all")

#%% plot 2 single avg
#plotter.perp_vs_p(p_start = 0.4, p_step= 0.01, trials = 50, nodes = 10, poly_order = 16, radius =1923)

#%% plot 2 single max

#plotter.perp_vs_p(p_start = 0.4, p_step= 0.01, trials = 50, nodes = 10, poly_order = 16, radius =1923, perp_type="max")
#%% plot 2 single min
#plotter.perp_vs_p(p_start = 0.4, p_step= 0.01, trials = 50, nodes = 10, poly_order = 16, radius =1923, perp_type="min")

#%% plot 2 all
#plotter.perp_vs_p(p_start = 0.4, p_step= 0.01, trials = 10, nodes = 10, poly_order = 16, radius =1923, perp_type = "all")

#max and min seem to give nonsense (min just gives 0 even though I've set geodesic = 0 ?)
#%% plot 3 metric path
#plotter.length_vs_p(p_start = 0.1, p_end = 2, p_step = 0.1, nodes = 50, trials = 10, radius = 1e78)

# =============================================================================
# both have intercepts at exactly 1 
#x-intercept for Shortest: [1.         0.90689235 0.79641628 0.70523713 0.61178326 0.56451401
#0.52595561]
#x-intercept for Longest: [9.18425967 3.37901698 2.44271483 1.63800233 1.32443848 1.11342973
 #1.        ]
# =============================================================================
