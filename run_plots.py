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

plotter.delta_vs_p(p_start = 0.01, p_end = 2, p_step = 0.1, nodes = 50, trials = 100, delta_type = "avg")