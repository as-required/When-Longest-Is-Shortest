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

#plotter.delta_vs_p(p_start = 0.01, p_end = 4, p_step = 0.1, nodes = 100, trials = 100, delta_type = "avg", radius = 1e100)

#p_start = 0.01, p_end = 4, p_step = 0.1, nodes = 100, trials = 100:
#Time taken: 1290.1173989772797
#avg delta_3.91 = 0.0850739936478105 so triangle inequality: True
#x-intercept for avg: [0.9998672]

"""
for p_start = 0.01, p_end = 4, p_step = 0.1, nodes = 10, trials = 50, delta_type = "avg", radius = 1e100

Time taken: 5.906562089920044
avg delta_3.91 = 0.1471470368848317 so triangle inequality: True
p-intercept: 0.999888851538669 +- 1.1859813574974518e-05 DOESN'T AGREE'
"""

"""
for p_start = 0.01, p_end = 4, p_step = 0.1, nodes = 100, trials = 100, delta_type = "avg", radius = 1e100
Time taken: 1331.1540732383728
avg delta_3.91 = 0.0850739936478105 so triangle inequality: True
p-intercept: 0.9998716468355326 +- 1.2042735034018568e-05
"""
#%% plot 2 single avg
plotter.perp_vs_p(p_start = 0.01, p_step= 0.01, trials = 100, nodes = 100, poly_order = 25, radius =1e100)

# for p_start = 0.01, p_step= 0.01, trials = 50, nodes = 10, poly_order = 25, radius =1e100
# 
"""Time taken: 25.0464289188385
The average difference in longest and shortest metric path deviations from the geodesic is: 0.2128978795213209
The minimum difference in longest and shortest metric path deviations from the geodesic is: 0.0
The maximum difference in longest and shortest metric path deviations from the geodesic is: 0.6386936385639627
x-intercept for avg: [-0.90856435  0.9974721   0.49830962  0.42101218]"""

"""
for p_start = 0.01, p_step= 0.01, trials = 20, nodes = 10, poly_order = 25, radius =1e100
used >0.001 for filtering roots
Time taken: 10.360100984573364
The average difference in longest and shortest metric path deviations from the geodesic is: 0.11663362187360189
The minimum difference in longest and shortest metric path deviations from the geodesic is: 0.0
The maximum difference in longest and shortest metric path deviations from the geodesic is: 0.417667145776281
p-intercept: 0.9995134950631189 +- 0.00027284965513214366 DOESN'T AGREEE'

"""
#%% plot 2 single max

#plotter.perp_vs_p(p_start = 0.4, p_step= 0.01, trials = 50, nodes = 10, poly_order = 16, radius =1923, perp_type="max")
#%% plot 2 single min
#plotter.perp_vs_p(p_start = 0.4, p_step= 0.01, trials = 50, nodes = 10, poly_order = 16, radius =1923, perp_type="min")

#%% plot 2 all
#plotter.perp_vs_p(p_start = 0.4, p_step= 0.01, trials = 10, nodes = 10, poly_order = 16, radius =1923, perp_type = "all")

#max and min seem to give nonsense (min just gives 0 even though I've set geodesic = 0 ?)
#%% plot 3 metric path
#plotter.length_vs_p(p_start = 0.9, p_end = 1.1, p_step = 0.005, nodes = 100, poly_order = 16, trials = 100, radius = 1e100)

# =============================================================================
# both have intercepts at exactly 1 
#x-intercept for Shortest: [1.         0.90689235 0.79641628 0.70523713 0.61178326 0.56451401
#0.52595561]
#x-intercept for Longest: [9.18425967 3.37901698 2.44271483 1.63800233 1.32443848 1.11342973
 #1.        ]
# =============================================================================
""" for p_start = 0.9, p_end = 1.1, p_step = 0.005, nodes = 20, poly_order = 16, trials = 10, radius = 1e100
Difference between geodesic and shortest metric path: 6.933881740533998e-05
Difference between geodesic and longest metric path: 0.08529098535171808
Difference between geodesic and shortest metric path: 6.933881740533998e-05
Difference between geodesic and longest metric path: 0.08529098535171808
x-intercept for Shortest: [-0.77409052 -0.64863987  1.          0.98631976  0.9476083   0.93557373
  0.91229288]
x-intercept for Longest: [-2.51935653  1.20055843  1.09757516  1.06211751  1.01320924  1.
  0.79057233]
/opt/anaconda3/lib/python3.8/site-packages/scipy/optimize/minpack.py:833: OptimizeWarning: Covariance of the parameters could not be estimated
  warnings.warn('Covariance of the parameters could not be estimated',"""

"""
after jackknife
for p_start = 0.9, p_end = 1.1, p_step = 0.005, nodes = 20, poly_order = 16, trials = 10, radius = 1e100
Time taken: 3.4676313400268555
Difference between geodesic and shortest metric path: 6.933881740533998e-05
Difference between geodesic and longest metric path: 0.08529098535171808
Difference between geodesic and shortest metric path: 6.933881740533998e-05
Difference between geodesic and longest metric path: 0.08529098535171808
/opt/anaconda3/lib/python3.8/site-packages/scipy/optimize/minpack.py:833: OptimizeWarning: Covariance of the parameters could not be estimated
  warnings.warn('Covariance of the parameters could not be estimated',
p-intercept for Shortest: 0.9998575588151902 +- 0.0001424412036686697
/opt/anaconda3/lib/python3.8/site-packages/scipy/optimize/minpack.py:833: OptimizeWarning: Covariance of the parameters could not be estimated
  warnings.warn('Covariance of the parameters could not be estimated',
p-intercept for Longest: 1.0001418565560798 +- 0.0001418565593946054
"""

"""
for plotter.length_vs_p(p_start = 0.9, p_end = 1.1, p_step = 0.005, nodes = 100, poly_order = 16, trials = 100, radius = 1e100)

Time taken: 641.8998870849609
Difference between geodesic and shortest metric path: 6.234673903815491e-07
Difference between geodesic and longest metric path: 0.10689144422187358
Difference between geodesic and shortest metric path: 6.234673903815491e-07
Difference between geodesic and longest metric path: 0.10689144422187358
/opt/anaconda3/lib/python3.8/site-packages/scipy/optimize/minpack.py:833: OptimizeWarning: Covariance of the parameters could not be estimated
  warnings.warn('Covariance of the parameters could not be estimated',
p-intercept for Shortest: 0.9998642842263675 +- 0.0001357134698645545
/opt/anaconda3/lib/python3.8/site-packages/scipy/optimize/minpack.py:833: OptimizeWarning: Covariance of the parameters could not be estimated
  warnings.warn('Covariance of the parameters could not be estimated',
p-intercept for Longest: 1.0001367206592697 +- 0.00013672066288506404

"""