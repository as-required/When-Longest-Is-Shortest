#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 19:53:55 2022

@author: ali and autumn
"""

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import networkx.algorithms.dag as dg
import networkx.algorithms.shortest_paths.unweighted as s_u
import networkx.algorithms.shortest_paths.weighted as s_w
import time
start_time = time.time()
np.random.seed(1) #noticed that if I put this within the class, each p will get the same nodes
# also using seed 1 as seed 0 gives node 2 right next to node 0,
# blocking its path to node 4 for the geodesic (doesn't really matter but doesn't look nice)
class Graph:
    
    def __init__(self, nodes = 10, p = 2, radius = 0.9, auto_radius = 1,\
                 geodesic = 0, nx_plot = 0, plt_plot = 0, show_weights = 0):
        """
        We have assumed D = 2 (the dimension)

        Parameters
        ----------
        nodes : TYPE, optional
            DESCRIPTION. The default is 10.
        p : TYPE, optional
            DESCRIPTION. The default is 2.
        auto_radius = sets radius to 5x the avg interparticle spacing
        seed: sets the random seed
        geodesic: adds geodesic with metric weight
        nx_plot : TYPE, optional
            DESCRIPTION. The default is 0.
        plt_plot : TYPE, optional
            DESCRIPTION. The default is 0.

        Returns
        -------
        None.

        """
        self._geodesic = geodesic
        self._show_weights = show_weights
        self._nodes = nodes
        self._p = p
        self._radius = radius
        self._start_time = time.time()
        
# =============================================================================
#         automate the radius
# =============================================================================
        if auto_radius == 1:
            rho = self._nodes # node density \rho = N/L^D
            # put the const as 5 as this seems to be the lowest const that connects for 5 nodes
            self._radius = 5 * (rho ** (-0.5) ) # avg internode spacing  r = const * \rho^{-1/D}
# =============================================================================
#       initialise the nodes
# =============================================================================
        self._x = []
        self._y = []
        self._r = [] # contains [[x_coord, y_coord],..]
        
        for i in range(0, self._nodes):
            x_i = np.random.uniform(0.0,1.0)
            y_i = np.random.uniform(0.0,1.0)
            self._x.append(x_i)
            self._y.append(y_i)

        # manually set first coords as (0.0, 0.0) and last as (1.0, 1.0)
        self._x[0] = 0.0
        self._x[-1] = 1.0
        self._y[0] = 0.0
        self._y[-1] = 1.0
        self._x.sort() # puts numbers in ascending order
        
        for i in range(len(self._x)):
            x_i = self._x[i]
            y_i = self._y[i]
            r_i = [x_i,y_i]
            self._r.append(r_i)
        
# =============================================================================
#         initialise the graph
# =============================================================================
        self._graph = nx.DiGraph()
        # lists to store ordered coords 
        # eg/ if i is before j, x_coords = [ [i_x, j_x],...  ]
        self._x_coords = []
        self._y_coords = []
        
        # initialise dictionary to store positions of the nodes
        self._pos = {}
        
        # populate graph with nodes to connect later
        self._graph.add_nodes_from(range(len(self._x)))
        
# =============================================================================
#         add edges to the graph
# =============================================================================
        self.add_edges()
        
# =============================================================================
#         plot nx graph
# =============================================================================
        if nx_plot == 1:
            self.draw_network(geodesic = self._geodesic,\
                              show_weights = self._show_weights)
            
# =============================================================================
#         plot plt graph
# =============================================================================
        if plt_plot == 1:
            self.draw_plt()
        
        print ("Time taken:", time.time() - start_time)


    def L_p(self, point_1, point_2):
        """
        
        I have verified that g2.L_p(g2._pos[0], g2._pos[g2._nodes-1]) = sqrt(2) as expected

        Parameters
        ----------
        point_1 : TYPE
            DESCRIPTION.
        point_2 : TYPE
            DESCRIPTION.

        Returns
        -------
        dist : TYPE
            DESCRIPTION.

        """
        
        sigma = (np.abs(point_1[0] - point_2[0]) ** self._p) + (np.abs(point_1[1] - point_2[1]) ** self._p)
        dist = sigma ** (1/self._p)
        
        return dist
        
        
    
    def add_edges(self):
        # adding position of node 0 to account for error where it has no value in the pos dict
        self._pos[0] = (0.0,0.0)
        # hard code in the last node
        self._pos[self._nodes - 1] = (1.0, 1.0)
        
        for i in range(len(self._x)): # i denotes the pair that pair j is being compared to
            for j in range(i+1, len(self._x)): # so now we don't need to do the x coord cube space check
        # iterates through the y coord of each pair
                # store coords of nodes i and j in pos (do this before cubespace criteria so that every node has a position)
                self._pos[i] = (self._r[i][0], self._r[i][1])
                self._pos[j] = (self._r[j][0], self._r[j][1])
        # define points for L_p calc
                point_1 = [self._r[i][0], self._r[i][1]]
                point_2 = [self._r[j][0], self._r[j][1]]
        # If it is greater than the current x AND y coord, join the points via a straight line
            #if r[j][0] > r[i][0] and r[j][1] > r[i][1]:
                # x coords already sorted so we only have to check y coords
# =============================================================================
# there is a problem here in that, if radius is too small, so no points are connected,
# there is a "Node has no position" error
# =============================================================================
                if self._r[j][1] > self._r[i][1] and \
                    (self.L_p(point_1, point_2) < self._radius): # I'm not sure if < r will work
                    # to add in L_p, put and statement where we choose a radius, measured by d_p
                    # gather the two x coords and two y coords, with i first to preserve causality
                    self._x_coords.append([self._r[i][0], self._r[j][0]]) 
                    self._y_coords.append([self._r[i][1], self._r[j][1]])
                    print()
                    # connect i and j on graph
                    edge_weight = self.L_p(self._pos[i], self._pos[j]) # adding weights
                    self._graph.add_edge(i,j,\
                                         weight = edge_weight)
                    
                    
            
    # if either j coord is less than the i coord, the cube space rule is violated
                else:
                    continue
                
    def longest_path(self):
        
        return dg.dag_longest_path(self._graph, weight = 1) #to override edge data
    
    def longest_path_length(self):
        """
        Counts EDGES
        
        but note the method = "djikstra" by default, I found an R package that has Minkowski
        and I'm trying to import it to python
        THIS IS FINE, THE NETWORK PATH IS UNAFFECTED AND THE METRIC PATH JUST RELIES ON THE WEIGHTS

        NOTE: this comes from the networkx.algorithms.dag library, which
        1) doesn't have anything for shortest paths
        2) uses a topological sort, which I don't think we want
        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return dg.dag_longest_path_length(self._graph, weight = 1) # need to specify weight = 1
        # as it automatically makes it weighted
    
    def longest_weighted_path(self):
        """
        Finds the longest metric path (which takes weights into account).
        
        """
        # dag_longest_path uses the edge weights by default
        return dg.dag_longest_path(self._graph)
    
    
    def longest_weighted_path_length(self):
        """
        Finds the length of the longest metric path (which takes weights into account).
        
        """
        # dag_longest_path uses the edge weights by default
        return dg.dag_longest_path_length(self._graph)
    
    def shortest_path(self):
        """
        display shortest path between start and end nodes as a list

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return nx.shortest_path(self._graph, source = list(self._graph.nodes)[0],\
                                target = list(self._graph.nodes)[-1])
    def shortest_path_length(self):
        """
        Counts EDGES
        
        but note the method = "djikstra" by default, I found an R package that has Minkowski
        and I'm trying to import it to python

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return nx.shortest_path_length(self._graph, source = list(self._graph.nodes)[0],\
                                target = list(self._graph.nodes)[-1])
    
    def shortest_weighted_path(self):
        
        return nx.shortest_path(self._graph, source = list(self._graph.nodes)[0],\
                                target = list(self._graph.nodes)[-1], weight = "weight")
    def shortest_weighted_path_length(self):
         
        return nx.shortest_path_length(self._graph, source = list(self._graph.nodes)[0],\
                                 target = list(self._graph.nodes)[-1], weight = "weight")
    
    def shortest_path_bidirectional(self):
        """
        This seems to be a more fundamental version of the shortest_path method
        which I initially first used. I don't see a reason to use this, so I made 
        a separate "shortest_path" method which seesm to give the same output regardless

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        start = list(self._graph.nodes)[0]
        end = list(self._graph.nodes)[-1]
        return s_u.bidirectional_shortest_path(self._graph,start, end)
                
    def draw_network(self, geodesic = 0, show_weights = 0):
        
        if geodesic == 1: #adds geodesic with L_p weight
            start = self._pos[0]
            end = self._pos[self._nodes -1]
            
            geo_dist = self.L_p(start, end)
            self._geo_dist = geo_dist
            self._graph.add_edge(0, self._nodes - 1, weight = geo_dist)
        edge_weights = nx.get_edge_attributes(self._graph,'weight')
        self._edge_weights = edge_weights
        
        fig, ax = plt.subplots()
        fig.set_size_inches(5, 5)
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
        ax.set_title(r"$p = {}$".format(self._p))
        #plt.axis('on') for some reason spyder crashes if using this
        
        if show_weights == 1:
            # round the edge weights to 2dp to display on fig
            edge_weights_2dp = edge_weights.copy()
            for i in edge_weights_2dp:
                edge_weights_2dp[i] = round(edge_weights_2dp[i], 2)
        
            # this only draws on the edge labels
            nx.draw_networkx_edge_labels(self._graph, self._pos,\
                                     edge_labels = edge_weights_2dp, ax = ax,\
                                         verticalalignment= "baseline")
        nx.draw_networkx(self._graph, self._pos, arrows = 1, ax = ax)
        plt.savefig("nx_graph_p={}.png".format(str(self._p)), dpi = 500, bbox_inches = "tight")
        plt.show()
        
    def draw_plt(self):
        
        plt.figure(figsize=(5,5))
        
        for i in range(len(self._x_coords)):
    
            plt.plot(self._x_coords[i], self._y_coords[i], marker = 'o')
        
        plt.title("r = {}".format(self._radius))
        plt.savefig("plt_graph", dpi = 500, bbox_inches = "tight")
        plt.show()   


    def triangle(self):
        # Removed "end_node", "start_node", and "intermediate node" arguments as cube space rule is now built in.
        """
        Verify the triangle inequality
        
        For p<1: order is imperative!
        Intermediate must be reachable from start by the cube space rule!
        need to do delta = d(start,intermediate) + d(intermediate, end) - d(start,end)
        [this doesn't affect p>1]

        """
        triangle_inequality_satisfied = False

        # Hard code in the positions of the initial and final nodes         
        self._pos[0] = (0.0,0.0)
        self._pos[self._nodes - 1] = (1.0, 1.0)
        
        # Apply cube space rule so that the correct nodes do not need to be manually added
        for i in range(len(self._x)):
            for j in range(i+1, len(self._x)):
                for k in range(i+2, len(self._x)):
                
                    # Store node positions
                    self._pos[i] = (self._r[i][0], self._r[i][1])
                    self._pos[j] = (self._r[j][0], self._r[j][1])
                    self._pos[k] = (self._r[k][0], self._r[k][1])
                    
                    # Define points for L_p calculation
                    point_1 = [self._r[i][0], self._r[i][1]]
                    point_2 = [self._r[j][0], self._r[j][1]]
                    point_3 = [self._r[k][0], self._r[k][1]]
                    
                    # x coordinates are already in ascending order, so only test y coordinates
                    # The following code is very messy - there must be a more concise way to do it?
                    
                    # Order i < j < k
                    if self._r[k][1] > self._r[j][1] > self._r[i][1] and \
                        (self.L_p(point_1, point_2) < self._radius) and \
                            (self.L_p(point_2, point_3) < self._radius):
                                
                                start_coords = self._pos[i]
                                intermediate_coords = self._pos[j]
                                end_coords = self._pos[k]
                    
                    
                    # Order i < k < j
                    elif self._r[j][1] > self._r[k][1] > self._r[i][1] and \
                        (self.L_p(point_1, point_3) < self._radius) and \
                            (self.L_p(point_3, point_2) < self._radius):
                                
                                start_coords = self._pos[i]
                                intermediate_coords = self._pos[k]
                                end_coords = self._pos[j]
                    
                    
                    # Order j < k < i
                    elif self._r[i][1] > self._r[k][1] > self._r[j][1] and \
                        (self.L_p(point_2, point_3) < self._radius) and \
                            (self.L_p(point_3, point_1) < self._radius):
                                
                                start_coords = self._pos[j]
                                intermediate_coords = self._pos[k]
                                end_coords = self._pos[i]
                                

                    # Order j < i < k
                    elif self._r[k][1] > self._r[i][1] > self._r[j][1] and \
                        (self.L_p(point_2, point_1) < self._radius) and \
                            (self.L_p(point_1, point_3) < self._radius):
                                
                                start_coords = self._pos[j]
                                intermediate_coords = self._pos[i]
                                end_coords = self._pos[k]
                                
                    
                    # Order k < i < j
                    elif self._r[j][1] > self._r[i][1] > self._r[k][1] and \
                        (self.L_p(point_3, point_1) < self._radius) and \
                            (self.L_p(point_1, point_2) < self._radius):
                                
                                start_coords = self._pos[k]
                                intermediate_coords = self._pos[i]
                                end_coords = self._pos[j]
                                
                       
                    # Order k < j < i
                    elif self._r[i][1] > self._r[j][1] > self._r[k][1] and \
                        (self.L_p(point_3, point_2) < self._radius) and \
                            (self.L_p(point_2, point_1) < self._radius):
                                
                                start_coords = self._pos[k]
                                intermediate_coords = self._pos[j]
                                end_coords = self._pos[i]
                        

                    # Modified Gromov delta: d(start,intermediate) + d(intermediate, end) - d(start,end)
                    delta = self.L_p(start_coords, intermediate_coords) + self.L_p(intermediate_coords, end_coords)\
                        - self.L_p(start_coords, end_coords)
        
                    if delta >= 0:
                        triangle_inequality_satisfied = True
            
                    print("delta_{} =".format(self._p), delta, "so triangle inequality:", triangle_inequality_satisfied)
                    return delta    
        
        
        # start_coords = self._pos[start_node]
        # intermediate_coords = self._pos[intermediate_node]
        # end_coords = self._pos[end_node]