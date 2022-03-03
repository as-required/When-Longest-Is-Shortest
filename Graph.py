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

np.random.seed(0)

class Graph:
    
    def __init__(self, nodes = 10, p = 2, radius = 0.9, nx_plot = 0, plt_plot = 0 ):
        """
        We have assumed D = 2 (the dimension)

        Parameters
        ----------
        nodes : TYPE, optional
            DESCRIPTION. The default is 10.
        p : TYPE, optional
            DESCRIPTION. The default is 2.
        nx_plot : TYPE, optional
            DESCRIPTION. The default is 0.
        plt_plot : TYPE, optional
            DESCRIPTION. The default is 0.

        Returns
        -------
        None.

        """
        
        self._nodes = nodes
        self._p = p
        self._radius = radius
        self._start_time = time.time()
        
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
            self.draw_network()
            
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
                    self._graph.add_edge(i,j)
                    
                    
            
    # if either j coord is less than the i coord, the cube space rule is violated
                else:
                    continue
                
    def longest_path(self):
        
        return dg.dag_longest_path(self._graph)
    
    def longest_path_length(self):
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
    
    def shortest_path_length(self):
        return len(self.shortest_path())
                
    def draw_network(self):
        
        fig, ax = plt.subplots()
        fig.set_size_inches(5, 5)
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
        ax.set_title(r"$p = {}$".format(self._p))
        #plt.axis('on') for some reason spyder crashes if using this
        nx.draw_networkx(self._graph, self._pos, arrows = True, ax = ax)
        plt.savefig("nx_graph_p={}.png".format(str(self._p)), dpi = 500, bbox_inches = "tight")
        plt.show()
        
    def draw_plt(self):
        
        plt.figure(figsize=(5,5))
        
        for i in range(len(self._x_coords)):
    
            plt.plot(self._x_coords[i], self._y_coords[i], marker = 'o')
        
        plt.title("r = {}".format(self._radius))
        plt.savefig("plt_graph", dpi = 500, bbox_inches = "tight")
        plt.show()   


    def triangle(self, n1 = 0, n2 = 1, n3 = 2):
        """
        Verify the triangle inequality

        """
        triangle_inequality_satisfied = False
        # Gromov delta: d(z,x) + d(z,y) - d(x,y) let n1 = x, n2 = y, n3 = z
        x = self._pos[n1]
        y = self._pos[n2]
        z = self._pos[n3]
        
        # I was being stupid here and calculated Euclidean distance
        #delta = ( ((self._pos[n3][0] - self._pos[n1][0]) ** 2 + (self._pos[n3][1] - self._pos[n1][1]) ** 2) ** 1/2) \
        #        + (( (self._pos[n3][0] - self._pos[n2][0]) ** 2 + (self._pos[n3][1] - self._pos[n2][1]) ** 2) ** 1/2) \
        #        - (( (self._pos[n1][0] - self._pos[n2][0]) ** 2 + (self._pos[n1][1] - self._pos[n2][1]) ** 2) ** 1/2)
        
        delta = self.L_p(z,x) + self.L_p(z,y) - self.L_p(x,y)
        
        if delta >= 0:
            triangle_inequality_satisfied = True
            
        print("delta_{} =".format(self._p), delta, "so triangle inequality:", triangle_inequality_satisfied)
        return delta        
    
    
    def geodesic(self):
        """
        Adds a geodesic (straight line) between (0,0) and (1,1).

        """
        start = list(self._pos)[0]
        end = list(self._pos)[-1]
        
        geo_dist = self.L_p(start, end)
        self._graph.add_edge(start, end, weight = geo_dist)
        plt.figure(figsize=(5,5))
        nx.draw_networkx(self._graph, self._pos, arrows = True)
        plt.show()
        
        
    def add_edge_weights(self):
        
        for i in range(len(self._x)):
            for j in range(i+1, len(self._x)):
        
                point_1 = [self._r[i][0], self._r[i][1]]
                point_2 = [self._r[j][0], self._r[j][1]]

                if self._r[j][1] > self._r[i][1] and \
                    (self.L_p(point_1, point_2) < self._radius):
                        
                    self._x_coords.append([self._r[i][0], self._r[j][0]]) 
                    self._y_coords.append([self._r[i][1], self._r[j][1]])
                    
                    self._pos[i] = (self._r[i][0], self._r[i][1])
                    self._pos[j] = (self._r[j][0], self._r[j][1])
                    
                    edge_weight = self.L_p(self._pos[i], self._pos[j])
        
                    self._graph.add_edge(i,j, weight=edge_weight)
        

    def longest_weighted_path(self):
        """
        Finds the longest metric path (which takes weights into account).
        
        """
        weights = []
        
        for i in range(len(self._pos)):
            for j in range(len(self._pos)-1):
                ij_weight = self.L_p(self._pos[i], self._pos[j])
                weights.append(ij_weight)
                
        return dg.dag_longest_path(self._graph, weight=weights)
    
    
    def longest_weighted_path_length(self):
        """
        Finds the length of the longest metric path (which takes weights into account).
        
        """
        weights = []
        
        for i in range(len(self._pos)):
            for j in range(len(self._pos)-1):
                ij_weight = self.L_p(self._pos[i], self._pos[j])
                weights.append(ij_weight)
                
        return dg.dag_longest_path_length(self._graph, weight=weights)
        