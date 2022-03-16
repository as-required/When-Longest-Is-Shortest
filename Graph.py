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
#np.random.seed(1) #noticed that if I put this within the class, each p will get the same nodes
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
        # adding position of node 0 to account for error where it has no value in the pos dict
        self._pos[0] = (0.0,0.0)
        # hard code in the last node
        self._pos[self._nodes - 1] = (1.0, 1.0)
        
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
    
    def cube_space(self, point_1, point_2):
        """
        Checks if point_1 can be connected to point_2 in cube space

        Returns
        -------
        bool

        """
        # need to check both x and y coord here to not have requirements on which one is point_1 and which is point_2
        if (point_2[0] > point_1[0]) and (point_2[1] > point_1[1]):
            return True
        
        else:
            return False
		
		
    def add_edges(self):
        """
        
        note we could tidy this up using the cube space method, but why fix what's not broken?
        Returns
        -------
        None.

        """
        for i in range(len(self._x)): # i denotes the pair that pair j is being compared to
            self._pos[i] = (self._r[i][0], self._r[i][1])
            point_1 = [self._r[i][0], self._r[i][1]]
            for j in range(i+1, len(self._x)): # so now we don't need to do the x coord cube space check
        # iterates through the y coord of each pair
                # store coords of nodes i and j in pos (do this before cubespace criteria so that every node has a position)
                self._pos[j] = (self._r[j][0], self._r[j][1])
        # define points for L_p calc
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
        # as it automatically makes it weighted, so we need to set all the weights to 1
    
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
        return dg.dag_longest_path_length(self._graph) # the method automatically uses the edge weights
    
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
        
        nx.draw_networkx(self._graph, self._pos, arrows = 1, ax = ax)
        
        if show_weights == 1:
            # round the edge weights to 2dp to display on fig
            edge_weights_2dp = edge_weights.copy()
            for i in edge_weights_2dp:
                edge_weights_2dp[i] = round(edge_weights_2dp[i], 2)
        
            # this only draws on the edge labels
            nx.draw_networkx_edge_labels(self._graph, self._pos,\
                                     edge_labels = edge_weights_2dp, ax = ax)
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
         
         checks: need delta_1 = 0 and delta_(p<1) < 0:
             works, tested for p=1,0.5,0.1
        """
        triangle_inequality_satisfied = False

        
        deltas = [] # list to store the deltas from the allowed triangles
        
        
# =============================================================================
#         Apply cube space rule so that the correct nodes do not need to be manually added
#         don't need to test x coord as nodes already ordered in ascending x coord numbers in pos dict
#         so i will always be start and k end
# =============================================================================
        for i in range(len(self._x)):
            point_i = [self._pos[i][0], self._pos[i][1]]
            for j in range(i+1, len(self._x)):
                point_j = [self._pos[j][0], self._pos[j][1]]
                # check if i can be connected to j before checking any k
                if self.cube_space(point_i, point_j) == False:
                    continue 
                for k in range(j+1, len(self._x)):
                # check i-k, j-k 
                    
                    point_k = [self._pos[k][0], self._pos[k][1]]
                    
                    # check if i-k can be connected and j-k can be connected
                    if self.cube_space(point_i, point_k) == True and self.cube_space(point_j, point_k) == True:
                        delta = self.L_p(point_i, point_j) + self.L_p(point_j, point_k)\
                            - self.L_p(point_i, point_k)
                        deltas.append(delta)
                        
                    # Metric is symmetric so only need to consider when each node is intermediate
                                

                    
        deltas = np.asarray(deltas) #resave deltas as an array to take mean
        avg_delta = np.mean(deltas)
        min_delta = np.min(deltas)
        max_delta = np.max(deltas)
        
        if round(avg_delta,17) >= 0: # avoid floating point errors by rounding to 17dp
        # would get e-19 deltas for p=1 triggering false triangle identity
            triangle_inequality_satisfied = True
            
        print("avg delta_{} =".format(self._p), round(avg_delta,17), "so triangle inequality:", triangle_inequality_satisfied)
        return avg_delta, min_delta, max_delta, deltas
    
	
    def perp_dist(self):
        
        """
        Finds the perpendicular distance, D, between a node at (v,w) and the line ax + by + c = 0.
        General formula: D = |av + bw + c| / (a**2 + b**2)**(1/2)
        This perpendicular distance is called the deviation.
        The line we are interested in is the geodesic, which is the straight line connecting (0,0) and (1,1).
        Hence, a = 1, b = -1 and c = 0.
        
        """        
        a = 1
        b = -1
        c = 0
        
        # List to store deviation of each node in the longest and shortest metric paths, respectively
        long_deviations = []
        short_deviations = []
        
        longest_path_nodes = self.longest_weighted_path()
        shortest_path_nodes = self.shortest_weighted_path()    
        
        # Calculate perpendicular distance of each node in the longest metric path from the geodesic
        for i in range(len(longest_path_nodes)):
            j =longest_path_nodes[i]
            v_i = self._r[j][0]     # x coordinate of node
            w_i = self._r[j][1]     # y coordinate of node
            
            D_i = (abs(a*v_i + b*w_i + c))/((a**2 + b**2)**(0.5))
            long_deviations.append(D_i)
        
        # Calculate perpendicular distance of each node in the shortest metric path from the geodesic
        for i in range(len(shortest_path_nodes)):
            j = shortest_path_nodes[i]
            v_i = self._r[j][0]     # x coordinate of node
            w_i = self._r[j][1]     # y coordinate of node
            
            D_i = (abs(a*v_i + b*w_i + c))/((a**2 + b**2)**(0.5))
            short_deviations.append(D_i)
            
        # Convert lists into arrays
        long_deviations = np.asarray(long_deviations)
        short_deviations = np.asarray(short_deviations)
        
        # Find mean, mininum and maximum for each array
        long_avg_dev = np.mean(long_deviations)
        long_min_dev = np.min(long_deviations)
        long_max_dev = np.max(long_deviations)
        short_avg_dev = np.mean(short_deviations)
        short_min_dev = np.min(short_deviations)
        short_max_dev = np.max(short_deviations)
        
        # Find difference between deviation of longest path and shortest path
        avg_diff = long_avg_dev - short_avg_dev
        min_diff = long_min_dev - short_min_dev
        max_diff = long_max_dev - short_max_dev
        
        # The longest path is a better approximation to the geodesic if diff < 0 (we expect this for p < 1)
        # The shortest path is a better approximation to the geodesic if diff > 0 (we expect this for p > 1)
        
        print('The deviation of each node from the geodesic for nodes in the longest metric path are:', long_deviations)
        print('The deviation of each node from the geodesic for nodes in the shortest metric path are:', short_deviations)
        print('The average difference in longest and shortest metric path deviations from the geodesic is:', avg_diff)
        print('The minimum difference in longest and shortest metric path deviations from the geodesic is:', min_diff)
        print('The maximum difference in longest and shortest metric path deviations from the geodesic is:', max_diff)
        
        return avg_diff, min_diff, max_diff
    
    
    def exclude_geo(self):
        
        """
        Makes a copy of the original graph and then removes the geodesic.
        Returns longest path and shortest path excluding the geodesic.
        
        """     
        # node names
        start_node = 0
        end_node = self._nodes - 1
        
        copied_graph = self._graph.copy()
        
        if self._graph.has_edge(start_node, end_node): # check if the geodesic already exists or not
            copied_graph.remove_edge(0, self._nodes-1) # MUST USE THE NAME OF THE NODE, NOT THE POS
        
        new_shortest_length = nx.shortest_path_length(copied_graph, source = start_node,\
                                 target = end_node, weight = "weight")
            
        new_longest_length = dg.dag_longest_path_length(copied_graph) # ONLY PUT weight =1 for network lengths!
        
        return new_shortest_length, new_longest_length
        
        
        
    def path_geo_diff(self):
        
        """
        Finds the difference between the L_p distance of the longest/shortest metric path and the geodesic.
        
        """
        
        start_node = self._pos[0]
        end_node = self._pos[self._nodes-1]
        
        #Lp_shortest = self.shortest_weighted_path_length()
        #Lp_longest = self.longest_weighted_path_length()
        
        Lp_shortest = self.exclude_geo()[0]
        Lp_longest = self.exclude_geo()[1]
        Lp_geo = self.L_p(start_node, end_node)

        diff_shortest = abs(Lp_shortest - Lp_geo)
        diff_longest = abs(Lp_longest - Lp_geo)
        
        print('Difference between geodesic and shortest metric path:', diff_shortest)
        print('Difference between geodesic and longest metric path:', diff_longest)
        
        return diff_shortest, diff_longest
    
    
    
         
