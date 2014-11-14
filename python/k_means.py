# -*- coding: utf-8 -*-
"""
Created on Tue May 27 15:28:39 2014

@author: Serdar Ozsezen (c)

License: BSD-3-Clause

k-means clustering example

"""

from numpy import array, mean
from numpy.linalg import norm
from matplotlib import pyplot
from copy import deepcopy as equals

class Point:
    """
    
    Describes the data points with their coordinates. 
    
    """
    def __init__(self,coords=None):
        if coords is None:
            self.coords = array([])
        else:
            self.coords = array(coords)
    
    def __repr__(self):
        return "Point("+str(self.coords)+")"
        
    def __str__(self):
        return "[" + str(self.coords[0]) + ", " + str(self.coords[1]) +"]" +\
        " Cluster: " + str(self.cluster_no)
    
    def __getitem__(self,key):
        return self.coords[key]

    def __setitem__(self,idx,val):
        self.coords[idx] = val
    
    def set_cluster(self,cluster_no):
        self.cluster_no = cluster_no
    
    def get_cluster(self):
        return self.cluster_no
    
class Centroid:
    """

    Describes the centroid points similarly to Point class.
    
    """
    def __init__(self,coords=None):
        if coords is None:
            self.coords = array([])
        else:
            self.coords = array(coords)
            
    def __eq__(self,other):
        return self.coords == other.coords

    def __repr__(self):
        return "Centroid("+str(self.coords)+")"
            
    def __str__(self):
        return "[" + str(self.coords[0]) + ", " + str(self.coords[1]) +"]"
    
    def __getitem__(self,key):
        return self.coords[key]

    def __setitem__(self,idx,val):
        self.coords[idx] = val

def distance(point1, point2):
    """

    Calculates the distance between points. 
    
    """
    
    return norm(point1.coords - point2.coords)
    
    
def assign_points(given_points,centroids):
    """

    Assigns coordinates to Point objects defined above along with centroids.  
    
    """
    point_objects = []
    for i in range(len(given_points)):
        point = Point(given_points[i])
        # Assign cluster numbers according to distances between centroid
        # and point
        if distance(centroids[0],point) > distance(centroids[1],point):
            point.set_cluster(1)
        else:
            point.set_cluster(0)
        point_objects.append(point)
    return point_objects
    
def calculate_centroids(data,centroid):
    """
    
    Calculates centroid of each cluster seperately. 
    
    """
    coordinates = []
    for i in range(len(centroid)):
        for point in data:
            if point.get_cluster() == i:
                coordinates.append(point.coords)
            cent = mean(coordinates,axis=0)
            centroid[i].coords = cent 
        coordinates = []
    return centroid
    
    
if __name__ == "__main__":
    
    points = [[0.1,0.2],[0.2,0.4],[0.3,0.2],[0.2,0.5],[0.4,1.0],[0.4,0.3],
          [1.0,1.0],[0.9,0.4],[1.0,0.3],[1.0,0.2],[1.0,0.0]]
    no_clusters = 2
    
    # Assign initial centroids to arbitrary coordinates    

    centroids =[Centroid([0.0,0.0]),Centroid([1.,1.])]
    print "Initial centroids: " + str(centroids)
    print
    # Assign centroids to points

    data = assign_points(points,centroids)
    print "Initial clusters:"
    for i in data: print i
    print "------------------------------"
    
    # Continue to calculate centroids and reassign points to accoding to new
    # centroid positions until centroid positions do not change anymore.

    iteration_no = 1
    while True:
        cent = [equals(centroids[0]),equals(centroids[1])]
        print "Iteration " + str(iteration_no) + "\n"
        centroids = calculate_centroids(data,centroids)
        print "Centroids: " + str(centroids) + "\n"
        data = assign_points(points,centroids)
        print "Clusters: "
        for i in data: print i
        print
        iteration_no += 1
        # If centroid positions are the same with the one calculated before,
        # stop calculating
        if (cent[0] == centroids[0]).all() and (cent[1] == centroids[1]).all():
            break
        
    # Draw plots for showing the clusters
    
    cluster1x = []; cluster1y = []
    cluster2x = []; cluster2y = []
    for point in data:
        if point.get_cluster() == 1:
            cluster1x.append(point[0])
            cluster1y.append(point[1])
        else:
            cluster2x.append(point[0])
            cluster2y.append(point[1])
    pyplot.scatter(cluster1x,cluster1y,c='r', marker='s',s=200)
    pyplot.scatter(cluster2x,cluster2y,c='b',s=200)
    pyplot.scatter(centroids[1][0],centroids[1][1],marker='*',c='r',s=400)
    pyplot.scatter(centroids[0][0],centroids[0][1],marker='*',c='b',s=400)
    pyplot.legend(("Cluster 1","Cluster 2","Cluster 1 centroid",
                   "Cluster 2 centroid"),loc=2)
    font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}
    pyplot.rc('font',**font)
    pyplot.xlim((-0.1,1.1))
    pyplot.ylim((-0.1,1.1))
    pyplot.grid()