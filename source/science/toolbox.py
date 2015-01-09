#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

"""
The toolbox module contains a number of functions that are used throughout the BETSE project.


Note the clip(subjectPolygon, clipPolygon) is a function defining
the Sutherland-Hodgman polygon clipping algorithm -- this is courtesy of Rosetta Code:
http://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#Python.
"""

# FIXME update the concave hull to be Sess' faster algorithm!!!

import numpy as np
import scipy.spatial as sps
import math

def flatten(ls_of_ls):
    """
    Flattens (i.e. un-nests) a nested python "list of lists".

    Parameters
    ----------
    ls_of_ls        a nested list of lists, as in: [[a,b,c],[d,e,f],[g,h,i]]


    Returns
    -------
    ls_flat        a flattened version of the input list, as in: [a,b,c,d,e,f,g,h,i]

    Notes
    -------
    Requires python nested lists of lists. Numpy arrays have their own tools for this.

    """

    ls_flat = []
    ind_map =[]
    for i, sublist in enumerate(ls_of_ls):
        for j, val in enumerate(sublist):
            ls_flat.append(val)
            ind_map.append([i,j])

    return ls_flat, ind_map



def area(p):

    """
    Calculates the area of an arbitrarily shaped polygon defined by a
    set of counter-clockwise oriented points in 2D.

    Parameters
    ----------
    p               xy list of polygon points


    Returns
    -------
    area            area of a polygon in square meters

    Notes
    -------
    The algorithm is an application of Green's theorem for the functions -y and x,
    exactly in the way a planimeter works.

    """

    return 0.5 * abs(sum(x0*y1 - x1*y0 for ((x0, y0), (x1, y1)) in zip(p, p[1:] + [p[0]])))



def clip(subjectPolygon, clipPolygon):  # This is the Sutherland-Hodgman polygon clipping algorithm

   def inside(p):

      return(cp2[0]-cp1[0])*(p[1]-cp1[1]) > (cp2[1]-cp1[1])*(p[0]-cp1[0])

   def computeIntersection():
      dc = [ cp1[0] - cp2[0], cp1[1] - cp2[1] ]
      dp = [ s[0] - e[0], s[1] - e[1] ]
      n1 = cp1[0] * cp2[1] - cp1[1] * cp2[0]
      n2 = s[0] * e[1] - s[1] * e[0]
      n3 = 1.0 / (dc[0] * dp[1] - dc[1] * dp[0])
      return [(n1*dp[0] - n2*dc[0]) * n3, (n1*dp[1] - n2*dc[1]) * n3]

   assert isinstance(subjectPolygon, list)
   assert isinstance(clipPolygon, list)
   assert len(subjectPolygon)
   assert len(clipPolygon)

   outputList = subjectPolygon
   cp1 = clipPolygon[-1]

   for clipVertex in clipPolygon:
      cp2 = clipVertex
      inputList = outputList
      outputList = []
      s = inputList[-1]

      for subjectVertex in inputList:
         e = subjectVertex
         if inside(e):
            if not inside(s):
               outputList.append(computeIntersection())
            outputList.append(e)
         elif inside(s):
            outputList.append(computeIntersection())
         s = e
      cp1 = cp2
   return(outputList)


def alpha_shape(points, alpha):
    """
    Calculate the alpha_shape of a cluster of points in 2D.

    Parameters
    ----------
    points              A numpy array listing [[x1,y1],[x2,y2]...] for a collection of 2D points

    alpha               The filtering parameter (gauges which triangles to remove from Delaunay triangulation)
                        Note, for this application alpha = 1/d_cell is ideal.

    Returns
    --------
    concave_hull    A list of the indices to vertices in the points structure which define the concave hull
                    (these are all of the points on the boundary).

    Notes
    --------
    Unlike the convex hull, the alpha shape method and concave hull work for complex, concave geometries.
    The result depends on the alpha parameter. A value of alpha = 1/d_cell gives suitable results.

    """

    tri = sps.Delaunay(points)
    tri_edges = []
    circum_r_list = []

    # loop over triangles:
    # ia, ib, ic = indices of corner points of the
    # triangle
    for ia, ib, ic in tri.vertices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]

        # Lengths of sides of triangle
        a = math.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
        b = math.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
        c = math.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)

        # Semiperimeter of triangle
        s = (a + b + c)/2.0

        # Area of triangle by Heron's formula
        area = math.sqrt(s*(s-a)*(s-b)*(s-c))

        if area > 0:
            circum_r = a*b*c/(4.0*area)

        if area == 0:
            circum_r = a*b*c/(4.0*1e-25)

        # circum_r = a*b*c/(4.0*area)
        # circum_r_list.append(circum_r)

        # Here's the radius filter:

        if circum_r < 1.0/alpha:
            #self.tri_edges_append([ia, ib])
            tri_edges.append([ia, ib])
            tri_edges.append([ib, ic])
            tri_edges.append([ia, ic])

    for i, edge in enumerate(tri_edges):  # First organize the list so that all [i,j] and [j,i] are equalized
        pt1 = edge[0]
        pt2 = edge[1]
        if pt1 > pt2:
            tri_edges[i]=[pt2,pt1]

    tri_edges.sort()

    tri_edges_len = len(tri_edges)
    concave_hull = []

    i = 1
    # Now step through to find and remove all duplicating entities and add them to a new list
    while i < tri_edges_len:
        j = i + 1
        tri_edge_i = tri_edges[i]
    
        while j < tri_edges_len and tri_edge_i == tri_edges[j]:
            j += 1
    
        if i == j - 1:
            concave_hull.append(tri_edge_i)
            i += 1
        else:
            i = j


    # Sess' improved concave hull finding algorithm.
    # self.exterior_edges = set()
    # self.interior_edges = set()
    #
    # def tri_edges_append(self, edge_list):
    #     edge_tuple = tuple(edge_list)
    #     if edge_tuple in self.exterior_edges:
    #         self.exterior_edges.remove(edge_tuple)
    #         self.interior_edges.add(edge_tuple)
    #     elif edge_tuple not in self.interior_edges:
    #         self.exterior_edges.add(edge_tuple)
    #
    # concave_hull = []
    # for edge_tuple in self.exterior_edges:
    #     concave_hull.append(list(edge_tuple))
    #
    # self.exterior_edges = None
    # self.interior_edges = None

    return concave_hull
