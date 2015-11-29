#!/usr/bin/env python3
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# This is a training module

import numpy as np
import scipy as sp
import scipy.spatial as sps
import matplotlib.pyplot as plt
from matplotlib.path import Path
import math
import time
from pprint import pprint
from matplotlib import collections  as col
import matplotlib.colors as colour

start_time = time.time()  # get a start value for timing the simulation

# define the basic class that holds variables
class NumVars(object):
    pass

# Define numerical constants for world set-up and simulation

const = NumVars()  # define the object that holds main parameters
const.wsx = 1000e-6  # the x-dimension of the world space [m]
const.wsy = 1000e-6  # the y-dimension of the world space [m]
const.rc = 5e-6  # radius of single cell
const.dc = const.rc * 2  # diameter of single cell
const.nx = int(const.wsx / const.dc)  # number of lattice sites in world x index
const.ny = int(const.wsy / const.dc)  # number of lattice sites in world y index
const.ac = 1e-6  # cell-cell separation for drawing
const.dc = const.rc * 2  # cell diameter
const.nl = 0.6  # noise level for the lattice
const.wsx = const.wsx + 5 * const.nl * const.dc  # readjust the world size for noise
const.wsy = const.wsy + 5 * const.nl * const.dc

# Step 1: define the set of 2d lattice points representing the simulation world

def makeSeeds():
    """
    lklklk

    Parameters
    ----------
    lklkllk

    Returns
    -------
    opopopop

    Notes
    -------
    Uses Numpy arrays

    """
    # first begin with linear vectors which are the "ticks" of the x and y dimensions
    x_v = np.linspace(0, (const.nx - 1) * (const.dc + const.ac), const.nx)  # create lattice vector x
    y_v = np.linspace(0, (const.ny - 1) * (const.dc + const.ac), const.ny)  # create lattice vector y

    # next define a 2d array of lattice points using the x- and y- vectors
    x_2d, y_2d = np.meshgrid(x_v, y_v)  # create 2D array of lattice points

    # define the 2d array that will store values of data applied to each cell (e.g. V, conc). For now call it 'z'.
    z_2d = np.zeros(x_2d.shape)

    # now create a matrix of points that will add a +/- deviation to each point centre
    x_rnd = const.nl * const.dc * (np.random.rand(const.ny, const.nx) - 0.5)  # create a mix of random deltas x dir
    y_rnd = const.nl * const.dc * (np.random.rand(const.ny, const.nx) - 0.5)  # create a mix of random deltas x dir

    # add the noise effect to the world point matrices and redefine the results
    x_2d = x_2d + x_rnd
    y_2d = y_2d + y_rnd

    # define a data structure that holds [x,y] coordinate points of each 2d grid-matrix entry (this is req'd for comp tissue)
    xypts = np.vstack((x_2d.ravel(), y_2d.ravel())).T

    z_xy = np.zeros(xypts.shape)   # create a similar data structure specific to cells data

    # define limits for the world space (used in plotting)
    wsxmin = -5 * const.nl * const.dc
    wsxmax = const.wsx
    wsymin = -5 * const.nl * const.dc
    wsymax = const.wsy

    # define geometric limits and centre for the cluster of points
    xmin = np.min(x_2d)
    xmax = np.max(x_2d)
    ymin = np.min(y_2d)
    ymax = np.max(y_2d)
    centx = (xmax - xmin) / 2
    centy = (ymax - ymin) / 2

# Step 2: Mask out points to define a cell cluster of a certain shape

# define points of a circle that will define shape of collective
cres = 50  # how many points desired in polygon
d_circ = xmax - xmin  # diameter of circle in x-direction  TODO try coding in an ellipse!
r_circ = d_circ / 2  # radius of circle
ind1 = np.linspace(0, 1, cres + 1)  # indices of angles defining circle points

angs = ind1 * 360 * (np.pi / 180)  # angles in radians defining circle points
circ_ptsx = r_circ * np.cos(angs) + centx  # points of the circle
circ_ptsy = r_circ * np.sin(angs) + centy  # points of the circle

crop_pts = np.vstack((circ_ptsx, circ_ptsy)).T  # reorganize points of the circle as [x,y] pairs

# turn the circle into a vector path structure so we can find the points inside of it
crop_path = Path(crop_pts, closed=True)

# create a boolean matrix mask which is 1 for points inside the circle and 0 for points outside
ws_mask_xy = crop_path.contains_points(xypts)  # create the mask for point inside the path
ws_mask = ws_mask_xy.reshape(x_2d.shape)  # reshape the boolean mask to correspond to the data grid
clust_x2d = np.ma.masked_array(x_2d, ~ws_mask)  # created a masked data structure of the x-grid
clust_y2d = np.ma.masked_array(y_2d, ~ws_mask)  # create a masked data structure of the y-grid
clust_z2d = np.ma.masked_array(z_2d,~ws_mask)  # create a masked array for the z-data structure (colour)


# create a data structure of [x,y] points that only contains the points of the cluster

clust_xy = np.array([0, 0])  # initialize the x,y points array to zero

for i in range(0, len(y_v)):  # indices to step through grid
    for j in range(0, len(x_v)):  # indices to step through grid

         if clust_x2d[i,j] is not np.ma.masked:  # if point is not masked (i.e. in the cell cluster)...
                aa=[x_2d[i,j],y_2d[i, j]]        # get the value of the x,y point by accessing x-grid and y-grid
                clust_xy = np.vstack((clust_xy,aa))  # augment the points list by adding in the new value

clust_xy = np.delete(clust_xy, 0, 0)    # delete the initialization value.


# Step 3: Computational Geometry data structures

# Delaunay Triangulation:
trimesh = sps.Delaunay(clust_xy)

# Convex Hull
bound_pts = sps.ConvexHull(clust_xy)

# Voronoi Dataset
cells = sps.Voronoi(clust_xy)

# Define a cropping boundary circle for the Voronoi diagram
ptp_bound = cells.points.ptp(axis=0)
vor_center = cells.points.mean(axis=0)

cres = 15  # how many points desired in polygon
d_circ = ptp_bound.max()  # diameter of circle in x-direction  TODO try coding in an ellipse!
r_circ = 1.08*(d_circ / 2)  # radius of circle
ind1 = np.linspace(0, 1, cres + 1)  # indices of angles defining circle points
angs = ind1 * 360 * (np.pi / 180)  # angles in radians defining circle points
circ_ptsx = r_circ * np.cos(angs) + centx  # points of the circle
circ_ptsy = r_circ * np.sin(angs) + centy  # points of the circle

crop_pts = np.vstack((circ_ptsx, circ_ptsy)).T  # reorganize points of the circle as [x,y] pairs
crop_path = Path(crop_pts, closed=True) # define a path for finding voronoi corners outside of region

crop_segs = []
plt.figure()
for i in range(0,len(crop_pts)-1):
    pt1 = crop_pts[i]
    pt2 = crop_pts[i+1]
    xpts=[pt1[0],pt2[0]]
    ypts=[pt1[1],pt2[1]]
    crop_segs.append([pt1,pt2])
    plt.plot(xpts,ypts)
crop_segs = np.asarray(crop_segs)
plt.axis('equal')
plt.show()

# rejiggering the Voronoi diagram: adding in missing -1 vertices to ridges and regions

i = -1   # index

ptp_bound = cells.points.ptp(axis=0)
vor_center = cells.points.mean(axis=0)

for pnt_indx, vor_edge in zip(cells.ridge_points, cells.ridge_vertices):
    vor_edge = np.asarray(vor_edge)

    i = i+1 # update the count-through index

    if np.any(vor_edge < 0):
            new_edge = vor_edge[vor_edge >= 0][0]  # find the ridge vertice that's not equal to -1!
            tang = cells.points[pnt_indx[1]] - cells.points[pnt_indx[0]]  # calculate the tangent of two points sharing ridge
            tang /= np.linalg.norm(tang)  # make the tangent a unit vector
            norml = np.array([-tang[1], tang[0]])  # calculate the normal of the two points sharing the ridge

            midpoint = cells.points[pnt_indx].mean(axis=0)   # calculate the midpoint between the two points of the ridge
            direction = np.sign(np.dot(midpoint - vor_center, norml)) * norml
            far_point = cells.vertices[new_edge] + direction * ptp_bound.max()

            vor_ind = cells.vertices.shape[0]    # get the current size of the voronoi vertices array, this will be the n+1 index

            cells.vertices = np.vstack((cells.vertices,far_point)) # add the new point to the vertices array
            cells.ridge_vertices[i] = [new_edge,vor_ind]  # add the new index to the ridge index vector at the right spot

            j=-1 # initialize another index
            for boo in cells.regions:    # step through each polygon region
                j = j+1  # update the index
                if -1 in boo and new_edge in boo:  # if the region has the edge of interest
                    a = boo.index(-1)              # find the index in the region that is -1
                    cells.regions[j][a] = vor_ind # add in the new voronoi vertice index to the appropriate region
                    verts = cells.vertices[boo]   # get the vertices for this region
                    region = np.asarray(boo)      # convert boo to a numpy array so it can be sorted
                    cent = verts.mean(axis=0)     # calculate the centre point
                    angles = np.arctan2(verts[:,1]-cent[1], verts[:,0] - cent[0])  # calculate the angle of each point
                    cells.regions[j] = region[np.argsort(angles)]   # sort the region indices counter-clockwise


plt.plot(cells.vertices[:,0],cells.vertices[:,1],'.')
plt.plot(crop_pts[:,0],crop_pts[:,1],'b')
plt.axis('equal')
plt.show()

# Sutherland-Hodgman polygon clipping function:
def clip(subjectPolygon, clipPolygon):
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

plt.figure()
# clipping the voronoi diagram to the crop circle
new_vor_verts = []
i=-1 # counting index

for poly_ind in cells.regions:  # step through each cell's polygonal regions...
    i = i+1                     # update the counting index
    cell_poly = cells.vertices[poly_ind]  # get the coordinates of the polygon vertices for the cell of interest
    cell_polya = cell_poly.tolist()  # convert data structures to python lists for cropping algorithm...
    crop_ptsa = crop_pts.tolist()

    if len(cell_polya)>=3:                     # if the polygon region has at least 3 vertices
        aa=clip(cell_polya,crop_ptsa)        # then send it to the clipping algorithm
        new_vor_verts.append(clip(cell_polya,crop_ptsa))     # append the points to the new vertices list
        plt.fill(*zip(*aa),alpha=0.6)                                  # plot the resulting polygon


#plt.plot(cells.points[:,0],cells.points[:,1],'ko')
plt.axis('equal')
plt.show()


mycmap = plt.get_cmap('jet')
# colorize
for edge in cells.ridge_vertices:
    if not -1 in edge:
        vedge = cells.vertices[edge]
        xedge= [vedge[0][0],vedge[1][0]]
        yedge = [vedge[0][1],vedge[1][1]]
        cscale=vedge[0][0]/200e-6
        colr=mycmap(cscale)
        plt.plot(xedge,yedge,color=colr)

plt.axis('equal')
plt.show()
