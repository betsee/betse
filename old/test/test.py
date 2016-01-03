#!/usr/bin/env python3
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# This is a training module

# Questions



import numpy as np
import scipy as sp
import scipy.spatial as sps
import matplotlib.pyplot as plt
from matplotlib.path import Path
import math
import time
from pprint import pprint


start_time = time.time()  # get a start value for timing the simulation

# define the basic class that holds variables
class Parameters(object):
    pass

# Define numerical constants for world set-up and simulation

const = Parameters()  # define the object that holds main parameters
const.wsx = 200e-6  # the x-dimension of the world space [m]
const.wsy = 200e-6  # the y-dimension of the world space [m]
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

for new_edge in range(0, len(y_v)):  # indices to step through grid
    for j in range(0, len(x_v)):  # indices to step through grid

         if clust_x2d[new_edge,j] is not np.ma.masked:  # if point is not masked (i.e. in the cell cluster)...
                aa=[x_2d[new_edge,j],y_2d[new_edge, j]]        # get the value of the x,y point by accessing x-grid and y-grid
                clust_xy = np.vstack((clust_xy,aa))  # augment the points list by adding in the new value

clust_xy = np.delete(clust_xy, 0, 0)    # delete the initialization value.


# Step 3: Computational Geometry data structures

# Delaunay Triangulation:
trimesh = sps.Delaunay(clust_xy)

# Convex Hull
bound_pts = sps.ConvexHull(clust_xy)

# Voronoi Dataset
cells = sps.Voronoi(clust_xy)    # TODO figure out how to do a bounded voronoi diagram

                                # idea for bounded Voronoi is to have a circle with a radius about 10 um larger than the masking circle
                                # which will have vertices where the unbounded voronoi edges intersect it....?

# Nearest Neighbour Algorithm Tree for points in cell cluster
tree = sps.KDTree(clust_xy)  # this tree is used to find nearest neighbours for a point in the cell cluster

tree_um = sps.KDTree(xypts)   # this tree is used to get indices of unmasked vector for mapping back to 2d grids

# function calculating area of a polygon

def area(poly):
    return 0.5 * abs(sum(x0*y1 - x1*y0
                         for ((x0, y0), (x1, y1)) in segments(poly)))

def segments(poly):
    return zip(poly, poly[1:] + [poly[0]])

# defining function to find intersection of two lines

def line_intersect(line1,line2):
    """
    Determine whether two 2D lines intersect
    and return their point of intersection

    Parameters
    ----------
    line1 : 2x2 array of line points as [[x1,y1],[x2,y2]]
    line2 : 2x2 array of line points as [[x3,y3],[x4,y4]]

    Returns
    -------
    int_ptn : [x,y] co-ordinates of the line intersection
    or None for no intersection

    Notes
    -------
    Uses Numpy arrays

    """
    x1 = line1[0][0]
    y1 = line1[0][1]
    x2 = line1[1][0]
    y2 = line1[1][1]

    x3 = line2[0][0]
    y3 = line2[0][1]
    x4 = line2[1][0]
    y4 = line2[1][1]

    m1 = (y2 - y1)/(x2 - x1)   # slope of first line
    m2 = (y4 - y3)/(x4 - x3)   # slope of second line

    b1 = y1 - m1*x1   # intercept of first line
    b2 = y3 - m2*x3   # intercept of second line

    deno = m1 - m2  # the denomenator of intersection; test if lines are parallel
    if deno == 0:
        int_ptn = None
    else:
        intx = (b2 - b1)/deno
        inty = (m1*b2 - m2*b1)/deno
        int_ptn = [intx,inty]

    return int_ptn

# Sutherland-Hodgman polygon clipping algorithm

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






# Rejiggering the Voronoi diagram -- for pretty plotting
plt.figure()

for vor_edge in cells.ridge_vertices:
        vor_edge = np.asarray(vor_edge)
        if np.all(vor_edge >= 0):
            plt.plot(cells.vertices[vor_edge,0], cells.vertices[vor_edge,1], 'k-')

ptp_bound = cells.points.ptp(axis=0)
vor_center = cells.points.mean(axis=0)

for pnt_indx, vor_edge in zip(cells.ridge_points, cells.ridge_vertices):
    vor_edge = np.asarray(vor_edge)
    if np.any(vor_edge < 0):
            new_edge = vor_edge[vor_edge >= 0][0]  # finite end Voronoi vertex

            tang = cells.points[pnt_indx[1]] - cells.points[pnt_indx[0]]  # tangent
            tang /= np.linalg.norm(tang)
            norml = np.array([-tang[1], tang[0]])  # normal

            midpoint = cells.points[pnt_indx].mean(axis=0)
            direction = np.sign(np.dot(midpoint - vor_center, norml)) * norml
            far_point = cells.vertices[new_edge] + direction * ptp_bound.max()

            plt.plot([cells.vertices[new_edge,0], far_point[0]],
                    [cells.vertices[new_edge,1], far_point[1]], 'k--')

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

# ****************************************************************************

# clipping the voronoi diagram to the crop circle

plt.figure()

for cell_centre, poly_ind in zip(cells.points, cells.regions):  # step through the cell and polygon simultaneously...

    cell_poly = cells.vertices[poly_ind]  # get the coordinates of the polygon vertices for the cell of interest
    cell_polya = cell_poly.tolist()
    crop_ptsa = crop_pts.tolist()

    if len(cell_polya):
        newpts = clip(cell_polya,crop_ptsa)
        plt.plot(*zip(*newpts))

plt.plot(cells.points[:,0],cells.points[:,1],'ko')
plt.axis('equal')
plt.show()


# plt.figure()
# # colorize
# for region in cell_regions:
#     polygon = cell_vertices[region]
#     plt.fill(*zip(*polygon), alpha=0.4)
#
# plt.plot(clust_xy[:,0], clust_xy[:,1], 'ko')
# plt.axis('equal')
# plt.title('Cell Network: Neighbour and Boundary Detection')
# plt.show()



plt.figure()
# colorize
for region in cells.regions:
    if not -1 in region:
        polygon = cells.vertices[region]
        plt.fill(*zip(*polygon), alpha=0.4)

plt.plot(clust_xy[:,0], clust_xy[:,1], 'ko')
plt.axis('equal')
plt.title('Cell Network: Neighbour and Boundary Detection')
plt.show()

# plt.figure()
# # colorize
# # for region in cells.regions:
# #     if not -1 in region:
# #         polygon = cells.vertices[region]
# #         plt.fill(*zip(*polygon), alpha=0.4)
# plt.plot(cells.vertices[:,0], cells.vertices[:,1], 'bo')
# plt.plot(clust_xy[:,0], clust_xy[:,1], 'ko')
# plt.axis('equal')
# plt.title('Cell Network: Neighbour and Boundary Detection')
# plt.show()


# Step 4: Do a computation and tag cell points with that "data" -- it will be visualized by colour

z_clust = np.zeros(clust_xy.shape)   # initialize the array that will hold data specific to the cell cluster seeds.
xi = clust_xy[:,0]  # define the xi and yi as vectors to keep things easy on the brain
yi = clust_xy[:,1]
z_clust[:,0] = 1e9*(xi**2) + 1e9*(yi**2)   # define z-data as some spatially conditional function -- later it will be nn operation


# will need a function to map data from points list to 2d-grid
# first work out the steps, later write the mapping function
_,z_ind = tree_um.query(clust_xy)
z_xy[z_ind] = z_clust   # map the masked data to the full vector
z_xyi = z_xy[:,0]   # get the column that holds the data
z_2d = z_xyi.reshape(x_2d.shape)   # reshape the matrix to a 2d-grid
clust_z2d = np.ma.masked_array(z_2d,~ws_mask)  # create a masked array for the z-data structure (colour)

# Visualize results

rs = 1e6

plt.figure(1)
plt.triplot(rs*clust_xy[:, 0], rs*clust_xy[:, 1], trimesh.simplices)
plt.plot(rs*clust_xy[bound_pts.vertices, 0], rs*clust_xy[bound_pts.vertices, 1], 'r', lw=2)  # FIXME want a closed path for the convex hull!!!
plt.scatter(rs*clust_x2d, rs*clust_y2d)
plt.title('Cell Network: Neighbour and Boundary Detection')
plt.axis('equal')
#plt.axis([wsxmin, wsxmax, wsymin, wsymax])
plt.show(block=False)


# THIS IS AN AMAZING KEY FIGURE!!!
# Don't need to go back and forth to the grid -- work with data cluster directly.
plt.figure(2)
plt.tripcolor(rs*clust_xy[:, 0], rs*clust_xy[:, 1], z_clust[:,0],shading='gouraud')
#plt.triplot(clust_xy[:, 0], clust_xy[:, 1], trimesh.simplices)
#plt.scatter(clust_x2d,clust_y2d,c=clust_z2d,s=100)
cb=plt.colorbar()
cb.set_label('Fake data values [Foo Units]')
plt.title('Cell Network: Smoothed Node Data')
plt.ylabel('Spatial y [um]')
plt.xlabel('Spatial x [um]')
plt.axis('equal')
#plt.axis([wsxmin, wsxmax, wsymin, wsymax])
plt.show(block=False)

plt.figure(3)
plt.tripcolor(rs*clust_xy[:, 0], rs*clust_xy[:, 1], z_clust[:,0],shading='gouraud')
plt.triplot(rs*clust_xy[:, 0], rs*clust_xy[:, 1], trimesh.simplices)
#plt.scatter(clust_x2d,clust_y2d,c=clust_z2d,s=100)
cb=plt.colorbar()
cb.set_label('Fake data values [Foo Units]')
plt.title('Cell Network: Connections with Smoothed Node Data')
plt.ylabel('Spatial dimension y [m]')
plt.xlabel('Spatial dimension x [m]')
plt.axis('equal')
#plt.axis([wsxmin, wsxmax, wsymin, wsymax])
plt.show(block=False)


plt.figure(4)
#plt.tripcolor(clust_xy[:, 0], clust_xy[:, 1], z_clust[:,0],shading='gouraud')
plt.triplot(rs*clust_xy[:, 0], rs*clust_xy[:, 1], trimesh.simplices)
plt.scatter(rs*clust_x2d,rs*clust_y2d,c=clust_z2d,s=100)
cb=plt.colorbar()
cb.set_label('Fake Data Values [Foo Units]')
plt.title('Cell Network: Data on Nodes')
plt.ylabel('Spatial y [um]')
plt.xlabel('Spatial x [um]')
plt.axis('equal')
#plt.axis([wsxmin, wsxmax, wsymin, wsymax])
plt.show(block=False)


# Voronoi diagram

# plotting edges of Voronoi with data on edges

# mycmap = plt.get_cmap('jet')
# # colorize
# for edge in cells.ridge_vertices:
#     if not -1 in edge:
#         vedge = cells.vertices[edge]
#         xedge= [vedge[0][0],vedge[1][0]]
#         yedge = [vedge[0][1],vedge[1][1]]
#         cscale=vedge[0][0]/200e-6
#         colr=mycmap(cscale)
#         plt.plot(xedge,yedge,color=colr)
#
# plt.axis('equal')
# plt.show()

#sps.voronoi_plot_2d(cells)
# plt.title('Voronoi Representation of Cell Network')
# plt.ylabel('Spatial y [m]')
# plt.xlabel('Spatial x [m]')
# plt.axis('equal')
# #plt.axis([wsxmin, wsxmax, wsymin, wsymax])
# plt.show(block=False)

print('The simulation took', time.time() - start_time, 'seconds to complete')

xa=np.linspace(0,100,1000)
ya=np.sin(xa)

plt.figure(5)
plt.plot(xa,ya)
plt.ylim([-1.5,1.5])
plt.title('Some Foo Data')
plt.ylabel('Foo Data [units]')
plt.xlabel('Distance from Foo [m]')
plt.show(block=False)

plt.show()

# This is how to plot a collection of line segments that are coloured to a value:
mycmap = plt.get_cmap('jet')  # need to call up a colormap
# colorize
for edge in cells.ridge_vertices:   # looping through each edge instance in the voronoi diagram
    if not -1 in edge:   # if it's not out of bounds
        vedge = cells.vertices[edge]         # get the array object listing the xy points of the edges
        xedge= [vedge[0][0],vedge[1][0]]
        yedge = [vedge[0][1],vedge[1][1]]
        cscale=vedge[0][0]/200e-6
        colr=mycmap(cscale)
        plt.plot(xedge,yedge,color=colr)

plt.axis('equal')
plt.show()


# Wastelands



#StreamplotExample
#
# Y, X = np.mgrid[-3:3:100j, -3:3:100j]
# U = -1 - X**2 + Y
# V = 1 + X - Y**2
# speed = np.sqrt(U*U + V*V)
#
# plt.streamplot(X, Y, U, V, color=U, linewidth=2, cmap=plt.cm.autumn)
# plt.colorbar()
#
# f, (ax1, ax2) = plt.subplots(ncols=2)
# ax1.streamplot(X, Y, U, V, density=[0.5, 1])
#
# lw = 5*speed/speed.max()
# ax2.streamplot(X, Y, U, V, density=0.6, color='k', linewidth=lw)
#
# plt.show(block=False)

# can also do:
#ws_mask=ws_mask.astype(int)  # defines the mask as 0 or 1
#clust_x=x_2d*ws_mask
#clust_y=y_2d*ws_mask

# indexing array element:
#vv=clust_x
#for i in range(10):
#    for j in range (10):
#        testmask=clust_x[i,j] is ma.masked
#        if testmask==False:
#           vv[i,j]=22

# Generally speaking, this is usually the best way to print the contents of an object: "pretty print" the variables of
# such object. ("Pretty print" is compsci slang for "print this thing in a well-formatted, readable manner.")
# pprint(vars(prms))

# This also works, of course! In general, it's best to avoid accessing object internals (that is, "__"-prefixed fields
# and methods) from outside that object. Consider calling the var() builtin, instead.
#parameters=prms.__dict__

#X = np.linspace(-10, 10, 250, endpoint=True)
#C = np.cos(X)

#plt.plot(X, C)
#plt.show()

# plot the scatter of points
# plt.figure(1)
# plt.scatter(x_2d, y_2d)
# plt.axis('equal')
# plt.axis([wsxmin, wsxmax, wsymin, wsymax])
# plt.title('Seed points for cell cluster')
# plt.xlabel('x-dimension [m]')
# plt.ylabel('y-dimension [m]')
# plt.show(block=False)
#
# plt.figure(2)
# plt.scatter(clust_x2d, clust_y2d)
# plt.plot(circ_ptsx, circ_ptsy, 'r', linewidth=5.0)
# plt.title('Seed points with cluster boundary, after masking')
# plt.axis('equal')
# plt.axis([wsxmin, wsxmax, wsymin, wsymax])
# plt.show(block = False)

#sps.voronoi_plot_2d(cells)
#plt.title('Voronoi Diagram of Cell Network')
#plt.axis('equal')
#plt.axis([wsxmin, wsxmax, wsymin, wsymax])
#plt.show(block=False)
#
# def voronoi_finite_polygons_2d(vor, radius=None):
#     """
#     Reconstruct infinite voronoi regions in a 2D diagram to finite
#     regions.
#
#     Parameters
#     ----------
#     vor : Voronoi
#         Input diagram
#     radius : float, optional
#         Distance to 'points at infinity'.
#
#     Returns
#     -------
#     regions : list of tuples
#         Indices of vertices in each revised Voronoi regions.
#     vertices : list of tuples
#         Coordinates for revised Voronoi vertices. Same as coordinates
#         of input vertices, with 'points at infinity' appended to the
#         end.
#
#     """
#
#     if vor.points.shape[1] != 2:
#         raise ValueError("Requires 2D input")
#
#     new_regions = []
#     new_vertices = vor.vertices.tolist()
#
#     center = vor.points.mean(axis=0)
#     if radius is None:
#         radius = vor.points.ptp().max()*2
#
#     # Construct a map containing all ridges for a given point
#     all_ridges = {}
#     for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
#         all_ridges.setdefault(p1, []).append((p2, v1, v2))
#         all_ridges.setdefault(p2, []).append((p1, v1, v2))
#
#     # Reconstruct infinite regions
#     for p1, region in enumerate(vor.point_region):
#         vertices = vor.regions[region]
#
#         if all([v >= 0 for v in vertices]):
#             # finite region
#             new_regions.append(vertices)
#             continue
#
#         # reconstruct a non-finite region
#         ridges = all_ridges[p1]
#         new_region = [v for v in vertices if v >= 0]
#
#         for p2, v1, v2 in ridges:
#             if v2 < 0:
#                 v1, v2 = v2, v1
#             if v1 >= 0:
#                 # finite ridge: already in the region
#                 continue
#
#             # Compute the missing endpoint of an infinite ridge
#
#             t = vor.points[p2] - vor.points[p1] # tangent
#             t /= np.linalg.norm(t)
#             n = np.array([-t[1], t[0]])  # normal
#
#             midpoint = vor.points[[p1, p2]].mean(axis=0)
#             direction = np.sign(np.dot(midpoint - center, n)) * n
#             far_point = vor.vertices[v2] + direction * radius
#
#             new_region.append(len(new_vertices))
#             new_vertices.append(far_point.tolist())
#
#         # sort region counterclockwise
#         vs = np.asarray([new_vertices[v] for v in new_region])
#         c = vs.mean(axis=0)
#         angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
#         new_region = np.array(new_region)[np.argsort(angles)]
#
#         # finish
#         new_regions.append(new_region.tolist())
#
#     return new_regions, np.asarray(new_vertices)
#
#
# cell_regions, cell_vertices = voronoi_finite_polygons_2d(cells, r_circ)
#
# print("--")
# print(cells.vertices)
# print("--")
# print(cells.regions)