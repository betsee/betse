#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME add in assertions for various methods
# FIXME redesign for non-equal x and y world dimensions
# FIXME method to create a hexagonal lattice base
# FIXME figure out how to scale each voronoi polygon to um instead of m dimensions when plotting
# FIXME need nx,ny (normal) and tx,ty (tangent) to each cell edge
# FIXME need boundary flags for cell_mids and ecm_verts
# FIXME documentation (especially alpha shape and concave hull) not up to date!
# FIXME allow user to specify their own set of points for clipping in points and voronoi clips (make circle funciton)
# FIXME update the concave hull to be Sess' faster algorithm

"""
The world module contains the class World, which holds
all data structures relating to the size of the environment,
the extent of the cell cluster, and the co-ordinates of cell
centre points.

The initialization method of the World class sets-up
and crops the cell cluster to an optional user-defined geometry input
as a set of points arranged in counter-clockwise order and
defining a closed polygon. Other methods define the cell centres of each
cell polygon, their area, and create nearest neighbour and edge matrices
for each cell.

"""

import numpy as np
import scipy.spatial as sps
from matplotlib.path import Path
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.collections import PolyCollection
import matplotlib.cm as cm
import math


class World(object):
    """
    The WorldSeeds object stores data structures relating to
    the geometric properties of the environmental grid and cell
    centre points.

    Parameters
    ----------
    WorldSeeds requires an instance of NumVars, see the Parameters module.
    Optional: crop_mask (default = None) a set of counter-clockwise
    arranged points defining a closed polygon to clip the cluster of
    cell seed points

    Fields
    -------
    self.x_v        linear numpy vectors as 'tick marks' of x and y
    self.y_v        axes

    self.x_2d       numpy 2d arrays of x or y gridded points
    self.y_2d       (irregular lattice)

    self.xmin, self.xmax      dimensions of world grid (after noise, before cropping)
    self.ymin, self.ymax

    self.centre     [x,y] coordinate of world lattice co-ords (after noise, before cropping)

    self.xypts      numpy array holding [x,y] points of irregular world grid

    self.clust_xy   numpy array of [x,y] points corresponding to only those of
                    the cropped cell cluster

    self.cluster_axis       maximum breadth and centre (as [x,y] point) of cropped cell cluster
    self.cluster_center

    self.ecm_verts     a nested python list containing Voronoi cell regions which specify [x,y] points of region
                        vertices for a clipped Voronoi diagram (arranged to cell index)

    self.cell_area     a list of areas of each Voronoi polygon (arranged to cell index)

    self.cell_centres    a list of [x,y] points defining the cell centre (arranged to cell index)

    self.cell_nn         a nested array of integer indices of each nearest neighbour for a particular cell (arranged
                        to cell index)

    self.bound_flags    a numpy array flagging cells on the envirnomental boundary with 1

    self.ecm_edges     a python nested listing of two point line segments of each membrane domain for each cell


    Methods
    -------
    makeSeeds()                             Create an irregular lattice of seed points in 2d space
    cropSeeds(crop_mask)                    Crop the points cluster to a polygonal shape (circle)
    makeVoronoi(vorclose = None)            Make and clip/close a Voronoi diagram from the seed points
    clip(subjectPolygon, clipPolygon)       The Sutherland-Hodgman polygon clipping algorithm
    area(p)                                 Find the area of a polygon defined by a list of [x,y] points
    vor_area()                              Returns the area of each polygon in the closed Voronoi diagram
    cell_index()                            Returns a list of [x,y] points defining the cell centres in order
    near_neigh()                            Calculate the nearest neighbour (nn) array for each cell
    cellEdges()                             List of membrane domains as two-point line segments for each cell
    plotPolyData(zdata = None,clrmap = None)  Plot cell polygons with data attached

    Notes
    -------
    Uses Numpy
    Uses Scipy spatial

    """

    def __init__(self,constants,crop_mask=None, vorclose=None):

        self.vorclose = vorclose
        self.crop_mask = crop_mask

        self.d_cell = constants.dc  # diameter of single cell
        self.nx = constants.nx   # number of lattice sites in world x index
        self.ny = constants.ny   # number of lattice sites in world y index
        self.ac = constants.ac  # cell-cell separation
        self.nl = constants.nl  # noise level for the lattice
        self.wsx = constants.wsx  # World size
        self.wsy = constants.wsy # World size
        self.search_d = constants.search_d  # distance to search for nearest neighbours (relative to d_cell)
        self.sf = constants.sf              # scale factor to take cell vertices in from extracellular space
        self.cell_sides = constants.cell_sides # minimum number of membrane domains per cell


        self.makeSeeds()    # Create the grid for the system (irregular)
        self.cropSeeds(self.crop_mask)      # Crop the grid to a geometric shape to define the cell cluster
        self.makeVoronoi(self.vorclose)    # Make, close, and clip the Voronoi diagram
        self.cell_index()            # Calculate the correct centre and index for each cell
        self.cellVerts()   # create individual cell polygon vertices
        self.vor_area()              # Calculate the area of each cell polygon
       # self.cell_index()            # Calculate the correct centre and index for each cell
        self.near_neigh()    # Calculate the nn array for each cell
        self.cellEdges()    # create a nested list of all membrane and ecm domains for each cell
        self.cellMids()    # calculate the midpoint of membrane domains and ecm segments
        self.boundTag()   # flag cells, ecm, and mem domains on the environmental boundary


    def makeSeeds(self):

        """
        makeSeeds returns an irregular scatter
        of points defined on a world space
        with dimensions wsx, wsy in [m].

        The amount of deviation from a square
        grid is specified by nl, defined from
        0 (perfect square grid) to 1 (full noise).


        Creates
        -------
        self.xypts      numpy array listing [x,y] of world seed points

        self.x_2d       numpy 2d arrays of x or y gridded points
        self.y_2d       (irregular lattice)

        self.x_v        linear numpy vectors as 'tick marks' of x and y
        self.y_v        axes

        self.xmin, self.xmax      dimensions of world grid (after noise, before cropping)
        self.ymin, self.ymax

        self.centre     [x,y] coordinate of world centre (after noise, before cropping)

        Notes
        -------
        Uses Numpy arrays

        """

        # first begin with linear vectors which are the "ticks" of the x and y dimensions
        self.x_v = np.linspace(0, (self.nx - 1) * (self.d_cell + self.ac), self.nx)  # create lattice vector x
        self.y_v = np.linspace(0, (self.ny - 1) * (self.d_cell + self.ac), self.ny)  # create lattice vector y

        # next define a 2d array of lattice points using the x- and y- vectors
        x_2d, y_2d = np.meshgrid(self.x_v, self.y_v)  # create 2D array of lattice points

        # now create a matrix of points that will add a +/- deviation to each point centre
        x_rnd = self.nl * self.d_cell * (np.random.rand(self.ny, self.nx) - 0.5)  # create a mix of random deltas x dir
        y_rnd = self.nl * self.d_cell * (np.random.rand(self.ny, self.nx) - 0.5)  # create a mix of random deltas x dir

        # add the noise effect to the world point matrices and redefine the results
        self.x_2d = x_2d + x_rnd
        self.y_2d = y_2d + y_rnd

        # define a data structure that holds [x,y] coordinate points of each 2d grid-matrix entry
        self.xypts = np.vstack((self.x_2d.ravel(), self.y_2d.ravel())).T

        # define geometric limits and centre for the cluster of points
        self.xmin = np.min(self.x_2d)
        self.xmax = np.max(self.x_2d)
        self.ymin = np.min(self.y_2d)
        self.ymax = np.max(self.y_2d)

        self.centre = self.xypts.mean(axis=0)

    def cropSeeds(self, crop_mask):

        """
        cropSeeds returns a geometrically
        cropped version of an irregular points scatter in 2D.

        The option crop_mask specifies the type of cropping where
        crop_mask=None gives no cropping and crop_mask='circle' crops
        to a circle with the diameter of the points scatter.

        Parameters
        ----------
        crop_mask          None = no cropping, 'circle'= crop to circle


        Creates
        -------
        self.clust_xy      an array listing [x,y] points of each cell seed
                            in the cropped cluster
        Notes
        -------
        Uses Numpy arrays

        """

        if crop_mask==None:  # if there's no crop-mask specified (default)
            self.clust_xy=self.xypts         # set cluster points to grid points

        elif crop_mask =='circle': # if 'circle' is specified:

            cres = 15  # how many points desired in polygon
            d_circ = self.xmax - self.xmin  # diameter of circle in x-direction
            r_circ = d_circ / 2  # radius of circle
            ind1 = np.linspace(0, 1, cres + 1)  # indices of angles defining circle points

            angs = ind1 * 360 * (np.pi / 180)  # angles in radians defining circle points
            circ_ptsx = r_circ * np.cos(angs) + self.centre[0]  # points of the circle
            circ_ptsy = r_circ * np.sin(angs) + self.centre[1]  # points of the circle

            crop_pts = np.vstack((circ_ptsx, circ_ptsy)).T  # reorganize points of the circle as [x,y] pairs

            crop_path = Path(crop_pts, closed=True)  # transform cropping points to a functional path
            # create a boolean matrix mask which is 1 for points inside the circle and 0 for points outside
            ws_mask_xy = crop_path.contains_points(self.xypts)  # create the mask for point inside the path
            ws_mask = ws_mask_xy.reshape(self.x_2d.shape)  # reshape the boolean mask to correspond to the data grid
            self.clust_x2d = np.ma.masked_array(self.x_2d, ~ws_mask)  # created a masked data structure of the x-grid
            self.clust_y2d = np.ma.masked_array(self.y_2d, ~ws_mask)  # create a masked data structure of the y-grid

            # finally, create a data structure of [x,y] points that only contains the points of the cluster

            self.clust_xy = np.array([0,0])  # initialize the x,y points array

            for new_edge in range(0, len(self.y_v)):  # indices to step through grid
                for j in range(0, len(self.x_v)):  # indices to step through grid
                    # if point is not masked (i.e. in the cell cluster)...
                     if self.clust_x2d[new_edge,j] is not np.ma.masked:
                          # get the value of the x,y point by accessing x-grid and y-grid
                            aa=[self.x_2d[new_edge,j],self.y_2d[new_edge, j]]
                           # augment the points list by adding in the new value
                            self.clust_xy = np.vstack((self.clust_xy,aa))

            self.clust_xy = np.delete(self.clust_xy, 0, 0)    # delete the initialization value.

    def makeVoronoi(self, vorclose = None):

        """
        makeVoronoi calculates the Voronoi diagram for an input
        of cell seed points.

        The option vorclose specifies the Voronoi diagram to be clipped (closed)
        to a polygon (circle) corresponding to the seed cluster maximum breadth.
        If vorclose=None, there is no cropping, while vorclose = 'circle' crops
        to a 15 point polygon (circle).

        Parameters
        ----------
        clust_xy            a numpy array listing [x,y] of seed points
        vorclose            None = no cropping, 'circle'= crop to circle

        Returns
        -------
        self.cells          a scipy spatial Voronoi diagram instance
        self.ecm_verts      nested python list specifying polygonal region
                            and vertices as [x,y] for each Voronoi cell in the
                            clipped/closed Voronoi diagram.

        Notes
        -------
        Uses Numpy arrays
        Uses Scipy spatial

        Uses clip(subjectPolygon, clipPolygon) a method defining
        the Sutherland-Hodgman polygon clipping algorithm -- courtesy of Rosetta Code:
        http://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#Python

        """

        self.vor = sps.Voronoi(self.clust_xy)

        self.cluster_axis = self.vor.points.ptp(axis=0)
        self.cluster_center = self.vor.points.mean(axis=0)

        # complete the Voronoi diagram by adding in undefined vertices to ridges and regions
        i = -1   # enumeration index

        for pnt_indx, vor_edge in zip(self.vor.ridge_points, self.vor.ridge_vertices):
            vor_edge = np.asarray(vor_edge)

            i = i+1 # update the count-through index

            if np.any(vor_edge < 0): # if either of the two ridge values are undefined (-1)

                # find the ridge vertice that's not equal to -1
                    new_edge = vor_edge[vor_edge >= 0][0]
                # calculate the tangent of two seed points sharing that ridge
                    tang = self.vor.points[pnt_indx[1]] - self.vor.points[pnt_indx[0]]
                    tang /= np.linalg.norm(tang)  # make the tangent a unit vector
                    norml = np.array([-tang[1], tang[0]])  # calculate the normal of the two points sharing the ridge

                    # calculate the midpoint between the two points of the ridge
                    midpoint = self.vor.points[pnt_indx].mean(axis=0)
                    # now there's enough information to calculate the missing direction and location of missing point
                    direction = np.sign(np.dot(midpoint - self.cluster_center, norml)) * norml
                    #far_point = self.vor.vertices[new_edge] + direction * self.cluster_axis.max()
                    far_point = self.vor.vertices[new_edge] + direction * self.d_cell

                    # get the current size of the voronoi vertices array, this will be the n+1 index after adding point
                    vor_ind = self.vor.vertices.shape[0]

                    self.vor.vertices = np.vstack((self.vor.vertices,far_point)) # add the new point to the vertices array
                    self.vor.ridge_vertices[i] = [new_edge,vor_ind]  # add the new index at the right spot

                    for j, region in enumerate(self.vor.regions):    # step through each polygon region

                        if len(region):

                            if -1 in region and new_edge in region:  # if the region has edge of interest...
                                a = region.index(-1)              # find index in the region that is undefined (-1)
                                self.vor.regions[j][a] = vor_ind # add in the new vertex index to the appropriate region

                            verts = self.vor.vertices[region]   # get the vertices for this region
                            region = np.asarray(region)      # convert region to a numpy array so it can be sorted
                            cent = verts.mean(axis=0)     # calculate the centre point
                            angles = np.arctan2(verts[:,1]-cent[1], verts[:,0] - cent[0])  # calculate point angles
                            #self.vor.regions[j] = region[np.argsort(angles)]   # sort indices counter-clockwise
                            sorted_region = region[np.argsort(angles)]   # sort indices counter-clockwise
                            sorted_region_b = sorted_region.tolist()
                            self.vor.regions[j] = sorted_region_b   # add sorted list to the regions structure


        # finally, clip the Voronoi diagram to polygon, if user-specified by vorclose option
        if vorclose==None:
            self.ecm_verts=[]
            for region in self.vor.regions:
                if len(region):
                    cell_poly = self.vor.vertices[region]
                    if len(cell_poly)>3:
                        self.ecm_verts.append(self.vor.vertices[region])


        elif vorclose=='circle':
            #cluster_axis = vor.points.ptp(axis=0)    # calculate the extent of the cell points
            #centx = vor.points.mean(axis=0)       # calculate the centre of the cell points

            cres = 15  # how many points desired in cropping polygon
            #d_circ = cluster_axis.max()  # diameter of cropping polygon
            d_circ = self.xmax - self.xmin
            r_circ = 1.01*(d_circ / 2)  # radius of cropping polygon
            ind1 = np.linspace(0, 1, cres + 1)  # indices of angles defining polygon points
            angs = ind1 * 360 * (np.pi / 180)  # angles in radians defining polygon points
            #circ_ptsx = r_circ * np.cos(angs) + centx[0]  # points of the polygon
            #circ_ptsy = r_circ * np.sin(angs) + centx[1]  # points of the polygon
            circ_ptsx = r_circ * np.cos(angs) + self.centre[0]  # points of the polygon
            circ_ptsy = r_circ * np.sin(angs) + self.centre[1]  # points of the polygon

            self.crop_pts = np.vstack((circ_ptsx, circ_ptsy)).T  # reorganize polygon points as [x,y] pairs
            crop_path = Path(self.crop_pts, closed=True)  # transform cropping points to a functional path

            crop_ptsa = self.crop_pts.tolist()   # a python list version to use with the clipping algorithm

            # Now clip the voronoi diagram to the cropping polygon
            self.ecm_verts = []

            for poly_ind in self.vor.regions:  # step through each cell's polygonal regions...

                if len(poly_ind) >= self.cell_sides: # check to make sure we're defining a polygon

                    cell_poly = self.vor.vertices[poly_ind]  # get the coordinates of the polygon vertices

                    inpath = crop_path.contains_points(cell_poly)    # get a boolean matrix

                    if inpath.all() == False:
                        pass

                    else:
                        cell_polya = cell_poly.tolist()  # convert data structures to python lists for cropping algorithm...
                        aa=self.clip(cell_polya,crop_ptsa)        # then send it to the clipping algorithm
                        if len(aa) >= self.cell_sides:                        # check to make sure result is still a polygon
                            self.ecm_verts.append(aa)     # append points to new region point list

    def clip(self, subjectPolygon, clipPolygon):  # This is the Sutherland-Hodgman polygon clipping algorithm

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

    def area(self, p):

        """
        Calculates the area of an arbitrary polygon defined by a
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

    def vor_area(self):

        """
        Calculates the area of each cell in a closed 2D Voronoi diagram.

        Parameters
        ----------
        ecm_verts               nested list of [x,y] points defining each polygon


        Returns
        -------
        self.v_area            area of all polygons of the Voronoi diagram in square meters

        Notes
        -------
        Uses area(p) function. The Voronoi diagram must be closed (no outer edges extending to infinity!!!)

        """
        self.cell_area = []
        for poly in self.cell_verts:
            self.cell_area.append(self.area(poly))

    def cell_index(self):

        """
        Calculate the cell centre for each voronoi polygon and return a list
        with an index consistent with all other data lists for the cell cluster.

        Parameters
        ----------
        vor_verts               nested list of [x,y] points defining each closed polygon in the
                                Voronoi diagram

        Returns
        -------
        self.cell_centres      [x,y] coordinate of the centre of each cell as a numpy array

        Notes
        -------
        The Voronoi diagram must be closed (no outer edges extending to infinity!!!)

        """

        self.cell_centres = np.array([0,0])

        for poly in self.ecm_verts:
            aa = np.asarray(poly)
            aa = np.mean(aa,axis=0)
            self.cell_centres = np.vstack((self.cell_centres,aa))

        self.cell_centres = np.delete(self.cell_centres, 0, 0)

    def near_neigh(self):

        """
        Calculate the nearest neighbours for each cell centre in the cluster and return a numpy
        array of nn indices with an index consistent with all other data lists for the cluster.

        Parameters
        ----------
        cell_centres            A numpy array listing the [x,y] co-ordinates of each cell in the cluster
        search_d                Value defining the search distance to find nearest neighbours (search_d > cell_d)
        cell_d                  The average diameter of a cell


        Returns
        -------
        self.cell_nn            A nested list defining the indices of all nearest neighbours to each cell
        self.con_segs           A nested list defining the two [x,y] points for all cell-neighbour line connections

        Notes
        -------
        Uses numpy arrays
        Uses scipy spatial KDTree search algorithm

        """

        cell_tree = sps.KDTree(self.cell_centres)
        self.cell_nn=cell_tree.query_ball_point(self.cell_centres,self.search_d*self.d_cell)

        # define a listing of all cell-neighbour line segments
        self.con_segs = []
        for centre, indices in zip(self.cell_centres,self.cell_nn):
            for pt in indices:
                self.con_segs.append([centre,self.cell_centres[pt]])

    def boundTag(self):

        """

        Flag cells that are on the boundary to the environment by calculating the convex hull
        for the cell centre points cluster.

        Parameters
        ----------
        points            A python list of the [x,y] co-ordinates of a 2D points cluster

        alpha             The alpha parameter used to calculate the concave hull -- should be equal
                            to the average spacing between points in the cluster, which is d_cell

        Returns
        -------
        bound_flags       A python list with 0 indicating cell not on boundary and 1 indicating boundary
                                cell
        Notes
        -------
        Uses numpy arrays
        Uses alpha_shape function to calculate the concave hull

        """
        # initialize the boundary flag structure
        self.bound_flag = []
        for cellvals in self.ecm_verts:
            hoo = []
            for pts in cellvals:
                hoo.append(0)
            self.bound_flag.append(hoo)

        # get a points cluster and a map back to the original indices
        self.ecm_flat,indmap_ecm = self.flatten(self.ecm_verts)
        self.ecm_flat = np.asarray(self.ecm_flat)  # convert to numpy array for plotting, etc

        self.con_hull = self.alpha_shape(self.ecm_flat, 1/self.d_cell)  # get the concave hull for the membrane midpoints

        #bmems = np.asarray(bmems)
        #bmems_f = np.unique(bmems)                    # instead of a list of edges, get list of unique points

        for indis in self.con_hull:
            for val in indis:
                org_ind = indmap_ecm[val]
                self.bound_flag[org_ind[0]][org_ind[1]]=1  # set the boundary flag in the original data format to 1

        # self.tri_edges = self.alpha_shape(self.clust_xy,1/self.d_cell)
        # self.con_hull = self.concave_hull(self.tri_edges)


    def cellEdges(self):
        """

        Flag cells that are on the boundary to the environment by calculating the convex hull
        for the cell centre points cluster.

        Parameters
        ----------
        reg_verts            A nested python list of the [x,y] co-ordinates of each cell polygon in the cluster

        Returns
        -------
        self.ecm_edges      A nested python list of the [x,y] point pairs defining line segments of each membrane
                            domain in a cell polygon. The list has segments arranged in a counterclockwise manner.
        self.cell_edges      A nested python list of the [x,y] point pairs defining line segments of each membrane
                            domain in a cell polygon. The list has segments arranged in a counterclockwise manner.

        """

        self.ecm_edges = []
        for poly in self.ecm_verts:
            edge =[]
            for i in range(0,len(poly)):
                edge.append([poly[i-1],poly[i]])

            self.ecm_edges.append(edge)

        self.cell_edges = []
        for poly in self.cell_verts:
            edge =[]
            for i in range(0,len(poly)):
                edge.append([poly[i-1],poly[i]])

            self.cell_edges.append(edge)

    def cellMids(self):

        self.mem_mids = []

        for edges in self.cell_edges:
            hoo = []
            for edge in edges:
                pt1 = edge[0]
                pt2 = edge[1]
                mid = (pt1 + pt2)/2
                hoo.append(mid)
            self.mem_mids.append(hoo)

        self.ecm_mids = np.array([0,0])
        for edges in self.ecm_edges:
            for edge in edges:
                pt1 = edge[0]
                pt2 = edge[1]
                pt1 = np.asarray(pt1)
                pt2= np.asarray(pt2)
                mid = (pt1 + pt2)/2
                self.ecm_mids = np.vstack((self.ecm_mids,mid))
        self.ecm_mids = np.delete(self.ecm_mids,0,0)

    def cellVerts(self):
        """
        Calculate the true vertices of each individual cell from the extracellular matrix (ecm) vertices
        of the closed & clipped Voronoi diagram.

        Parameters
        ----------
        sf                   The scale factor by which the vertices are dilated (must be less than 1.0!)

        Creates
        -------
        self.cell_verts      A nested python list of the [x,y] point pairs defining vertices of each individual cell
                            polygon. The points of each polygon are arranged in a counterclockwise manner.

        Notes
        -------


        """
        self.cell_verts = []

        for centre,poly in zip(self.cell_centres,self.ecm_verts):
            pt_scale = []
            for vert in poly:
                pt_zero = vert - centre
                pt_scale.append(self.sf*pt_zero + centre)
            self.cell_verts.append(pt_scale)

    def plotPolyData(self,zdata = None,clrmap = None):
        """
        Assigns color-data to each polygon in a 2D Voronoi diagram and returns a plot instance (fig, axes)

        Parameters
        ----------
        vor_verts              Nested list of [x,y] points defining each polygon. May be ecm_verts or
                               cell_verts

        zdata                  A data array with each scalar entry corresponding to a polygon entry in
                               vor_verts. If not specified the default is z=1. If 'random'
                               is specified the method creates random vales from 0 to 1..

        clrmap                 The colormap to use for plotting. Must be specified as cm.mapname. A list of
                               available mapnames is supplied at
                               http://matplotlib.org/examples/color/colormaps_reference.html
                               Default is cm.rainbow. Good options are cm.coolwarm, cm.Blues, cm.jet


        Returns
        -------
        fig, ax                Matplotlib figure and axes instances for the plot.

        Notes
        -------
        Uses matplotlib.collections PolyCollection, matplotlib.cm, matplotlib.pyplot and numpy arrays
        Computationally slow -- not recommended for large collectives (500 x 500 um max)
        """
        if zdata == None:  # if user doesn't supply data
            z = np.ones(len(self.cell_verts)) # create flat data for plotting

        elif zdata == 'random':  # if user doesn't supply data
            z = np.random.random(len(self.cell_verts)) # create some random data for plotting

        else:
            z = zdata

        fig, ax = plt.subplots()    # define the figure and axes instances

        # Make the polygon collection and add it to the plot.
        if clrmap == None:
            clrmap = cm.rainbow

        coll = PolyCollection(self.cell_verts, array=z, cmap=clrmap, edgecolors='none')
        ax.add_collection(coll)
        ax.axis('equal')

        # Add a colorbar for the PolyCollection
        if zdata != None:
            fig.colorbar(coll, ax=ax)

        ax.autoscale_view()


        return fig,ax

    def plotCellData(self,zdata=None,clrmap=None,edgeOverlay = None,pointOverlay=None):
        """
        The work-horse of pre-defined plotting methods, this method assigns color-data to each node in cell_centres
        and interpolates data to generate a smooth surface plot. The method returns a plot instance (fig, axes)

        Parameters
        ----------
        zdata                  A data array with each scalar entry corresponding to a point in
                               cell_centres. If not specified the default is z=1. If 'random'
                               is specified the method creates random vales from 0 to 1..

        clrmap                 The colormap to use for plotting. Must be specified as cm.mapname. A list of
                               available mapnames is supplied at
                               http://matplotlib.org/examples/color/colormaps_reference.html
                               Default is cm.rainbow. Good options are cm.coolwarm, cm.Blues, cm.jet

        edgeOverlay             This option allows the user to specify whether or not they want cell edges overlayed.
                                Default is False, set to True to use.

        pointOverlay            This option allows user to specify whether or not they want cell_centre points plotted
                                Default is False, set to True to use.


        Returns
        -------
        fig, ax                Matplotlib figure and axes instances for the plot.

        Notes
        -------
        Uses matplotlib.pyplot and numpy arrays
        With edgeOverlay and pointOverlay == None, this is computationally fast and *is* recommended for plotting data
        on large collectives.


        """
        if zdata == None:  # if user doesn't supply data
            z = np.ones(len(self.cell_centres)) # create flat data for plotting

        elif zdata == 'random':  # if user doesn't supply data
            z = np.random.random(len(self.cell_centres)) # create some random data for plotting

        else:
            z = zdata

        if clrmap == None:
            clrmap = cm.rainbow

        fig, ax = plt.subplots()    # define the figure and axes instances

        sc = 1e6

        triplt = ax.tripcolor(self.cell_centres[:, 0], self.cell_centres[:, 1], z,shading='gouraud', cmap=clrmap)
        ax.axis('equal')

        # Add a colorbar for the z-data
        if zdata != None:
            fig.colorbar(triplt, ax=ax)

        if pointOverlay == True:
            ax.plot(self.cell_centres[:,0],self.cell_centres[:,1],'k.',alpha=0.5)

        if edgeOverlay == True:
            for poly in self.cell_edges:
                for edge in poly:
                    edge = np.asarray(edge)
                    ax.plot(edge[:,0],edge[:,1],color='k',alpha=0.5)

        ax.autoscale_view()


        return fig, ax

    def plotVertData(self,vor_verts,zdata=None,clrmap=None,edgeOverlay = None,pointOverlay=None):
        """
        The work-horse of pre-defined plotting methods, this method assigns color-data to each node in cell_verts,
        ecm_verts, cell_mids, or ecm_mids data structures and interpolates data to generate a smooth surface plot.
        The method returns a plot instance (fig, axes)

        Parameters
        ----------
        vor_verts              An instance of cell_verts, ecm_verts, cell_mids, or ecm_mids

        zdata                  A data array with each scalar entry corresponding to a point in
                               cell_centres. If not specified the default is z=1. If 'random'
                               is specified the method creates random vales from 0 to 1..

        clrmap                 The colormap to use for plotting. Must be specified as cm.mapname. A list of
                               available mapnames is supplied at
                               http://matplotlib.org/examples/color/colormaps_reference.html
                               Default is cm.rainbow. Good options are cm.coolwarm, cm.Blues, cm.jet

        edgeOverlay             This option allows the user to specify whether or not they want cell edges overlayed.
                                Default is False, set to True to use.

        pointOverlay            This option allows user to specify whether or not they want cell_centre points plotted
                                Default is False, set to True to use.


        Returns
        -------
        fig, ax                Matplotlib figure and axes instances for the plot.

        Notes
        -------
        Uses matplotlib.pyplot and numpy arrays
        With edgeOverlay and pointOverlay == None, this is computationally fast and *is* recommended for
        plotting data on large collectives
        """

        vor_verts_flat,_ = self.flatten(vor_verts)

        if zdata == None:  # if user doesn't supply data
            z = np.ones(len(vor_verts_flat)) # create flat data for plotting

        elif zdata == 'random':  # if user doesn't supply data
            z = np.random.random(len(vor_verts_flat)) # create some random data for plotting

        else:
            z = zdata

        if clrmap == None:
            clrmap = cm.rainbow

        fig, ax = plt.subplots()    # define the figure and axes instances

        sc = 1e6

        triplt = ax.tripcolor(vor_verts_flat[:, 0], vor_verts_flat[:, 1], z,shading='gouraud', cmap=clrmap)
        ax.axis('equal')

        # Add a colorbar for the z-data
        if zdata != None:
            fig.colorbar(triplt, ax=ax)

        if pointOverlay == True:
            ax.plot(self.cell_centres[:,0],self.cell_centres[:,1],'k.',alpha=0.5)

        if edgeOverlay == True:
            for poly in self.cell_edges:
                for edge in poly:
                    edge = np.asarray(edge)
                    ax.plot(edge[:,0],edge[:,1],color='k',alpha=0.5)


        ax.autoscale_view()

        return fig, ax

    def plotMemData(self,zdata=None,clrmap=None):
        """

        Assigns color-data to edges in a 2D Voronoi diagram and returns a plot instance (fig, axes)

        Parameters
        ----------
        zdata                  A data array with each scalar entry corresponding to a polygon entry in
                               vor_verts. If not specified the default is z=1. If 'random'
                               is specified the method creates random vales from 0 to 1..

        clrmap                 The colormap to use for plotting. Must be specified as cm.mapname. A list of
                               available mapnames is supplied at
                               http://matplotlib.org/examples/color/colormaps_reference.html
                               Default is cm.rainbow. Good options are cm.coolwarm, cm.Blues, cm.jet


        Returns
        -------
        fig, ax                Matplotlib figure and axes instances for the plot.

        Notes
        -------
        Uses matplotlib.collections LineCollection, matplotlib.cm, matplotlib.pyplot and numpy arrays
        Computationally slow -- not recommended for large collectives (500 x 500 um max)

        """
        fig, ax = plt.subplots()

        # Make a line collection for each cell and add it to the plot.
        for cell in self.cell_edges:
            if zdata == None:
                z = np.ones(len(cell))
            elif zdata == 'random':
                z = np.random.random(len(cell))
            else:
                z = zdata

            if clrmap == None:
                clrmap = cm.rainbow

            coll = LineCollection(cell, array=z, cmap=clrmap)
            ax.add_collection(coll)

        ax.axis('equal')

        # Add a colorbar for the Line Collection
        if zdata != None:
            fig.colorbar(coll, ax=ax)

        ax.autoscale_view()
        ax.axis('equal')

        return fig, ax

    def plotConnectionData(self,zdata=None,clrmap=None):
        """
        Assigns color-data to connections between a cell and its nearest neighbours and returns plot instance

        Parameters
        ----------

        zdata                  A data array with each scalar entry corresponding to a polygon entry in
                               vor_verts. If not specified the default is z=1. If 'random'
                               is specified the method creates random vales from 0 to 1..

        clrmap                 The colormap to use for plotting. Must be specified as cm.mapname. A list of
                               available mapnames is supplied at
                               http://matplotlib.org/examples/color/colormaps_reference.html
                               Default is cm.rainbow. Good options are cm.coolwarm, cm.Blues, cm.jet


        Returns
        -------
        fig, ax                Matplotlib figure and axes instances for the plot.

        Notes
        -------
        Uses matplotlib.collections LineCollection, matplotlib.cm, matplotlib.pyplot and numpy arrays
        Computationally slow -- not recommended for large collectives (500 x 500 um max)

        """
        fig, ax = plt.subplots()

        if zdata == None:
            z = np.ones(len(self.con_segs))

        elif zdata == 'random':
            z = np.random.random(len(self.con_segs))

        else:
            z = zdata

        if clrmap == None:
            clrmap = cm.rainbow

         # Make a line collection for each cell and add it to the plot.

        coll = LineCollection(self.con_segs, array=z, cmap=clrmap)
        ax.add_collection(coll)

        # Plot the cell centres
        ax.plot(self.cell_centres[:,0],self.cell_centres[:,1],'ko')

        ax.axis('equal')

        # Add a colorbar for the Line Collection
        if zdata != None:
            fig.colorbar(coll, ax=ax)

        ax.autoscale_view()


        return fig, ax

    def plotBoundCells(self):
        """

        :param cell_verts:
        :param cdata:
        :return:
        """
        fig, ax = plt.subplots()

        for flagset,cellset in zip(self.bound_flag,self.ecm_verts):
            for flag, cell, in zip(flagset,cellset):
                if flag == 0:
                    ax.plot(cell[0],cell[1],'ko')
                if flag == 1:
                    ax.plot(cell[0],cell[1],'ro')

        # ax.plot(self.clust_xy[:,0],self.clust_xy[:,1],'ko')
        #
        # for inds in self.con_hull:
        #     point = self.clust_xy[inds]
        #     ax.plot(point[:,0],point[:,1],'ro')



        ax.axis('equal')

        ax.autoscale_view()


        return fig, ax

    def flatten(self,ls_of_ls):
        # ls_flat = [val for sublist in ls_of_ls for val in sublist]
        # ls_flat = np.asarray(ls_flat)

        ls_flat = []
        ind_map =[]
        for i, sublist in enumerate(ls_of_ls):
            for j, val in enumerate(sublist):
                ls_flat.append(val)
                ind_map.append([i,j])

        return ls_flat, ind_map

    def alpha_shape(self,points, alpha):
        """
        Calculate the alpha_shape of a cluster of points in 2D.

        Parameters
        ----------
        points              A numpy array listing [[x1,y1],[x2,y2]...] for a collection of 2D points

        alpha               The filtering parameter (which triangles to remove from Delaunay triangulation)
                            Note, for this application alpha = d_cell

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

        return concave_hull














