#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME method to create a hexagonal lattice base
# FIXME figure out how to scale each voronoi polygon to um instead of m dimensions when plotting
# FIXME allow user to specify their own set of points for clipping in points and voronoi clips (make circle function)


"""
This module contains the class World, which holds
all data structures relating to the size of the environment,
the extent of the cell cluster, the co-ordinates of cell
centre points, and all kinds of data relating to individual cell properties.

The initialization method of the World class sets-up
and crops the cell cluster to an optional user-defined geometry input
as a set of points arranged in counter-clockwise order and
defining a closed polygon. Other methods define the cell centres of each
cell polygon, their area, and create nearest neighbour and edge matrices
for each cell. Finally, a suite of methods facilitate adding data (as colour)
to the various geometrical aspects of the cell cluster and return plot objects
that can be integrated into the QT (i.e. PySide) Gui.

"""

import numpy as np
import scipy.spatial as sps
from matplotlib.path import Path
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.collections import PolyCollection
import matplotlib.cm as cm
import math
from science import toolbox as tb


class World(object):
    """
    The WorldSeeds object creates and stores data structures relating to
    the geometric properties of the environmental grid and cell
    centre points and provides functions to facilitate data plotting on
    the geometric structures (cell areas, membranes, ect)

    Parameters
    ----------
    WorldSeeds requires an instance of NumVars, see the Parameters module.
    Optional: crop_mask (default = None) a set of counter-clockwise
    arranged points defining a closed polygon to clip the cluster of
    cell seed points.

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

    self.ecm_verts     a nested python list containing Voronoi cell regions which specify [x,y] points of region
                        vertices for a clipped Voronoi diagram (arranged to cell index)

    self.cell_verts     a nested python list specifying coordinates of verticies for each unique cell in the cluster

    self.cell_area     a list of areas of each cell (arranged to consistent cell index)

    self.cell_centres    a list of [x,y] points defining the cell centre (arranged to cell index)

    self.cell_nn         a nested array of integer indices of each nearest neighbour for a particular cell (arranged
                        to cell index)

    self.con_segs       a list of two [x,y] points defining line segments that join cells in the cluster

    self.ecm_edges_i     a python list of two point line segments of the ecm

    self.cell_edges    a python nested list of line segments defining each membrane domain of a cell (cc arrangement)

    self.mem_mids       a nest python list of [x,y] coordinates defining the midpoint of each membrane for each cell

    self.ecm_mids_i       a python list defining the [x,y] coordinates of the midpoint of each ecm segment

    self.bflags_mems    a list flagging cell membrane domains on the envirnomental boundary with 1

    self.bflags_ecm     a list flagging ecm domains on the environmental boundary with 1

    self.cell_vects     a numpy array specifying [x,y,nx,ny,tx,ty] specifying the normal and tangent to each membrane
                        domain of a cell. Normals point into the cell when positive.

    self.ecm_vects      a numpy array specifying [x,y,tx,ty] the tangent to each unique ecm segment


    Methods
    -------
    makeSeeds()                             Create an irregular lattice of seed points in 2d space
    cropSeeds(crop_mask)                    Crop the points cluster to a polygonal shape (circle)
    makeVoronoi(vorclose = None)            Make and clip/close a Voronoi diagram from the seed points
    vor_area()                              Returns the area of each polygon in the closed Voronoi diagram
    cell_index()                            Returns a list of [x,y] points defining the cell centres in order
    near_neigh()                            Calculate the nearest neighbour (nn) array for each cell
    boundTag(points)                        Creates index-matched boolean lists identifying elements on environ bound
    cellEdges()                             List of membrane domains as two-point line segments for each cell
    cellMids()                              Creates coordinate lists of membrane and ecm segment midpoints
    cellVerts()                             Creates unique vertices defining each cell from the Voronoi diagram
    cellVects()                             Creates normals and tangents to each cell membrane

    Plotting methods:

    plotPolyData(zdata = None,clrmap = None)                                   Plot cell polygons with data as colour
    plotCellData(zdata=None,clrmap=None,edgeOverlay = None,pointOverlay=None)  Plot smoothed cell-centre data as colour
    plotVertData(vor_verts,zdata=None,clrmap=None,edgeOverlay=None,pointOverlay=None)       Plot smoothed nested data
    plotMemData(zdata=None,clrmap=None)                                      Plot membrane domains with data as colour
    plotConnectionData(zdata=None,clrmap=None)                      Plot GJ connections with data as colour
    plotBoundCells()                                    Plot points flagged as existing on the environmental boundary
    plotVects()                                  Plot unit vectors corresponding to cell membrane and ecm


    Notes
    -------
    Uses Numpy
    Uses Scipy spatial
    Uses BETSE-specific toolbox

    """

    def __init__(self,constants,crop_mask=None, vorclose=None):
        # Extract the constants from the input object:
        self.vorclose = vorclose   # whether or not to close the voronoi
        self.crop_mask = crop_mask # whether or not to clip the cluster

        self.d_cell = constants.dc  # diameter of single cell
        self.nx = constants.nx   # number of lattice sites in world x index
        self.ny = constants.ny   # number of lattice sites in world y index
        self.ac = constants.ac  # cell-cell separation
        self.nl = constants.nl  # noise level for the lattice
        self.wsx = constants.wsx  # World size
        self.wsy = constants.wsy # World size
        self.search_d = constants.search_d  # distance to search for nearest neighbours (relative to d_cell)
        self.sf = constants.scale_cell              # scale factor to take cell vertices in from extracellular space
        self.cell_sides = constants.cell_sides # minimum number of membrane domains per cell
        self.sa = constants.scale_alpha        # amount to scale (1/d_cell) in boundary search algorithm (alpha_shape)

        #self.makeWorld()

    def makeWorld(self):

        """
        Call internal methods to set up the cell cluster.

        """

        self.makeSeeds()    # Create the grid for the system (irregular)
        self.cropSeeds(self.crop_mask)      # Crop the grid to a geometric shape to define the cell cluster
        self.makeVoronoi(self.vorclose)    # Make, close, and clip the Voronoi diagram
        self.ecm_verts_flat, self.indmap_ecm = tb.flatten(self.ecm_verts)   #  a flat list and ind map to ecm verts
        self.cell_index()            # Calculate the correct centre and index for each cell
        self.cellVerts()   # create individual cell polygon vertices
        self.vor_area()              # Calculate the area of each cell polygon
        self.near_neigh()    # Calculate the nn array for each cell
        self.cellEdges()    # create a nested list of all membrane and ecm domains for each cell
        self.cellMids()    # calculate the midpoint of membrane domains and ecm segments
        self.bflags_mems = self.boundTag(self.mem_mids)   # flag mem domains on the environmental boundary
        self.bflags_ecm = self.boundTag(self.ecm_verts)   # flag ecm domains on the environmental boundary
        self.cellVects()          # calculate the normals and tangents to each membrane domain of each cell

        self.cell_number = self.cell_centres.shape[0]

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

        self.xmin, self.xmax      dimensions of world grid (after noise)
        self.ymin, self.ymax

        self.centre     [x,y] coordinate of world centre (after noise)

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
        Uses Numpy arrays.

        Important: bug found that if points are cropped first, and then a Voronoi is created and closed,
        the result is sporadically totally messed up. Therefore, the points are not being pre-cropped and
        crop_mask is always = None in the instancing of this World class!

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
        vorclose            None = no cropping, 'circle'= crop to circle

        Creates
        -------
        self.ecm_verts      nested python list specifying polygonal region
                            and vertices as [x,y] for each Voronoi cell in the
                            clipped/closed Voronoi diagram. Arranged as: [ [ [a,b],[c,d],[e,f]],[[g,h],[i,j],[k,l] ] ]
                            These represent the vertices of the extracellular space around each cell.

        Notes
        -------
        Uses Numpy arrays
        Uses Scipy spatial

        """

        vor = sps.Voronoi(self.clust_xy)

        cluster_center = vor.points.mean(axis=0)

        # complete the Voronoi diagram by adding in undefined vertices to ridges and regions
        i = -1   # enumeration index

        for pnt_indx, vor_edge in zip(vor.ridge_points, vor.ridge_vertices):
            vor_edge = np.asarray(vor_edge)

            i = i+1 # update the count-through index

            if np.any(vor_edge < 0): # if either of the two ridge values are undefined (-1)

                # find the ridge vertice that's not equal to -1
                    new_edge = vor_edge[vor_edge >= 0][0]
                # calculate the tangent of two seed points sharing that ridge
                    tang = vor.points[pnt_indx[1]] - vor.points[pnt_indx[0]]
                    tang /= np.linalg.norm(tang)  # make the tangent a unit vector
                    norml = np.array([-tang[1], tang[0]])  # calculate the normal of the two points sharing the ridge

                    # calculate the midpoint between the two points of the ridge
                    midpoint = vor.points[pnt_indx].mean(axis=0)
                    # now there's enough information to calculate the missing direction and location of missing point
                    direction = np.sign(np.dot(midpoint - cluster_center, norml)) * norml
                    #far_point = self.vor.vertices[new_edge] + direction * self.cluster_axis.max()
                    far_point = vor.vertices[new_edge] + direction * self.d_cell

                    # get the current size of the voronoi vertices array, this will be the n+1 index after adding point
                    vor_ind = vor.vertices.shape[0]

                    vor.vertices = np.vstack((vor.vertices,far_point)) # add the new point to the vertices array
                    vor.ridge_vertices[i] = [new_edge,vor_ind]  # add the new index at the right spot

                    for j, region in enumerate(vor.regions):    # step through each polygon region

                        if len(region):

                            if -1 in region and new_edge in region:  # if the region has edge of interest...
                                a = region.index(-1)              # find index in the region that is undefined (-1)
                                vor.regions[j][a] = vor_ind # add in the new vertex index to the appropriate region

                            verts = vor.vertices[region]   # get the vertices for this region
                            region = np.asarray(region)      # convert region to a numpy array so it can be sorted
                            cent = verts.mean(axis=0)     # calculate the centre point
                            angles = np.arctan2(verts[:,1]-cent[1], verts[:,0] - cent[0])  # calculate point angles
                            #self.vor.regions[j] = region[np.argsort(angles)]   # sort indices counter-clockwise
                            sorted_region = region[np.argsort(angles)]   # sort indices counter-clockwise
                            sorted_region_b = sorted_region.tolist()
                            vor.regions[j] = sorted_region_b   # add sorted list to the regions structure


        # finally, clip the Voronoi diagram to polygon, if user-specified by vorclose option
        if vorclose==None:
            self.ecm_verts=[]
            for region in vor.regions:
                if len(region):
                    cell_poly = vor.vertices[region]
                    if len(cell_poly)>3:
                        self.ecm_verts.append(vor.vertices[region])


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

            for poly_ind in vor.regions:  # step through each cell's polygonal regions...

                if len(poly_ind) >= self.cell_sides: # check to make sure we're defining a polygon

                    cell_poly = vor.vertices[poly_ind]  # get the coordinates of the polygon vertices

                    inpath = crop_path.contains_points(cell_poly)    # get a boolean matrix

                    if inpath.all() == False:
                        pass

                    else:
                        cell_polya = cell_poly.tolist()  # convert data structures to python lists for cropping algorithm...
                        aa=tb.clip(cell_polya,crop_ptsa)        # then send it to the clipping algorithm
                        if len(aa) >= self.cell_sides:                        # check to make sure result is still a polygon
                            self.ecm_verts.append(aa)     # append points to new region point list

    def vor_area(self):

        """
        Calculates the area of each cell in a closed 2D Voronoi diagram.

        Returns
        -------
        self.cell_area            area of all polygons of the Voronoi diagram in square meters

        Notes
        -------
        Uses area(p) function.

        """
        self.cell_area = []
        for poly in self.cell_verts:
            self.cell_area.append(tb.area(poly))

    def cell_index(self):

        """
        Calculate the cell centre for each voronoi polygon and return a list
        with an index consistent with all other data lists for the cell cluster.


        Creates
        -------
        self.cell_centres      [x,y] coordinate of the centre of each cell as a numpy array

        Notes
        -------
        After the Voronoi diagram has been created,closed, and clipped, this method is required to
        create an ordering of cells that is consistent with the Voronoi polygons, membrane domains, and ecm polygons
        and segments.


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

        Creates
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
        len_ind = []  # initialize a list that will hold number of nns to a cell

        for centre, indices in zip(self.cell_centres,self.cell_nn):
            len_ind.append(len(indices))
            for pt in indices:
                self.con_segs.append([centre,self.cell_centres[pt]])

        self.average_nn = sum(len_ind)/len(len_ind)

    def boundTag(self,points):

        """

        Flag elements that are on the boundary to the environment by calculating the convex hull
        for a points cluster.

        Parameters
        ----------
        points          A nested set of [x,y] points as a list-of-lists


        Returns
        -------
        bflags       A python list with 0 indicating point not on boundary and 1 indicating boundary

        Notes
        -------
        Uses numpy arrays
        Uses alpha_shape function to calculate the concave hull
        Requires a nested input such as self.mem_mids or self.ecm_verts

        """
        # initialize the boundary flag structure
        bflags = []
        for cellvals in points:
            hoo = []
            for pts in cellvals:
                hoo.append(0)
            bflags.append(hoo)

        # get a points cluster and a map back to the original indices
        points_flat, indmap = tb.flatten(points)
        points_flat = np.asarray(points_flat)  # convert to numpy array for plotting, etc

        con_hull = tb.alpha_shape(points_flat, self.sa/self.d_cell)  # get the concave hull for the membrane midpoints

        for inds in con_hull:
            for val in inds:
                org_ind = indmap[val]
                bflags[org_ind[0]][org_ind[1]]=1  # set the boundary flag in the original data format to 1

        return bflags

    def cellEdges(self):

        """

        Calculate line segments corresponding to a cell's unique membrane domains or shared extracellular edge.
        These are arranged in a manner consistent with the cell index.

        Returns
        -------
        self.ecm_edges_i      A python list of the indices to point pairs in self.ecm_verts_flat, which define unique
                            line segments of the ecm upon calling self.ecm_verts_flat[self.ecm_edges_i].

        self.cell_edges      A nested python list of the [x,y] point pairs defining line segments of each membrane
                            domain in a cell polygon. The list has segments arranged in a counterclockwise manner.

        """

        ecm_edge_ind = []     # this will hold the indices to the self.ecm_verts_flat [x,y] points
        ecm_edges_unique = set()   # a set that will contain only the unique indices so no edges are repeated.

        for poly in self.ecm_verts:  # for every polygon defined in the self.ecm_verts data structure
            for i in range(0,len(poly)):   # for every vertice defining the polygon...
                edge_pt1 = poly[i-1]    # first point of line segment
                edge_pt2 = poly[i]      # second point of line segment
                edge_ind1 = self.ecm_verts_flat.index(edge_pt1)   # get the indices of the [x,y] points
                edge_ind2 = self.ecm_verts_flat.index(edge_pt2)
                ecm_edge_ind.append([edge_ind1,edge_ind2])      # append the indices to the list

        for i, edge in enumerate(ecm_edge_ind):  # re-jigger the list to make sure that all [a,b] == [b,a] elements
            pt1 = edge[0]
            pt2 = edge[1]
            if pt1 > pt2:
                ecm_edge_ind[i]=[pt2,pt1]

        ecm_edge_ind.sort()

        for edge in ecm_edge_ind:    # for each edge indices [a,b] convert into a tuple (a,b) and add to set
            edge_tuple = tuple(edge)
            ecm_edges_unique.add(edge_tuple)  # the set makes sure only one copy of (a,b) is stored...

        self.ecm_edges_i = []
        for edge_tuple in ecm_edges_unique:
            self.ecm_edges_i.append(list(edge_tuple))   # reconvert everything into a usable python list defining edges
                                                        # note: to go to the cell-specific edges use the indmap_ecm

        # now get edges which correspond to the unique membrane domains of each cell (these are automatically unique)
        self.cell_edges = []
        for poly in self.cell_verts:
            edge =[]
            for i in range(0,len(poly)):
                edge.append([poly[i-1],poly[i]])
            self.cell_edges.append(edge)

    def cellMids(self):

        """

        The midpoint between line segments defining the ecm and the cell membrane domains.

        Returns
        -------
        self.mem_mids       A nested python list of the [x,y] points defining the midpoint of each membrain domain
                            around a cell. Arranged counterclockwise rotation.
        self.ecm_mids_i      A nested python list of the [x,y] points defining the midpoint of each ecm segment around a
                            cell.

        """

        self.mem_mids = []
        self.ecm_mids_i = []

        for celledges in self.cell_edges:
            hoo = []
            for celledge in celledges:
                cellpt1 = celledge[0]
                cellpt2 = celledge[1]
                cellpt1 = np.asarray(cellpt1)
                cellpt2 = np.asarray(cellpt2)
                mid = (cellpt1 + cellpt2)/2
                hoo.append(mid)
            self.mem_mids.append(hoo)

        # ecm_mids_i are a bit more complex. First go through and calculate all midpoints from the cell-specific
        # self.ecm_verts data structure points. Then, flatten the resulting list-of-lists and get an indices map.
        # Finally, find the unique midpoints and store the result as indices to the flattened list of points.

        mids_list_nested = []

        for poly in self.ecm_verts:
            mids_list = []
            for i in range(0,len(poly)):
                edge_pt1 = poly[i-1]    # first point of line segment
                edge_pt2 = poly[i]      # second point of line segment
                edge_pt1 = np.asarray(edge_pt1)
                edge_pt2 = np.asarray(edge_pt2)
                midpoint = (edge_pt1 + edge_pt2)/2   # find the midpoint...
                midpoint = midpoint.tolist()   # convert it back to a python list...
                mids_list.append(midpoint)    # append the midpoint [x,y] to a list...
            mids_list_nested.append(mids_list)  # append the previous list to get cell-organized data

        self.ecm_mids_flat, self.indmap_mids = tb.flatten(mids_list_nested)   # get the flattened version and ind map

        mids_unique = set()   # a set that will contain only the unique indices so no mids are repeated.

        for poly in mids_list_nested:  # step through the nested list of points (which contain redundancies)
            for pt in poly:
                m_index = self.ecm_mids_flat.index(pt)  # get the index of the flattened midpoints list
                mids_unique.add(m_index)  # add it to the set (will only be added once if it's duplicate)

        for el in mids_unique:
            self.ecm_mids_i.append(el)     # get a python list of the unique indices for the ecm midpoints
                                            # note: can get cell-organized indices by using self.indmap_mids

    def cellVerts(self):
        """
        Calculate the true vertices of each individual cell from the extracellular matrix (ecm) vertices
        of the closed & clipped Voronoi diagram.

        Creates
        -------
        self.cell_verts      A nested python list of the [x,y] point pairs defining vertices of each individual cell
                            polygon. The points of each polygon are arranged in a counterclockwise manner.

        Notes
        -------
        The Voronoi diagram returns a connected graph. For this simulation, each cell needs unique vertices and edges.
        This method takes the vertices of the original diagram and scales them in to make unique cells.


        """
        self.cell_verts = []

        for centre,poly in zip(self.cell_centres,self.ecm_verts):
            pt_scale = []
            for vert in poly:
                pt_zero = vert - centre
                pt_scale.append(self.sf*pt_zero + centre)
            self.cell_verts.append(pt_scale)

    def cellVects(self):
        """
        Calculates unit vectors that are normal and tangent to each cell's membrane domain.
        Calculates unit vectors that are tangent to each ecm line segment.

        """

        cv_x=[]
        cv_y=[]
        cv_nx=[]
        cv_ny=[]
        cv_tx=[]
        cv_ty=[]

        ev_x=[]
        ev_y=[]
        ev_tx=[]
        ev_ty=[]

        for pnt_list in self.cell_edges:
            for pt1,pt2 in pnt_list:
                tang_a = pt2 - pt1
                tang = tang_a/np.linalg.norm(tang_a)
                normal = np.array([-tang[1],tang[0]])
                midpoint = (pt1 +pt2)/2
                cv_x.append(midpoint[0])
                cv_y.append(midpoint[1])
                cv_nx.append(normal[0])
                cv_ny.append(normal[1])
                cv_tx.append(tang[0])
                cv_ty.append(tang[1])

        self.cell_vects = np.array([cv_x,cv_y,cv_nx,cv_ny,cv_tx,cv_ty]).T

        for ind_list in self.ecm_edges_i:
            pt1 = self.ecm_verts_flat[ind_list[0]]
            pt2 = self.ecm_verts_flat[ind_list[1]]
            pt1 = np.asarray(pt1)
            pt2 = np.asarray(pt2)
            tang_a = pt2 - pt1
            tang = tang_a/np.linalg.norm(tang_a)
            midpoint = (pt1 +pt2)/2
            ev_x.append(midpoint[0])
            ev_y.append(midpoint[1])
            ev_tx.append(tang[0])
            ev_ty.append(tang[1])

        self.ecm_vects = np.array([ev_x,ev_y,ev_tx,ev_ty]).T

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

        ax.autoscale_view(tight=True)


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

        ax.autoscale_view(tight=True)


        return fig, ax

    def plotVertData(self,vor_verts,zdata=None,clrmap=None,edgeOverlay = None,pointOverlay=None):
        """
        The work-horse of pre-defined plotting methods, this method assigns color-data to each node in cell_verts,
        ecm_verts, cell_mids, or ecm_mids_i data structures and interpolates data to generate a smooth surface plot.
        The method returns a plot instance (fig, axes)

        Parameters
        ----------
        vor_verts              An instance of cell_verts, ecm_verts, cell_mids, or ecm_mids_i

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

        vor_verts_flat, _ = tb.flatten(vor_verts)

        vor_verts_flat = np.asarray(vor_verts_flat)

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


        ax.autoscale_view(tight=True)

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

        ax.axis('equal')
        ax.autoscale_view(tight=True)

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
        ax.plot(self.cell_centres[:,0],self.cell_centres[:,1],'k.')

        ax.axis('equal')

        # Add a colorbar for the Line Collection
        if zdata != None:
            fig.colorbar(coll, ax=ax)

        ax.autoscale_view(tight=True)


        return fig, ax

    def plotBoundCells(self, points, bflags):
        """
        Plot elements tagged on the boundary as red points.

        Parameters
        ----------
        points          A nested array of points corresponding to the bflags data structure

        bflags          A nested array of boolean flags indicating boundary tagging

        Returns
        -------
        fig, ax         Matplotlib plotting objects

        Note
        ------
        This particular plot is extremely slow -- intended for cross-checking purposes only!

        """
        fig, ax = plt.subplots()

        for flagset,cellset in zip(bflags,points):
            for flag, cell, in zip(flagset,cellset):
                if flag == 0:
                    ax.plot(cell[0],cell[1],'ko')
                if flag == 1:
                    ax.plot(cell[0],cell[1],'ro')

        ax.axis('equal')

        ax.autoscale_view(tight=True)

        return fig, ax

    def plotVects(self):
        """
        This function plots all unit vectors in the tissue system as a cross-check.
        Normals to cell membranes are shown as red arrows.
        Tangents to cell membranes are black arrows.
        Tangents to ecm edges are shown as green arrows.
        Cell membrane edges are drawn as blue lines.

        To plot streamline and vector plots with data use the pyplot quiver and streamplot functions, respectively.

        """
        fig, ax = plt.subplots()

        ax.quiver(self.cell_vects[:,0],self.cell_vects[:,1],self.cell_vects[:,4],self.cell_vects[:,5],color='k')
        ax.quiver(self.cell_vects[:,0],self.cell_vects[:,1],self.cell_vects[:,2],self.cell_vects[:,3],color='r')
        ax.quiver(self.ecm_vects[:,0],self.ecm_vects[:,1],self.ecm_vects[:,2],self.ecm_vects[:,3],color='g')

        for cell in self.cell_edges:
            coll = LineCollection(cell)
            ax.add_collection(coll)

        ax.axis('equal')

        ax.autoscale_view(tight=True)

        return fig, ax














