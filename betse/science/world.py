#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.



# FIXME create a few options for neat seed points: hexagonal or radial-spiral array


"""
This module contains the class World, which holds
all data structures relating to the size of the environment,
the extent of the cell cluster, the co-ordinates of cell
centre points, and all kinds of data relating to individual cell properties.

The initialization method of the World class sets-up
and crops the cell cluster to an optional user-defined geometry input
(a set of points arranged in counter-clockwise order and
defining a closed polygon). Other methods define the cell centres of each
cell polygon, their volume, and create cell-cell gap junctions (GJs) and membrane domains
for each cell. Finally, a suite of methods facilitate adding data (as colour)
to the various geometrical aspects of the cell cluster and return plot objects
that can be integrated into the QT (i.e. PySide) Gui.
"""


import numpy as np
import scipy.spatial as sps
import copy
import math
from betse.science import toolbox as tb
from betse.science.bitmapper import Bitmapper
import os, os.path
from betse.science import filehandling as fh
from betse.util.io import loggers

class World(object):
    """
    The World object creates and stores data structures relating to
    the geometric properties of the environmental grid and cell
    centre points and provides functions to facilitate data plotting on
    the geometric structures (cell areas, membranes, ect)

    Parameters
    ----------
    constants                           World requires an instance of NumVars, see the Parameters module.

    vorclose (default = None)           a set of counter-clockwise arranged points defining a closed
                                        polygon to clip the cluster of cells.

    worldtype (default = None)          'full' creates a complex world with individual membrane domains, extracellular
                                        matrix points, boundary flags, and normal and tangent vectors to each membrane
                                        domain and ecm edge, in addition to cell-cell GJ connections.

                                        'basic' creates a simple world with cell-cell GJ connections.

    Fields
    -------
    self.xmin, self.xmax      dimensions of world grid (after noise, before cropping)
    self.ymin, self.ymax

    self.centre     [x,y] coordinate of world lattice co-ords (after noise, before cropping)

    self.cluster_centre     x,y coordinate of cluster centre

    self.xypts      numpy array holding [x,y] points of irregular world grid

    self.ecm_verts_unique   a numpy array listing [x,y] coordinates of ecm vertices (flattened from ecm_verts)

    self.ecm_edges_i    a list of [n,m] indices into ecm_verts_unique which return two points defining a unique ecm
                        segment. Note the index of ecm_edges_i is *the* ecm index: ecm_i

    self.ecm_mids       a python list of midpoint [x,y] for each unique ecm segment (in ecm_i order)

    self.ecm_vol     a python list of the length of each ecm segment (in ecm_i order)

    self.ecm_vects      a numpy array of [x, y, tx, ty] for each ecm segment (in ecm_i order)

    self.ecmMatrix     a 2D numpy array allowing ecm fluxes to be propperly applied to ecm<-->ecm distributions

    self.mem_to_ecm    a numpy array mapping ecm spaces (ecm_i) to the appropriate membrane

    self.ecm_UpdateMatrix    a matrix updating ecm space concentrations for cell <---> ecm fluxes

    self.cell_UpdateMatrix   a matrix updating cell space concentrations for cell <----> ecm fluxes

    self.cell2ecm_map   a nested numpy array returns the k-ecm indices given the cell [cell_i] and membrane [mem_j] inds

    self.bflags_ecm     a python list of indices to ecm vertices on the env bound (ordered to ecm_verts_unique)

    self.bmask_ecm      a python list of boolean flags to ecm verts on the env bound = 1 (ordered to ecm_verts_unique)

    self.cell_verts     a nested python list specifying [x,y] of verts for each unique cell (arranged to cell_i)

    self.cell_vol     a list of volumes of each cell (arranged to cell_i)  [m3]

    self.cell_sa     a list of total surface area of each cell (arranged to cell_i) [m2]

    self.cell_centres    a numpy array of [x,y] points defining the cell centre (arranged to cell_i)

    self.cell_nn         a nested array of integer indices of each nearest neighbour for a particular cell (arranged
                        to cell_i)

    self.cell_number    a single value reporting the total number of cells in the cluster

    self.average_nn     a single value reporting the average number of nearest neighbour connections per cell

    self.gap_jun_i      a list of index pairs [a,b] to self.cell_i points defining unique cell-cell GJ connections
                        arranged to gap junction index: gj_i

    self.gj_vects       a numpy array of [x,y,tx,ty] defining tangent vectors to each unique gj (arranged to gj_i)

    self.cell2GJ_map    a nested list of indices to gj_i given a particular cell_i

    self.mem_to_cells   a numpy array allowing cell data to be mapped to respective membrane domains

    self.mem_edges     nested python list of segments defining each membrane domain of a cell (arranged cc to cell_i)

    self.mem_length       the length of each membrane domain

    self.mem_sa         the surface area of each membrane domain (pseudo 3D)

    self.mem_mids       nested python list of [x,y] coordinates defining midpoint of each membrane (arranged to cell_i)

    self.mem_mids_flat  numpy array of [x,y] coordinates as flattened version of mem_mids

    self.mem_vects_flat     a numpy array specifying [x,y,nx,ny,tx,ty] specifying the normal and tangent to each membrane
                        domain of a cell. Normals point into the cell when positive.

    self.cell_i         a python list of indices to cell data arrays (cell_i)

    self.ecm_i          a python list of indices to ecm data arrays (ecm_i)

    self.gj_i           a python list of indices to gj data arrays (gj_i)

    self.mem_i          a python list of indices to membrane data arrays (mem_i)

    self.gjMatrix       a matrix allowing a quantity, such as flux, to be properly distributed to cells in network



    Methods
    -------
    makeWorld()                             Create a cell cluster for simulation
    makeSeeds()                             Create an irregular lattice of seed points in 2d space
    cropSeeds(crop_mask =None)              Crop the points cluster to a polygonal shape (circle)
    makeVoronoi(vorclose = None)            Make and clip/close a Voronoi diagram from the seed points
    vor_area()                              Returns the area of each polygon in the closed Voronoi diagram
    cell_index()                            Returns a list of [x,y] points defining the cell centres in order
    near_neigh()                            Calculate the nearest neighbour (nn) array for each cell
    clean_ecm()                             Open ecm cells at the environmental boundary and reformulate data structs
    boundTag(points)                        Creates index-matched boolean lists identifying elements on environ bound
    cellVerts()                             Copy & scale in points from the ecm matrix to create unique polygonal cells
    cellGeo()                               Creates midpoints, lengths, volumes, normal and tangent vectors + more
    cleanUp()                               After computations, null unimportant fields to free up memory


    Notes
    -------
    Uses Numpy
    Uses Scipy spatial
    Uses BETSE-specific toolbox

    """

    def __init__(self, p, worldtype = 'basic'):
        # Extract the constants from the input object:
        # self.vorclose = vorclose   # whether or not to close the voronoi
        # self.crop_mask = crop_mask # whether or not to clip the cluster
        self.worldtype = worldtype # the complexity of cluster to create
        self.fileInit(p)

        self.um = 1e6    # multiplication factor to convert m to um
        self.do_once_cavity = True  # boolean to ensure any cavity preparation happens only once...
        self.do_once_cuts = True      #boolean to ensure cut lines are only prepared once...

        # initialize some parameters that may or may not be used...
        self.cavity_inds = []
        self.true_env_inds = []
        self.cavity_volume = []

    def fileInit(self,p):

        """
        Initializes file saving and loading directory as the betse cach.
        For now, automatically assigns file names, but later, will allow
        user-specified file names.

        """

        # Make the BETSE-specific cache directory if not found.
        betse_cache_dir = os.path.expanduser(p.init_path)
        os.makedirs(betse_cache_dir, exist_ok=True)

        # Define data paths for saving an initialization and simulation run:
        self.savedWorld = os.path.join(betse_cache_dir, p.world_filename)

    def makeWorld(self,p):

        """
        Call internal methods to set up the cell cluster.

        """

        if self.worldtype is None or self.worldtype == 'full':
            self.makeSeeds(p)    # Create the grid for the system (irregular)
            # self.cropSeeds(self.crop_mask,p)      # Crop the grid to a geometric shape to define the cell cluster
            self.makeVoronoi(p)    # Make, close, and clip the Voronoi diagram
            self.cell_index(p)            # Calculate the correct centre and index for each cell
            self.cellVerts(p)   # create individual cell polygon vertices
            self.bflags_ecm,self.bmask_ecm = self.boundTag(self.ecm_verts_unique,p,alpha=1.5)   # flag ecm domains on the env bound
            self.cellGeo(p,close_ecm='yes') # calculate volumes, surface areas, membrane domains, ecm segments and unit vectors
            self.bflags_ecm,_ = self.boundTag(self.ecm_mids,p,alpha=1.2)   # flag ecm domains on the env bound
            self.bflags_cells,_ = self.boundTag(self.cell_centres,p,alpha=0.8)  # flag cell centres on the env bound
            self.near_neigh(p)    # Calculate the nn array for each cell
            self.make_env_points(p)  # get the environmental interaction points for each boundary ecm
            self.cleanUp(p)       # Free up memory...

        elif self.worldtype == 'basic':
            self.makeSeeds(p)    # Create the grid for the system (irregular)
            # self.cropSeeds(self.crop_mask,p)      # Crop the grid to a geometric shape to define the cell cluster
            self.makeVoronoi(p)    # Make, close, and clip the Voronoi diagram
            self.cell_index(p)            # Calculate the correct centre and index for each cell
            self.cellVerts(p)   # create individual cell polygon vertices and membrane specific data structures
            self.bflags_cells,_ = self.boundTag(self.cell_centres,p)
            self.near_neigh(p)    # Calculate the nn array for each cell
            self.cleanUp(p)      # Free up memory...

    def makeSeeds(self,p):

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

        self.xmin, self.xmax      dimensions of world grid (after noise)
        self.ymin, self.ymax

        self.centre     [x,y] coordinate of world centre (after noise)

        Notes
        -------
        Uses Numpy arrays

        """

        # first begin with linear vectors which are the "ticks" of the x and y dimensions
        self.x_v = np.linspace(0, (p.nx - 1) * (p.d_cell + p.ac), p.nx)  # create lattice vector x
        self.y_v = np.linspace(0, (p.ny - 1) * (p.d_cell + p.ac), p.ny)  # create lattice vector y

        # next define a 2d array of lattice points using the x- and y- vectors
        self.x_2d, self.y_2d = np.meshgrid(self.x_v, self.y_v)  # create 2D array of lattice points

        # now create a matrix of points that will add a +/- deviation to each point centre
        x_rnd = p.nl * p.d_cell * (np.random.rand(p.ny, p.nx) - 0.5)  # create a mix of random deltas x dir
        y_rnd = p.nl * p.d_cell * (np.random.rand(p.ny, p.nx) - 0.5)  # create a mix of random deltas x dir

        # add the noise effect to the world point matrices and redefine the results
        x_2d = self.x_2d + x_rnd
        y_2d = self.y_2d + y_rnd

        # define a data structure that holds [x,y] coordinate points of each 2d grid-matrix entry
        self.xypts = np.vstack((x_2d.ravel(), y_2d.ravel())).T

        # define geometric limits and centre for the cluster of points
        self.xmin = np.min(self.xypts[:,0])
        self.xmax = np.max(self.xypts[:,0])
        self.ymin = np.min(self.xypts[:,1])
        self.ymax = np.max(self.xypts[:,1])

        self.centre = self.xypts.mean(axis=0)

        self.clust_xy = self.xypts

    def makeVoronoi(self, p):

        """
        Calculates, closes and clips the Voronoi diagram to cell seed points.

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

        loggers.log_info('Creating Voronoi geometry... ')
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
                    far_point = vor.vertices[new_edge] + direction * p.d_cell

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


        # self.ecm_verts_unique = vor.vertices

        self.ecm_verts = []

        # finally, clip the Voronoi diagram to polygon defined by clipping bitmap or the default circle:

        # load the bitmap used to clip the cell cluster and create a clipping function
        loggers.log_info('Clipping Voronoi geometry to cluster shape... ')
        bitmasker = Bitmapper(p,'clipping',self.xmin, self.xmax,self.ymin,self.ymax)

        for poly_ind in vor.regions: # step through the regions of the voronoi diagram

            if len(poly_ind) >= p.cell_sides:

                cell_poly = vor.vertices[poly_ind]
                point_check = np.zeros(len(cell_poly))

                for i, pnt in enumerate(cell_poly):

                    point_val = bitmasker.clipping_function(pnt[0],pnt[1])

                    if point_val != 0.0:

                        point_check[i] = 1.0

                if point_check.all() == 1.0:  # if all of the region's point are in the clipping func range

                    cell_polya = cell_poly.tolist()
                    self.ecm_verts.append(cell_polya)

        self.cluster_mask = bitmasker.clippingMatrix
        self.msize = bitmasker.msize

        # next redefine the set of unique vertex points from ecm_verts arrangement
        ecm_verts_flat,_,_ = tb.flatten(self.ecm_verts)

        ecm_verts_set = set()

        for vert in ecm_verts_flat:
            ptx = vert[0]
            pty = vert[1]
            ecm_verts_set.add((ptx,pty))

        self.ecm_verts_unique = [list(verts) for verts in list(ecm_verts_set)]

        # Finally, re-do indicies for ecm polygons in terms of unique vertices list
        # self.ecm_verts_unique = self.ecm_verts_unique.tolist()   # first convert to list to use indexing function

        self.ecm_polyinds = []    # define a new field to hold the indices of polygons in terms of unique vertices

        for poly in self.ecm_verts:
            verthold = []
            for vert in poly:
                ind = self.ecm_verts_unique.index(vert)
                verthold.append(ind)

            self.ecm_polyinds.append(verthold)

        self.ecm_verts_unique = np.asarray(self.ecm_verts_unique)  # convert to numpy array

        # #--------------------remove small edges---------------------------------------------------
        #
        loggers.log_info('Cleaning Voronoi geometry... ')
        perm_cut = 2*math.pi*p.rc*p.merge_cut_off # the threshhold edge length

        ecm_verts_2 = []

        for poly in self.ecm_verts: # step through each closed and clipped region of the Voronoi
            hold_verts = []

            if len(poly)<4:

                for vert in poly:
                    hold_verts.append(vert)  # if the region has less than 5 verts, we can only append the points

            elif len(poly) >= 4:      # if the region has greater than or equal to 5 vertices, then proceed

                for i,vert in enumerate(poly):

                    xo = poly[i-1][0]
                    yo = poly[i-1][1]
                    x1 = vert[0]
                    y1 = vert[1]

                    length = math.sqrt((x1-xo)**2 + (y1 -yo)**2)

                    if length > perm_cut:
                        hold_verts.append(poly[i-1])

            hold_verts = np.asarray(hold_verts)

            ecm_verts_2.append(hold_verts)

        self.ecm_verts = ecm_verts_2

        # next redefine the set of unique vertex points from ecm_verts arrangement
        ecm_verts_flat,_,_ = tb.flatten(self.ecm_verts)

        ecm_verts_set = set()

        for vert in ecm_verts_flat:
            ptx = vert[0]
            pty = vert[1]
            ecm_verts_set.add((ptx,pty))

        self.ecm_verts_unique = [list(verts) for verts in list(ecm_verts_set)]

        # Finally, re-do indicies for ecm polygons in terms of unique vertices list
        # self.ecm_verts_unique = self.ecm_verts_unique.tolist()   # first convert to list to use indexing function

        self.ecm_polyinds = []    # define a new field to hold the indices of polygons in terms of unique vertices

        for poly in self.ecm_verts:
            verthold = []
            for vert in poly:
                vert = list(vert)
                ind = self.ecm_verts_unique.index(vert)
                verthold.append(ind)

            self.ecm_polyinds.append(verthold)

        self.ecm_verts_unique = np.asarray(self.ecm_verts_unique)  # convert to numpy array

        # ensure every point in the regions are in order:

        for j, region in enumerate(self.ecm_polyinds):    # step through each polygon region

            verts = self.ecm_verts_unique[region]   # get the vertices for this region
            region = np.asarray(region)      # convert region to a numpy array so it can be sorted
            cent = verts.mean(axis=0)     # calculate the centre point
            angles = np.arctan2(verts[:,1]-cent[1], verts[:,0] - cent[0])  # calculate point angles
            #self.vor.regions[j] = region[np.argsort(angles)]   # sort indices counter-clockwise
            sorted_region = region[np.argsort(angles)]   # sort indices counter-clockwise
            sorted_region_b = sorted_region.tolist()
            self.ecm_polyinds[j] = sorted_region_b   # add sorted list to the regions structure

        # Go back yet again and redo the ecm verts with the organized polyinds
        self.ecm_verts = []
        for i, inds in enumerate(self.ecm_polyinds):
            verts = self.ecm_verts_unique[inds]
            self.ecm_verts.append(verts)

        #-----------------------------------------------------------------------------------------

        # now find the unique vertices used in the cell structure

        self.ecm_polyinds = np.asarray(self.ecm_polyinds)

    def vor_area(self,p):

        """
        Calculates the area of each cell in a closed 2D Voronoi diagram, and multiplying by height, returns cell volume

        Returns
        -------
        self.cell_vol            stores volume of each cell polygon of the Voronoi diagram in cubic meters

        Notes
        -------
        Uses area(p) function.

        """

        loggers.log_info('Calculating area of each cell... ')
        self.cell_vol = []
        for poly in self.cell_verts:
            self.cell_vol.append(p.cell_height*tb.area(poly))

    def cell_index(self,p):

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

    def near_neigh(self,p):

        """
        Calculate the nearest neighbours for each cell centre in the cluster and return a numpy
        array of nn indices with an index consistent with all other data lists for the cluster.

        Creates
        -------
        self.cell_nn            A nested list defining the indices of all nearest neighbours to each cell
        self.gap_jun_i          A list of index pairs to self.cell_centres, each pair defining a unique cell-cell GJ
        self.cell2GJ_map        Returns a list of indices to gap junctions for each cell index

        Notes
        -------
        Uses numpy arrays
        Uses scipy spatial KDTree search algorithm

        """

        loggers.log_info('Creating gap junctions... ')

        cell_tree = sps.KDTree(self.cell_centres)
        self.cell_nn=cell_tree.query_ball_point(self.cell_centres,p.search_d*p.d_cell)

        len_ind = []  # initialize a list that will hold number of nns to a cell

        for indices in self.cell_nn:
            len_ind.append(len(indices) -1)  # minus one because query cell is included in each nn list

        self.average_nn = (sum(len_ind)/len(len_ind))

        GJs = set()
        for cell1_ind, nn_inds in enumerate(self.cell_nn):
            for cell2_ind in nn_inds:
                if cell1_ind == cell2_ind:
                    pass
                elif cell1_ind < cell2_ind:
                    indpair = (cell1_ind,cell2_ind)
                    GJs.add(indpair)
                elif cell1_ind > cell2_ind:
                    indpair = (cell2_ind, cell1_ind)
                    GJs.add(indpair)

        self.gap_jun_i = []

        gv_x = []
        gv_y = []
        gv_tx = []
        gv_ty = []

        for val in GJs:
            vallist = list(val)
            pt1 = self.cell_centres[vallist[0]]
            pt2 = self.cell_centres[vallist[1]]
            pt1 = np.asarray(pt1)
            pt2 = np.asarray(pt2)
            mid = (pt1 + pt2)/2       # midpoint calculation
            tang_a = pt2 - pt1       # tangent
            tang = tang_a/np.linalg.norm(tang_a)
            gv_x.append(mid[0])
            gv_y.append(mid[1])
            gv_tx.append(tang[0])
            gv_ty.append(tang[1])
            self.gap_jun_i.append(vallist)

        self.gj_vects = np.array([gv_x,gv_y,gv_tx,gv_ty]).T

        # self.cell2GJ_map = []
        #
        # for i, neighs in enumerate(self.cell_nn):
        #     holdgj = []
        #     for j, ns in enumerate(neighs):
        #         if i < ns:
        #             gj_ind = self.gap_jun_i.index([i,ns])
        #             holdgj.append(gj_ind)
        #         if i > ns:
        #             gj_ind = self.gap_jun_i.index([ns,i])
        #             holdgj.append(gj_ind)
        #     self.cell2GJ_map.append(holdgj)

        self.gap_jun_i = np.asarray(self.gap_jun_i)

        # calculate lengths of gap junctions
        seg1=self.cell_centres[self.gap_jun_i[:,0]]
        seg2=self.cell_centres[self.gap_jun_i[:,1]]
        nn_diff_gj = (seg2 - seg1)**2
        self.gj_len = np.sqrt(nn_diff_gj[:,0] + nn_diff_gj[:,1])

        if self.worldtype == 'full':
            # repeat process for ecms
            loggers.log_info('Creating ecm junctions... ')
            ecm_tree = sps.KDTree(self.ecm_mids)
            nn_ecm = list(ecm_tree.query(self.ecm_mids,k=5))[1]

            nn_ecm_refined = []
            for i,ind_matches in enumerate(nn_ecm):
                aa = np.delete(ind_matches,0)   # get rid of the self indice
                nn_ecm_refined.append(aa)
            nn_ecm_refined = np.asarray(nn_ecm_refined)

            ecm_nn_set = set()

            for ecm_ind1, ecm_inds in enumerate(nn_ecm_refined):

                for ecm_ind2 in ecm_inds:
                    if ecm_ind1 == ecm_ind2:
                        pass
                    elif ecm_ind1 < ecm_ind2:
                        indpair = ecm_ind1,ecm_ind2
                        ecm_nn_set.add(indpair)
                    elif ecm_ind1 > ecm_ind2:
                        indpair = ecm_ind2,ecm_ind1
                        ecm_nn_set.add(indpair)

            self.ecm_nn_i = []
            for val in ecm_nn_set:
                vallist=list(val)
                self.ecm_nn_i.append(vallist)
            self.ecm_nn_i = np.asarray(self.ecm_nn_i)

            # calculate lengths of ecm-ecm junctions
            seg1=self.ecm_mids[self.ecm_nn_i[:,0]]
            seg2=self.ecm_mids[self.ecm_nn_i[:,1]]
            nn_diff_ecm = (seg2 - seg1)**2
            self.len_ecm_junc = np.sqrt(nn_diff_ecm[:,0] + nn_diff_ecm[:,1])

            ec_x = []
            ec_y = []
            ec_tx = []
            ec_ty = []

            # next calculate tangent unit vectors between ecm and its neighbour:
            for ind_pair in self.ecm_nn_i:
                ind1 = ind_pair[0]
                ind2 = ind_pair[1]
                pt1 = self.ecm_mids[ind1]
                pt2 = self.ecm_mids[ind2]
                mid = (pt1 + pt2)/2
                tang_a = (pt2 - pt1)
                tang = tang_a/np.linalg.norm(tang_a)
                ec_x.append(mid[0])
                ec_y.append(mid[1])
                ec_tx.append(tang[0])
                ec_ty.append(tang[1])

            self.ecm_vects = np.array([ec_x,ec_y,ec_tx,ec_ty]).T

    def make_env_points(self,p):

        """
        Defines points in contact with the environment, which are extrapolated from
        extracellular space boundary points.

        """

        loggers.log_info('Creating environmental points... ')

        delta = p.cell_space

        ecm_bound_points = self.ecm_mids[self.bflags_ecm]
        ecm_bound_norm = self.ecm_seg_vects[self.bflags_ecm,4:6]
        self.env_points = ecm_bound_points + delta*ecm_bound_norm

    def clean_ecm(self,p,clean=None):

        """
        Calculates ecm points on the environmental boundary using the alpha-shape concave hull method,
        deletes these points, and updates data structures referring to the ecm vertices.

        Note: if clean = 'no' this function returns data structures without deleting ecm verticies.

        """

        # get the concave hull for ecm vertices on the outer boundary of the cluster

        if clean is None or clean == 'yes':

            con_hull = tb.alpha_shape(self.ecm_verts_unique, p.scale_alpha/p.d_cell)
            con_hull = np.asarray(con_hull)

            boundverts = np.unique(con_hull)    # get the value of unique indices from segments

             # Re-do indicies for ecm polygons in terms of revised unique vertices list

            ecm_polyinds2 =[]

            b=set(boundverts)

            for i, polyinds in enumerate(self.ecm_polyinds):
                a = set(polyinds)
                ans = list(a & b)
                if len(ans)>0:
                    for val in ans:
                        eind = polyinds.index(val)
                        self.ecm_polyinds[i][eind] = False   # flag value to be removed

            for i, poly in enumerate(self.ecm_polyinds):
                holdval = []
                for j, val in enumerate(poly):
                    if val != False:
                        holdval.append(val)
                ecm_polyinds2.append(holdval)

            self.ecm_polyinds =[]
            self.ecm_polyinds = copy.deepcopy(ecm_polyinds2)

        elif clean == 'no':
            pass

        # Now re-do the ecm_verts

        self.ecm_verts = []

        for poly in self.ecm_polyinds:
            verts = self.ecm_verts_unique[poly]  # [x,y] coordinates of polygon
            self.ecm_verts.append(verts)

        if clean is None or clean == 'yes':
            self.ecm_verts_unique = np.delete(self.ecm_verts_unique,boundverts,0)   # delete indices from ecm verts list

        self.ecm_verts_unique = self.ecm_verts_unique.tolist()   # first convert to list to use indexing function

        self.ecm_polyinds = []    # define a new field to hold the indices of polygons in terms of unique vertices

        for poly in self.ecm_verts:
            verthold = []
            for vert in poly:
                vert =vert.tolist()
                ind = self.ecm_verts_unique.index(vert)
                verthold.append(ind)

            self.ecm_polyinds.append(verthold)

        self.ecm_verts_unique = np.asarray(self.ecm_verts_unique)  # convert back to numpy array
        # now find the unique vertices used in the cell structure

        polyinds_flat,_,_ = tb.flatten(self.ecm_polyinds)
        polyinds_flat = np.asarray(polyinds_flat)
        polyinds_unique = np.unique(polyinds_flat)

        self.ecm_verts_unique = self.ecm_verts_unique[polyinds_unique]

    def boundTag(self,points,p,alpha=1.0):

        """

        Flag elements that are on the boundary to the environment by calculating the convex hull
        for a points cluster.

        Parameters
        ----------
        points          A numpy array of [x,y] points. This may be ecm_verts_unique, cell_centres, or mem_mids_flat.


        Returns
        -------
        bflags       A python list of indices of points that are on the boundary
        bmask        A numpy array of boolean flags of points that are on the boundary (order indexed to points_Flat)

        Notes
        -------
        Uses numpy arrays
        Uses alpha_shape function to calculate the concave hull
        Requires a nested input such as self.mem_mids or self.ecm_verts

        """

        loggers.log_info('Tagging environmental boundary points... ')

        con_hull = tb.alpha_shape(points, alpha/p.d_cell)  # get the concave hull for the membrane midpoints
        con_hull = np.asarray(con_hull)

        bflags = np.unique(con_hull)    # get the value of unique indices from segments

        bmask = np.array((points[:,0]))
        bmask[:] = 0
        bmask[bflags] = 1

        return bflags, bmask

    def cellVerts(self,p):
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

        loggers.log_info('Creating cell vertices and membrane transit vectors... ')

        for centre,poly in zip(self.cell_centres,self.ecm_verts):
            pt_scale = []
            for vert in poly:
                pt_zero = vert - centre
                pt_scale.append(p.scale_cell*pt_zero + centre)
            self.cell_verts.append(np.asarray(pt_scale))

        self.cell_verts = np.asarray(self.cell_verts)

        self.cell_vol = []   # storage for cell volumes
        self.cell_sa = []    # whole cell surface areas
        self.cell_area = []

        self.mem_edges = []  # storage for membrane edge points
        self.mem_length = []   # storage for membrane surface area values
        self.mem_mids = []   # storage for membrane midpoints

        # storage for various vector properties of membrane
        cv_x=[]
        cv_y=[]
        cv_nx=[]
        cv_ny=[]
        cv_tx=[]
        cv_ty=[]

        perim = 2*math.pi*p.rc*p.cell_height    # area of perimeter of cell (general value)

        for polyc in self.cell_verts:
            # First calculate individual cell volumes from cell vertices:
            poly = [x for x in reversed(polyc)]
            self.cell_vol.append(p.cell_height*tb.area(poly))
            # self.cell_sa.append(2*tb.area(poly))   # surface area of whole cell [m2]
            # self.cell_area.append(tb.area(poly))
            # Next calculate individual membrane domains, midpoints, and vectors:
            edge = []
            mps = []
            surfa = []

            for i in range(0,len(polyc)):
                pt1 = polyc[i-1]
                pt2 = polyc[i]
                pt1 = np.asarray(pt1)
                pt2 = np.asarray(pt2)
                edge.append([pt1,pt2])
                mid = (pt1 + pt2)/2       # midpoint calculation
                mps.append(mid)

                lgth = np.sqrt((pt2[0] - pt1[0])**2 + (pt2[1]-pt1[1])**2)  # length of membrane domain
                sa = lgth*p.cell_height    # surface area of membrane
                surfa.append(lgth)

                tang_a = pt2 - pt1       # tangent
                tang = tang_a/np.linalg.norm(tang_a)
                # normal = np.array([-tang[1],tang[0]])
                normal = np.array([tang[1],-tang[0]])
                cv_x.append(mid[0])
                cv_y.append(mid[1])
                cv_nx.append(normal[0])
                cv_ny.append(normal[1])
                cv_tx.append(tang[0])
                cv_ty.append(tang[1])

            self.mem_edges.append(edge)
            self.mem_mids.append(mps)
            self.mem_length.append(surfa)

        self.mem_vects_flat = np.array([cv_x,cv_y,cv_nx,cv_ny,cv_tx,cv_ty]).T

        self.mem_mids_flat, self.indmap_mem, _ = tb.flatten(self.mem_mids)
        self.mem_mids_flat = np.asarray(self.mem_mids_flat)  # convert the data structure to an array

        # self.cell_sa = np.asarray(self.cell_sa)

    def cellGeo(self,p,close_ecm=None):
        """
        Calculates a number of geometric properties relating to cells, membrane domains, and ecm segments.

        Parameters
        ----------
        close_ecm           Default = None. If close_ecm = 'yes' then ecm segments around boundary cells remain
                            unbroken.

        Creates
        --------
        self.cell_vol
        self.mem_edges
        self.mem_length
        self.mem_mids
        self.mem_vects_flat
        self.ecm_vol
        self.ecm_edges_i
        self.ecm_mids
        self.ecm_vects
        self.cell2ecm_map
        """

        self.cell_vol = []   # storage for cell volumes
        self.cell_sa = []    # whole cell surface areas
        self.cell_area = []

        self.mem_edges = []  # storage for membrane edge points
        self.mem_length = []   # storage for membrane surface area values
        self.mem_mids = []   # storage for membrane midpoints

        # storage for various vector properties of membrane
        cv_x=[]
        cv_y=[]
        cv_nx=[]
        cv_ny=[]
        cv_tx=[]
        cv_ty=[]

        loggers.log_info('Creating membrane transit vectors... ')

        perim = 2*math.pi*p.rc*p.cell_height    # area of perimeter of cell (general value)

        for polyc in self.cell_verts:
            # First calculate individual cell volumes from cell vertices:
            poly = [x for x in reversed(polyc)]
            self.cell_vol.append(p.cell_height*tb.area(poly))
            # self.cell_sa.append(2*tb.area(poly))   # surface area of whole cell [m2]
            # self.cell_area.append(tb.area(poly))
            # Next calculate individual membrane domains, midpoints, and vectors:
            edge = []
            mps = []
            surfa = []

            for i in range(0,len(polyc)):
                pt1 = polyc[i-1]
                pt2 = polyc[i]
                pt1 = np.asarray(pt1)
                pt2 = np.asarray(pt2)
                edge.append([pt1,pt2])
                mid = (pt1 + pt2)/2       # midpoint calculation
                mps.append(mid)

                lgth = np.sqrt((pt2[0] - pt1[0])**2 + (pt2[1]-pt1[1])**2)  # length of membrane domain
                sa = lgth*p.cell_height    # surface area of membrane
                surfa.append(lgth)

                tang_a = pt2 - pt1       # tangent
                tang = tang_a/np.linalg.norm(tang_a)
                # normal = np.array([-tang[1],tang[0]])
                normal = np.array([tang[1],-tang[0]])
                cv_x.append(mid[0])
                cv_y.append(mid[1])
                cv_nx.append(normal[0])
                cv_ny.append(normal[1])
                cv_tx.append(tang[0])
                cv_ty.append(tang[1])

            self.mem_edges.append(edge)
            self.mem_mids.append(mps)
            self.mem_length.append(surfa)

        self.mem_vects_flat = np.array([cv_x,cv_y,cv_nx,cv_ny,cv_tx,cv_ty]).T
        self.mem_mids_flat, self.indmap_mem, self.rindmap_mem = tb.flatten(self.mem_mids)
        self.mem_mids_flat = np.asarray(self.mem_mids_flat)  # convert the data structure to an array

        # self.cell_sa = np.asarray(self.cell_sa)

        # Extracellular matrix specific data

        loggers.log_info('Creating extracellular space vectors... ')

        ecm_edge_ind = set()     # this will hold the unique index pairs to the self.ecm_verts_unique [x,y] points

        ecmverts_list = self.ecm_verts_unique.tolist()

        for poly in self.ecm_verts:  # for every polygon defined in the self.ecm_verts data structure

            for i in range(0,len(poly)):   # for every vertex defining the polygon...

                edge_pt1 = poly[i-1].tolist()   # first point of line segment
                edge_pt2 = poly[i].tolist()     # second point of line segment
                edge_ind1 = ecmverts_list.index(edge_pt1)   # get the indices of the [x,y] points in the flat array
                edge_ind2 = ecmverts_list.index(edge_pt2)
                ind1_flag = self.bmask_ecm[edge_ind1]   # get the boolean boundary flag value of point 1
                ind2_flag = self.bmask_ecm[edge_ind2]   # get the boolean boundary flag of point 2

                if close_ecm == 'yes':  # if we want a closed ecm then artificially mark both cells as off bounds
                    ind1_flag = 0
                    ind2_flag = 0

                # in the case that both ecm points are not on the non-boundary (but one may be):
                if (ind1_flag ==0 and ind2_flag == 0) or (ind1_flag ==0 and ind2_flag == 1) or (ind1_flag ==1 and ind2_flag == 0):

                    if edge_ind1 == edge_ind2: # if the indices are equal for some reason, it's not an edge so pass
                        pass

                    if edge_ind1 < edge_ind2:
                        edgepair = (edge_ind1,edge_ind2)
                        ecm_edge_ind.add(edgepair)      # append the indices to the list

                    if edge_ind1 > edge_ind2:
                        edgepair = (edge_ind2, edge_ind1)
                        ecm_edge_ind.add(edgepair)        # append the indices to the list

        self.ecm_edges_i = []   # reconvert everything into a usable python list defining edges
        for edge_tuple in ecm_edge_ind:
            self.ecm_edges_i.append(list(edge_tuple))

        self.cell2ecm_map = []

        len_unique_edges = len(self.ecm_edges_i)

        self.ecm_mids = [0]*len_unique_edges
        self.ecm_vol = [0]*len_unique_edges

        ev_x=[0]*len_unique_edges
        ev_y=[0]*len_unique_edges
        ev_tx=[0]*len_unique_edges
        ev_ty=[0]*len_unique_edges
        ev_nx=[0]*len_unique_edges
        ev_ny=[0]*len_unique_edges


        for i, poly in enumerate(self.ecm_verts):
            holdinds = []
            for j in range(0, len(poly)):
                edge_pt1 = poly[j-1].tolist()    # first point of line segment
                edge_pt2 = poly[j].tolist()      # second point of line segment
                edge_ind1 = ecmverts_list.index(edge_pt1)   # get the indices of the [x,y] points in the flat array
                edge_ind2 = ecmverts_list.index(edge_pt2)
                ind1_flag = self.bmask_ecm[edge_ind1]   # get the boolean boundary flag value of point 1
                ind2_flag = self.bmask_ecm[edge_ind2]   # get the boolean boundary flag of point 2

                if close_ecm == 'yes':  # if we want a closed ecm then artificially mark both cells as off bounds
                    ind1_flag = 0
                    ind2_flag = 0

                # if both of the points are not on the boundary it's a connector
                if (ind1_flag ==0 and ind2_flag == 0) or (ind1_flag == 1 and ind2_flag == 0) or (ind1_flag ==0 and ind2_flag == 1):

                    if edge_ind1 == edge_ind2: # if the indices are equal, it's not an edge so pass
                        pass

                    if edge_ind1 < edge_ind2:
                        mapval = self.ecm_edges_i.index([edge_ind1,edge_ind2])   # get the index to the unique ecm edge
                        holdinds.append(mapval)
                        pnt1 = self.ecm_verts_unique[edge_ind1]
                        pnt2 = self.ecm_verts_unique[edge_ind2]
                        midpoint = (pnt1 + pnt2)/2   # find the midpoint...
                        lgth = np.sqrt((pnt2[0] - pnt1[0])**2 + (pnt2[1]-pnt1[1])**2)  # length of membrane domain
                        vol = lgth*p.cell_height*p.cell_space
                        self.ecm_mids[mapval] = midpoint  # add the midpoint to its list, keeping the same ordering
                        self.ecm_vol[mapval] = vol
                        tang_a = pnt2 - pnt1
                        tang = tang_a/np.linalg.norm(tang_a)
                        ev_x[mapval] = midpoint[0]
                        ev_y[mapval] = midpoint[1]
                        ev_tx[mapval] = tang[0]
                        ev_ty[mapval] = tang[1]
                        ev_nx[mapval] = tang[1]
                        ev_ny[mapval] = -tang[0]

                    if edge_ind2 < edge_ind1:
                        mapval = self.ecm_edges_i.index([edge_ind2,edge_ind1])
                        holdinds.append(mapval)
                        pnt1 = self.ecm_verts_unique[edge_ind2]
                        pnt2 = self.ecm_verts_unique[edge_ind1]
                        midpoint = (pnt1 + pnt2)/2   # find the midpoint...
                        lgth = np.sqrt((pnt2[0] - pnt1[0])**2 + (pnt2[1]-pnt1[1])**2)  # length of membrane domain
                        vol = lgth*p.cell_height*p.cell_space
                        self.ecm_mids[mapval] = midpoint  # add the midpoint to its list, keeping the same ordering
                        self.ecm_vol[mapval] = vol
                        tang_a = pnt2 - pnt1
                        tang = tang_a/np.linalg.norm(tang_a)
                        ev_x[mapval] = midpoint[0]
                        ev_y[mapval] = midpoint[1]
                        ev_tx[mapval] = tang[0]
                        ev_ty[mapval] = tang[1]
                        ev_nx[mapval] = -tang[1]
                        ev_ny[mapval] = tang[0]

            self.cell2ecm_map.append(holdinds)

        self.ecm_seg_vects = np.array([ev_x,ev_y,ev_tx,ev_ty,ev_nx,ev_ny]).T
        self.ecm_mids = np.array(self.ecm_mids)
        self.ecm_edges_i = np.asarray(self.ecm_edges_i)
        self.ecm_vol = np.asarray(self.ecm_vol)

    def cleanUp(self,p):

        """
        Nulls unused data structures to free up memory.
        Creates index data structures for unique cell, gap junction and ecm segments (used in simulation during
        randomization of progression through the data structure).

        """

        self.cell_i = [x for x in range(0,len(self.cell_centres))]
        self.gj_i = [x for x in range(0,len(self.gap_jun_i))]

        self.mem_i = [x for x in range(0,len(self.mem_mids_flat))]

        self.cell_vol = np.asarray(self.cell_vol)

        self.R = np.sqrt(self.cell_vol/(math.pi*p.cell_height))    # effective radius of each cell

        # self.cell_sa = 2*math.pi*self.R*p.cell_height

        self.mem_length,_,_ = tb.flatten(self.mem_length)

        self.mem_length = np.asarray(self.mem_length)

        self.mem_sa = self.mem_length*p.cell_height


        loggers.log_info('Creating computational matrices... ')

        # calculating centre, min, max of cluster after all modifications
        self.clust_centre = np.mean(self.cell_centres)
        self.clust_x_max = np.max(self.cell_centres[:,0])
        self.clust_x_min = np.min(self.cell_centres[:,0])
        self.clust_y_max = np.max(self.cell_centres[:,1])
        self.clust_y_min = np.min(self.cell_centres[:,1])

        # calculating matrix for gap junction flux calculation between cells
        self.gjMatrix = np.zeros((len(self.gj_i),len(self.cell_i)))
        for igj, pair in enumerate(self.gap_jun_i):
            ci = pair[0]
            cj = pair[1]
            self.gjMatrix[igj,ci] = -1
            self.gjMatrix[igj,cj] = 1

        # define map allowing a dispatch from cell index to each respective membrane
        self.indmap_mem = np.asarray(self.indmap_mem)

        self.mem_to_cells = self.indmap_mem[self.mem_i][:,0]   # gives cell index for each mem_i index placeholder

        # data structures for plotting current streamlines
        self.xpts_Igj = np.hstack((self.gj_vects[:,0],self.mem_vects_flat[:,0]))
        self.ypts_Igj = np.hstack((self.gj_vects[:,1],self.mem_vects_flat[:,1]))
        self.nx_Igj = np.hstack((self.gj_vects[:,2],self.mem_vects_flat[:,2]))
        self.ny_Igj = np.hstack((self.gj_vects[:,3],self.mem_vects_flat[:,3]))

        # structures for plotting interpolated data on cell centres:
        xgrid = np.linspace(self.xmin,self.xmax,self.msize)
        ygrid = np.linspace(self.ymin,self.ymax,self.msize)
        # xgrid = np.linspace(self.clust_x_min - p.clip,self.clust_x_max + p.clip,self.msize)
        # ygrid = np.linspace(self.clust_y_min - p.clip,self.clust_y_max + p.clip,self.msize)
        self.Xgrid, self.Ygrid = np.meshgrid(xgrid,ygrid)

        # structures for interpolating data to calculate clean gradients and laplacians:

        self.x_lin = np.linspace(self.xmin,self.xmax,p.grid_size)
        self.y_lin = np.linspace(self.ymin,self.ymax,p.grid_size)

        self.X_cells, self.Y_cells,self.dx_cells, self.dy_cells = tb.makegrid(self.cell_centres[:,0],
            self.cell_centres[:,1],p.grid_size,self)

        self.X_gj, self.Y_gj, self.dx_gj, self.dy_gj = tb.makegrid(self.gj_vects[:,0],self.gj_vects[:,1],p.grid_size,self)

        # compute mapping between cell and gj:
        self.cell_to_gj =[[] for x in range(0,len(self.cell_i))]

        for i, inds in enumerate(self.gap_jun_i):
            ind1 = inds[0]
            ind2 = inds[1]
            self.cell_to_gj[ind1].append(i)
            self.cell_to_gj[ind2].append(i)

        self.cell_to_gj = np.asarray(self.cell_to_gj)

        self.cell_to_mems = []   # construct a mapping giving membrane index for each cell_i

        for cell_index in self.cell_i:

            index2mems = list(*(self.mem_to_cells == cell_index).nonzero())
            self.cell_to_mems.append(index2mems)

        self.cell_to_mems = np.asarray(self.cell_to_mems)

        # cell surface area:
        self.cell_sa = []
        for grp in self.cell_to_mems:
            cell_sa = sum(self.mem_sa[grp])
            self.cell_sa.append(cell_sa)

        self.cell_sa = np.asarray(self.cell_sa)

        # define matrix for updating cells with fluxes from membranes:
        if self.worldtype == 'full':

            self.ecm_i = [x for x in range(0,len(self.ecm_edges_i))]

            self.env_i = [x for x in range(0,len(self.env_points))]

            # calculate the mapping between ecm indices and membrane indices:
            cell_tree = sps.KDTree(self.ecm_mids)
            matches = cell_tree.query(self.mem_mids_flat)
            self.mem_to_ecm = list(matches)[1]

            ecm_n = len(self.ecm_i)
            cell_n = len(self.cell_i)
            mem_n = len(self.mem_i)

            self.cell_UpdateMatrix = np.zeros((mem_n,cell_n))
            self.ecm_UpdateMatrix = np.zeros((mem_n,ecm_n))

            for i, cell_index in enumerate(self.mem_to_cells):
                self.cell_UpdateMatrix[i,cell_index] =1

            for i, ecm_index in enumerate(self.mem_to_ecm):
                self.ecm_UpdateMatrix[i, ecm_index] = 1

            self.ecmMatrix = np.zeros((len(self.ecm_nn_i),len(self.ecm_i)))
            for iecm, pair in enumerate(self.ecm_nn_i):
                ci = pair[0]
                cj = pair[1]
                self.ecmMatrix[iecm,ci] = -1
                self.ecmMatrix[iecm,cj] = 1

            # create a flattened version of cell_verts that will serve as membrane verts:
            self.mem_verts,_,_ = tb.flatten(self.cell_verts)
            self.mem_verts = np.asarray(self.mem_verts)

            cellVertTree = sps.KDTree(self.mem_verts)

            # create a map from flattened mem midpoints to corresponding edge points of mem segment:
            self.index_to_mem_verts = []
            for cell_nest in self.mem_edges:
                for mem_points in cell_nest:
                    pt_ind1 = list(cellVertTree.query(mem_points[0]))[1]
                    pt_ind2 = list(cellVertTree.query(mem_points[1]))[1]
                    self.index_to_mem_verts.append([pt_ind1,pt_ind2])
            self.index_to_mem_verts = np.asarray(self.index_to_mem_verts)

            # create a matrix that will map and interpolate data on mem mids to the mem verts
            # it will work as data on verts = dot( data on mids, matrixMap2Verts ):

            self.matrixMap2Verts = np.zeros((len(self.mem_mids_flat),len(self.mem_verts)))
            for i, indices in enumerate(self.index_to_mem_verts):
                self.matrixMap2Verts[i,indices[0]]=1/2
                self.matrixMap2Verts[i,indices[1]]=1/2


            # create a mapping from each vert to each membrane segment, mem_seg_i:

            self.mem_seg_i = []

            self.mem_edges_flat, _, _ = tb.flatten(self.mem_edges)
            self.mem_edges_flat = np.asarray(self.mem_edges_flat)

            vertTree = sps.KDTree(self.mem_verts)

            for seg in self.mem_edges_flat:
                pt1 = seg[0]
                pt2 = seg[1]
                seg_ind1 = vertTree.query(pt1)[1]
                seg_ind2 = vertTree.query(pt2)[1]
                self.mem_seg_i.append([seg_ind1,seg_ind2])

            self.mem_seg_i = np.asarray(self.mem_seg_i)  # pairs two indices to mem_verts defining line segment

            # now to go from membrane vert data to mid data by calculating the pseudo-inverse:
            self.matrixMap2Mids = np.linalg.pinv(self.matrixMap2Verts)

             # calculating matrix for membrane flux calculation between connected vertices:
            self.memMatrix = np.zeros((len(self.mem_seg_i),len(self.mem_i)))
            for igj, pair in enumerate(self.mem_seg_i):
                ci = pair[0]
                cj = pair[1]
                self.memMatrix[igj,ci] = -1
                self.memMatrix[igj,cj] = 1

            # Create a map from cell to ecm space
            self.cell_to_ecm = []

            for i in self.cell_i:

                inds_mtc = (self.mem_to_cells ==i).nonzero()
                ecm_inds = self.mem_to_ecm[inds_mtc]
                self.cell_to_ecm.append(ecm_inds)

            self.cell_to_ecm = np.asarray(self.cell_to_ecm)
            self.bcell_to_ecm = self.cell_to_ecm[self.bflags_cells]
            self.bcell_to_ecm,_,_ = tb.flatten(self.bcell_to_ecm)

            # structures for plotting interpolated data and streamlines:
            self.plot_xy = np.vstack((self.mem_mids_flat,self.mem_verts))

            self.xpts_Iecm = self.ecm_vects[:,0]
            self.ypts_Iecm = self.ecm_vects[:,1]
            self.nx_Iecm = self.ecm_vects[:,2]
            self.ny_Iecm = self.ecm_vects[:,3]

            self.nx_Ienv = self.ecm_seg_vects[:,4][self.bflags_ecm]
            self.ny_Ienv = self.ecm_seg_vects[:,5][self.bflags_ecm]

             # structures for interpolating data to calculate clean gradients and laplacians:

            self.X_ecm, self.Y_ecm, self.dx_ecm, self.dy_ecm = tb.makegrid(self.ecm_mids[:,0],self.ecm_mids[:,1],
                p.grid_size, self)

            self.X_ej, self.Y_ej, self.dx_ej, self.dy_ej = tb.makegrid(self.ecm_vects[:,0],self.ecm_vects[:,1],
                p.grid_size, self)

            loggers.log_info('Cleaning up unnecessary data structures... ')

            # self.indmap_mem = None
            # self.rindmap_mem = None
            # # self.ecm_verts = None
            # self.cell_area = None
            # self.cell2ecm_map = None
            # self.ecm_polyinds = None
            # self.ecm_verts_unique = None
            # self.cell2GJ_map = None


        self.cell_number = self.cell_centres.shape[0]
        self.sim_ECM = p.sim_ECM

        self.clust_xy = None
        self.cell_nn = None

        self.mem_mids = np.asarray(self.mem_mids)

        loggers.log_info('Cell cluster creation complete!')

    def redo_gj(self,dyna,p,savecells =True):

        # profile_names = list(p.tissue_profiles.keys())  # names of each tissue profile...
        profile_names = dyna.tissue_profile_names
        new_gj_nn_i = []
        new_gj_vects = []
        removal_flags = np.zeros(len(self.gj_i))

        for name in profile_names:

            cell_targets = dyna.cell_target_inds[name]   # get the cell target inds for this tissue
            insular_flag = p.tissue_profiles[name]['insular gj']

            # step through gj's and find cases where connection is split between cells in different tissues:
            if insular_flag == True:

                for i, ind_pair in enumerate(self.gap_jun_i):

                    ind_a = ind_pair[0]
                    ind_b = ind_pair[1]
                    check_a = len(list(*(cell_targets == ind_a).nonzero()))
                    check_b = len(list(*(cell_targets == ind_b).nonzero()))

                    if check_a != check_b:
                        removal_flags[i] = 1.0

        for i, (flag, ind_pair) in enumerate(zip(removal_flags,self.gap_jun_i)):

            gj_vects = self.gj_vects[i,:]

            if flag == 0:
                new_gj_nn_i.append(ind_pair)
                new_gj_vects.append(gj_vects)

        self.gap_jun_i = np.asarray(new_gj_nn_i)
        self.gj_vects = np.asarray(new_gj_vects)

        # remake gap junction properties based on new configuration:
        self.gj_i = [x for x in range(0,len(self.gap_jun_i))]

        self.gjMatrix = np.zeros((len(self.gj_i),len(self.cell_i)))
        for igj, pair in enumerate(self.gap_jun_i):
            ci = pair[0]
            cj = pair[1]
            self.gjMatrix[igj,ci] = -1
            self.gjMatrix[igj,cj] = 1

        # recompute mapping between cell and gj:
        self.cell_to_gj =[[] for x in range(0,len(self.cell_i))]

        for i, inds in enumerate(self.gap_jun_i):
            ind1 = inds[0]
            ind2 = inds[1]
            self.cell_to_gj[ind1].append(i)
            self.cell_to_gj[ind2].append(i)

        self.cell_to_gj = np.asarray(self.cell_to_gj)

        if savecells == True:

            # save the cell cluster
            loggers.log_info('Saving the cell cluster... ')

            # celf = copy.deepcopy(self)
            datadump = [self,p]
            fh.saveSim(self.savedWorld,datadump)
            message = 'Cell cluster saved to' + ' ' + self.savedWorld
            loggers.log_info(message)

    def redoCells(self,p):


        loggers.log_info('Re-creating computational matrices... ')

        # calculating matrix for gap junction flux calculation between cells
        self.gjMatrix = np.zeros((len(self.gj_i),len(self.cell_i)))
        for igj, pair in enumerate(self.gap_jun_i):
            ci = pair[0]
            cj = pair[1]
            self.gjMatrix[igj,ci] = -1
            self.gjMatrix[igj,cj] = 1

        # define map allowing a dispatch from cell index to each respective membrane
        self.mem_mids_flat, self.indmap_mem, self.rindmap_mem = tb.flatten(self.mem_mids)
        self.indmap_mem = np.asarray(self.indmap_mem)

        self.mem_to_cells = self.indmap_mem[self.mem_i][:,0]   # gives cell index for each mem_i index placeholder

        # data structures for plotting current streamlines
        self.xpts_Igj = np.hstack((self.gj_vects[:,0],self.mem_vects_flat[:,0]))
        self.ypts_Igj = np.hstack((self.gj_vects[:,1],self.mem_vects_flat[:,1]))
        self.nx_Igj = np.hstack((self.gj_vects[:,2],self.mem_vects_flat[:,2]))
        self.ny_Igj = np.hstack((self.gj_vects[:,3],self.mem_vects_flat[:,3]))

        # structures for plotting interpolated data on cell centres:
        xgrid = np.linspace(self.xmin,self.xmax,self.msize)
        ygrid = np.linspace(self.ymin,self.ymax,self.msize)
        self.Xgrid, self.Ygrid = np.meshgrid(xgrid,ygrid)

        # compute mapping between cell and gj:
        self.cell_to_gj =[[] for x in range(0,len(self.cell_i))]

        for i, inds in enumerate(self.gap_jun_i):
            ind1 = inds[0]
            ind2 = inds[1]
            self.cell_to_gj[ind1].append(i)
            self.cell_to_gj[ind2].append(i)

        self.cell_to_gj = np.asarray(self.cell_to_gj)

        self.cell_to_mems = []   # construct a mapping giving membrane index for each cell_i

        for cell_index in self.cell_i:

            index2mems = list(*(self.mem_to_cells == cell_index).nonzero())
            self.cell_to_mems.append(index2mems)

        self.cell_to_mems = np.asarray(self.cell_to_mems)

        # define matrix for updating cells with fluxes from membranes:
        if self.worldtype == 'full':

            # self.ecm_i = [x for x in range(0,len(self.ecm_edges_i))]
            #
            # self.env_i = [x for x in range(0,len(self.env_points))]

            # calculate the mapping between ecm indices and membrane indices:
            cell_tree = sps.KDTree(self.ecm_mids)
            matches = cell_tree.query(self.mem_mids_flat)
            self.mem_to_ecm = list(matches)[1]

            ecm_n = len(self.ecm_i)
            cell_n = len(self.cell_i)
            mem_n = len(self.mem_i)

            self.cell_UpdateMatrix = np.zeros((mem_n,cell_n))
            self.ecm_UpdateMatrix = np.zeros((mem_n,ecm_n))

            for i, cell_index in enumerate(self.mem_to_cells):
                self.cell_UpdateMatrix[i,cell_index] =1

            for i, ecm_index in enumerate(self.mem_to_ecm):
                self.ecm_UpdateMatrix[i, ecm_index] = 1

            self.ecmMatrix = np.zeros((len(self.ecm_nn_i),len(self.ecm_i)))
            for iecm, pair in enumerate(self.ecm_nn_i):
                ci = pair[0]
                cj = pair[1]
                self.ecmMatrix[iecm,ci] = -1
                self.ecmMatrix[iecm,cj] = 1

            # create a flattened version of cell_verts that will serve as membrane verts:
            self.mem_verts,_,_ = tb.flatten(self.cell_verts)
            self.mem_verts = np.asarray(self.mem_verts)

            cellVertTree = sps.KDTree(self.mem_verts)

            # create a map from flattened mem midpoints to corresponding edge points of mem segment:
            self.index_to_mem_verts = []
            for cell_nest in self.mem_edges:
                for mem_points in cell_nest:
                    pt_ind1 = list(cellVertTree.query(mem_points[0]))[1]
                    pt_ind2 = list(cellVertTree.query(mem_points[1]))[1]
                    self.index_to_mem_verts.append([pt_ind1,pt_ind2])
            self.index_to_mem_verts = np.asarray(self.index_to_mem_verts)

            # create a matrix that will map and interpolate data on mem mids to the mem verts
            # it will work as data on verts = dot( data on mids, matrixMap2Verts ):

            self.matrixMap2Verts = np.zeros((len(self.mem_mids_flat),len(self.mem_verts)))
            for i, indices in enumerate(self.index_to_mem_verts):
                self.matrixMap2Verts[i,indices[0]]=1/2
                self.matrixMap2Verts[i,indices[1]]=1/2

            # Create a map from cell to ecm space
            self.cell_to_ecm = []

            for i in self.cell_i:

                inds_mtc = (self.mem_to_cells ==i).nonzero()
                ecm_inds = self.mem_to_ecm[inds_mtc]
                self.cell_to_ecm.append(ecm_inds)

            self.cell_to_ecm = np.asarray(self.cell_to_ecm)
            self.bcell_to_ecm = self.cell_to_ecm[self.bflags_cells]
            self.bcell_to_ecm,_,_ = tb.flatten(self.bcell_to_ecm)

            # structures for plotting interpolated data and streamlines:
            self.plot_xy = np.vstack((self.mem_mids_flat,self.mem_verts))

            self.xpts_Iecm = self.ecm_vects[:,0]
            self.ypts_Iecm = self.ecm_vects[:,1]
            self.nx_Iecm = self.ecm_vects[:,2]
            self.ny_Iecm = self.ecm_vects[:,3]

            self.nx_Ienv = self.ecm_seg_vects[:,4][self.bflags_ecm]
            self.ny_Ienv = self.ecm_seg_vects[:,5][self.bflags_ecm]

        self.cell_number = self.cell_centres.shape[0]
        self.sim_ECM = p.sim_ECM

        loggers.log_info('Cell cluster re-creation complete!')





