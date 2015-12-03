#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME create a few options for neat seed points: hexagonal or radial-spiral array

"""
This module contains the class `World`, which holds
all data structures relating to the size of the environment,
the extent of the cell cluster, the co-ordinates of cell
centre points, and all kinds of data relating to individual cell properties.

The initialization method of the `World` class sets-up
and crops the cell cluster to an optional user-defined geometry input
(a set of points arranged in counter-clockwise order and
defining a closed polygon). Other methods define the cell centres of each
cell polygon, their volume, and create cell-cell gap junctions (GJs) and membrane domains
for each cell.
"""

import math
import os
import os.path

import numpy as np
import scipy.spatial as sps
from scipy import interpolate as interp
from scipy import ndimage

from betse.science import toolbox as tb
from betse.science import finitediff as fd
from betse.science.tissue.bitmapper import BitMapper
from betse.science import filehandling as fh
from betse.util.io import loggers


class World(object):
    """
    The World object creates and stores data structures relating to
    the geometric properties of the environmental grid and cell cluster, constructs
    matrices allowing for direct computation of gradients and laplacians on
    cell and environmental points, and provides functions to facilitate data plotting on
    the geometric structures (cell areas, membranes, ect)

    Parameters
    ----------
    constants                           World requires an instance of NumVars, see the Parameters module.

    worldtype (default = None)          'full' creates a complex world with extracellular
                                        matrix points, in addition to cell-cell GJ connections.

                                        'basic' creates a simple world with cell-cell GJ connections.

    Methods
    -------
    makeWorld()                       Create a cell cluster for simulation
    fileInit()                        Create directories for file saving
    makeSeeds()                       Create an irregular lattice of seed points in 2d space
    makeVoronoi()                     Make and clip/close a Voronoi diagram from the seed points
    cell_index()                      Returns a list of [x,y] points defining the cell centres in order
    near_neigh()                      Calculate the nearest neighbour (nn) array for each cell (make gap junctions)
    boundTag(points)                  Creates index-matched boolean lists identifying elements on environ bound
    cellVerts()                       Copy & scale in points from the ecm matrix to create unique polygonal cells
    cleanUp()                         Accessory calculations, including creating of computational matrices
    makeECM()                         Make the Marker and Cell (MACs) grid for extracellular calculations
    environment()                     Calculate details for the extracellular calculations, including mappings
    graphLaplacian()                  Creates an abstract discrete Laplacian for the irregular Voronoi-based cell grid
    """

    def __init__(self, p, worldtype = 'basic'):
        # Extract the constants from the input object:
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
        Initializes file saving and loading directory as the BETSE cache, which is
        automatically assigned from the user-specified path in the configuration file.
        """

        # Make the BETSE-specific cache directory if not found.
        betse_cache_dir = os.path.expanduser(p.init_path)
        os.makedirs(betse_cache_dir, exist_ok=True)

        # Define data paths for saving an initialization and simulation run:
        self.savedWorld = os.path.join(betse_cache_dir, p.world_filename)

    def makeWorld(self,p):

        """
        Calls internal methods to set up the cell cluster.

        """

        if self.worldtype is None or self.worldtype == 'full':
            self.makeSeeds(p)    # Create the grid for the system (irregular)
            self.makeVoronoi(p)    # Make, close, and clip the Voronoi diagram
            self.cell_index(p)            # Calculate the correct centre and index for each cell
            self.cellVerts(p)   # create individual cell polygon vertices
            self.bflags_mems,_ = self.boundTag(self.mem_mids_flat,p,alpha=0.8)  # flag membranes on the cluster bound
            self.bflags_cells,_ = self.boundTag(self.cell_centres,p,alpha=1.0)  # flag membranes on the cluster bound
            self.near_neigh(p)    # Calculate the nn array for each cell
            self.cleanUp(p)       # Free up memory...
            self.voronoiGrid(p)
            self.makeECM(p)       # create the ecm grid
            self.environment(p)   # define features of the ecm grid
            self.grid_len =len(self.xypts)
            self.memWork(p)

            if p.deformation is True:

                self.deformationMatrix(p)


        elif self.worldtype == 'basic':
            self.makeSeeds(p)    # Create the grid for the system (irregular)
            self.makeVoronoi(p)    # Make, close, and clip the Voronoi diagram
            self.cell_index(p)            # Calculate the correct centre and index for each cell
            self.cellVerts(p)   # create individual cell polygon vertices and membrane specific data structures
            self.bflags_mems,_ = self.boundTag(self.mem_mids_flat,p,alpha=0.8)  # flag membranes on the cluster bound
            self.bflags_cells,_ = self.boundTag(self.cell_centres,p,alpha=1.0)  # flag membranes on the cluster bound
            self.near_neigh(p)    # Calculate the nn array for each cell
            self.cleanUp(p)      # Free up memory...
            self.voronoiGrid(p)
            self.makeECM(p)       # create the ecm grid
            self.environment(p)   # features of the environment, without Poisson solvers...
            self.memWork(p)

            if p.deformation is True:

                self.deformationMatrix(p)

    def makeSeeds(self,p):

        """
        Returns an irregular scatter
        of points defined on a world space
        with dimensions supplied by p.wsx in [m].

        The amount of deviation from a square
        grid is specified by p.nl, defined from
        0 (perfect square grid) to 1 (full noise).

        Parameters
        -----------
        p                   An instance of the Parameters object.

        Creates
        -----------
        self.xmin, self.xmax, self.ymin, self.ymax      Min/max points of global world space
        self.centre                                     Centre of global world space
        self.clust_xy                                   List of x,y coordinates of seed points

        Notes
        -------
        Uses Numpy arrays

        """

        # first begin with linear vectors which are the "ticks" of the x and y dimensions
        x_v = np.linspace(0, (p.nx - 1) * (p.d_cell + p.ac), p.nx)  # create lattice vector x
        y_v = np.linspace(0, (p.ny - 1) * (p.d_cell + p.ac), p.ny)  # create lattice vector y

        # next define a 2d array of lattice points using the x- and y- vectors
        x_2d, y_2d = np.meshgrid(x_v, y_v)  # create 2D array of lattice points

        # now create a matrix of points that will add a +/- deviation to each point centre
        x_rnd = p.nl * p.d_cell * (np.random.rand(p.ny, p.nx) - 0.5)  # create a mix of random deltas x dir
        y_rnd = p.nl * p.d_cell * (np.random.rand(p.ny, p.nx) - 0.5)  # create a mix of random deltas x dir

        # add the noise effect to the world point matrices and redefine the results
        x_2d = x_2d + x_rnd
        y_2d = y_2d + y_rnd

        # define a data structure that holds [x,y] coordinate points of each 2d grid-matrix entry
        xypts = np.vstack((x_2d.ravel(), y_2d.ravel())).T

        # define geometric limits and centre for the cluster of points
        self.xmin = np.min(xypts[:,0])
        self.xmax = np.max(xypts[:,0])
        self.ymin = np.min(xypts[:,1])
        self.ymax = np.max(xypts[:,1])

        self.centre = xypts.mean(axis=0)

        self.clust_xy = xypts

    def makeVoronoi(self, p):

        """
        Calculates the Voronoi diagram from cell seed points.

        The Voronoi diagram is then closed at the global (square) boundaries of the world.

        Finally, cells of the Voronoi diagram are removed to define a cluster shape.

        Parameters
        ----------
        p                   An instance of the Parameters object.

        Creates
        ---------
        self.ecm_verts              x,y points of Voronoi cell vertices (nulled at world creation endpoint)
        self.ecm_verts_unique       x,y points of unique Voronoi cell vertices (nulled at world creation endpoint)
        self.ecm_polyinds           indices into self.ecm_verts_unique list defining each Voronoi polygon
                                    (nulled at world creation endpoint)
        self.cluster_mask           Matrix of booleans defining masked shape of cell cluster
        self.msize                  Size of bitmap (side pixel number)


        Notes
        -------
        Uses Numpy arrays
        Uses Scipy spatial

        """

        loggers.log_info('Creating Voronoi geometry... ')

        # define the Voronoi diagram from the seed points:
        vor = sps.Voronoi(self.clust_xy)

        vor.vertices = np.round(vor.vertices,6)  

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


        self.ecm_verts = [] # voronoi verts of clipped cluster

        self.voronoi_verts = []  # keeps track of all voronoi cells, even those not in cluster

        # Clip the Voronoi diagram to polygon defined by clipping bitmap or the default circle:

        # Load the bitmap used to clip the cell cluster and create a clipping function:
        loggers.log_info('Clipping Voronoi geometry to cluster shape...')
        self.bitmasker = BitMapper(
            p.clipping_bitmap_matcher,
            self.xmin, self.xmax, self.ymin, self.ymax)

        for poly_ind in vor.regions: # step through the regions of the voronoi diagram
            if len(poly_ind) >= p.cell_sides:
                cell_poly = vor.vertices[poly_ind]
                point_check = np.zeros(len(cell_poly))

                for i, pnt in enumerate(cell_poly):
                    point_val = self.bitmasker.clipping_function(pnt[0],pnt[1])

                    if point_val != 0.0:
                        point_check[i] = 1.0

                cell_polya = cell_poly.tolist()
                self.voronoi_verts.append(cell_polya)

                if point_check.all() == 1.0:  # if all of the region's point are in the clipping func range
                    self.ecm_verts.append(cell_polya)

        self.cluster_mask = self.bitmasker.clipping_matrix
        self.msize = self.bitmasker.msize

        # next redefine the set of unique vertex points from ecm_verts arrangement:
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

        # Finally, re-do indicies for ecm polygons in terms of unique vertices list:
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

        # Convert ecm_polyinds into a Numpy ndarray:
        self.ecm_polyinds = np.asarray(self.ecm_polyinds)

    def cell_index(self,p):

        """
        Calculate the cell centre for each Voronoi polygon and return a list
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
        # self.voronoi_centres = []

        for poly in self.ecm_verts:
            aa = np.asarray(poly)
            aa = np.mean(aa,axis=0)
            self.cell_centres = np.vstack((self.cell_centres,aa))

        self.cell_centres = np.delete(self.cell_centres, 0, 0)

    def cell_index_quick(self):

        cell_centres_x = np.dot(self.cell_cents_M,self.ecm_verts_unique[:,0])
        cell_centres_y = np.dot(self.cell_cents_M,self.ecm_verts_unique[:,1])

        self.cell_centres = np.column_stack((cell_centres_x,cell_centres_y))

    def near_neigh(self,p):

        """
        Calculate the nearest neighbours for each cell centre in the cluster and return a numpy
        array of nn indices with an index consistent with all other data lists for the cluster.

        Creates
        -------
        self.cell_nn            Indices of all nearest neighbours to each cell (ordered to self.cell_i)
        self.num_nn             Number of nearest neighbours for each cell (ordered to self.cell_i)
        self.average_nn         Average number of nearest neighbours for entire cluster
        self.nn_i               Non-unique list of index pairs to cells, each pair defining a cell-cell GJ
        self.nn_len             Length of each GJ [m]
        self.nn_vects           Normal and tangent vectors to each gj

        Notes
        -------
        Uses numpy arrays
        Uses scipy spatial KDTree search algorithm

        """

        loggers.log_info('Creating gap junctions... ')

        cell_tree = sps.KDTree(self.cell_centres)
        self.cell_nn=cell_tree.query_ball_point(self.cell_centres,p.search_d*p.d_cell)

        # first need to go through and modify cell nn list to get rid of any empty returns:
        temp_cell_nn = []
        for inds in self.cell_nn:
            if len(inds) != 0:
                temp_cell_nn.append(inds)

        self.cell_nn = temp_cell_nn

        self.num_nn = []  # initialize a list that will hold number of nns to a cell

        for indices in self.cell_nn:
            self.num_nn.append(len(indices) -1)  # minus one because query cell is included in each nn list

        self.average_nn = (sum(self.num_nn)/len(self.num_nn))

        #--------------------------------------------------------------------------------------------
        # first -- need to re-do cell nearest neighbours list to remove the self value for each entry:
        nn_list = []

        for i, inds in enumerate(self.cell_nn):
            sublist = []
            for j in inds:
                if i!= j:
                    sublist.append(j)
            nn_list.append(sublist)
        nn_list = np.asarray(nn_list)
        self.cell_nn = nn_list

        # next, calculate a list of non-unique nn pairs for each cell, in addition to nn midpoints and unit vectors:
        self.nn_i = []

        nn_x = []
        nn_y = []
        nn_tx = []
        nn_ty = []

        self.nn_len = []

        for i, inds in enumerate(self.cell_nn):

            for j in inds:

                pt1 = self.cell_centres[i]
                pt2 = self.cell_centres[j]

                mid = (pt1 + pt2)/2       # midpoint calculation
                tang_a = pt2 - pt1       # tangent
                tang = tang_a/np.linalg.norm(tang_a)
                nn_x.append(mid[0])
                nn_y.append(mid[1])
                nn_tx.append(tang[0])
                nn_ty.append(tang[1])

                length = np.sqrt(tang_a[0]**2 + tang_a[1]**2)

                self.nn_len.append(length)

                self.nn_i.append([i,j])

        self.nn_i = np.asarray(self.nn_i)

        self.nn_vects = np.array([nn_x,nn_y,nn_tx,nn_ty]).T

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

        if p.deformation is False:
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

        The BETSE cell grid has each cell defined by unique vertices, which are scaled in from the ecm points.

        Creates
        -------
        self.cell_verts      A nested python list of the [x,y] point pairs defining vertices of each individual cell
                            polygon. The points of each polygon are arranged in a counterclockwise manner.
        self.cell_vol       Volume of each cell [m3]
        self.cell_sa        Whole cell surface area

        self.mem_edges          membrane edge points [x,y] coordinates
        self.mem_length         membrane surface area values [m2]
        self.mem_mids           membrane edge midpoints [x,y] coordinates nested to self.cell_i
        self.mem_mids_flat      unraveled list of membrane edge midpoints [x,y] arranged to self.mem_i
        self.mem_vects_flat     list of normal and tangent vectors (non nested) arranged to self.mem_i


        Notes
        -------
        The Voronoi diagram returns a connected graph. For this simulation, each cell needs unique vertices and edges.
        This method takes the vertices of the original diagram and scales them in to make unique cells.


        """
        self.cell_verts = []
        # self.all_voronoi_verts = []

        if p.deformation is False:

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

        for polyc in self.cell_verts:
            # First calculate individual cell volumes from cell vertices:
            poly = [x for x in reversed(polyc)]
            self.cell_vol.append(p.cell_height*tb.area(poly))



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
                mps.append(mid.tolist())

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

    def cellVerts_quick(self,p):

        # subtract the appropriate cell centre from each ecm vertex:
        ecm_verts_extended = np.dot(self.cell_verts_M,self.ecm_verts_unique)
        cell_centred = ecm_verts_extended - self.cell_centres[self.mem_to_cells]

        self.cell_verts_unique = cell_centred*p.scale_cell + self.cell_centres[self.mem_to_cells]


        # repackage the cell verts into a nested array structure (used for plotting, unfortunately...)

        self.cell_verts = [] # null the original cell verts data structure...

        for i in range(0,len(self.cell_to_mems)):

            vert_nest = self.cell_verts_unique[self.cell_to_mems[i]]

            self.cell_verts.append(vert_nest)

        self.cell_verts = np.asarray(self.cell_verts)   # Voila! Deformed cell_verts!


        # calculate cell area using a matrix version of the Shoelace formula:
        A_flat = 0.5*(self.cell_verts_unique[self.even_vert_inds,0]*self.cell_verts_unique[self.special_vert_inds,1] -
                self.cell_verts_unique[self.even_vert_inds,1]*self.cell_verts_unique[self.special_vert_inds,0])

        self.cell_vol = np.abs(np.dot(self.M_sum_mems,A_flat)*p.cell_height)

        # calculate membrane properties:
        # membrane edges:
        self.mem_edges_flat[:,0,:] = self.cell_verts_unique[self.even_vert_inds]
        self.mem_edges_flat[:,1,:] = self.cell_verts_unique[self.special_vert_inds]

        # membrane mids:
        self.mem_mids_flat = (self.cell_verts_unique[self.even_vert_inds] +
                              self.cell_verts_unique[self.special_vert_inds])/2

        # membrane length and surface area:
        self.mem_length = np.sqrt((self.cell_verts_unique[self.even_vert_inds,0] -
                                        self.cell_verts_unique[self.special_vert_inds,0])**2 +
                                        (self.cell_verts_unique[self.even_vert_inds,1] -
                                         self.cell_verts_unique[self.special_vert_inds,1])**2)

        self.mem_sa = p.cell_height*self.mem_length

        # membrane tangent and normal vectors:
        x1 = self.cell_verts_unique[self.special_vert_inds][:,0]
        y1 = self.cell_verts_unique[self.special_vert_inds][:,1]
        x2 = self.cell_verts_unique[self.even_vert_inds][:,0]
        y2 = self.cell_verts_unique[self.even_vert_inds][:,1]

        tang_x = x2 - x1
        tang_y = y2 - y1
        norm_tang = np.sqrt(tang_x**2 + tang_y**2)
        tang_ax = tang_x/norm_tang
        tang_ay = tang_y/norm_tang

        cv_x = self.mem_mids_flat[:,0]
        cv_y = self.mem_mids_flat[:,1]
        cv_nx = tang_ay
        cv_ny = -tang_ax
        cv_tx = tang_ax
        cv_ty = tang_ay

        self.mem_vects_flat = np.column_stack((cv_x,cv_y,cv_nx,cv_ny,cv_tx,cv_ty))

        # redo boundary tagging:
        self.bflags_mems,_ = self.boundTag(self.mem_mids_flat,p,alpha=0.8)  # flag membranes on the cluster bound
        self.bflags_cells,_ = self.boundTag(self.cell_centres,p,alpha=1.0)  # flag membranes on the cluster bound

    def makeECM(self,p):

        """
        makeECM returns an regular scatter
        of points defined on the world space
        with dimensions wsx, wsy in [m].


        Creates
        -------
        self.delta          spacing between points of square environmental grid
        self.grid_obj       an instance of FiniteDiffSolver, with MACs grid defined within
        self.X, self.Y      major X and Y grid (cell centre points) of MACs grid
        self.map_ij2k       mapping between (i,j) indices of 2D grids and the k indice of the unravelled grid
                            (for grid-cell centres)
        self.index_k        unravelled grid k-index
        self.ecm_sa         surface area of the extracellular space grid-cell face
        self.ecm_vol        volume of the extracellular space grid-cell
        self.xypts          numpy array listing [x,y] of world seed points

        self.xmin, self.xmax      dimensions of world grid
        self.ymin, self.ymax

        self.centre     [x,y] coordinate of world centre

        Notes
        -------
        Uses Numpy arrays

        """

        # base parameter definitions
        # self.delta = p.d_cell*(2/3)
        self.delta = (p.wsx/p.grid_size) # spacing between grid points -- approximately 2/3 of one cell

        self.grid_obj = fd.FiniteDiffSolver()

        self.grid_obj.cell_grid(self.delta,self.xmin,self.xmax,self.ymin,self.ymax)

        self.X = self.grid_obj.cents_X
        self.Y = self.grid_obj.cents_Y

        self.xypts = self.grid_obj.xy_cents
        self.map_ij2k = self.grid_obj.map_ij2k_cents

        # linear k index:
        self.index_k = [x for x in range(0,len(self.xypts))]

        # properties of ecm spaces:
        self.ecm_sa = self.delta*p.cell_height # surface area of ecm space in direction of cell flux
        self.ecm_vol = (p.cell_height*self.delta**2)*np.ones(len(self.xypts))  # volume of ecm space
        self.ecm_r = ((3/4)*(self.ecm_vol/math.pi))**(1/3)  # virtual radius of ecm space

    def environment(self,p):

        """
        Defines conditions for points in contact with the global environment at the outer boundary
        of the square world.

        Notes
        -------

        This is how to access all elements in a main-grid format from the k-vector:
        Z[cells.map_ij2k[:,0],cells.map_ij2k[:,1]]

        or access elements in a main-grid format to a subset of the k-vector:

        accessing bottom boundary:
        Z[cells.map_ij2k[bBot_k][:,0],cells.map_ij2k[bBot_k][:,1]]

        accessing all cells:
        Z[cells.map_ij2k[map_cell2ecm][:,0],cells.map_ij2k[map_cell2ecm][:,1]] = 10.0

        """

        loggers.log_info('Setting global environmental conditions... ')

        # first obtain a structure to map to total xypts vector index:
        self.points_tree = sps.KDTree(self.xypts)

        # define a mapping between a cell and its ecm space in the full list of xy points for the world:
        self.map_cell2ecm = list(self.points_tree.query(self.cell_centres))[1]
        self.map_mem2ecm = list(self.points_tree.query(self.mem_mids_flat,k=1))[1]

        # get a list of all membranes for boundary cells:
        all_bound_mem_inds = self.cell_to_mems[self.bflags_cells]
        all_bound_mem_inds, _ ,_ = tb.flatten(all_bound_mem_inds)

        self.ecm_bound_k = self.map_mem2ecm[self.bflags_mems]  # k indices to xypts for ecms on cluster boundary

        self.ecm_allbound_k = self.map_mem2ecm[all_bound_mem_inds]

        # self.all_clust_pts = np.vstack((self.cell_centres,self.mem_mids_flat))

        # get a list of k indices to the four exterior (global) boundaries of the rectangular world:
        bBot_x = self.X[0,:]
        bTop_x = self.X[-1,:]
        bL_x = self.X[:,0]
        bR_x = self.X[:,-1]

        bBot_y = self.Y[0,:]
        bTop_y = self.Y[-1,:]
        bL_y = self.Y[:,0]
        bR_y = self.Y[:,-1]

        bBot_pts = np.column_stack((bBot_x, bBot_y))
        bTop_pts = np.column_stack((bTop_x, bTop_y))
        bL_pts = np.column_stack((bL_x, bL_y))
        bR_pts = np.column_stack((bR_x, bR_y))

        self.bBot_k = list(self.points_tree.query(bBot_pts))[1]
        self.bTop_k = list(self.points_tree.query(bTop_pts))[1]
        self.bL_k = list(self.points_tree.query(bL_pts))[1]
        self.bR_k = list(self.points_tree.query(bR_pts))[1]

        #-----------------------------------------------------------
        # calculate the cluster masking matrix
        self.make_maskM(p)
        #-------------------------------------------------------------------------
        if p.sim_ECM is True:

            # Create a matrix to update ecm from mem fluxes
            self.ecm_UpdateMatrix = np.zeros((len(self.mem_i),len(self.xypts)))

            for i, ecm_index in enumerate(self.map_mem2ecm):
                self.ecm_UpdateMatrix[i,ecm_index] = 1

            loggers.log_info('Creating environmental Poisson solver for voltage...')
            self.lapENV, self.lapENVinv = self.grid_obj.makeLaplacian()
            self.lapENV = None   # get rid of the non-inverse matrix as it only hogs memory...

            loggers.log_info('Creating environmental Poisson solver for pressure...')
            bdic = {'N':'flux','S':'flux','E':'flux','W':'flux'}
            self.lapENV_P, self.lapENV_P_inv = self.grid_obj.makeLaplacian(bound=bdic)

            self.lapENV_P = None # get rid of the non-inverse matrix as it only hogs memory...

    def short_environment(self,p):

        if self.points_tree is None:

            # first obtain a structure to map to total xypts vector index:
            self.points_tree = sps.KDTree(self.xypts)


        # define a mapping between a cell and its ecm space in the full list of xy points for the world:
        self.map_cell2ecm = list(self.points_tree.query(self.cell_centres))[1]
        self.map_mem2ecm = list(self.points_tree.query(self.mem_mids_flat,k=1))[1]

        # get a list of all membranes for boundary cells:
        all_bound_mem_inds = self.cell_to_mems[self.bflags_cells]
        all_bound_mem_inds, _ ,_ = tb.flatten(all_bound_mem_inds)

        self.ecm_bound_k = self.map_mem2ecm[self.bflags_mems]  # k indices to xypts for ecms on cluster boundary

        self.ecm_allbound_k = self.map_mem2ecm[all_bound_mem_inds]

        if p.sim_ECM is True:

            # Create a matrix to update ecm from mem fluxes
            self.ecm_UpdateMatrix = np.zeros((len(self.mem_i),len(self.xypts)))

            for i, ecm_index in enumerate(self.map_mem2ecm):
                self.ecm_UpdateMatrix[i,ecm_index] = 1

    def graphLaplacian(self,p):
        '''
        Defines an abstract Laplacian that is used to solve Poisson's equation on the
        irregular grid of the cell cluster.

        Parameters
        ----------
        p               An instance of the Parameters object

        Creates
        ----------
        self.lapGJ, self.lapGJ

        '''

        # now define a Laplacian matrix for this cell collection
        self.lapGJ = np.zeros((len(self.cell_i,), len(self.cell_i)))

        nn_inds = self.nn_i.tolist()

        for i, inds in enumerate(self.cell_nn):

            idiag = 0
            ave_mem = self.av_mem_sa[i]
            vol = self.cell_vol[i]

            alpha = ave_mem/vol

            for j in inds:

                nn_index = nn_inds.index([i,j])
                L = self.nn_len[nn_index]

                idiag = idiag - (1/L)*alpha
                self.lapGJ[i,j] = (1/L)*alpha

            self.lapGJ[i,i] = idiag

        self.lapGJinv = np.linalg.pinv(self.lapGJ)

        # null out the original matrix to save memory:
        self.lapGJ = None

    def cleanUp(self,p):

        """
        Nulls unused data structures to free up memory.
        Creates index data structures for unique cell, gap junction and ecm segments (used in simulation during
        randomization of progression through the data structure).

        """

        self.cell_i = [x for x in range(0,len(self.cell_centres))]
        self.nn_index = [x for x in range(0,len(self.nn_i))]

        self.mem_i = [x for x in range(0,len(self.mem_mids_flat))]

        self.cell_vol = np.asarray(self.cell_vol)

        self.R = ((3/4)*(self.cell_vol/math.pi))**(1/3)    # effective radius of each cell

        self.mem_length,_,_ = tb.flatten(self.mem_length)

        self.mem_length = np.asarray(self.mem_length)

        self.mem_sa = self.mem_length*p.cell_height

        loggers.log_info('Creating computational matrices for discrete transfers... ')

        # define map allowing a dispatch from cell index to each respective membrane
        self.indmap_mem = np.asarray(self.indmap_mem)

        self.mem_to_cells = self.indmap_mem[self.mem_i][:,0]   # gives cell index for each mem_i index placeholder

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

        # get an average cell_sa for each membrane domain (used in graph Laplacian calculation):
        self.num_nn = np.asarray(self.num_nn)
        nn_zero = (self.num_nn == 0).nonzero()
        self.num_nn[nn_zero] = 1

        self.av_mem_sa = self.cell_sa/self.num_nn

        self.ave_sa_all = np.mean(self.cell_sa)/np.mean(self.num_nn)

        self.mem_edges_flat, _, _ = tb.flatten(self.mem_edges)
        self.mem_edges_flat = np.asarray(self.mem_edges_flat)

        # structures for plotting interpolated data and streamlines:
        # create a flattened version of cell_verts that will serve as membrane verts:
        self.mem_verts,_,_ = tb.flatten(self.cell_verts)
        self.mem_verts = np.asarray(self.mem_verts)

        self.plot_xy = np.vstack((self.mem_mids_flat,self.mem_verts))

        # do gj stuff as we need it for later:
        self.gj_stuff(p)

        # define matrix for updating cells with fluxes from membranes:
        if self.worldtype == 'full':

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


            self.mem_edges_flat, _, _ = tb.flatten(self.mem_edges)
            self.mem_edges_flat = np.asarray(self.mem_edges_flat)

            self.cell_UpdateMatrix = np.zeros((len(self.mem_i),len(self.cell_i)))

            for i, cell_index in enumerate(self.mem_to_cells):
                self.cell_UpdateMatrix[i,cell_index] = 1

        # matrix for summing property on membranes for each cell and a count of number of mems per cell:
        self.M_sum_mems = np.zeros((len(self.cell_i),len(self.mem_i)))
        self.num_mems = []

        for i, inds in enumerate(self.cell_to_mems):
            n = 0
            for j in inds:
                self.M_sum_mems[i,j] = 1
                n = n+1

            self.num_mems.append(n)

        self.num_mems = np.asarray(self.num_mems)

        self.mem_DivM_inv = np.linalg.pinv(self.M_sum_mems)

        # if studying lateral movement of pumps and channels in membrane,
        # create a matrix that will take a continuous gradient for a value on a cell membrane:
        if p.sim_eosmosis is True:
            self.gradMem = np.zeros((len(self.mem_i),len(self.mem_i)))

            for i, inds in enumerate(self.cell_to_mems):

                inds = np.asarray(inds)

                inds_p1 = np.roll(inds,1)
                inds_n1 = np.roll(inds,-1)
                inds_o = np.roll(inds,0)

                dist = self.mem_mids_flat[inds_p1] - self.mem_mids_flat[inds_n1]
                len_mem = np.sqrt(dist[:,0]**2 + dist[:,1]**2)
                dist_sign = np.sign(self.mem_mids_flat[inds_p1] - self.mem_mids_flat[inds_n1])

                tangx = (self.mem_vects_flat[inds_p1,4] + self.mem_vects_flat[inds_n1,4])/2
                tangy = (self.mem_vects_flat[inds_p1,5] + self.mem_vects_flat[inds_n1,5])/2

                self.gradMem[inds_o,inds_p1] = (1*(tangx/dist_sign[:,0]) + 1*(tangy/dist_sign[:,1]))/len_mem
                self.gradMem[inds_o,inds_n1] = (-1*(tangx/dist_sign[:,0]) - 1*(tangy/dist_sign[:,1]))/len_mem

        #---------------------------------------------------------------------------

        self.cell_number = self.cell_centres.shape[0]
        self.sim_ECM = p.sim_ECM

        self.mem_mids = np.asarray(self.mem_mids)

        # get rid of fields that aren't required any more:
        self.clust_xy = None
        self.ecm_polyinds = None

    def short_cleanUp(self,p):

        self.cell_vol = np.asarray(self.cell_vol)

        self.R = ((3/4)*(self.cell_vol/math.pi))**(1/3)    # effective radius of each cell

        self.mem_length,_,_ = tb.flatten(self.mem_length)

        self.mem_length = np.asarray(self.mem_length)

        self.mem_sa = self.mem_length*p.cell_height

        # cell surface area:
        self.cell_sa = []
        for grp in self.cell_to_mems:
            cell_sa = sum(self.mem_sa[grp])
            self.cell_sa.append(cell_sa)

        self.cell_sa = np.asarray(self.cell_sa)

        # get an average cell_sa for each membrane domain (used in graph Laplacian calculation):
        self.num_nn = np.asarray(self.num_nn)
        nn_zero = (self.num_nn == 0).nonzero()
        self.num_nn[nn_zero] = 1

        self.av_mem_sa = self.cell_sa/self.num_nn

        self.ave_sa_all = np.mean(self.cell_sa)/np.mean(self.num_nn)

        self.mem_verts = self.cell_verts_unique[:]

        self.mem_edges_flat, _, _ = tb.flatten(self.mem_edges)
        self.mem_edges_flat = np.asarray(self.mem_edges_flat)

        self.plot_xy = np.vstack((self.mem_mids_flat,self.cell_verts_unique))

        self.quick_maskM(p)
        # self.make_maskM(p)


        # if studying lateral movement of pumps and channels in membrane,
        # create a matrix that will take a continuous gradient for a value on a cell membrane:
        if p.sim_eosmosis is True:
            #
            # self.mem_edges_flat, _, _ = tb.flatten(self.mem_edges)
            # self.mem_edges_flat = np.asarray(self.mem_edges_flat)

            self.gradMem = np.zeros((len(self.mem_i),len(self.mem_i)))

            for i, inds in enumerate(self.cell_to_mems):

                inds = np.asarray(inds)

                inds_p1 = np.roll(inds,1)
                inds_n1 = np.roll(inds,-1)
                inds_o = np.roll(inds,0)

                dist = self.mem_mids_flat[inds_p1] - self.mem_mids_flat[inds_n1]
                len_mem = np.sqrt(dist[:,0]**2 + dist[:,1]**2)
                dist_sign = np.sign(self.mem_mids_flat[inds_p1] - self.mem_mids_flat[inds_n1])

                tangx = (self.mem_vects_flat[inds_p1,4] + self.mem_vects_flat[inds_n1,4])/2
                tangy = (self.mem_vects_flat[inds_p1,5] + self.mem_vects_flat[inds_n1,5])/2

                self.gradMem[inds_o,inds_p1] = (1*(tangx/dist_sign[:,0]) + 1*(tangy/dist_sign[:,1]))/len_mem
                self.gradMem[inds_o,inds_n1] = (-1*(tangx/dist_sign[:,0]) - 1*(tangy/dist_sign[:,1]))/len_mem

        # redo boundary tagging:
        self.bflags_mems,_ = self.boundTag(self.mem_mids_flat,p,alpha=0.8)  # flag membranes on the cluster bound
        self.bflags_cells,_ = self.boundTag(self.cell_centres,p,alpha=1.0)  # flag membranes on the cluster bound

    def redo_gj(self,dyna,p,savecells =True):

        # profile_names = list(p.tissue_profiles.keys())  # names of each tissue profile...
        profile_names = dyna.tissue_profile_names

        flag_cell_nn = [ [] for x in range(0,len(self.cell_i))]

        for name in profile_names:

            cell_targets = dyna.cell_target_inds[name]   # get the cell target inds for this tissue
            insular_flag = p.tissue_profiles[name]['insular gj']

            # step through gj's and find cases where connection is split between cells in different tissues:
            if insular_flag is True:

                for i, inds in enumerate(self.cell_nn):

                    # see if the cell is in the tissue region:
                    check_a = int(i in cell_targets)

                    # check and see if the cell's neighbours are in the tissue region:
                    for ind_b in inds:

                        check_b = int(ind_b in cell_targets)

                        if check_a != check_b: # if both cells are not in or out of the region:
                            flag_cell_nn[i].append(ind_b)

        new_cell_nn = [ [] for x in range(0,len(self.cell_i))]

        for i, flags in enumerate(flag_cell_nn):

            neighs = self.cell_nn[i]  # get a list of all current neighbours to this cell

            for cell in neighs:

                if cell not in flags:  # if the cell is not in the removal list

                    new_cell_nn[i].append(cell)

        self.cell_nn = np.asarray(new_cell_nn)
        # next, re-calc list of non-unique nn pairs for each cell
        self.nn_i = []

        nn_x = []
        nn_y = []
        nn_tx = []
        nn_ty = []

        self.nn_len = []

        for i, inds in enumerate(self.cell_nn):

            for j in inds:

                pt1 = self.cell_centres[i]
                pt2 = self.cell_centres[j]

                mid = (pt1 + pt2)/2       # midpoint calculation
                tang_a = pt2 - pt1       # tangent
                tang = tang_a/np.linalg.norm(tang_a)
                nn_x.append(mid[0])
                nn_y.append(mid[1])
                nn_tx.append(tang[0])
                nn_ty.append(tang[1])

                length = np.sqrt(tang_a[0]**2 + tang_a[1]**2)

                self.nn_len.append(length)

                self.nn_i.append([i,j])

        self.nn_i = np.asarray(self.nn_i)

        self.nn_vects = np.array([nn_x,nn_y,nn_tx,nn_ty]).T

        # now that basics are done, do the remaining calculations for gap junctions:
        self.gj_stuff(p)

    def gj_stuff(self,p):

         # remake gap junction properties based on new configuration:
        self.nn_index = [x for x in range(0,len(self.nn_i))]

        #------------------------------------------------------------

        # compute mapping between cell and nn with outwards vectors:
        self.cell_to_nn =[[] for x in range(0,len(self.cell_i))]

        # mapping between cell and all nn with duplicates:
        self.cell_to_nn_full = [[] for x in range(0,len(self.cell_i))]

        nn_inds = self.nn_i.tolist()

        for i, inds in enumerate(self.cell_nn):

            for j in inds:
                nn_index = nn_inds.index([i,j])
                self.cell_to_nn[i].append(nn_index)

                self.cell_to_nn_full[i].append(nn_index)
                self.cell_to_nn_full[j].append(nn_index)

        self.cell_to_nn = np.asarray(self.cell_to_nn)
        self.cell_to_nn_full = np.asarray(self.cell_to_nn_full)


        # recalculate matrix for gj divergence of the flux calculation:
        self.gjMatrix = np.zeros((len(self.cell_centres), len(self.nn_i)))

        for igj, pair in enumerate(self.nn_i):

            ci = pair[0]

            sa_i = self.av_mem_sa[ci]

            vol_i = self.cell_vol[ci]

            self.gjMatrix[ci,igj] = -1*(sa_i/vol_i)

         # matrix for averaging values on gap junctions to each cell:

        self.gj2cellMatrix = np.zeros((len(self.cell_i),len(self.nn_i)))

        for i, inds in enumerate(self.cell_to_nn):
            ave_fact = len(inds)
            for j in inds:
                self.gj2cellMatrix[i,j] = 1/ave_fact

        # matrix for summing values on gap junctions to each cell:
        self.gj2cellSum = np.zeros((len(self.cell_i),len(self.nn_i)))

        for i, inds in enumerate(self.cell_to_nn):
            for j in inds:
                self.gj2cellSum[i,j] = 1


        # the nnAveMatrix will take a property defined from two cells onto a single gap junction and average
        # the property to provide one unique result:

        if p.gj_flux_sensitive is True:# if the user desires flux sensitive gj, construct the very large nnAveMatrix:
            self.nnAveMatrix = np.zeros((len(self.nn_i),len(self.nn_i)))

            nn_list = self.nn_i.tolist()

            for i, (pt1, pt2) in enumerate(self.nn_i):
                # find the index of the duplicate point in the gj matrix:
                nn_dupe = nn_list.index([pt2,pt1])

                self.nnAveMatrix[i,i] = 1/2
                self.nnAveMatrix[i,nn_dupe] = 1/2

    def recalc_gj_vects(self,p):
        """
        Recalculate nearest neighbour (gap junction)
        vectors.

        Used in deformation sequence.
        """

        mid = (self.cell_centres[self.nn_i[:,0]] + self.cell_centres[self.nn_i[:,1]])/2

        tang_ax = self.cell_centres[self.nn_i[:,1]][:,0] - self.cell_centres[self.nn_i[:,0]][:,0]
        tang_ay = self.cell_centres[self.nn_i[:,1]][:,1] - self.cell_centres[self.nn_i[:,0]][:,1]

        tang_a = np.sqrt(tang_ax**2 + tang_ay**2)

        tang_x = tang_ax/tang_a
        tang_y = tang_ay/tang_a

        self.nn_vects = np.column_stack((mid[:,0],mid[:,1],tang_x,tang_y))
        self.nn_len = np.sqrt(tang_ax**2 + tang_ay**2)

        # recalculate matrix for gj divergence of the flux calculation:
        self.gjMatrix = np.zeros((len(self.cell_centres), len(self.nn_i)))

        for igj, pair in enumerate(self.nn_i):

            ci = pair[0]

            sa_i = self.av_mem_sa[ci]

            vol_i = self.cell_vol[ci]

            self.gjMatrix[ci,igj] = -1*(sa_i/vol_i)

    def save_cluster(self,p,savecells = True):
        '''
        Saves the cell cluster using a python pickle.

        Parameters
        ----------
        p               Instance of the Parameters object
        savecells       Boolean indicating whether the cluster should be saved (that's kind of dumb if you're calling
                        the function anyway!)

        '''

        if savecells is True:

            self.points_tree = None

            # save the cell cluster
            loggers.log_info('Saving the cell cluster... ')

            datadump = [self,p]
            fh.saveSim(self.savedWorld,datadump)
            message = 'Cell cluster saved to' + ' ' + self.savedWorld
            loggers.log_info(message)

    def gaussMatrix(self,p):
        """
        Defines a matrix that can calculate voltages within cells and the extracellular spaces
        with each precisely defined using the Voronoi-based cluster lattice.

        The construction of this matrix is based on Gauss' law for the electric flux in both the
        enclosed cell or extracellular space.

        """

        cell_stack = np.hstack((self.cell_i,self.mem_i))
        self.a = 0
        self.b = len(self.cell_i) +1
        self.c = self.b -1
        self.d = len(self.mem_i) + self.b

        VMatrix = np.zeros((len(cell_stack),len(cell_stack)))

        ave_mem = np.mean(self.mem_sa)
        ave_vol = np.mean(self.cell_vol)

        ecm_vol = p.cell_space*ave_mem*p.cell_height

        term_cell = ave_mem*(1/(p.tm))
        term_ecm = ave_mem*(1/(p.tm))

        for cell_i in self.cell_i:

            mem_i_set = self.cell_to_mems[cell_i]
            mem_sum = len(mem_i_set)

            VMatrix[cell_i,cell_i] = mem_sum*term_cell

            for mem_i in mem_i_set:

                flags = list((cell_stack == mem_i).nonzero())[0]

                if len(flags) == 2:

                    stack_ind = flags[1]

                else:

                    stack_ind = flags[0]

                VMatrix[cell_i,stack_ind] = -1*term_cell

        for j, mem_i in enumerate(cell_stack[self.c:self.d]):

            cell_i = self.mem_to_cells[mem_i]

            VMatrix[j,cell_i] = VMatrix[j,cell_i] -1*term_ecm
            VMatrix[j,j] = 2*term_ecm

        self.VMatrix_inv = np.linalg.pinv(VMatrix)

    def memWork(self,p):
        """
        Find opposing membranes and those on the boundary.
        For use in deformation calculations.

        """

        self.mem_distance = p.cell_space + 2*p.tm

        sc = (p.rc/2.5)*(p.scale_cell)
        memTree = sps.KDTree(self.mem_mids_flat)

        mem_nn_o = memTree.query_ball_point(self.mem_mids_flat,sc)
        mem_nn = np.zeros((len(self.mem_i),2),dtype=np.int16)
        mem_bound = []
        self.mem_tx = np.zeros(len(self.mem_i))
        self.mem_ty = np.zeros(len(self.mem_i))

        for i, ind_pair in enumerate(mem_nn_o):


            if len(ind_pair) == 1:

                mem_bound.append(*ind_pair)

                mem_nn[i,:] = [i, ind_pair[0]]


            elif len(ind_pair) == 2:

                mem_nn[i,:] = ind_pair

                ta = (self.mem_mids_flat[ind_pair[1]] - self.mem_mids_flat[ind_pair[0]])
                tang = ta/np.linalg.norm(ta)
                self.mem_tx[i] = tang[0]
                self.mem_ty[i] = tang[1]

            elif len(ind_pair) > 2:

        #         index_of_i = ind_pair.index(i)
                i_n = [self.mem_vects_flat[i,2],self.mem_vects_flat[i,3]]

                for j in ind_pair:
                    a = [self.mem_vects_flat[j,2],self.mem_vects_flat[j,3]]
                    ia = round(np.dot(i_n,a),1)

                    if ia == -1.0:

                        mem_nn[i,:] = [i,j]

                        ta = (self.mem_mids_flat[j] - self.mem_mids_flat[i])
                        tang = ta/np.linalg.norm(ta)
                        self.mem_tx[i] = tang[0]
                        self.mem_ty[i] = tang[1]

        self.mem_nn = np.asarray(mem_nn)
        self.mem_bound = np.asarray(mem_bound)

    def voronoiGrid(self,p):

        """
        Creates a set of unique, flat points
        corresponding to cluster and 'ghost' points
        of cell centres, membrane mids, and ecm spaces.

        """

        voronoi_grid = set()

        for verts in self.voronoi_verts:

            for v in verts:
                voronoi_grid.add((v[0],v[1]))


        voronoi_grid = [list(x) for x in voronoi_grid]
        self.voronoi_grid = np.asarray(voronoi_grid)

        # vertTree = sps.KDTree(self.voronoi_grid)
        # self.map_voronoi2ecm = list(vertTree.query(self.ecm_verts_unique))[1]

    def deformationMatrix(self,p):

        """
        Calculate a matrix used to map strain deformation
        from membrane midpoints to the cell vertices.

        """

        loggers.log_info('Creating computational matrices for mechanical deformation... ')
        # calculate some quantities used only in deformation sequences:
        # calculate the chords from cell centre to membrane midpoint
        chord_mag = []
        for i, mids in enumerate(self.mem_mids):
            cent = self.cell_centres[i]

            chords = mids - cent
            chord_m = np.sqrt(chords[:,0]**2 + chords[:,1]**2)
            chord_mag.append(chord_m)

        self.chord_mag, _ , _ = tb.flatten(chord_mag)
        self.chord_mag = np.asarray(self.chord_mag)

        self.cell_verts_unique, _, _ = tb.flatten(self.cell_verts)
        self.cell_verts_unique = np.asarray(self.cell_verts_unique)

        #--------------------------------------------------------------------

        # build search trees for the flattened cell verts and mem mids:
        memTree = sps.KDTree(self.mem_mids_flat)
        vertTree = sps.KDTree(self.cell_verts_unique)

        # create the deformation matrix, which will apply strain at mem mids to the vertices
        # (np.dot(cell_verts_unique,strain)):
        self.deforM = np.zeros((len(self.cell_verts_unique),len(self.mem_mids_flat)))

        for i_cell, cell_verts in enumerate(self.cell_verts):

            mem_mids = self.mem_mids[i_cell]

            if len(mem_mids) == len(cell_verts):

                seq_i_verts = np.arange(0,len(mem_mids))
                seq_ip_mem = np.roll(seq_i_verts,0)
                seq_im_mem = np.roll(seq_i_verts,-1)

                mem_mids = np.asarray(mem_mids)

                vert_points = cell_verts[seq_i_verts]
                mem_points_p = mem_mids[seq_ip_mem]
                mem_points_m = mem_mids[seq_im_mem]

                # find these points in the flattened vectors:
                vert_inds = list(vertTree.query(vert_points))[1]
                mem_inds_p = list(memTree.query(mem_points_p))[1]
                mem_inds_m = list(memTree.query(mem_points_m))[1]

                self.deforM[vert_inds,mem_inds_p] = 1
                self.deforM[vert_inds,mem_inds_m] = 1

        # next build a matrix that will merge the right ecm verts to
        # maintain the same number of unique ecm vertices:
        ecm_verts_flat, map_a, map_b = tb.flatten(self.ecm_verts)
        ecm_verts_flat = np.asarray(ecm_verts_flat)

        # for a set of flattened ecm vertices, find the mapping that will keep originally unique verts merged
        # also, know how to repackage the ecm verts
        ecmTree = sps.KDTree(ecm_verts_flat)
        dist_uniqueECM = list(ecmTree.query(self.ecm_verts_unique,k=2))[0]
        inds_uniqueECM = list(ecmTree.query(self.ecm_verts_unique,k=2))[1]

        # we are going to use these indices to build a matrix that will merge the right vertices of the ecm_verts_flat array,
        # keeping the same number of unique ECM points as the original
        # to use this the syntax is:  ecm_uniques = np.dot(ecm_unique_M,ecm_verts_flat)

        self.ecm_unique_M = np.zeros((len(self.ecm_verts_unique),len(ecm_verts_flat)))

        for i, ind_pair in enumerate(inds_uniqueECM):

            # determine if the second find is a match or a neighbour:
            dist = dist_uniqueECM[i][1]

            if dist == 0.0:  # if distance is indeed zero, then this point has duplicates in the flat list:
                # in the matrix math, the two points of the flat array will be averaged to give the unique list point:
                self.ecm_unique_M[i,ind_pair[0]] = 0.5
                self.ecm_unique_M[i,ind_pair[1]] = 0.5

            else: # this point is unique in the flat and unique lists. Create an identity condition:
                self.ecm_unique_M[i,ind_pair[0]] = 1.0

        # Finally, build a list of inds (self.ecmInds) to map between unique and flattened ecm_verts vectors:
        ecm_verts_flat, map_a, map_b = tb.flatten(self.ecm_verts)
        ecm_verts_flat = np.asarray(ecm_verts_flat)

        ecmTree = sps.KDTree(self.ecm_verts_unique)

        self.ecmInds = list(ecmTree.query(ecm_verts_flat))[1]

        #----matrices for calculating trans-membrane Laplacian and inverse Laplacian:

        self.mem_LapM = np.zeros((len(self.cell_i), len(self.cell_i)))



        for cell_i in self.cell_i:

            # get the set of membrane indices for the cell
            mem_inds = self.cell_to_mems[cell_i]

            diag_multi = self.num_mems[cell_i]  # number of nearest neighbours

            # diagonal element will be the negative of the number of neighbours:
            self.mem_LapM[cell_i,cell_i] = -diag_multi*(1/self.mem_distance)*(self.cell_sa[cell_i])

            for mem_i in mem_inds:

                # find out which membrane the mem_i is partnered to:
                mem_partners = self.mem_nn[mem_i]

                if mem_partners[0] == mem_partners[1]: # then we know we're on a boundary
                    # we know the membrane belongs to this cell, already set so do nothing
                    pass

                elif mem_partners[0] == mem_i: # otherwise, if the first partner is the index from the query cell

                    # find out which cell the partner belongs to:
                    cell_j = self.mem_to_cells[mem_partners[1]]
                    self.mem_LapM[cell_i, cell_j] = (1/self.mem_distance)*(self.mem_sa[mem_i])
                    # self.mem_LapM[cell_j, cell_i] = (1/self.mem_distance)*(self.mem_sa[mem_i])


                elif mem_partners[1] == mem_i:

                    cell_j = self.mem_to_cells[mem_partners[0]]
                    self.mem_LapM[cell_i, cell_j] = (1/self.mem_distance)*(self.mem_sa[mem_i])
                    # self.mem_LapM[cell_j, cell_i] = (1/self.mem_distance)*(self.mem_sa[mem_i])




        self.mem_LapM_inv = np.linalg.pinv(self.mem_LapM)  # take the inverse to solve poisson equation

        # matrices for re-calculating cell verts quickly from new ecm verts:
        ecmTree = sps.KDTree(self.ecm_verts_unique)
        cellTree = sps.KDTree(self.cell_verts_unique)

        self.cell_verts_M = np.zeros((len(self.cell_verts_unique),len(self.ecm_verts_unique)))
        self.cell_cents_M = np.zeros((len(self.cell_i),len(self.ecm_verts_unique)))
        self.ecm_to_cell = np.zeros(len(self.ecm_verts_unique),dtype=np.int64) # maps cell centre to ecm index

        for cell_i, verts in enumerate(self.ecm_verts):

            verts_cll = self.cell_verts[cell_i]
            # get inds to ecm_verts_unique for these verts:
            ecm_verts = list(ecmTree.query(verts))[1]
            cell_verts = list(cellTree.query(verts_cll))[1]

            self.cell_verts_M[cell_verts,ecm_verts] = 1

            num_verts = len(verts)

            for v in ecm_verts:

                self.cell_cents_M[cell_i,v] = 1/num_verts

                self.ecm_to_cell[v] = cell_i

        # speedy way to calculate cell area using matrices:
        self.even_vert_inds = []
        self.special_vert_inds = []

        for verts in self.cell_verts:

            flat_verts = list(cellTree.query(verts))[1]
            flat_verts_rollup = np.roll(flat_verts,1)

            for i in range(len(flat_verts)):
                self.even_vert_inds.append(flat_verts[i])
                self.special_vert_inds.append(flat_verts_rollup[i])

    def make_maskM(self,p):
        """
        Create structures for plotting interpolated data on cell centres
        and differentiating between the cell cluster and environment.

        """

        voronoiTree = sps.KDTree(self.voronoi_grid)
        self.map_voronoi2ecm = list(voronoiTree.query(self.ecm_verts_unique))[1]

        self.voronoi_mask = np.zeros(len(self.voronoi_grid))
        self.voronoi_mask[self.map_voronoi2ecm]=1

        xv = np.linspace(self.xmin,self.xmax,p.plot_grid_size)
        yv = np.linspace(self.xmin,self.xmax,p.plot_grid_size)

        X,Y = np.meshgrid(xv,yv)

        self.Xgrid = X
        self.Ygrid = Y

        self.maskM = interp.griddata((self.voronoi_grid[:,0],self.voronoi_grid[:,1]),self.voronoi_mask,(self.Xgrid,self.Ygrid),
                             method='linear',fill_value=0)

        self.maskM = ndimage.filters.gaussian_filter(self.maskM, 2, mode='nearest')
        self.maskM = np.round(self.maskM,0)


        maskECM = interp.griddata((X.ravel(),Y.ravel()),self.maskM.ravel(), (self.X, self.Y), method='linear',fill_value=0)
        maskECM = ndimage.filters.gaussian_filter(maskECM, 2, mode='nearest')
        maskECM = np.round(maskECM,0)

        self.inds_env = list(*(maskECM.ravel() == 0).nonzero())

    def quick_maskM(self,p):

        self.maskM = interp.griddata((self.voronoi_grid[:,0],self.voronoi_grid[:,1]),self.voronoi_mask,
                                     (self.Xgrid,self.Ygrid),
                                        method='linear',fill_value=0)

        self.maskM = ndimage.filters.gaussian_filter(self.maskM, 2, mode='nearest')
        self.maskM = np.round(self.maskM,0)


        maskECM = interp.griddata((self.Xgrid.ravel(),self.Ygrid.ravel()),self.maskM.ravel(), (self.X, self.Y),
                                  method='linear',fill_value=0)
        maskECM = ndimage.filters.gaussian_filter(maskECM, 2, mode='nearest')
        maskECM = np.round(maskECM,0)

        self.inds_env = list(*(maskECM.ravel() == 0).nonzero())


#-----------WASTELANDS-------------------------------------------------------------------------------------------------

