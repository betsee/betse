#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME allow user to specify their own set of points for clipping in points and voronoi clips (make circle function)

# FIXME create a few options for neat seed points: hexagonal or radial-spiral array
# FIXME allow user to save and load a world

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
from matplotlib.path import Path
import copy
import math
from betse.science import toolbox as tb
from betse.science.parameters import params as p

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

    self.xypts      numpy array holding [x,y] points of irregular world grid

    self.ecm_verts     a nested python list containing Voronoi cell regions which specify [x,y] points of region
                        vertices for a clipped Voronoi diagram (arranged to cell index)

    self.ecm_verts_unique   a numpy array listing [x,y] coordinates of ecm vertices (flattened from ecm_verts)

    self.indmap_ecm     a python list returning [n,m] indices into ecm_verts for each index k of ecm_verts_unique

    self.rindmap_ecm    a python list-of-lists returning indek k of ecm_verts at [n][m] indices

    self.ecm_edges_i    a list of [n,m] indices into ecm_verts_unique which return two points defining a unique ecm
                        segment. Note the index of ecm_edges_i is *the* ecm index: ecm_i

    self.ecm_mids       a python list of midpoint [x,y] for each unique ecm segment (in ecm_i order)

    self.ecm_vols     a python list of the length of each ecm segment (in ecm_i order)

    self.ecm_vects      a numpy array of [x, y, tx, ty] for each ecm segment (in ecm_i order)

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

    self.mem_edges     nested python list of segments defining each membrane domain of a cell (arranged cc to cell_i)

    self.mem_length       the length of each membrane domain

    self.mem_mids       nested python list of [x,y] coordinates defining midpoint of each membrane (arranged to cell_i)

    self.mem_mids_flat  numpy array of [x,y] coordinates as flattened version of mem_mids

    self.indmap_mem     returns indices to mem_mids from index in mem_mids_flat

    self.rindmap_mem    returns incice to mem_mids_flat from nested indices corresponding to mem_mids

    self.bflags_mems    a list of membrane midpoints on the boundary (indices to mem_mids_flat)

    self.bmask_mems     a list of booleans with midpoints on boundary = 1 (arranged to mem_mids_flat)

    self.mem_vects_flat     a numpy array specifying [x,y,nx,ny,tx,ty] specifying the normal and tangent to each membrane
                        domain of a cell. Normals point into the cell when positive.

    self.cell_i         a python list of indices to cell data arrays (cell_i)

    self.ecm_i          a python list of indices to ecm data arrays (ecm_i)

    self.gj_i           a python list of indices to gj data arrays (gj_i)

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

    def __init__(self,crop_mask=None, vorclose=None, worldtype = None):
        # Extract the constants from the input object:
        self.vorclose = vorclose   # whether or not to close the voronoi
        self.crop_mask = crop_mask # whether or not to clip the cluster
        self.worldtype = worldtype # the complexity of cluster to create

        self.um = 1e6    # multiplication factor to convert m to um

    def makeWorld(self):

        """
        Call internal methods to set up the cell cluster.

        """

        if self.worldtype == None or self.worldtype == 'full':
            self.makeSeeds()    # Create the grid for the system (irregular)
            self.cropSeeds(self.crop_mask)      # Crop the grid to a geometric shape to define the cell cluster
            self.makeVoronoi(self.vorclose)    # Make, close, and clip the Voronoi diagram
            self.cell_index()            # Calculate the correct centre and index for each cell
            self.cellVerts()   # create individual cell polygon vertices
            self.clean_ecm(clean='no')  # pop ecm vertices around the outer cell membranes
            self.bflags_ecm,self.bmask_ecm = self.boundTag(self.ecm_verts_unique)   # flag ecm domains on the env bound
            self.cellGeo(close_ecm='yes') # calculate volumes, surface areas, membrane domains, ecm segments and unit vectors
            self.mem_mids_flat, self.indmap_mem, self.rindmap_mem = tb.flatten(self.mem_mids)
            self.mem_mids_flat = np.asarray(self.mem_mids_flat)  # convert the data structure to an array
            self.bflags_mems,self.bmask_mems = self.boundTag(self.mem_mids_flat)   # flag mem domains on the env bound
            self.near_neigh()    # Calculate the nn array for each cell
            self.cleanUp()       # Free up memory...

        elif self.worldtype == 'basic':
            self.makeSeeds()    # Create the grid for the system (irregular)
            self.cropSeeds(self.crop_mask)      # Crop the grid to a geometric shape to define the cell cluster
            self.makeVoronoi(self.vorclose)    # Make, close, and clip the Voronoi diagram
            self.cell_index()            # Calculate the correct centre and index for each cell
            self.cellVerts()   # create individual cell polygon vertices
            self.vor_area()              # Calculate the area of each cell polygon
            self.near_neigh()    # Calculate the nn array for each cell
            self.cleanUp()      # Free up memory...

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
        x_2d, y_2d = np.meshgrid(self.x_v, self.y_v)  # create 2D array of lattice points

        # now create a matrix of points that will add a +/- deviation to each point centre
        x_rnd = p.nl * p.d_cell * (np.random.rand(p.ny, p.nx) - 0.5)  # create a mix of random deltas x dir
        y_rnd = p.nl * p.d_cell * (np.random.rand(p.ny, p.nx) - 0.5)  # create a mix of random deltas x dir

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

        self.ecm_verts_unique = vor.vertices


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

                if len(poly_ind) >= p.cell_sides: # check to make sure we're defining a polygon
                    cell_poly = vor.vertices[poly_ind]  # get the coordinates of the polygon vertices
                    inpath = crop_path.contains_points(cell_poly)    # get a boolean matrix

                    if inpath.all() == False:  # if all of the polygon's points are outside of the crop path, ignore it
                        pass

                    else:
                        cell_polya = cell_poly.tolist()  # convert data structures to python lists for cropping algorithm...
                        aa=tb.clip(cell_polya,crop_ptsa)        # then send it to the clipping algorithm
                        if len(aa) >= p.cell_sides:                        # check to make sure result is still a polygon
                            self.ecm_verts.append(aa)     # append points to new region point list

        # Finally, re-do indicies for ecm polygons in terms of unique vertices list
        self.ecm_verts_unique = self.ecm_verts_unique.tolist()   # first convert to list to use indexing function

        self.ecm_polyinds = []    # define a new field to hold the indices of polygons in terms of unique vertices

        for poly in self.ecm_verts:
            verthold = []
            for vert in poly:
                ind = self.ecm_verts_unique.index(vert)
                verthold.append(ind)

            self.ecm_polyinds.append(verthold)

        self.ecm_verts_unique = np.asarray(self.ecm_verts_unique)  # convert back to numpy array

        # now find the unique vertices used in the cell structure

        #self.ecm_polyinds = np.asarray(self.ecm_polyinds)

    def vor_area(self):

        """
        Calculates the area of each cell in a closed 2D Voronoi diagram, and multiplying by height, returns cell volume

        Returns
        -------
        self.cell_vol            stores volume of each cell polygon of the Voronoi diagram in cubic meters

        Notes
        -------
        Uses area(p) function.

        """
        self.cell_vol = []
        for poly in self.cell_verts:
            self.cell_vol.append(p.cell_height*tb.area(poly))

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
        self.gap_jun_i          A list of index pairs to self.cell_centres, each pair defining a unique cell-cell GJ
        self.cell2GJ_map        Returns a list of indices to gap junctions for each cell index

        Notes
        -------
        Uses numpy arrays
        Uses scipy spatial KDTree search algorithm

        """

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

        self.cell2GJ_map = []

        for i, neighs in enumerate(self.cell_nn):
            holdgj = []
            for j, ns in enumerate(neighs):
                if i < ns:
                    gj_ind = self.gap_jun_i.index([i,ns])
                    holdgj.append(gj_ind)
                if i > ns:
                    gj_ind = self.gap_jun_i.index([ns,i])
                    holdgj.append(gj_ind)
            self.cell2GJ_map.append(holdgj)

        self.gap_jun_i = np.asarray(self.gap_jun_i)

    def clean_ecm(self,clean=None):

        """
        Calculates ecm points on the environmental boundary using the alpha-shape concave hull method,
        deletes these points, and updates data structures referring to the ecm vertices.

        Note: if clean = 'no' this function returns data structures without deleting ecm verticies.

        """

        # get the concave hull for ecm vertices on the outer boundary of the cluster

        if clean == None or clean == 'yes':

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

        if clean == None or clean == 'yes':
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

    def boundTag(self,points):

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

        con_hull = tb.alpha_shape(points, p.scale_alpha/p.d_cell)  # get the concave hull for the membrane midpoints
        con_hull = np.asarray(con_hull)

        bflags = np.unique(con_hull)    # get the value of unique indices from segments

        bmask = np.array((points[:,0]))
        bmask[:] = 0
        bmask[bflags] = 1

        return bflags, bmask

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
                pt_scale.append(p.scale_cell*pt_zero + centre)
            self.cell_verts.append(np.asarray(pt_scale))

        self.cell_verts = np.asarray(self.cell_verts)

    def cellGeo(self,close_ecm=None):
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
        self.ecm_vols
        self.ecm_edges_i
        self.ecm_mids
        self.ecm_vects
        self.cell2ecm_map


        """

        self.cell_vol = []   # storage for cell volumes
        self.cell_sa = []    # whole cell surface areas

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

        for poly in self.cell_verts:
            # First calculate individual cell volumes from cell vertices:
            self.cell_vol.append(p.cell_height*tb.area(poly))
            self.cell_sa.append(2*tb.area(poly) + perim)   # surface area of whole cell [m2]
            # Next calculate individual membrane domains, midpoints, and vectors:
            edge = []
            mps = []
            surfa = []

            for i in range(0,len(poly)):
                pt1 = poly[i-1]
                pt2 = poly[i]
                pt1 = np.asarray(pt1)
                pt2 = np.asarray(pt2)
                edge.append([pt1,pt2])
                mid = (pt1 + pt2)/2       # midpoint calculation
                mps.append(mid)

                lgth = np.sqrt((pt2[0] - pt1[0])**2 + (pt2[1]-pt1[1])**2)  # length of membrane domain
                sa = lgth*p.cell_height    # surface area
                surfa.append(lgth)

                tang_a = pt2 - pt1       # tangent
                tang = tang_a/np.linalg.norm(tang_a)
                normal = np.array([-tang[1],tang[0]])
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
        self.cell_sa = np.asarray(self.cell_sa)

        # Extracellular matrix specific data

        ecm_edge_ind = set()     # this will hold the unique index pairs to the self.ecm_verts_unique [x,y] points

        ecmverts_list = self.ecm_verts_unique.tolist()

        for poly in self.ecm_verts:  # for every polygon defined in the self.ecm_verts data structure

            for i in range(0,len(poly)):   # for every vertex defining the polygon...

                edge_pt1 = poly[i-1].tolist()    # first point of line segment
                edge_pt2 = poly[i].tolist()      # second point of line segment
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
        self.ecm_vols = [0]*len_unique_edges

        ev_x=[0]*len_unique_edges
        ev_y=[0]*len_unique_edges
        ev_tx=[0]*len_unique_edges
        ev_ty=[0]*len_unique_edges

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
                        self.ecm_vols[mapval] = vol
                        tang_a = pnt2 - pnt1
                        tang = tang_a/np.linalg.norm(tang_a)
                        ev_x[mapval] = midpoint[0]
                        ev_y[mapval] = midpoint[1]
                        ev_tx[mapval] = tang[0]
                        ev_ty[mapval] = tang[1]

                    if edge_ind2 < edge_ind1:
                        mapval = self.ecm_edges_i.index([edge_ind2,edge_ind1])
                        holdinds.append(mapval)
                        pnt1 = self.ecm_verts_unique[edge_ind2]
                        pnt2 = self.ecm_verts_unique[edge_ind1]
                        midpoint = (pnt1 + pnt2)/2   # find the midpoint...
                        lgth = np.sqrt((pnt2[0] - pnt1[0])**2 + (pnt2[1]-pnt1[1])**2)  # length of membrane domain
                        vol = lgth*p.cell_height*p.cell_space
                        self.ecm_mids[mapval] = midpoint  # add the midpoint to its list, keeping the same ordering
                        self.ecm_vols[mapval] = vol
                        tang_a = pnt2 - pnt1
                        tang = tang_a/np.linalg.norm(tang_a)
                        ev_x[mapval] = midpoint[0]
                        ev_y[mapval] = midpoint[1]
                        ev_tx[mapval] = tang[0]
                        ev_ty[mapval] = tang[1]

            self.cell2ecm_map.append(holdinds)

        self.ecm_vects = np.array([ev_x,ev_y,ev_tx,ev_ty]).T
        self.ecm_mids = np.array(self.ecm_mids)
        self.ecm_edges_i = np.asarray(self.ecm_edges_i)
        self.ecm_vols = np.asarray(self.ecm_vols)

    def cleanUp(self):

        """
        Nulls unused data structures to free up memory.
        Creates index data structures for unique cell, gap junction and ecm segments (used in simulation during
        randomization of progression through the data structure).

        """

        self.x_v = None
        self.y_v = None
        self.x_2d = None
        self.y_2d = None
        self.clust_xy = None

        self.cell_i = [x for x in range(0,len(self.cell_centres))]
        self.ecm_i = [x for x in range(0,len(self.ecm_edges_i))]
        self.gj_i = [x for x in range(0,len(self.gap_jun_i))]

        self.cell_vol = np.asarray(self.cell_vol)

        self.gjMatrix = np.zeros((len(self.gj_i),len(self.cell_i)))
        for igj, pair in enumerate(self.gap_jun_i):
            ci = pair[0]
            cj = pair[1]
            self.gjMatrix[igj,ci] = -1
            self.gjMatrix[igj,cj] = 1

