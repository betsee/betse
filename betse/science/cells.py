#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS                            }....................
import math
import os
import os.path
import numpy as np
import scipy.spatial as sps
from numpy import ndarray
from scipy import interpolate as interp
from scipy import ndimage
from betse.exceptions import BetseSequenceException, BetseSimConfigException
from betse.lib.numpy import arrays
from betse.science import filehandling as fh
from betse.science.math import finitediff as fd
from betse.science.math import toolbox as tb
from betse.science.tissue.bitmapper import BitMapper
from betse.util.io.log import logs
from betse.util.type.call.memoizers import property_cached
from betse.util.type.types import type_check, SequenceTypes

# ....................{ CLASSES                            }....................
# FIXME create a new option for seed points: Fibonacci radial-spiral array
class Cells(object):
    '''
    High-level tissue simulation object encapsulating the cell population.

    Specifically, this object:

    * Creates and stores data structures relating to the geometric properties
      of the cell cluster grid and (optional) environmental computational grid.
    * Provides functions to facilitate data association and plotting on the
      geometric structures (e.g., cell volumes, membrane surface area, etc).
    * Constructs computational matrices used in simulation for direct
      computation of gradients, divergence, and Laplacians on cell, membrane
      and environmental points (as required).

    Methods
    -------
    fileInit()                        Create directories for file saving and loading
    seed()                       Create a cell cluster for simulation purposes
    deformWorld(p, ecm_verts)         Modifies Cells object under a deformation
    makeSeeds()                       Create a 2D random scatter of points which will serve as cell centres
    makeVoronoi()                     Make and clip/close a Voronoi diagram from the seed points
    cell_index()                      Returns a list of [x,y] points defining the cell centres in order
    cellVerts()                       Copy & scale in points from the ecm matrix to create unique polygonal cells
    quickVerts()                      Reformulates an exisitng cell world
    cellMatrices()                    Calculates for matrices for BETSE
    cell_vols()                       Calculates cell patch volumes
    near_neigh()                      Calculate the nearest neighbour (nn) array for each cell (make gap junctions)
    voronoiGrid()
    makeECM()                         Make the Marker and Cell (MACs) grid for extracellular calculations
    environment()                     Calculate details for the extracellular calculations, including mappings
    graphLaplacian()                  Creates an abstract discrete Laplacian for the irregular Voronoi-based cell grid
    redo_gj()                         Create gap junction connection network after assessing tissue profile requests

    Attributes (Cell)
    ----------
    cell_centres : ndarray
        Two-dimensional Numpy array of the coordinates of the center points of
        all cells, whose:
        . First dimension indexes cells, whose length is the number of cells.
        . Second dimension indexes the coordinates of the center point of the
          current cell, whose length is unconditionally guaranteed to be 2
          *and* whose:
          . First element is the X coordinate of the current cell center.
          . Second element is the Y coordinate of the current cell center.
    cell_i : list
        One-dimensional list indexing each cell such that each element is that
        cell's index (i.e., ``[0, 1, ..., n-2, n-1]`` for the number of cells
        ``n``), required for efficient Numpy slicing.
    cell_verts : ndarray
        Three-dimensional Numpy array of the coordinates of the vertices of all
        cells, whose:
        . First dimension indexes each cell such that each element is a
          matplotlib-compatible **polygon patch** (i.e., a two-dimensional
          Numpy array of all vertex coordinates defining the current cell's
          polygon), suitable for passing as is to the
          :meth:`matplotlib.patches.Polygon.__init__` method.
        . Second dimension indexes each vertex of the current cell (in
          counterclockwise order).
        . Third dimension indexes each coordinate of the current cell vertex,
          whose length is guaranteed to be 2 *and* whose:
          . First element is the X coordinate of the current cell vertex.
          . Second element is the Y coordinate of the current cell vertex.
    R : ndarray
        One-dimensional Numpy array indexing each cell such that each element is
        the radius of that cell.

    Attributes (Cell Membrane)
    ----------
    cell_to_mems : ndarray
        Two-dimensional Numpy array of the indices of all membranes for all
        cells, whose:
        . First dimension indexes each cell.
        . Second dimension indexes each membrane of the current cell such that
          each element is the index of the corresponding membrane in
          membrane-centric arrays (e.g., :attr:`mem_mids_flat`).
        For example:
        * ``cell_to_mems[0][0]`` is the index of the first cell's first membrane
          -- which is *always* 0.
        * ``cell_to_mems[-1][-1]`` is the index of the last cell's last membrane
          -- which is *always* ``mem_i[-1]``.
    mem_to_cells : ndarray
        One-dimensional Numpy array of length the number of cell membranes such
        that each element is the index of the cell containing the membrane
        indexed by that element. For example:
        * ``mem_to_cells[0]`` is the index of the cell containing the first
          membrane.
        * ``mem_to_cells[-1]`` is the index of the cell containing the last
          membrane.
    mem_i : list
        One-dimensional list of length the number of cell membranes such that
        each element is that cell membrane's index (i.e.,
        ``[0, 1, ..., m-2, m-1]`` for the number of cell membranes ``m``),
        required for efficient Numpy slicing.
    mem_mids_flat : ndarray
        Two-dimensional Numpy array of the coordinates of the midpoints of all
        cell membranes, whose:
        . First dimension indexes each cell membrane.
        . Second dimension indexes each coordinate of the midpoint of the
          current cell membrane, whose length is guaranteed to be 2 *and*
          whose:
          . First element is the X coordinate of the current membrane midpoint.
          . Second element is the Y coordinate of the current membrane
            midpoint.
    mem_vects_flat : ndarray
        Two-dimensional Numpy array of the coordinates of various vectors of
        all cell membranes, whose:
        . First dimension indexes each cell membrane.
        . Second dimension indexes each coordinate of a vector describing the
          current cell membrane, whose length is guaranteed to be 6 *and* whose:
          . First element is the X coordinate of the current membrane midpoint,
            equivalent to :attr:`mem_mids_flat[mem_i,0]` where ``mem_i`` is the
            index of the current membrane.
          . Second element is the Y coordinate of the current membrane
            midpoint, equivalent to :attr:`mem_mids_flat[mem_i,1]`.
          . Third element is the X coordinate of the normal unit vector
            orthogonal to the tangent unit vector of the current membrane
            (defined by the corresponding element of the fifth and sixth
            elements of this array).
          . Fourth element is the Y coordinate of this normal unit vector.
          . Fifth element is the X coordinate of the tangent unit vector
            parallel to the vector implied by the pair of coordinates defining
            the current membrane.
          . Sixth element is the Y coordinate of this tangent unit vector.
    num_mems : ndarray
        One-dimensional Numpy array indexing each cell such that each element
        is the number of cell membranes contained by the current cell.
    M_sum_mems : ndarray
        Numpy matrix (i.e., two-dimensional array) of size ``m x n``, where:
        * `m` is the total number of cells.
        * `n` is the total number of cell membranes.
        For each cell ``i`` and membrane ``j``, element ``M_sum_mems[i, j]`` is:
        * 0 if this cell does *not* contain this membrane. Since most cells do
          *not* contain most membranes, most entries of this matrix are zero,
          implying this matrix to typically (but *not* necessarily) be sparse.
        * 1 if this cell contains this membrane.
        The dot product of this matrix by a Numpy vector (i.e., one-dimensional
        array) of size ``n`` containing cell membrane-specific data yields
        another Numpy vector of size ``m`` containing cell-specific data
        totalized for each cell over all membranes this cell contains, where
        ``m`` and ``n`` are as defined above.

    Attributes (Cell Membrane Vertices)
    ----------
    index_to_mem_verts : ndarray
        Two-dimensional Numpy array of the indices of the vertices of all cell
        membranes in the :attr:`mem_verts` array, whose:
        . First dimension indexes each cell membrane.
        . Second dimension indexes the index of each vertex defining the
          current membrane, whose length is unconditionally guaranteed to be 2
          _and_ whose:
          . First element is the index of the first vertex defining this
            membrane in the :attr:`mem_verts` subarray, guaranteed to be
            counterclockwise from the second vertex defining this membrane.
          . Second element is the index of the second vertex defining this
            membrane in the :attr:`mem_verts` subarray, guaranteed to be
            clockwise from the first vertex defining this membrane.
    mem_verts : ndarray
        Two-dimensional Numpy array of the coordinates of the vertices of all
        cell membranes, whose:
        . First dimension indexes cell membrane vertices, whose length is the
          number of cell membrane vertices in this cluster. This length is
          strictly greater than the number of cell membranes in this cluster.
          Since this array contains _no_ duplicate vertices, this length is
          strictly less than twice the number of cell membranes.
        . Second dimension indexes the coordinates of the current membrane
          vertex, whose length is unconditionally guaranteed to be 2 _and_
          whose:
          . First element is the X coordinate of the current membrane vertex.
          . Second element is the Y coordinate of the current membrane vertex.

    Attributes (Cell Gap Junctions)
    ----------
    cell_nn_i : ndarray
        Two-dimensional Numpy array of the indices of the pairs of adjacent
        cells comprising all gap junctions connecting one cell to another in
        this cluster, whose:
        . First dimension indexes cell membranes, whose length is the number of
          cell membranes in this cluster.
        . Second dimension indexes the indices of the pair of cells defining
          the gap junction connecting the current membrane to another membrane,
          whose length is unconditionally guaranteed to be 2 _and_ whose:
          . First element is the index of the cell containing the current
            membrane.
          . Second element is the index of an adjacent cell containing the
            adjacent membrane to which the current membrane is connected.
    gj_len : float
        Uniform length (in meters) of each gap junction, equivalent to the
        uniform distance between each pair of neighbouring cells.
    mem_nn : ndarray
        Two-dimensional Numpy array of the indices of the pairs of adjacent
        cell membranes comprising all gap junctions connecting one membrane to
        another in this cluster, whose:
        . First dimension indexes cell membranes, whose length is the number of
          cell membranes in this cluster.
        . Second dimension indexes the indices of the pair of cell membranes
          defining the gap junction connecting the current membrane to another
          membrane, whose length is unconditionally guaranteed to be 2 _and_
          whose:
          . First element is the index of the current membrane.
          . Second element is:
            * If the current membrane is adjacent to no other membrane (e.g.,
              due to being situated at either the periphery of this cluster
              _or_ a discontiguous hole in this cluster), the index of the
              current membrane again. A solitary membrane thus connects to
              itself with a self-referential gap junction "loop."
            * If the current membrane is adjacent to exactly one other membrane
              (which is the common case), the index of that membrane.
            * If the current membrane is adjacent to two or more other
              membranes, the index of the membrane to which the current
              membrane is most adjacent (excluding itself).
        Hence, _all_ membranes are guaranteed to participate in exactly one
        gap junction connection.
    nn_i : ndarray
        One-dimensional Numpy array of length the number of cell membranes such
        that each element is the index of the nearest neighbouring cell
        membrane of the cell membrane indexed by that element, comprising the
        gap junction connecting these membranes. Each element thus maps each
        membrane to its nearest partner. For example:
        * `nn_i[0]` is the index of the membrane nearest to the first membrane.
        * `nn_i[-1]` is the index of the membrane nearest to the last membrane.
        Indexing a Numpy array of length the number of cell membranes providing
        data spatially situated at these membranes by this array yields another
        array of the same length providing the corresponding data spatially
        situated at their adjacent membranes.

    Attributes (Cell Surface Area)
    ----------
    cell_sa : ndarray
        One-dimensional Numpy array of length the number of cells such that each
        element is the surface area of the cell indexed by that element,
        defined as the summation of the surface areas of all membranes
        comprising that cell.
    mem_sa : ndarray
        One-dimensional Numpy array of length the number of cell membranes such
        that each element is the surface area of the cell membrane indexed by
        that element.

    Attributes (Extracellular Grid)
    ----------
    The following attributes are defined *only* if the :meth:`makeECM` has been
    called, implying the current simulation to enable extracellular spaces. If
    that method has *not* yet been called, these attributes remain undefined.

    delta : NumericTypes
        Distance in meters between between each extracellular grid point *and*
        between each extracellular grid space, uniformally applied in both the
        X and Y dimensions.
    grid_obj : fd.FiniteDiffSolver
        Finite difference solver defining the extracellular grid.
    index_k : list
        List of the indices of all extracellular grid spaces, equivalent to
        ``[0, 1, ..., len(self.xypts)-2, len(self.xypts)-1]``.
    map_ij2k : ndarray
        See the :attr:`fd.FiniteDiffSolver.map_ij2k_cents` array for further
        details.
    X : ndarray
        Two-dimensional Numpy array of the X coordinates of the centres of all
        extracellular grid spaces. See the :attr:`fd.FiniteDiffSolver.cents_X`
        array for further details.
    Y : ndarray
        Two-dimensional Numpy array of the Y coordinates of the centres of all
        extracellular grid spaces. See the :attr:`fd.FiniteDiffSolver.cents_Y`
        array for further details.
    xypts : ndarray
        Two-dimensional Numpy array of the Cartesian coordinates of the centres
        of all extracellular grid spaces. See the
        :attr:`fd.FiniteDiffSolver.xy_cents` array for further details.

    Attributes (Voronoi Diagram)
    ----------
    cell_to_grid : ndarray
        One-dimensional Numpy array indexing each cell such that each element
        is the index of the Voronoi region whose vertices most closely
        spatially align with those of that cell. For example:
        * `cell_to_grid[0]` is the index of the region most closely spatially
          aligned with the first cell.
        * `cell_to_grid[-1]` is the index of the region most closely spatially
          aligned with the last cell.
        Assigning a Numpy array of length the number of regions indexed by this
        array from a Numpy array of length the number of cells maps the
        cell-specific data defined by the latter into region-specific data
        (e.g., `regions_data[cells.cell_to_grid] = cells_centre_data`).
    voronoi_centres : ndarray
        Two-dimensional Numpy array of the coordinates of the center points of
        all polygonal regions in the Voronoi diagram producing this cell
        cluster, whose:
        . First dimension indexes regions in this diagram, whose length is the
          total number of regions in this diagram.
        . Second dimension indexes the coordinates of the center point of the
          current region, whose length is unconditionally guaranteed to be 2
          _and_ whose:
          . First element is the X coordinate of the current region center.
          . Second element is the Y coordinate of the current region center.
    voronoi_grid : ndarray
        Two-dimensional Numpy array of the vertex coordinates of all
        polygonal regions in the Voronoi diagram producing this cell cluster,
        whose:
        . First dimension indexes the vertices of all regions in this diagram,
          whose length is the total number of vertices in this entire diagram.
        . Second dimension indexes the coordinates of the current vertex, whose
          length is unconditionally guaranteed to be 2 _and_ whose:
          . First element is the X coordinate of the current region vertex.
          . Second element is the Y coordinate of the current region vertex.
        This array is equivalent to the :attr:`voronoi_verts` array flattened
        over the first dimension of that array.
    voronoi_verts : ndarray
        Three-dimensional Numpy array of the vertex coordinates of all
        polygonal regions in the Voronoi diagram producing this cell cluster,
        whose dimensions are structured as those of the :attr:`cell_verts`
        attribute (replacing "cell" with "Voronoi region").
    '''

    def __init__(self, p):
        '''
        Initialize this cell cluster.
        '''

        # Extract the constants from the input object.
        self.fileInit(p)

    def fileInit(self, p):
        """
        Initializes file saving and loading directory as the BETSE cache, which
        is automatically assigned from the user-specified path in the
        configuration file.
        """

        # Make the BETSE-specific cache directory if not found.
        betse_cache_dir = os.path.expanduser(p.init_path)
        os.makedirs(betse_cache_dir, exist_ok=True)

        # Define data paths for saving an initialization and simulation run:
        self.savedWorld = os.path.join(betse_cache_dir, p.world_filename)

    def makeWorld(self, p):
        """
        Calls internal methods to set up the cell cluster.
        """

        logs.log_info('Creating Voronoi geometry... ')
        self.makeSeeds(p, seed_type=p.lattice_type)  # Create the grid for the system (irregular)
        self.makeVoronoi(p)  # Make, close, and clip the Voronoi diagram

        logs.log_info('Cleaning Voronoi geometry... ')
        self.cleanVoronoi(p)  # clean the Voronoi diagram of empty data structures

        logs.log_info('Defining cell-specific geometric properties... ')
        self.cell_index(p)  # Calculate the correct centre and index for each cell
        self.cellVerts(p)  # create individual cell polygon vertices

        logs.log_info('Creating computational matrices for cell-cell transfers... ')
        self.cellMatrices(p)  # creates a variety of matrices used in routine cells calculations
        # self.intra_updater(p)    # creates matrix used for finite volume integration on cell patch

        self.cell_vols(p)   # calculate the volume of cell and its internal regions
        # self.cellDivM(p)    # create matrix to invert divergence
        # self.memLaplacian()   # creates an inverse matrix to calculate voltage on individual membranes

        logs.log_info('Creating gap junctions... ')
        self.mem_processing(p)  # calculates membrane nearest neighbours, ecm interaction, boundary tags, etc
        self.near_neigh(p)  # Calculate the nn array for each cell
        self.voronoiGrid(p)

        logs.log_info('Setting global environmental conditions... ')
        self.makeECM(p)  # create the ecm grid
        self.environment(p)  # define features of the ecm grid
        self.grid_len = len(self.xypts)
        self.make_maskM(p)
        # self.env_weighting(p)

        if p.sim_ECM is True:

            self.sim_ECM = True

            # logs.log_info("Creating Maxwell Capacitance Matrix voltage solver for cell cluster...")
            # self.maxwellCapMatrix(p)  # create Maxwell Capacitance Matrix solver for voltages

            logs.log_info('Creating environmental Poisson solver for voltage...')
            bdic = {'N': 'value', 'S': 'value', 'E': 'value', 'W': 'value'}
            self.lapENV, self.lapENVinv = self.grid_obj.makeLaplacian(bound=bdic)
            self.lapENV = None  # get rid of the non-inverse matrix as it only hogs memory...

            # logs.log_info('Creating environmental Poisson solver for currents...')
            # bdic = {'N': 'flux', 'S': 'flux', 'E': 'flux', 'W': 'flux'}
            # self.lapENV_P, self.lapENV_P_inv = self.grid_obj.makeLaplacian(bound=bdic)
            #
            # self.lapENV_P = None  # get rid of the non-inverse matrix as it only hogs memory...

            # logs.log_info('Creating environmental Screened Poisson Equation solver...')
            # self.lapENVScreen, self.lapENVScreen_inv = self.grid_obj.makeScreenedLaplacian(ko = 1.0e9)

            # logs.log_info('Creating finite volume grid integrator...')
            # self.gridInt = self.grid_obj.makeIntegrator()

        else:

            self.sim_ECM = False

        # set all Laplacian matrices to None fields to allow for flexible creation
        # Laplacians and inverses on the cell grid (two boundary conditions sets)
        self.lapGJinv = None
        self.lapGJ_P_inv = None
        self.lapGJ = None
        self.lapGJ_P = None

        # Lapalcian inverses on the env grid
        # self.lapENVinv = None

        # other matrices
        self.M_sum_mem_to_ecm = None   # used for deformation
        self.gradMem = None  # used for electroosmosis

        self.cell_i = np.asarray(self.cell_i) # we need this to be an array for advanced indexing & assignments

    def deformWorld(self,p, ecm_verts):
        """
        Runs necessary methods to recalculate essential world
        data structures after a mechanical deformation.

        Note: each cell is assumed to be incompressible by default. Therefore, cell
        volumes and total surface area are not updated in a deformation routine.

        """

        # begin by creating new cell centres from the passed Voronoi patch vertices:

        self.cell_centres = np.array([0,0])

        for poly in ecm_verts:
            aa = np.asarray(poly)
            aa = np.mean(aa,axis=0)
            self.cell_centres = np.vstack((self.cell_centres,aa))

        self.cell_centres = np.delete(self.cell_centres, 0, 0)

        #---------------------------------------------------

        #FIXME: Uh oh! We're in Don't Repeat Yourself (DRY) trouble here. All
        #or most of the following code appears to be duplicated verbatim from
        #the cellVerts() method. Since this code is non-trivial and hence
        #fragile, this is bad... Would it be feasible to simply call
        #cellVerts() here? If not, how can we unify this logic? Max axolotl!

        # calculate basic properties such as volume, surface area, normals, etc for the cell array
        self.cell_verts = []

        for centre, poly in zip(self.cell_centres, ecm_verts):
            pt_scale = []
            for vert in poly:
                pt_zero = vert - centre
                pt_scale.append(p.scale_cell * pt_zero + centre)
            self.cell_verts.append(np.asarray(pt_scale))

        self.cell_verts = np.asarray(self.cell_verts)

        mem_edges = []  # storage for membrane edge points
        mem_length = []  # storage for membrane surface area values
        mem_mids = []  # storage for membrane midpoints

        # storage for various vector properties of membrane
        cv_x = []
        cv_y = []
        cv_nx = []
        cv_ny = []
        cv_tx = []
        cv_ty = []

        for polyc in self.cell_verts:

            #FIXME: Wouldn't the "mps" list of membrane midpoints be more
            #efficiently calculatable via a single
            #"mps = np.mean(mps, axis=0)" assignment rather than iteratively
            #constructed via looping?

            # Calculate individual membrane domains, midpoints, and vectors:
            edge = []
            mps = []
            surfa = []

            for i in range(0, len(polyc)):
                pt1 = polyc[i - 1]
                pt2 = polyc[i]
                pt1 = np.asarray(pt1)
                pt2 = np.asarray(pt2)
                edge.append([pt1, pt2])
                mid = (pt1 + pt2) / 2  # midpoint calculation
                mps.append(mid.tolist())

                lgth = np.sqrt((pt2[0] - pt1[0]) ** 2 + (pt2[1] - pt1[1]) ** 2)  # length of membrane domain

                # sa = lgth * p.cell_height  # surface area of membrane
                surfa.append(lgth)

                tang_a = pt2 - pt1  # tangent
                tang = tang_a / np.linalg.norm(tang_a)
                # normal = np.array([-tang[1],tang[0]])
                normal = np.array([tang[1], -tang[0]])
                cv_x.append(mid[0])
                cv_y.append(mid[1])
                cv_nx.append(normal[0])
                cv_ny.append(normal[1])
                cv_tx.append(tang[0])
                cv_ty.append(tang[1])

            mem_edges.append(edge)
            mem_mids.append(mps)
            mem_length.append(surfa)

        self.mem_vects_flat = np.array([cv_x, cv_y, cv_nx, cv_ny, cv_tx, cv_ty]).T

        # ---post processing and calculating peripheral structures-----------------------------------------------------

        self.mem_mids_flat, indmap_mem, _ = tb.flatten(mem_mids)
        self.mem_mids_flat = np.asarray(self.mem_mids_flat)  # convert the data structure to an array

        # Finish up by creating indices vectors and converting to Numpy arrays where needed:

        #FIXME: Reduce the following two assignments to simply:
        #
        #    self.cell_i = np.asarray(range(len(self.cell_centres)))
        #    self.mem_i =  np.asarray(range(len(self.mem_mids_flat)))
        #
        #No need for intermediate list comprehensions. Fire-beast flagons!

        self.cell_i = [x for x in range(0, len(self.cell_centres))]
        self.mem_i =  [x for x in range(0, len(self.mem_mids_flat))]

        # convert mem_length into a flat vector
        mem_length, _, _ = tb.flatten(mem_length)
        mem_length = np.asarray(mem_length)

        self.mem_sa = mem_length * p.cell_height

        self.mem_edges_flat, _, _ = tb.flatten(mem_edges)
        self.mem_edges_flat = np.asarray(self.mem_edges_flat)

        # create a flattened version of cell_verts that will serve as membrane verts:
        self.mem_verts, _, _ = tb.flatten(self.cell_verts)
        self.mem_verts = np.asarray(self.mem_verts)

        # structures for plotting interpolated data and streamlines:
        self.plot_xy = np.vstack((self.mem_mids_flat, self.mem_verts))

        # cell surface area:
        self.cell_sa = []
        for grp in self.cell_to_mems:
            cell_sa = sum(self.mem_sa[grp])
            self.cell_sa.append(cell_sa)

        self.cell_sa = np.asarray(self.cell_sa)

        #------------------------------------------------------
        # next obtain the set of *unique* vertex points from the total ecm_verts arrangement:
        ecm_verts_flat,_,_ = tb.flatten(ecm_verts)

        ecm_verts_set = set()

        for vert in ecm_verts_flat:
            ptx = vert[0]
            pty = vert[1]
            ecm_verts_set.add((ptx,pty))

        ecm_verts_unique = [list(verts) for verts in list(ecm_verts_set)]

        ecm_verts_unique = np.asarray(ecm_verts_unique)  # convert to numpy array

        # redo world boundaries used in plotting, if necessary:
        xmin = np.min(ecm_verts_unique[:,0])
        xmax = np.max(ecm_verts_unique[:,0])
        ymin = np.min(ecm_verts_unique[:,1])
        ymax = np.max(ecm_verts_unique[:,1])

        if xmin < self.xmin:
            self.xmin = xmin

        if xmax > self.xmax:
            self.xmax = xmax

        if ymin < self.ymin:
            self.ymin = ymin

        if ymax > self.ymax:
            self.ymax = ymax

    def makeSeeds(self, p: 'Parameters', seed_type = 'rect') -> None:

        '''
        Creates the irregular scatter lattice of seed points defined on a 2D
        world space, with dimensions supplied by `p.wsx` in [m].

        Specifically, this method defines the following object attributes:

        * `xmin`, `xmax`, `ymin`, `ymax`, the minimum and maximum points of the
          global world space.
        * `centre`, the centre of the global world space.
        * `clust_xy`, the list of `(x, y)` Cartesian coordinates of seed
          points.

        The amount of deviation from a square grid is specified by `p.nl`,
        defined from 0 (perfect square grid) to 1 (full noise).

        Parameters
        -----------
        p : Parameters
            Current simulation configuration.
        seed_type     'hex' produces a hexagonal lattice seed, while 'rect' produces a
                       rectangular grid. Both are perturbed by the lattice noise option p.nl

        '''

        if seed_type == 'rect': # prepare a standard rectangular grid of points

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

            xypts = np.vstack((x_2d.ravel(), y_2d.ravel())).T

        elif seed_type == 'hex': # prepare a hexagonal grid of points

            # recalculate the number of lattice sites

            # r = (p.d_cell/2)*np.sqrt(3)
            r = p.d_cell/2

            # blue and red simply refer to the different seed point arrays used to
            # create the hexagonal lattice.

            a_hex = np.sqrt(3) * r  # distance flat radius (y-axis hex radius)
            b_hex = 3 * r  # x-spacing between midpoint of blue and red grids

            p.wsy = p.wsx

            blue_y_start = a_hex / 2
            blue_y_end = p.wsx - a_hex / 2
            blue_x_start = 0
            blue_x_end = p.wsy

            blue_space_y = a_hex
            blue_space_x = b_hex

            blue_x_num = int((blue_x_end - blue_x_start) / blue_space_x)
            blue_y_num = int((blue_y_end - blue_y_start) / blue_space_y)

            red_y_start = 0
            red_y_end = p.wsx
            red_x_start = b_hex / 2
            red_x_end = p.wsy - b_hex / 2

            red_space_y = a_hex
            red_space_x = b_hex

            red_x_num = int((red_x_end - red_x_start) / red_space_x)
            red_y_num = int((red_y_end - red_y_start) / red_space_y)

            # first begin with linear vectors which are the "ticks" of the x and y dimensions
            x_blue = np.linspace(blue_x_start, blue_x_end, blue_x_num)  # create lattice vector x
            y_blue = np.linspace(blue_y_start, blue_y_end, blue_y_num)  # create lattice vector y

            # create the staggered points for alternative hex-cell centres:
            x_red = np.linspace(red_x_start, red_x_end, red_x_num)  # create lattice vector x
            y_red = np.linspace(red_y_start, red_y_end, red_y_num)  # create lattice vector y

            # next define a 2d array of lattice points using the x- and y- vectors
            x_2d_blue, y_2d_blue = np.meshgrid(x_blue, y_blue)  # create 2D array of lattice points
            x_2d_red, y_2d_red = np.meshgrid(x_red, y_red)  # create 2D array of lattice points

            # interleave the arrays:
            x_2d = np.hstack((x_2d_blue.ravel(), x_2d_red.ravel()))
            y_2d = np.hstack((y_2d_blue.ravel(), y_2d_red.ravel()))

            # now create a matrix of points that will add a +/- deviation to each point centre
            x_rnd = p.nl * p.d_cell * (np.random.rand(len(x_2d)) - 0.5)  # create a mix of random deltas x dir
            y_rnd = p.nl * p.d_cell * (np.random.rand(len(y_2d)) - 0.5)  # create a mix of random deltas x dir

            # # add the noise effect to the world point matrices and redefine the results
            x_2d = x_2d + x_rnd
            y_2d = y_2d + y_rnd

            xypts = np.vstack((x_2d, y_2d)).T

        else:

            raise BetseSimConfigException("You have chosen a world option lattice type that is "
                                           "not available. Options are 'hex' and 'rect'."
                                           " Please check your config file and try again.")

        # define a data structure that holds [x,y] coordinate points of each 2d grid-matrix entry


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

        The Voronoi diagram is then closed at the square boundaries of the
        global environment. Finally, the cells of the Voronoi diagram are
        removed to define a cluster shape.

        Parameters
        ----------
        p : Parameters
            An instance of the Parameters object.

        Creates
        ---------
        self.ecm_verts              x,y points of Voronoi cell vertices (nulled at world creation endpoint)
        self.ecm_verts_unique       x,y points of unique Voronoi cell vertices (nulled at world creation endpoint)
        self.ecm_polyinds           indices into self.ecm_verts_unique list defining each Voronoi polygon
                                    (nulled at world creation endpoint)
        self.cluster_mask           Matrix of booleans defining masked shape of cell cluster
        self.msize                  Size of bitmap (side pixel number)
        """


        # define the Voronoi diagram from the seed points:
        vor = sps.Voronoi(self.clust_xy)

        # round the x,y values of the vertices so that duplicates aren't formed when we use search algorithms later:
        vor.vertices = np.round(vor.vertices,6)

        # calculate the centre of the diagram
        cluster_center = vor.points.mean(axis=0)

        # complete the Voronoi diagram by adding in undefined vertices to ridges and regions
        i = -1   # enumeration index

        for pnt_indx, vor_edge in zip(vor.ridge_points, vor.ridge_vertices):
            vor_edge = np.asarray(vor_edge)

            i = i+1 # update the count-through index

            if np.any(vor_edge == -1): # if either of the two ridge values are undefined (-1)

                # find the ridge vertice that's not equal to -1
                    new_edge = vor_edge[vor_edge != -1][0]
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

        # Clip the Voronoi cluster to the shape of the clipping bitmap -------------------------------------------------
        self.ecm_verts = [] # voronoi verts of clipped cluster

        self.voronoi_verts = []  # track all voronoi cells, even those not in cluster (used as grid for masking)

        # Load the bitmap used to clip the cell cluster and create a clipping function:
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
                    self.ecm_verts.append(cell_polya) # This makes a jagged boundary of cells

                # if point_check.any() == 1.0: # if any of the region's points are in the clipping func range
                #     self.ecm_verts.append(cell_polya)   # this makes a more solid boundary of cells

        self.cluster_mask = self.bitmasker.clipping_matrix  # keep track of cluster mask and its size
        self.msize = self.bitmasker.msize

        # next obtain the set of *unique* vertex points from the total ecm_verts arrangement:
        ecm_verts_flat,_,_ = tb.flatten(self.ecm_verts)

        ecm_verts_set = set()

        for vert in ecm_verts_flat:
            ptx = vert[0]
            pty = vert[1]
            ecm_verts_set.add((ptx,pty))

        #FIXME: The following three lines safely reduce to the one-liner:
        #
        #    # Action, Jackson!
        #    self.ecm_verts_unique = np.asarray(list(ecm_verts_set))
        self.ecm_verts_unique = [
            list(ecm_verts) for ecm_verts in list(ecm_verts_set)]
        self.ecm_verts_unique = np.asarray(self.ecm_verts_unique)  # convert to numpy array

    def cleanVoronoi(self,p):
        """
        Removes empty data structures from the Voronoi object
        and defines a set of unique ecm vertices.

        """
        # first need to go through and make sure each patch has unique vertex points

        #-----clipping out redundant points--------------------------------------------------------

        # perm_cut = 2 * math.pi * p.rc * p.merge_cut_off  # the threshhold edge length
        #
        ecm_verts_2 = []
        #
        for poly in self.ecm_verts:  # step through each closed and clipped region of the Voronoi
            hold_verts = []

            for i, vert in enumerate(poly):

                xo = poly[i - 1][0]
                yo = poly[i - 1][1]
                x1 = vert[0]
                y1 = vert[1]

                length = math.sqrt((x1 - xo) ** 2 + (y1 - yo) ** 2)

                if length > 0.0:
                    hold_verts.append(poly[i - 1])

                # hold_verts = np.asarray(hold_verts)

            ecm_verts_2.append(hold_verts)

        self.ecm_verts = ecm_verts_2

        #--------------------------------------------------------------------------------------------------------

        # next redefine the set of unique vertex points from ecm_verts arrangement
        ecm_verts_flat, _, _ = tb.flatten(self.ecm_verts)

        ecm_verts_set = set()

        for vert in ecm_verts_flat:
            ptx = vert[0]
            pty = vert[1]
            ecm_verts_set.add((ptx, pty))

        self.ecm_verts_unique = [list(verts) for verts in list(ecm_verts_set)]

        self.ecm_verts_unique = np.asarray(self.ecm_verts_unique)  # convert to numpy array

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

        for poly in self.ecm_verts:
            aa = np.asarray(poly)
            aa = np.mean(aa,axis=0)
            self.cell_centres = np.vstack((self.cell_centres,aa))

        self.cell_centres = np.delete(self.cell_centres, 0, 0)

    def cellVerts(self,p):
        """
        Calculate the true vertices of each individual cell from the extracellular matrix (ecm) vertices
        of the closed & clipped Voronoi diagram.

        The BETSE cell grid has each cell defined by unique vertices, which are scaled in from the ecm points.

        Creates
        -------
        self.cell_verts      A nested python list of the [x,y] point pairs defining vertices of each individual cell
                            polygon. The points of each polygon are arranged in a counterclockwise manner.
        self.mem_verts          Flat list of self.cell_verts
        self.cell_vol       Volume of each cell [m3]
        self.cell_sa        Whole cell surface area

        self.mem_edges_flat          membrane edge points [x,y] coordinates
        self.mem_length         membrane surface area values [m2]
        self.mem_mids_flat      unraveled list of membrane edge midpoints [x,y] arranged to self.mem_i
        self.mem_vects_flat     list of normal and tangent vectors (non nested) arranged to self.mem_i


        Notes
        -------
        The Voronoi diagram returns a connected graph. For this simulation, each cell needs unique vertices and edges.
        This method takes the vertices of the original diagram and scales them in to make unique cells.


        """

        self.gj_len = p.cell_space      # distance between gap junction (as "pipe length")

        # calculate basic properties such as volume, surface area, normals, etc for the cell array
        self.cell_verts = []

        for centre,poly in zip(self.cell_centres,self.ecm_verts):
            pt_scale = []
            for vert in poly:
                pt_zero = vert - centre
                pt_scale.append(p.scale_cell*pt_zero + centre)
            self.cell_verts.append(np.asarray(pt_scale))

        self.cell_verts = np.asarray(self.cell_verts)

        # self.cell_vol = []   # storage for cell volumes
        self.cell_sa = []    # whole cell surface areas
        self.cell_area = []

        mem_edges = []  # storage for membrane edge points
        mem_length = []   # storage for membrane surface area values
        mem_mids = []   # storage for membrane midpoints

        # storage for various vector properties of membrane
        cv_x=[]
        cv_y=[]
        cv_nx=[]
        cv_ny=[]
        cv_tx=[]
        cv_ty=[]

        for polyc in self.cell_verts:
            # # First calculate individual cell volumes from cell vertices:
            # poly = [x for x in reversed(polyc)]
            # self.cell_vol.append(p.cell_height*tb.area(poly))

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
                # sa = lgth*p.cell_height    # surface area of membrane
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

            mem_edges.append(edge)
            mem_mids.append(mps)
            mem_length.append(surfa)

        #FIXME: For readability, it would be great if we could extract the
        #last four columns of this array into two new arrays with
        #human-readable names resembling the "mem_mids_flat" array: e.g.,
        #
        #* "self.mem_norms_flat", providing the normal membrane unit vectors.
        #* "self.mem_tangs_flat", providing the tangent membrane unit vectors.
        #
        #Currently, we reference these columns with non-human-readable magic
        #numbers like "self.mem_vects_flat[:,3]", which is fairly hard to
        #mentally parse when perusing the code. Calm qualms in an oceanic quay!
        #FIXME: The first two columns of this array are exact duplicates of the
        #first (and only) two columns of the "mem_mids_flat" array, defined
        #below. Since the "mem_mids_flat" array is more human-readable than
        #this array, that array should probably be preferred everywhere for
        #obtaining the coordinates of membrane midpoints, in which case the
        #first two columns of this array (i.e., "cv_x" and "cv_y") should
        #probably be removed entirely from this array. Idle Ides of March!
        self.mem_vects_flat = np.array([cv_x,cv_y,cv_nx,cv_ny,cv_tx,cv_ty]).T

        #---post processing and calculating peripheral structures-----------------------------------------------------

        self.mem_mids = np.asarray(mem_mids)

        self.mem_mids_flat, indmap_mem, _ = tb.flatten(mem_mids)
        self.mem_mids_flat = np.asarray(self.mem_mids_flat)  # convert the data structure to an array

        # Finish up by creating indices vectors and converting to Numpy arrays where needed:

        self.cell_i = [x for x in range(0,len(self.cell_centres))]

        self.mem_i = [x for x in range(0,len(self.mem_mids_flat))]

        # convert mem_length into a flat vector
        mem_length,_,_ = tb.flatten(mem_length)
        mem_length = np.asarray(mem_length)

        self.mem_sa = mem_length*p.cell_height

        self.mem_edges_flat, _, _ = tb.flatten(mem_edges)
        self.mem_edges_flat = np.asarray(self.mem_edges_flat)

        # create a flattened version of cell_verts that will serve as membrane verts:
        self.mem_verts,_,_ = tb.flatten(self.cell_verts)
        self.mem_verts = np.asarray(self.mem_verts)

        # structures for plotting interpolated data and streamlines:
        self.plot_xy = np.vstack((self.mem_mids_flat,self.mem_verts))

        # define map allowing a dispatch from cell index to each respective membrane -------------------------------
        indmap_mem = np.asarray(indmap_mem)

        #FIXME: This assignment is safely reducible to:
        #    self.mem_to_cells = indmap_mem[:,0]
        #Indexing by "[self.mem_i]" is superfluous here. Lively jumping snakes!
        self.mem_to_cells = indmap_mem[self.mem_i][:,0]   # gives cell index for each mem_i index placeholder

        # construct a mapping giving membrane index for each cell_i------------------------------------------------
        self.cell_to_mems = []

        for cell_index in self.cell_i:

            #FIXME: This assignment is safely reducible to:
            #    index2mems = (self.mem_to_cells == cell_index).nonzero()[0]
            #Unpacking the contents of this array into a list is unnecessary.
            #Striving leaping laser-sharks!
            index2mems = list(*(self.mem_to_cells == cell_index).nonzero())
            self.cell_to_mems.append(index2mems)

        self.cell_to_mems = np.asarray(self.cell_to_mems)

        #----------------------------------------------------------
        # cell surface area:
        self.cell_sa = []
        for grp in self.cell_to_mems:
            cell_sa = sum(self.mem_sa[grp])
            self.cell_sa.append(cell_sa)

        self.cell_sa = np.asarray(self.cell_sa)

        #----------------------------------------------------------------------
        # Construct an array indexing vertices of the membrane vertices array.

        cellVertTree = sps.KDTree(self.mem_verts)

        self.index_to_mem_verts = []
        for cell_nest in mem_edges:
            for mem_points in cell_nest:
                _, pt_ind1 = cellVertTree.query(mem_points[0])
                _, pt_ind2 = cellVertTree.query(mem_points[1])
                self.index_to_mem_verts.append([pt_ind1, pt_ind2])
        self.index_to_mem_verts = np.asarray(self.index_to_mem_verts)

        # create radial vectors for each cell, defined from their centre to each membrane midpoint
        self.rads = self.mem_mids_flat - self.cell_centres[self.mem_to_cells]

        # magnitude (R) and unit vectors (n) of rads:
        self.R_rads = np.sqrt(self.rads[:,0]**2 + self.rads[:,1]**2)
        self.nx_rads = self.rads[:,0]/self.R_rads
        self.ny_rads = self.rads[:,1]/self.R_rads

    def quickVerts(self, p):

        # calculate basic properties such as volume, surface area, normals, etc for the cell array
        self.cell_verts = []

        for centre,poly in zip(self.cell_centres,self.ecm_verts):
            pt_scale = []
            for vert in poly:
                pt_zero = vert - centre
                pt_scale.append(p.scale_cell*pt_zero + centre)
            self.cell_verts.append(np.asarray(pt_scale))

        self.cell_verts = np.asarray(self.cell_verts)

        mem_edges = []  # storage for membrane edge points
        mem_length = []   # storage for membrane surface area values
        mem_mids = []   # storage for membrane midpoints

        # storage for various vector properties of membrane
        cv_x=[]
        cv_y=[]
        cv_nx=[]
        cv_ny=[]
        cv_tx=[]
        cv_ty=[]

        for polyc in self.cell_verts:
            # # First calculate individual cell volumes from cell vertices:
            # poly = [x for x in reversed(polyc)]
            # self.cell_vol.append(p.cell_height*tb.area(poly))

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
                # sa = lgth*p.cell_height    # surface area of membrane
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

            mem_edges.append(edge)
            mem_mids.append(mps)
            mem_length.append(surfa)

        self.mem_vects_flat = np.array([cv_x,cv_y,cv_nx,cv_ny,cv_tx,cv_ty]).T

        #---post processing and calculating peripheral structures-----------------------------------------------------

        self.mem_mids_flat, indmap_mem, _ = tb.flatten(mem_mids)
        self.mem_mids_flat = np.asarray(self.mem_mids_flat)  # convert the data structure to an array

        # Finish up by creating indices vectors and converting to Numpy arrays where needed:

        self.cell_i = [x for x in range(0,len(self.cell_centres))]

        self.mem_i = [x for x in range(0,len(self.mem_mids_flat))]

        # convert mem_length into a flat vector
        mem_length,_,_ = tb.flatten(mem_length)
        mem_length = np.asarray(mem_length)

        self.mem_sa = mem_length*p.cell_height

        self.mem_edges_flat, _, _ = tb.flatten(mem_edges)
        self.mem_edges_flat = np.asarray(self.mem_edges_flat)

        # create a flattened version of cell_verts that will serve as membrane verts:
        self.mem_verts,_,_ = tb.flatten(self.cell_verts)
        self.mem_verts = np.asarray(self.mem_verts)

        # structures for plotting interpolated data and streamlines:
        self.plot_xy = np.vstack((self.mem_mids_flat,self.mem_verts))

        # cell surface area:
        self.cell_sa = []
        for grp in self.cell_to_mems:
            cell_sa = sum(self.mem_sa[grp])
            self.cell_sa.append(cell_sa)

        self.cell_sa = np.asarray(self.cell_sa)

    def cellMatrices(self, p):
        """
        Creates the main matrices used in routine calculations on the Cell Grid.

        """

        #----------------MATRIX CALCULATIONs----------------------------------------

        # create a matrix that will map and interpolate data on mem mids to the mem verts -----------------------------
        # it will work as data on verts = dot( data on mids, matrixMap2Verts ):
        # self.matrixMap2Verts = np.zeros(
        #     (len(self.mem_mids_flat), len(self.mem_verts)))
        #
        # for i, indices in enumerate(self.index_to_mem_verts):
        #     self.matrixMap2Verts[i, indices[0]] = 1/2
        #     self.matrixMap2Verts[i, indices[1]] = 1/2

        # matrix for summing property on membranes for each cell and a count of number of mems per cell:---------------
        self.M_sum_mems = np.zeros((len(self.cell_i),len(self.mem_i)))
        self.num_mems = []

        for i, inds in enumerate(self.cell_to_mems):
            n = 0
            for j in inds:
                self.M_sum_mems[i,j] = 1
                n = n+1

            self.num_mems.append(n)

        # self.cellDivM(p)

        self.num_mems = np.asarray(self.num_mems)  # number of membranes per cell

        self.mem_distance = p.cell_space + 2*p.tm # distance between two adjacent intracellluar spaces

        self.cell_number = self.cell_centres.shape[0]

        # matrix for calculating gradients around the cell circumference:
        self.gradTheta = np.zeros((len(self.mem_i), len(self.mem_i)))

        for cell_i, mem_i in enumerate(self.cell_to_mems):

            mem_io = np.roll(mem_i, 1)
            li = self.mem_mids_flat[mem_i] - self.mem_mids_flat[mem_io]

            lm = np.sqrt(li[:, 0] ** 2 + li[:, 1] ** 2)

            self.gradTheta[mem_i, mem_i] = 1 / lm
            self.gradTheta[mem_i, mem_io] = -1 / lm


    def memLaplacian(self):

        # matrix for computing divergence of a property defined on a membrane of each cell patch:
        lapGJmem = np.zeros((len(self.mem_i), len(self.mem_i)))

        for cell_i, mem_i in enumerate(self.cell_to_mems):

            num_mems = self.num_mems[cell_i]

            for nj, j in enumerate(mem_i):

                mem_j = np.roll(mem_i, -1 - nj)

                memjj = mem_j[0:-1]

                lapGJmem[j, j] = ((num_mems - 1)/num_mems)*(self.mem_sa[j]/self.mem_vol[j])
                lapGJmem[j, memjj] = -(1/num_mems)*(self.mem_sa[memjj]/self.mem_vol[memjj])

        self.lapGJmem_inv = np.linalg.pinv(lapGJmem)

    def cell_vols(self,p):
        """
        Each cell is divided into a spider-web like volume with an outer
        "pie-box" region volume associated with each membrane, and a region
        around the central point (centroid region).

        This helper-functino calculates volumes of cell, membrane "pie-boxes" and centroids,
        as well as associated parameters of cell chords, which are the distance between the
        cell center and each membrane midpoint.

        """

        # calculate cell chords
        # self.chords = []
        # for i, memMid in enumerate(self.mem_mids_flat):
        #     # get cell index for the membrane:
        #     celli = self.mem_to_cells[i]
        #     dist = memMid - self.cell_centres[celli]
        #     chrd = np.sqrt(dist[0] ** 2 + dist[1] ** 2)
        #     self.chords.append(chrd)
        #
        # self.chords = np.asarray(self.chords)

        # calculate the volume of each membrane's outer "pie-box":

        self.mem_vol = (1 / 2) * self.R_rads * self.mem_sa

        # calaculate cell volume by suming up the large pies:
        self.cell_vol = np.dot(self.M_sum_mems, self.mem_vol)

        self.R = ((3 / 4) * (self.cell_vol / math.pi)) ** (1 / 3)  # effective radius of each cell

        # Finally, create an easy divergence inverse term for cells
        self.diviterm = (self.cell_vol / self.cell_sa)

    def mem_processing(self,p):
        """
        Finds membrane-specific nearest neighbour pairs, interaction with
        the extracellular data points (of voronoi grid) and membranes on the boundary.
        """

        #-- find nearest neighbour cell-cell junctions via adjacent membranes-------------------------------------------

        sc = (p.rc/2.4)*(p.scale_cell)  # threshhold for searching nearest-neighbour membranes
        memTree = sps.KDTree(self.mem_mids_flat)

        mem_nn_o = memTree.query_ball_point(self.mem_mids_flat,sc)
        mem_nn = [[] for x in self.mem_i]
        mem_bound = []

        for i, ind_pair in enumerate(mem_nn_o):
            if len(ind_pair) == 1:
                mem_bound.append(i)
                mem_nn[i].append(i)
                mem_nn[i].append(i)

            elif len(ind_pair) == 2:
                mem_nn[i].append(ind_pair[0])
                mem_nn[i].append(ind_pair[1])

            #FIXME: It'd be great if we could document exactly what and how
            #this algorithm is doing. Are we searching multiple possible
            #neighboring membranes for the nearest neighboring of the current
            #membrane to find the membranes participating in this gap junction?
            elif len(ind_pair) > 2:
                i_n = [self.mem_vects_flat[i,2],self.mem_vects_flat[i,3]]

                for j in ind_pair:
                    a = [self.mem_vects_flat[j,2],self.mem_vects_flat[j,3]]
                    ia = round(np.dot(i_n,a),1)

                    if ia == -1.0:
                        mem_nn[i] = []
                        mem_nn[i].append(i)
                        mem_nn[i].append(j)

                    else:  # in rare cases, tag as self instead of leaving a blank spot:
                        mem_nn[i] =[]
                        mem_nn[i].append(i)
                        mem_nn[i].append(i)

        #---------------------------------------------------------------------------------------------------------------

        self.mem_nn = np.asarray(mem_nn)

        # Tag membranes and cells on the outer boundary of the cell cluster---------------------------------------------
        self.bflags_mems = np.asarray(mem_bound)

         # get the boundary cells associated with these membranes:
        self.bflags_cells = []

        for mem_i in self.bflags_mems:

            cell_i = self.mem_to_cells[mem_i]

            self.bflags_cells.append(cell_i)

        self.bflags_cells = np.asarray(self.bflags_cells)

        # midpoints of extracellular matrix lattice:
        # calculate midpoints of each ecm (voronoi cell) segment:
        ecm_mids = set()

        for verts in self.ecm_verts:
            for i in range(0,len(verts)):
                        pt1 = verts[i-1]
                        pt2 = verts[i]
                        pt1 = np.asarray(pt1)
                        pt2 = np.asarray(pt2)
                        mid = (pt1 + pt2)/2       # midpoint calculation
                        mx = mid[0]
                        my = mid[1]
                        ecm_mids.add((mx,my))

        ecm_mids = list(ecm_mids)
        self.ecm_mids = np.asarray(ecm_mids)

    def near_neigh(self,p):

        """
        Calculate the nearest neighbours for each cell centre in the cluster and return a numpy
        array of nn indices with an index consistent with all other data lists for the cluster.

        Creates
        -------
        self.cell_nn            Indices of all nearest neighbours to each cell (ordered to self.cell_i)
        self.num_nn             Number of nearest neighbours for each cell (ordered to self.cell_i)
        self.average_nn         Average number of nearest neighbours for entire cluster
        self.nn_i                           Non-unique list of index pairs to cells, each pair defining a cell-cell GJ
        self.gj_len                         Length of each GJ [m]
        self.nn_tx, self.nn_ty              Tangent vector coordinates to each gj

        Notes
        -------
        Uses numpy arrays
        Uses scipy spatial KDTree search algorithm

        """

        self.nn_i = [] # gives the partnering membrane index at the vectors' index
        self.cell_nn_i = [[] for x in self.mem_i] # stores the two connecting cell indices at a shared membrane

        for i, (mem_i,mem_j) in enumerate(self.mem_nn):

            if mem_i == mem_j:  # we're on a boundary cell

                self.nn_i.append(i)
                cell_i = self.mem_to_cells[i]
                self.cell_nn_i[mem_i].append(cell_i)
                self.cell_nn_i[mem_i].append(cell_i)

            elif i == mem_i and i != mem_j:

                self.nn_i.append(mem_j)
                cell_i = self.mem_to_cells[mem_i]
                cell_j = self.mem_to_cells[mem_j]

                self.cell_nn_i[i].append(cell_i)
                self.cell_nn_i[i].append(cell_j)

            elif i == mem_j and i != mem_i:

                self.nn_i.append(mem_i)
                cell_i = self.mem_to_cells[mem_j]
                cell_j = self.mem_to_cells[mem_i]

                self.cell_nn_i[i].append(cell_i)
                self.cell_nn_i[i].append(cell_j)

            else:
                logs.log_info("WARNING: entry not placed in seed nearest neighbour construction. "
                                 "Results may not be accurate.")

        self.nn_i = np.asarray(self.nn_i)
        self.cell_nn_i = np.asarray(self.cell_nn_i)

        # Next find the nearest neighbour set for each cell:
        self.cell_nn = []
        for cell_i, mem_i_set in enumerate(self.cell_to_mems):

            cell_neigh_set = []

            for mem_i in mem_i_set:

                mem_j = self.nn_i[mem_i]  # find the partner to this membrane...

                if mem_j == mem_i:  # if the indices are equal, we're on a neighborless boundary cell
                    pass

                else:

                    cell_j = self.mem_to_cells[mem_j]

                    if cell_i != cell_j:  # cross-check that values are not the same
                        cell_neigh_set.append(cell_j)

            self.cell_nn.append(cell_neigh_set)

        self.num_nn = []  # initialize a list that will hold number of nns to a cell

        for indices in self.cell_nn:
            self.num_nn.append(len(indices))

        self.average_nn = (sum(self.num_nn)/len(self.num_nn))

        self.num_nn = np.asarray(self.num_nn)

        self.cell_nn = np.asarray(self.cell_nn)

        # nearest neighbours to the boundary cells:
        nn_bound = self.cell_nn[self.bflags_cells]
        nn_bound, _,_ = tb.flatten(nn_bound)

        self.nn_bound = []
        for ind in nn_bound:  # take out the shared values:

            if ind not in self.bflags_cells:
                self.nn_bound.append(ind)

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
        self.delta = ((self.xmax-self.xmin)/p.grid_size) # spacing between grid points

        self.grid_obj = fd.FiniteDiffSolver()

        self.grid_obj.cell_grid(self.delta,self.xmin,self.xmax,self.ymin,self.ymax)

        self.X = self.grid_obj.cents_X
        self.Y = self.grid_obj.cents_Y

        self.xypts = self.grid_obj.xy_cents
        self.map_ij2k = self.grid_obj.map_ij2k_cents

        #FIXME: The following complex assignment can be reduced to simply:
        #
        #    self.index_k = range(len(self.xypts))
        #
        #Neat-o, eh?

        # linear k index:
        self.index_k = [x for x in range(0,len(self.xypts))]

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


        # global properties of ecm spaces:
        self.ecm_sa = 4*self.delta*p.cell_height # surface area of ecm space in direction of cell flux
        # self.ecm_vol = (p.cell_height*self.delta**2)*np.ones(len(self.xypts))  # volume of ecm space
        self.ecm_vol = (p.cell_height*self.delta**2)  # volume of ecm space

        # ratio of mean cell volume to ecm square volume (gives approx num cells per ecm square)
        self.ratio_cell2ecm = self.ecm_vol/self.cell_vol.mean()

        #-------------------------

        # first obtain a structure to map to total xypts vector index:
        self.points_tree = sps.KDTree(self.xypts)

        # define a mapping between a cell and its ecm space in the full list of xy points for the world:
        _, self.map_cell2ecm = self.points_tree.query(self.cell_centres)
        _, self.map_mem2ecm  = self.points_tree.query(self.mem_mids_flat, k=1)

        # get a list of all membranes for boundary cells:
        all_bound_mem_inds = self.cell_to_mems[self.bflags_cells]
        all_bound_mem_inds, _ ,_ = tb.flatten(all_bound_mem_inds)

        # need these to obtain cluster membrane values from the ECM perspective, or it won't write to the array!
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

        _, self.bBot_k = self.points_tree.query(bBot_pts)
        _, self.bTop_k = self.points_tree.query(bTop_pts)
        _, self.bL_k = self.points_tree.query(bL_pts)
        _, self.bR_k = self.points_tree.query(bR_pts)

        # get a mapping specifying which mem mids an ecm space interacts with:
        self.map_ecm2mem = [[] for ind in self.xypts]

        for ind_mem, ind_ecm in enumerate(self.map_mem2ecm):
            self.map_ecm2mem[ind_ecm].append(ind_mem)

        # next, find out the total set of ecm spaces that interact with membranes
        # and develop the "weight-paint" functions:
        self.envInds_inClust = []
        self.memSa_per_envSquare = np.zeros(len(self.xypts))
        self.mems_per_envSquare = np.zeros(len(self.xypts))

        # create an array to hold the "true" extracellullar volume,
        # and populate it initially with the environmental square volume:
        self.true_ecm_vol = np.ones(len(self.xypts))*self.ecm_vol

        for ind_ecm, lst in enumerate(self.map_ecm2mem):

            if len(lst) > 0:
                self.envInds_inClust.append(ind_ecm)
                sas = np.sum(self.mem_sa[lst])
                self.memSa_per_envSquare[ind_ecm] = sas
                self.mems_per_envSquare[ind_ecm] = len(lst)

                self.true_ecm_vol[ind_ecm] = sas*p.cell_space*(1/2)

        self.envInds_inClust = np.asarray(self.envInds_inClust)

        # correction coefficient for converting from cell to env divergences:
        self.cell2env_corrF = (self.cell_vol / self.true_ecm_vol[self.map_cell2ecm]) * (self.ecm_sa / self.cell_sa)

        # calculate indices to tag TJ at boundary
        neigh_to_bcells, _, _ = tb.flatten(self.cell_nn[self.bflags_cells])
        all_bound_mem_inds_o = self.cell_to_mems[self.bflags_cells]
        interior_bound_mem_inds_o = self.cell_to_mems[neigh_to_bcells]
        interior_bound_mem_inds_o, _, _ = tb.flatten(interior_bound_mem_inds_o)
        all_bound_mem_inds_o, _, _ = tb.flatten(all_bound_mem_inds_o)

        self.all_bound_mem_inds = self.map_mem2ecm[all_bound_mem_inds_o]
        self.interior_bound_mem_inds = self.map_mem2ecm[interior_bound_mem_inds_o]
        self.inds_outmem = self.map_mem2ecm[self.bflags_mems]
        self.ecm_inds_bound_cell = self.map_cell2ecm[self.bflags_cells]

    def graphLaplacian(self,p):
        '''
        Defines an abstract inverse Laplacian that is used to solve Poisson's equation on the
        irregular Voronoi grid of the cell cluster.

        Parameters
        ----------
        p               An instance of the Parameters object



        Creates
        ----------
        self.lapGJinv          Solver for Poisson equation with Dirchlet (zero value) boundary
        self.lapGJ_P_inv       Solver for Poisson equation with Neumann (zero gradient) boundary

        '''

        # zero-value fixed boundary version (Dirchlet condition)
        lapGJ = np.zeros((len(self.cell_i,), len(self.cell_i)))
        # zero-gradient, free boundary version (Neumann condition)
        lapGJ_P = np.zeros((len(self.cell_i,), len(self.cell_i)))

        cell_nn_pairs = self.cell_nn_i.tolist()

        for cell_i, cell_inds in enumerate(self.cell_nn):

            vol = self.cell_vol[cell_i]

            for cell_j in cell_inds:

                # get the distance between the cell centres of the pair:
                lx = self.cell_centres[cell_j, 0] - self.cell_centres[cell_i, 0]
                ly = self.cell_centres[cell_j, 1] - self.cell_centres[cell_i, 1]
                len_ij = np.sqrt(lx ** 2 + ly ** 2)

                # find the shared membrane index for the pair:
                mem_ij = cell_nn_pairs.index([cell_i, cell_j])

                mem_sa = self.mem_sa[mem_ij]

                # correction constant as membrane normals and segments connecting centers aren't quite parallel:
                norm_cor = 0.5*(self.nn_tx[mem_ij]*self.mem_vects_flat[mem_ij, 2] +
                                self.nn_ty[mem_ij]*self.mem_vects_flat[mem_ij, 3])
                # norm_cor = 0.5

                lapGJ[cell_i, cell_i] = lapGJ[cell_i, cell_i] - (1 / (len_ij)) * (mem_sa / vol) * norm_cor

                lapGJ_P[cell_i, cell_i] = lapGJ_P[cell_i, cell_i] - (1 / (len_ij)) * (mem_sa / vol) * norm_cor
                lapGJ_P[cell_i, cell_j] = lapGJ_P[cell_i, cell_j] + (1 / (len_ij)) * (mem_sa / vol) * norm_cor

                if cell_j not in self.bflags_cells:
                    lapGJ[cell_i, cell_j] = lapGJ[cell_i, cell_j] + (1 / (len_ij)) * (mem_sa / vol) * norm_cor

            # deal with boundary values:
            if cell_i in self.bflags_cells:

                lapGJ[cell_i,cell_i] = 0

        self.lapGJinv = np.linalg.pinv(lapGJ)
        self.lapGJ_P_inv = np.linalg.pinv(lapGJ_P)

        if p.td_deform is True:
            # if time0dependent deformation is selected, also save the direct Laplacian operator:
            self.lapGJ = lapGJ
            self.lapGJ_P = lapGJ_P

        # weighting function for the voronoi lattice:
        self.geom_weight = np.dot(self.M_sum_mems, self.mem_sa / self.mem_vol) * p.cell_height

    def cellDivM(self, p):

        """
        Defines a matrix, self.divCell_inv, which
        performs the inverse of divergence operation
        for membranes on individual cells of the cluster.

        If we consider a scalar property Phi, which is
        defined at points immediately inside and outside
        of cell membranes, for the Poisson equation:

         Del^2 Phi = rho

        self.divMem_inv returns Del Phi which is
        the gradient of Phi with respect to Phi's
        values inside and outside of the cell at each
        membrane.

        (only normal components to the membrane
        are considered possible)

        """

        # matrix for computing divergence of a property defined on a membrane of each cell patch:
        divCell = np.zeros((len(self.cell_i), len(self.mem_i)))

        for i, inds in enumerate(self.cell_to_mems):

            # get the volume for the cell:
            cell_vol = self.cell_vol[i]

            for j in inds:
                # self.M_sum_mems[i,j] = 1

                # get the membrane surface area for this membrane patch:
                mem_sa = self.mem_sa[j]

                # set the divergence term:
                divCell[i, j] = mem_sa/cell_vol


        # calculate the inverse of the divergence matrix:
        self.divCell_inv = np.linalg.pinv(divCell)

    def maxwellCapMatrix(self, p):
        """
        This method defines the Maxwell Capacitance matrix
        for the collection of cells with their structured
        Voronoi lattice.

        Each cell and respective extracellular space are
        considered to be conductors separated by the insulating
        region of the cell membrane.

        Each cell interacts with its N nearest ecm spaces.
        Likewise, each ecm space interacts with two neighbouring cells,
        or if the ecm space is on an external boundary, with one
        neighbouring cell.

        The Maxwell Capacitance matrix is created by solving for the
        charge in each cell or ecm space, given the voltages of the space
        and the capacitive connections between spaces. The matrix is
        then inverted, so that we can use it to solve for voltages
        knowing charges.

        """

        envLen = self.envInds_inClust

        data_length = len(self.mem_i) + len(envLen)
        # define ranges within the total data length where we can
        # work with cell centres or ecm mids specifically:
        self.mem_range_a = 0
        self.mem_range_b = len(self.mem_i)
        self.env_range_a = self.mem_range_b
        self.env_range_b = len(envLen) + len(self.mem_i)

        M_max_cap = np.zeros((data_length, data_length))

        # first do cells -- where index of Maxwell vector is equal to cell index
        for mem_i in range(self.mem_range_a, self.mem_range_b):
            cm_mem = p.cm * self.mem_sa[mem_i]  # capacitance for cell diagonal term
            cs_mem = p.electrolyte_screening * self.mem_sa[mem_i]  # calculate self-capacitance

            # get the raw env spaces for each membrane in Max Cap matrix index format:
            raw_env_i = self.map_mem2ecm[mem_i]

            # get the index of the environnmental square:
            env_i_o = int((self.envInds_inClust == raw_env_i).nonzero()[0])

            # convert the index to the Max Cap indexing format:
            env_i = env_i_o + len(self.mem_i)

            # set the diagonal element for cells:
            M_max_cap[mem_i, mem_i] = cm_mem + cs_mem  # membrane plus self-capacitance

            # set the off-diagonal elements for cells:
            M_max_cap[mem_i, env_i] = - cm_mem

        # next do ecm spaces -- index of maxwell vector equal to ecm index - len(cell_i)
        for env_i in range(self.env_range_a, self.env_range_b):
            env_i_o = env_i - len(self.mem_i)  # get the ecm index wrt to the envInds_inClust

            raw_env_i = self.envInds_inClust[env_i_o]  # get the env index wrt to the Env Grid

            mem_set = self.map_ecm2mem[raw_env_i]  # get the set of mem inds corresponding to each ecm space

            cm_set = p.cm * self.mem_sa[mem_set]  # get the capacitance of the individual membranes
            cs_set = p.electrolyte_screening * self.mem_sa[mem_set]  # get the self-capacitance of the ecm spaces

            M_max_cap[env_i, env_i] = cm_set.sum() + cs_set.sum()  # plus self capacitance
            M_max_cap[env_i, mem_set] = -cm_set

        # get the inverse of the matrix:
        self.M_max_cap_inv = np.linalg.pinv(M_max_cap)
        # self.M_max_cap = M_max_cap

    def redo_gj(self,dyna,p,savecells =True):

        # profile_names = list(p.profiles.keys())  # names of each tissue profile...
        profile_names = dyna.tissue_profile_names

        flag_cell_nn = [ [] for x in range(0,len(self.cell_i))]

        for name in profile_names:

            cell_targets = dyna.cell_target_inds[name]   # get the cell target inds for this tissue
            insular_flag = p.profiles[name]['insular gj']

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

        # Redo the number and average nearest neighbours per cell:
        self.num_nn = []  # initialize a list that will hold number of nns to a cell

        for indices in self.cell_nn:
            self.num_nn.append(len(indices))

        self.average_nn = (sum(self.num_nn)/len(self.num_nn))

        self.num_nn = np.asarray(self.num_nn)

        for cell_i, nn_cell_i_set in enumerate(self.cell_nn):

            mem_i_set = self.cell_to_mems[cell_i]  # get all the membranes for this cell

            for mem_i in mem_i_set:

                mem_j = self.nn_i[mem_i]  # get the current neighbour mem and cell...

                if mem_j != mem_i:  # if we're not on a boundary

                    cell_j = self.mem_to_cells[mem_j]

                    if cell_j not in nn_cell_i_set:  # if the partner cell is no longer listed as a nn...

                        #...then set both the membrane and cell neighbour spot to "self":
                        self.nn_i[mem_i] = mem_i
                        self.cell_nn_i[mem_i] = [cell_i,cell_i]

        # calculate gap junction vectors
        self.calc_gj_vects(p)

    def calc_gj_vects(self,p):

        """
        Recalculate nearest neighbour (gap junction)
        vectors.

        Used in deformation sequence.
        """

        self.nn_mids = []

        self.nn_tx = []  # tangent vector to gap junction (through neighboring cell centres)
        self.nn_ty = []


        self.nn_len = []  # distance between neighbouring cell centres

        self.nn_edges = [[] for x in self.mem_i]  # line segment between neighbouring cell centres

        for mem_i, mem_j in enumerate(self.nn_i):

            cell_i, cell_j = self.cell_nn_i[mem_i]

            # calculate vectors for the pairing:
            pt1_mem = self.mem_mids_flat[mem_i]
            pt2_mem = self.mem_mids_flat[mem_j]

            pt1_cell = self.cell_centres[cell_i]
            pt2_cell = self.cell_centres[cell_j]

            tang_o = pt2_mem - pt1_mem

            tang_x_o = tang_o[0]
            tang_y_o = tang_o[1]

            tang_mag = np.sqrt(tang_x_o**2 + tang_y_o**2)

            if tang_mag == 0.0:
                tang_x = 0
                tang_y = 0

            else:

                tang_x = tang_x_o/tang_mag
                tang_y = tang_y_o/tang_mag

            mid = (pt1_mem + pt2_mem)/2
            self.nn_mids.append(mid)

            # calculate length
            len_o = pt2_cell - pt1_cell

            len_xo  = len_o[0]
            len_yo = len_o[1]

            len_mag = np.sqrt(len_xo**2 + len_yo**2)

            if len_mag == 0.0:

                self.nn_len.append(-1)

            else:

                self.nn_len.append(len_mag)

            self.nn_tx.append(tang_x)
            self.nn_ty.append(tang_y)

            self.nn_edges[mem_i].append(pt1_cell)
            self.nn_edges[mem_i].append(pt2_cell)

        self.nn_mids = np.asarray(self.nn_mids)

        self.nn_tx = np.asarray(self.nn_tx)
        self.nn_ty = np.asarray(self.nn_ty)
        self.nn_len = np.asarray(self.nn_len)
        self.nn_edges = np.asarray(self.nn_edges)

        self.cell_nn_tx = []
        self.cell_nn_ty = []

        for cell_i, cell_j in self.cell_nn_i:

            pt1 = self.cell_centres[cell_i]
            pt2 = self.cell_centres[cell_j]

            tang_o = pt2 - pt1
            norm_tang = np.sqrt(tang_o[0]**2 + tang_o[1]**2)

            if norm_tang != 0:
                tang = tang_o/norm_tang

            else:
                norm_tang = 1
                tang = tang_o/norm_tang
                tang[0] = 0
                tang[1] = 0

            self.cell_nn_tx.append(tang[0])
            self.cell_nn_ty.append(tang[1])

        self.cell_nn_tx = np.asarray(self.cell_nn_tx)
        self.cell_nn_ty = np.asarray(self.cell_nn_ty)

        # mapping between gap junction index and cell:
        self.cell_to_nn_full = [[] for x in range(len(self.cell_i))]

        for i, (cell_i, cell_j) in enumerate(self.cell_nn_i):

            if cell_i != cell_j:  # if it's not a boundary membrane...

                self.cell_to_nn_full[cell_i].append(i)
                self.cell_to_nn_full[cell_j].append(i)

        self.cell_to_nn_full = np.asarray(self.cell_to_nn_full)

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

            for key, valu in vars(p).items():
                if type(valu) == interp.interp1d or callable(valu):
                    setattr(p,key,None)

            # save the cell cluster
            logs.log_info('Saving the cell cluster... ')

            datadump = [self,p]
            fh.saveSim(self.savedWorld,datadump)
            message = 'Cell cluster saved to' + ' ' + self.savedWorld
            logs.log_info(message)

    def voronoiGrid(self,p):
        """
        Creates a set of unique, flat points corresponding to cells in the
        cluster in addition to "ghost" points of cell centres present in the
        original Voronoi diagram but removed due to cluster shape.
        """

        # first process voronoi_verts to clip out structures larger than the desired size:
        voronoi_verts = []

        for verts in self.voronoi_verts:
            a = np.asarray(verts)

            inds_highx = (a[:,0] <= self.xmax).nonzero()
            inds_lowx  = (a[:,0] >= self.xmin).nonzero()
            inds_highy = (a[:,1] <= self.ymax).nonzero()
            inds_lowy  = (a[:,1] >= self.ymin).nonzero()

            if (
                len(inds_highx[0]) == len(a[:,0]) and
                len(inds_lowx [0]) == len(a[:,0]) and
                len(inds_highy[0]) == len(a[:,0]) and
                len(inds_lowy[0] ) == len(a[:,0])
            ):
                voronoi_verts.append(verts)

        self.voronoi_verts = np.asarray(voronoi_verts)

        #----------------------------------------------
        voronoi_grid = set()

        for verts in self.voronoi_verts:
            for v in verts:
                voronoi_grid.add((v[0],v[1]))

        voronoi_grid = [list(x) for x in voronoi_grid]
        self.voronoi_grid = np.asarray(voronoi_grid)

        # Create cell centres for the whole voronoi grid:
        self.voronoi_centres = np.array([0,0])
        # self.voronoi_centres = []

        for poly in self.voronoi_verts:
            aa = np.asarray(poly)

            # Center point of this Voronoi region, defined as a 2-element array
            # whose first and second elements are the X and Y coordinates of
            # this center point.
            aa = np.mean(aa,axis=0)

            #FIXME: Inefficient. Consider optimizing by redefining
            #"self.voronoi_centres = []", appending to that list here, and then
            #converting that list to a proper Numpy array below.

            # Append this center point to this array of these points.
            self.voronoi_centres = np.vstack((self.voronoi_centres,aa))

        self.voronoi_centres = np.delete(self.voronoi_centres, 0, 0)

        # define a mapping between the voronoi cell centres and the cluster cell centres:
        vertTree = sps.KDTree(self.voronoi_centres)
        _, self.cell_to_grid = vertTree.query(self.cell_centres)

    def make_maskM(self,p):
        """
        Create structures for plotting interpolated data on cell centres
        and differentiating between the cell cluster and environment.

        """

        voronoiTree = sps.KDTree(self.voronoi_grid)
        _, self.map_voronoi2ecm = voronoiTree.query(self.ecm_verts_unique)

        self.voronoi_mask = np.zeros(len(self.voronoi_grid))
        self.voronoi_mask[self.map_voronoi2ecm]=1

        xv = np.linspace(self.xmin,self.xmax,p.plot_grid_size)
        yv = np.linspace(self.xmin,self.xmax,p.plot_grid_size)

        X,Y = np.meshgrid(xv,yv)

        self.Xgrid = X
        self.Ygrid = Y

        self.maskM = interp.griddata((self.voronoi_grid[:,0],self.voronoi_grid[:,1]),
            self.voronoi_mask,(self.Xgrid,self.Ygrid),
                             method='linear',fill_value=0)

        self.maskM = ndimage.filters.gaussian_filter(self.maskM, 1, mode='nearest')
        self.maskM = np.round(self.maskM,0)

        self.maskECM = interp.griddata((X.ravel(),Y.ravel()),self.maskM.ravel(), (self.X, self.Y), method='linear',fill_value=0)
        self.maskECM = ndimage.filters.gaussian_filter(self.maskECM, 1, mode='nearest')
        self.maskECM = np.round(self.maskECM,0)

        self.inds_env = list(*(self.maskECM.ravel() == 0).nonzero())
        self.inds_clust = list(*(self.maskECM.ravel() == 1).nonzero())

    def quick_maskM(self,p):

        self.maskM = interp.griddata((self.voronoi_grid[:,0],self.voronoi_grid[:,1]),self.voronoi_mask,
                                     (self.Xgrid,self.Ygrid),method='linear',fill_value=0)

        self.maskM = ndimage.filters.gaussian_filter(self.maskM, 2, mode='nearest')
        self.maskM = np.round(self.maskM,0)


        maskECM = interp.griddata((self.Xgrid.ravel(),self.Ygrid.ravel()),self.maskM.ravel(), (self.X, self.Y),
                                  method='linear',fill_value=0)
        maskECM = ndimage.filters.gaussian_filter(maskECM, 2, mode='nearest')
        maskECM = np.round(maskECM,0)

        self.inds_env = list(*(maskECM.ravel() == 0).nonzero())

        self.inds_clust = list(*(self.maskECM.ravel() == 1).nonzero())

    def intra_updater(self,p):
        """
        Calculates a matrix that takes values on membrane midpoints,
        interpolates them to cell vertices, and adds them all together as half
        of a finite volume integration for the pie-box regions. The other half
        will come from the centroid region value.
        """

        self.M_int_mems = np.zeros((len(self.mem_i), len(self.mem_i)))

        for i, inds in enumerate(self.cell_to_mems):

            # get the set of indices for the cell:
            inds = np.asarray(inds)

            inds_p1 = np.roll(inds, 1)
            inds_o = np.roll(inds, 0)
            inds_n1 = np.roll(inds, -1)

            self.M_int_mems[inds_o, inds_o] = (1/3)
            self.M_int_mems[inds_o, inds_p1] = (1/12)
            self.M_int_mems[inds_o, inds_n1] = (1/12)

    def interp_to_mem(self,f, interp_method = 'linear'):
        """
        Interpolates a parameter defined on cell centres to the membrane
        midpoints.

        Parameters
        -----------
        f                A parameter defined on cell centres
        interp_method    Interpolation to use with scipy gridddata ('nearest', 'linear', 'cubic')
        Returns
        -----------
        f_mem
        Interpolation from cell centres to membrane midpoints
        """
        # interpolate f to mems:
        f_mem = interp.griddata((self.cell_centres[:,0],self.cell_centres[:,1]),f,
                                (self.mem_mids_flat[:,0],self.mem_mids_flat[:,1]),
                                fill_value = 0, method=interp_method)

        return f_mem

    def integrator(self, f, fmem) -> tuple:
        """
        Finite volume integrator for the irregular Voronoi cell grid.

        Interpolates a parameter defined on cell centres to membrane midpoints
        and then uses a centre-midpoint interpolation scheme to calculate the
        working 2D integral (volume independent).

        Parameters
        -----------
        f                  A parameter defined on cell centres
        fmem               Same parameter defined on cell membranes

        Returns
        -----------
        fcent          Finite volume interpolation integral over each cell grid (volume independent)
        fmemi          Interpolation of parameter between adjacent membranes
        """

        # average the parameter between adjacent membranes:
        fmemi = (fmem[self.nn_i] + fmem[self.mem_i])/2

        # average the values at the cell centre point:
        fcent = (1/2)*(f + (np.dot(self.M_sum_mems, fmemi)/self.num_mems))

        return fcent, fmemi

    def curl(self, Fx, Fy, phi_z) -> tuple:
        """
        Curl of a vector field defined on cell centres.

        Parameters
        -----------
        Fx       x-component of a 2D vector field
        Fy:      y-component of a 2D vector field
        phi_z:   z- component of a vector field with only a z-component

        Returns
        ----------
        curl_x      Components of curl.
        curl_y
        curl_z
        """

        if type(phi_z) == int:
            # x and y components of the fields:
            grad_F_cell_x = (Fx[self.cell_nn_i[:,1]] - Fx[self.cell_nn_i[:,0]])/(self.nn_len)

            dFx_dy = grad_F_cell_x*self.cell_nn_ty

            grad_F_cell_y = (Fy[self.cell_nn_i[:,1]] - Fy[self.cell_nn_i[:,0]])/(self.nn_len)

            dFy_dx = grad_F_cell_y*self.cell_nn_tx

            curlF_o = dFy_dx - dFx_dy

            curl_z = np.dot(self.M_sum_mems,curlF_o)/self.num_mems

            curl_x = 0
            curl_y = 0

        elif type(Fx) == int and type(Fy) == int:
            # the z-component of the field:
            grad_phi = (phi_z[self.cell_nn_i[:,1]] - phi_z[self.cell_nn_i[:,0]])/(self.nn_len)

            dphi_dx_o = grad_phi*self.cell_nn_tx
            dphi_dy_o = grad_phi*self.cell_nn_ty

            curl_phi_x_o = dphi_dy_o
            curl_phi_y_o = -dphi_dx_o

            curl_x = np.dot(self.M_sum_mems,curl_phi_x_o)/self.num_mems
            curl_y = np.dot(self.M_sum_mems,curl_phi_y_o)/self.num_mems

            curl_z = 0

        else:
            raise BetseSimConfigException("Input to cells.curl not defined properly."
                                           "It takes (Fx = 0, Fy=0, phi)  or "
                                           "(Fx,Fy,phi=0). Also, the 0 must"
                                           "be an integer, not 0.0.")

        return curl_x, curl_y, curl_z

    def zero_div_cell(self, Fn, rho=0.0, bc = 0.0, open_bounds=True):
        """
        Calculates a divergence free field on the cell cluster.

        Parameters
        -----------
        Fn              Force field normal to cell membranes
        rho             Places (such as boundary) with divergence (defined on cell centers)
        bc              Values at the boundary (on bflags_mems)
        cells           An instance of BETSE cells
        open_bounds     whether or not the bounds can support outward flux

        Returns
        -------
        Fn              Divergence corrected normal component of field
        F_cell_x        x-axis component of force field defined on cell centres
        F_cell_y        y-axis component of force field defined on cell centres
        """

        # calculate divergence as the sum of this vector x each surface area, divided by cell volume:
        div_F = (np.dot(self.M_sum_mems, Fn * self.mem_sa) / self.cell_vol)

        fxo = Fn*self.nn_tx
        fyo = Fn*self.nn_ty

        if open_bounds is True:
            Phi = np.dot(self.lapGJinv, div_F + rho)

        else:
            Phi = np.dot(self.lapGJ_P_inv, div_F + rho)


        gPhi = (Phi[self.cell_nn_i[:, 1]] - Phi[self.cell_nn_i[:, 0]]) / (self.nn_len)

        gPx = gPhi*self.nn_tx
        gPy = gPhi*self.nn_ty

        # make the field divergence-free:
        Fx = fxo - gPx
        Fy = fyo - gPy

        Fn = Fx*self.cell_vects_flat[:,2] + Fy*self.cell_vects_flat[:,3]

        # # assign the boundary condition:
        # Fn[self.bflags_mems] = bc

        # Fx = Fn * self.mem_vects_flat[:, 2]
        # Fy = Fn * self.mem_vects_flat[:, 3]


        # calculate the net displacement of cell centres under the applied force under incompressible conditions:
        F_cell_x = np.dot(self.M_sum_mems, Fx) / self.num_mems
        F_cell_y = np.dot(self.M_sum_mems, Fy) / self.num_mems

        return Fn, F_cell_x, F_cell_y

    def HH_cells(self, Fx, Fy, bounds_closed = True, rot_only = False):
        """
        This awesomely magic algorithm calculates a Helmholtz-Hodge
        decomposition to return the divergence-free component of a vector field.

        Thank you, Bhatia et al (2013): "The Helmholtz-Hodge Decompotion: a
        Survey"!

        Parameters
        ------------
        Fx:                # x-coordinate of field on cell centers
        Fy:                # y-coordinate of field on cell centers
        bounds_closed:     # cluster boundary moves flux out, or is sealed
        rot_only           # return only the divergence-free field components
                           # (i.e. rotation only)

        Returns
        ------------
        AA, Ax, Ay     The vector potential AA, as well as divergence-free
                       field components Ax and Ay

        BB, Bx, By     The scalar potential BB, as well as the curl-free
                       field components Bx and By
        """

        nx = self.mem_vects_flat[:, 2]
        ny = self.mem_vects_flat[:, 3]

        _, _, curlF = self.curl(Fx, Fy, 0)

        if bounds_closed is True:

            AA = np.dot(self.lapGJinv, -curlF)

        else:

            AA = np.dot(self.lapGJ_P_inv, -curlF)

        Ax, Ay, _ = self.curl(0, 0, AA)

        # -----------------

        if rot_only is False:

            Fxm = self.interp_to_mem(Fx)
            Fym = self.interp_to_mem(Fy)

            Fn = Fxm * nx + Fym * ny

            divF = np.dot(self.M_sum_mems, Fn * self.mem_sa) / self.cell_vol

            BB = np.dot(self.lapGJ_P_inv, divF)

            gBB = (BB[self.cell_nn_i[:, 1]] - BB[self.cell_nn_i[:, 0]]) / (self.nn_len)

            Bxm = gBB * nx
            Bym = gBB * ny

            Bx = np.dot(self.M_sum_mems, Bxm) / self.num_mems
            By = np.dot(self.M_sum_mems, Bym) / self.num_mems

        else:
            BB = 0
            Bx = 0
            By = 0

        return AA, Ax, Ay, BB, Bx, By

    def deform_tools(self,p):

        # Data structures specific for deformation option------------------------------

        logs.log_info('Creating computational tools for mechanical deformation... ')

        # create a data structure that will allow us to repackage ecm_verts and re-build the
        # cells world after deforming ecm_verts_unique:

        # first get the search-points tree:
        ecmTree = sps.KDTree(self.ecm_verts_unique)

        self.inds2ecmVerts = []

        for verts in self.ecm_verts:

            sublist = []

            for v in verts:
                _, ind = ecmTree.query(v)
                sublist.append(ind)
            self.inds2ecmVerts.append(sublist)

        self.inds2ecmVerts = np.asarray(self.inds2ecmVerts)

    def eosmo_tools(self,p):

        # if studying lateral movement of pumps and channels in membrane,
        # create a matrix that will take a continuous gradient for a value on a cell membrane:
        # returns gradient tangent to cell membrane

        self.gradMem = np.zeros((len(self.mem_i),len(self.mem_i)))

        for i, inds in enumerate(self.cell_to_mems):

            inds = np.asarray(inds)

            inds_p1 = np.roll(inds,1)
            inds_o = np.roll(inds,0)

            dist = self.mem_mids_flat[inds_p1] - self.mem_mids_flat[inds_o]
            len_mem = np.sqrt(dist[:,0]**2 + dist[:,1]**2)

            self.gradMem[inds_o,inds_p1] = 1/len_mem.mean()
            self.gradMem[inds_o,inds_o] = -1/len_mem.mean()

    #..........{ PROPERTIES ~ membrane                   }.....................
    @property_cached
    def membranes_normal_unit_x(self) -> ndarray:
        '''
        One-dimensional Numpy array indexing each cell membrane such that each
        element is the X coordinate of the normal unit vector orthogonal to the
        tangent unit vector for this membrane.
        '''

        return self.mem_vects_flat[:, 2]


    @property_cached
    def membranes_normal_unit_y(self) -> ndarray:
        '''
        One-dimensional Numpy array indexing each cell membrane such that each
        element is the Y coordinate of the normal unit vector orthogonal to the
        tangent unit vector for this membrane.
        '''

        return self.mem_vects_flat[:, 3]

    #..........{ PROPERTIES ~ mappers                    }.....................
    #FIXME: For readability, rename to membranes_midpoint_to_vertices().
    @property_cached
    def matrixMap2Verts(self) -> ndarray:
        '''
        Numpy matrix (i.e., two-dimensional array) of size ``m x n``, where:

        * ``m`` is the total number of cell membranes.
        * ``n`` is the total number of cell membrane vertices.

        For each membrane ``i`` and membrane vertex ``j``, element
        ``matrixMap2Verts[i, j]`` is:

        * 0 if this vertex is *not* one of the two vertices defining this
          membrane. Since most vertices do *not* define most membranes, most
          entries of this matrix are zero, implying this matrix to typically
          (but *not* necessarily) be sparse.
        * 0.5 if this vertex is one of the two vertices defining this membrane,
          thus averaging membrane data defined at membrane midpoints over the
          vertex pairs defining these membranes.

        Usage
        -----------
        The dot product of a Numpy vector (i.e., one-dimensional array) of size
        ``m`` containing membrane-specific data by this matrix yields another
        Numpy vector of size ``n`` containing membrane vertex-specific data
        interpolated from these membranes over these vertices, where ``m`` and
        ``n`` are as defined above. Technically, the dot product of a vector by a
        matrix is undefined. To facilitate this otherwise invalid operation,
        Numpy implicitly converts:

        * The input vector of size ``m`` into a matrix of size ``1 x m``.
        * The output matrix of size ``n x 1`` into an output vector of size
          ``n``.

        This matrix is cached *only* on the first access of this property.
        '''

        # Zero this matrix to the expected dimensions.
        matrixMap2Verts = np.zeros(
            (len(self.mem_mids_flat), len(self.mem_verts)))

        # For the indices of each membrane and the two vertices terminating
        # that membrane, interpolate arbitrary data defined at the former over
        # the latter.
        for cell_membrane_index, cell_membrane_vertices_index in enumerate(
            self.index_to_mem_verts):
            matrixMap2Verts[
                cell_membrane_index, cell_membrane_vertices_index[0]] = 1/2
            matrixMap2Verts[
                cell_membrane_index, cell_membrane_vertices_index[1]] = 1/2

        # Cache this matrix.
        return matrixMap2Verts


    #FIXME: Eventually we want to switch this up. This data structure should
    #replace "self.M_sum_mems" everywhere; after doing so, "self.M_sum_mems"
    #should be removed.
    @property_cached
    def membranes_midpoint_to_cells_centre(self) -> ndarray:
        '''
        Numpy matrix (i.e., two-dimensional array) of size ``m x n``, where:

        * ``m`` is the total number of cell membranes.
        * ``n`` is the total number of cells.

        For each membrane ``i`` and cell ``j``, element
        ``mems_midpoint_to_cells_centre[i, j]`` is:

        * 0 if this cell does *not* contain this membrane. Since most cells do
          *not* contain most membranes, most entries of this matrix are zero,
          implying this matrix to typically (but *not* necessarily) be sparse.
        * ``1/k`` if this cell contains this membrane, where ``k`` is the number
          of membranes this cell contains.

        Usage
        -----------
        The dot product of a Numpy vector (i.e., one-dimensional array) of size
        ``m`` containing arbitrary data spatially situated at cell membrane
        midpoints by this matrix by yields another Numpy vector of size ``n``
        containing the same data resituated at cell centres, where ``m`` and
        ``n`` are as defined above.

        This matrix is cached *only* on the first access of this property.
        '''

        # Dismantled, this is:
        #
        # * "self.M_sum_mems.T", the transpose of the "M_sum_mems" matrix. Since
        #   this matrix is of size m x n for m the number of cells and n is the
        #   number of membranes, this transpose is a matrix of size n x m. Each
        #   element of this transpose is either:
        #   * 0 if this cell does *NOT* contain this membrane.
        #   * 1 if this cell contains this membrane.
        # * "... / self.num_mems", normalizing each cell membrane element of
        #   this matrix by the number of membranes in that cell. Since
        #   "num_mems" is a row vector of length m whose elements are the number
        #   of membranes in that cell, each column of this transpose is divided
        #   by the corresponding element of this row vector.
        return self.M_sum_mems.T / self.num_mems

    #..........{ MAPPERS                                 }.....................
    @type_check
    def map_membranes_midpoint_to_cells_centre(
        self, membranes_midpoint_data: SequenceTypes) -> ndarray:
        """
        Interpolate a one- to two-dimensional sequence of arbitrary data
        spatially situated at cell membrane midpoints into a Numpy array or
        matrix of the same data resituated at cell centres.

        Each element of the output array is the average of the passed membrane
        data over all membranes of the corresponding cell, thus interpolating
        from fine-grained data defined at membrane midpoints to coarse-grained
        data defined at cell centers.

        Parameters
        -----------
        membranes_midpoint_data : SequenceTypes
            Either:
            * One-dimensional sequence indexing each cell membrane, such that
              each element is arbitrary data spatially situated at the midpoint
              of that cell membrane.
            * Two-dimensional sequence whose:
              * First dimension is of arbitrary length, typically indexing each
                simulation time step.
              * Second dimension indexes each cell membrane, such that each
                element is arbitrary data spatially situated at the midpoint of
                that cell membrane.

        Returns
        -----------
        ndarray
            If the passed array is:
            * One-dimensional, the returned array is also one-dimensional of
              length the number of cells, indexing the same data as that of the
              passed array resituated at cell centres.
            * Two-dimensional, the returned array is also two-dimensional
              whose:
              * First dimension is the same length as that of the first
                dimension of the passed array.
              * Second dimension has length the number of cells, indexing the
                same data as that of the passed array resituated at cell
                centres.
        """

        # Numpy array converted from the passed sequence.
        membranes_midpoint_data = arrays.from_sequence(membranes_midpoint_data)

        # If the last dimension of this array does *NOT* index all cell
        # membranes and hence is *NOT* spatially situated at cell membrane
        # midpoints, raise an exception.
        if membranes_midpoint_data.shape[-1] != self.mem_mids_flat.shape[0]:
            raise BetseSequenceException(
                'Input array not spatially situated at cell membrane midpoints '
                '(i.e., last array dimension length {} not {}).'.format(
                    membranes_midpoint_data.shape[-1],
                    self.mem_mids_flat.shape[0]))

        # If this array is three- or more-dimensional, raise an exception.
        if len(membranes_midpoint_data.shape) >= 3:
            raise BetseSequenceException(
                'Input array of dimensionality {} neither one- nor '
                'two-dimensional.'.format(len(membranes_midpoint_data.shape)))

        # Map this array from cell membrane midpoints onto cell centres. By
        # design, this efficiently supports both one- and two-dimensional input
        # arrays as is.
        return np.dot(
            membranes_midpoint_data, self.membranes_midpoint_to_cells_centre)


    def map_cells_centre_to_membranes_midpoint(
        self,
        cells_centre_data: SequenceTypes,
        interp_method: str = 'linear',
    ) -> ndarray:
        """
        Interpolate a sequence of arbitrary data defined on cell centres
        into a Numpy array of the same data defined on cell membrane midpoints.

        Parameters
        -----------
        cells_centre_data : SequenceTypes
            One-dimensional sequence of length the number of cells such that
            each element is arbitrary cell data spatially situated at the
            center of the cell with that 0-based index.
        interp_method : str
            Interpolation type to pass to the
            :func:`scipy.interpolate.gridddata` function (e.g., `nearest`,
            `linear`, `cubic`).

        Returns
        -----------
        ndarray
            One-dimensional Numpy array of length the number of cell membranes
            such that each element is arbitrary cell membrane data spatially
            situated at the midpoint of the membrane with that 0-based index.
        """

        # Numpy array converted from the passed sequence.
        cells_centre_data = arrays.from_sequence(cells_centre_data)

        # If the last dimension of this array does *NOT* index all cells and
        # hence is *NOT* spatially situated at cell centres, raise an exception.
        if cells_centre_data.shape[-1] != self.cell_centres.shape[0]:
            raise BetseSequenceException(
                'Input array not spatially situated at cell centres '
                '(i.e., last array dimension length {} not {}).'.format(
                    cells_centre_data.shape[-1], self.cell_centres.shape[0]))

        # If this array is three- or more-dimensional, raise an exception.
        if len(cells_centre_data.shape) >= 3:
            raise BetseSequenceException(
                'Input array of dimensionality {} neither one- nor '
                'two-dimensional.'.format(len(cells_centre_data.shape)))

        # 2-tuple of the X and Y coordinates of all cell centres.
        cell_centres_points = (self.cell_centres[:,0], self.cell_centres[:,1])

        # 2-tuple of the X and Y coordinates of all cell membrane midpoints.
        cell_membrane_midpoints = (
            self.mem_mids_flat[:,0], self.mem_mids_flat[:,1])

        # If this array is one-dimensional, map this array from cell centres
        # onto cell membrane midpoints with a single interpolation call.
        if len(cells_centre_data.shape) == 1:
            return interp.griddata(
                cell_centres_points,
                cells_centre_data,
                cell_membrane_midpoints,
                fill_value=0,
                method=interp_method,
            )
        # Else, this array is two-dimensional. Since SciPy requires that the
        # second parameter to the interp.griddata() function be one-dimensional,
        # map this array from cell centres onto cell membrane midpoints with
        # multiple interpolation calls. While inefficient, SciPy appears to
        # provide no alternative optimization.
        else:
            # Two-dimensional list aggregating all interpolation results.
            cells_centre_data_interpolated = []

            # For each one-dimensional array in this two-dimensional array...
            for cells_centre_data_one in cells_centre_data:
                # One-dimensional array interpolated from this array.
                cells_centre_data_one_interpolated = interp.griddata(
                    cell_centres_points,
                    cells_centre_data_one,
                    cell_membrane_midpoints,
                    fill_value=0,
                    method=interp_method,
                )

                # Append this array to this list.
                cells_centre_data_interpolated.append(
                    cells_centre_data_one_interpolated)

            # Return this list converted back to an array.
            return arrays.from_sequence(cells_centre_data_interpolated)


    #FIXME: To reduce code duplication:
    #
    #    # Globally replace all instances of this...
    #    np.dot(cells.M_sum_mems, some_array) / cells.num_mems
    #
    #    # ...with this.
    #    cells.calculate_divergence(some_array)

        # ...which should internally do something like this:
    #    np.dot(some_array, cells.M_sum_mems.T) / cells.num_mems
    #
    #After doing so, consider replacing the "M_sum_mems" array with its
    #tranpspose, which is a # the more general-purpose and efficient array.
    #Unlike the former, the latter is applicable to both Numpy arrays and
    #matrices.
