#!/usr/bin/env python3
# ....................{ LICENSE                           }....................
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS                           }....................
import math
import numpy as np
from numpy import ndarray
from scipy import interpolate as interp
from scipy import ndimage
from scipy.spatial import Voronoi, cKDTree
from betse.exceptions import BetseSequenceException, BetseSimConfException
from betse.science import filehandling as fh
from betse.science.enum.enumconf import CellLatticeType
from betse.science.math import finitediff as fd
from betse.science.math import toolbox as tb
from betse.util.math.geometry.polygon.geopolyconvex import clip_counterclockwise
from betse.util.math.geometry.polygon.geopoly import orient_counterclockwise, is_convex
from betse.science.phase.phasecls import SimPhase
from betse.util.io.log import logs
from betse.util.type.decorator.decmemo import property_cached
from xml.dom import minidom
from betse.util.type.types import (
    type_check, NumericOrSequenceTypes, SequenceTypes)
from betse.lib.numpy import nparray
from betse.util.type.text import regexes
from betse.science.tissue.picker.tispickimage import TissuePickerImage
from betse.science.math.mesh import DECMesh

# ....................{ CLASSES                           }....................
#FIXME: Create a new option for seed points: Fibonacci radial-spiral array.
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
    cellVerts()                       Copy & scale in points from the ecm matrix to create unique polygonal cells
    quickVerts()                      Reformulates an exisitng cell world
    cellMatrices()                    Calculates for matrices for BETSE
    cell_vols()                       Calculates cell patch volumes
    near_neigh()                      Calculate the nearest neighbour (nn) array for each cell (make gap junctions)
    voronoiGrid()
    makeECM()                         Make the Marker and Cell (MACs) grid for extracellular calculations
    environment()                     Calculate details for the extracellular calculations, including mappings
    graphLaplacian()                  Creates an abstract discrete Laplacian for the irregular Voronoi-based cell grid

    Attributes (Cell)
    ----------
    cell_centres : ndarray
        Two-dimensional Numpy array of the coordinates of the center points of
        all cells, whose:

        #. First dimension indexes cells, whose length is the number of cells.
        #. Second dimension indexes the coordinates of the center point of the
           current cell, whose length is unconditionally guaranteed to be 2
           *and* whose:

           #. First item is the X coordinate of the current cell center.
           #. Second item is the Y coordinate of the current cell center.
    cell_i : list
        One-dimensional list indexing each cell such that each item is that
        cell's index (i.e., ``[0, 1, ..., n-2, n-1]`` for the number of cells
        ``n``), required for efficient Numpy slicing.
    cell_verts : ndarray
        Three-dimensional Numpy array of the coordinates of the vertices of all
        cells, whose:

        #. First dimension indexes each cell such that each item is a
           matplotlib-compatible **polygon patch** (i.e., a two-dimensional
           Numpy array of all vertex coordinates defining the current cell's
           polygon), suitable for passing as is to the
           :meth:`matplotlib.patches.Polygon.__init__` method.
        #. Second dimension indexes each vertex of the current cell (in
           counterclockwise order).
        #. Third dimension indexes each coordinate of the current cell vertex,
           whose length is guaranteed to be 2 *and* whose:

           #. First item is the X coordinate of the current cell vertex.
           #. Second item is the Y coordinate of the current cell vertex.
    R : ndarray
        One-dimensional Numpy array indexing each cell such that each item is
        the radius of that cell.

    Attributes (Cell Membrane)
    ----------
    cell_to_mems : ndarray
        Two-dimensional Numpy array of the indices of all membranes for all
        cells, whose:

        #. First dimension indexes each cell.
        #. Second dimension indexes each membrane of the current cell such that
           each item is the index of the corresponding membrane in
           membrane-centric arrays (e.g., :attr:`mem_mids_flat`).

        For example:

        * ``cell_to_mems[0][0]`` is the index of the first cell's first
          membrane (guaranteed to be 0).
        * ``cell_to_mems[-1][-1]`` is the index of the last cell's last
          membrane (guaranteed to be ``mem_i[-1]``).
    mem_to_cells : ndarray
        One-dimensional Numpy array of length the number of cell membranes such
        that each item is the index of the cell containing the membrane
        indexed by that item. For example:

        * ``mem_to_cells[0]`` is the index of the cell containing the first
          membrane.
        * ``mem_to_cells[-1]`` is the index of the cell containing the last
          membrane.
    mem_i : list
        One-dimensional list of length the number of cell membranes such that
        each item is that cell membrane's index (i.e.,
        ``[0, 1, ..., m-2, m-1]`` for the number of cell membranes ``m``),
        required for efficient Numpy slicing.
    mem_mids_flat : ndarray
        Two-dimensional Numpy array of the coordinates of the midpoints of all
        cell membranes, whose:

        #. First dimension indexes each cell membrane.
        #. Second dimension indexes each coordinate of the midpoint of the
           current cell membrane, whose length is guaranteed to be 2 and whose:

           #. First item is the X coordinate of the current membrane midpoint.
           #. Second item is the Y coordinate of the current membrane
              midpoint.
    mem_vects_flat : ndarray
        Two-dimensional Numpy array of the coordinates of various vectors of
        all cell membranes, whose:

        #. First dimension indexes each cell membrane.
        #. Second dimension indexes each coordinate of a vector describing the
           current cell membrane, whose length is guaranteed to be 6 and whose:

           #. First item is the X coordinate of the current membrane midpoint,
              equivalent to :attr:`mem_mids_flat[mem_i,0]` where ``mem_i`` is
              the index of the current membrane.
           #. Second item is the Y coordinate of the current membrane
              midpoint, equivalent to :attr:`mem_mids_flat[mem_i,1]`.
           #. Third item is the X coordinate of the normal unit vector
              orthogonal to the tangent unit vector of the current membrane
              (defined by the corresponding item of the fifth and sixth
              elements of this array).
           #. Fourth item is the Y coordinate of this normal unit vector.
           #. Fifth item is the X coordinate of the tangent unit vector
              parallel to the vector implied by the pair of coordinates
              defining the current membrane.
           #. Sixth item is the Y coordinate of this tangent unit vector.
    num_mems : ndarray
        One-dimensional Numpy array indexing each cell such that each item
        is the number of cell membranes contained by the current cell.
    M_sum_mems : ndarray
        Numpy matrix (i.e., two-dimensional array) of size ``m x n``, where:

        * `m` is the total number of cells.
        * `n` is the total number of cell membranes.

        For each cell ``i`` and membrane ``j``, item ``M_sum_mems[i, j]`` is:

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

        #. First dimension indexes each cell membrane.
        #. Second dimension indexes the index of each vertex defining the
           current membrane, whose length is unconditionally guaranteed to be 2
           *and* whose:

           #. First item is the index of the first vertex defining this
              membrane in the :attr:`mem_verts` subarray, guaranteed to be
              counterclockwise from the second vertex defining this membrane.
           #. Second item is the index of the second vertex defining this
              membrane in the :attr:`mem_verts` subarray, guaranteed to be
              clockwise from the first vertex defining this membrane.
    mem_verts : ndarray
        Two-dimensional Numpy array of the coordinates of the vertices of all
        cell membranes, whose:

        #. First dimension indexes cell membrane vertices, whose length is the
           number of cell membrane vertices in this cluster. This length is
           strictly greater than the number of cell membranes in this cluster.
           Since this array contains *no* duplicate vertices, this length is
           strictly less than twice the number of cell membranes.
        #. Second dimension indexes the coordinates of the current membrane
           vertex, whose length is unconditionally guaranteed to be 2 *and*
           whose:

           #. First item is the X coordinate of the current membrane vertex.
           #. Second item is the Y coordinate of the current membrane vertex.

    Attributes (Cell Gap Junctions)
    ----------
    cell_nn_i : ndarray
        Two-dimensional Numpy array of the indices of the pairs of adjacent
        cells comprising all gap junctions connecting one cell to another in
        this cluster, whose:

        #. First dimension indexes cell membranes, whose length is the number
           of cell membranes in this cluster.
        #. Second dimension indexes the indices of the pair of cells defining
           the gap junction connecting the current membrane to another
           membrane, whose length is unconditionally guaranteed to be 2 *and*
           whose:

           #. First item is the index of the cell containing the current
              membrane.
           #. Second item is the index of an adjacent cell containing the
              adjacent membrane to which the current membrane is connected.
    gj_len : float
        Uniform length (in meters) of each gap junction, equivalent to the
        uniform distance between each pair of neighbouring cells.
    mem_nn : ndarray
        Two-dimensional Numpy array of the indices of the pairs of adjacent
        cell membranes comprising all gap junctions connecting one membrane to
        another in this cluster, whose:

        #. First dimension indexes cell membranes, whose length is the number
           of cell membranes in this cluster.
        #. Second dimension indexes the indices of the pair of cell membranes
           defining the gap junction connecting the current membrane to another
           membrane, whose length is unconditionally guaranteed to be 2 *and*
           whose:

           #. First item is the index of the current membrane.
           #. Second item is:

              * If the current membrane is adjacent to no other membrane (e.g.,
                due to being situated at either the periphery of this cluster
                *or* a discontiguous hole in this cluster), the index of the
                current membrane again. A solitary membrane thus connects to
                itself with a self-referential gap junction "loop."
              * If the current membrane is adjacent to exactly one other
                membrane (which is the common case), the index of that membrane.
              * If the current membrane is adjacent to two or more other
                membranes, the index of the membrane to which the current
                membrane is most adjacent (excluding itself).

        Hence, *all* membranes are guaranteed to participate in exactly one
        gap junction connection.
    nn_i : ndarray
        One-dimensional Numpy array of length the number of cell membranes such
        that each item is the index of the nearest neighbouring cell
        membrane of the cell membrane indexed by that item, comprising the
        gap junction connecting these membranes. Each item thus maps each
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
        One-dimensional Numpy array of length the number of cells such that
        each item is the surface area of the cell indexed by that item, defined
        as the summation of the surface areas of all membranes comprising that
        cell.
    mem_sa : ndarray
        One-dimensional Numpy array of length the number of cell membranes such
        that each item is the surface area of the cell membrane indexed by
        that item.

    Attributes (Extracellular Grid)
    ----------
    These attributes are defined *only* if the :meth:`makeECM` method has been
    called, implying the current simulation to enable extracellular spaces. If
    that method has *not* yet been called, these attributes remain undefined.

    delta : NumericSimpleTypes
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
    points_tree : scipy.spatial.cKDTree
        Kd-tree on the :attr:`xypts` array, enabling efficient mapping from
        arbitrary Cartesian coordinates to their nearest extracellular grid
        spaces.

    Attributes (Voronoi Diagram)
    ----------
    cell_to_grid : ndarray
        One-dimensional Numpy array indexing each cell such that each item
        is the index of the Voronoi region whose vertices most closely
        spatially align with those of that cell. For example:

        * ``cell_to_grid[0]`` is the index of the region most closely spatially
          aligned with the first cell.
        * ``cell_to_grid[-1]`` is the index of the region most closely
          spatially aligned with the last cell.

        Assigning a Numpy array of length the number of regions indexed by this
        array from a Numpy array of length the number of cells maps the
        cell-specific data defined by the latter into region-specific data
        (e.g., ``regions_data[cells.cell_to_grid] = cells_centre_data``).
    voronoi_centres : ndarray
        Two-dimensional Numpy array of the coordinates of the center points of
        all polygonal regions in the Voronoi diagram producing this cell
        cluster, whose:

        #. First dimension indexes regions in this diagram, whose length is the
           total number of regions in this diagram.
        #. Second dimension indexes the coordinates of the center point of the
           current region, whose length is unconditionally guaranteed to be 2
           *and* whose:

           #. First item is the X coordinate of the current region center.
           #. Second item is the Y coordinate of the current region center.
    voronoi_grid : ndarray
        Two-dimensional Numpy array of the vertex coordinates of all
        polygonal regions in the Voronoi diagram producing this cell cluster,
        whose:

        #. First dimension indexes the vertices of all regions in this diagram,
           whose length is the total number of vertices in this entire diagram.
        #. Second dimension indexes the coordinates of the current vertex, whose
           length is unconditionally guaranteed to be 2 *and* whose:

           #. First item is the X coordinate of the current region vertex.
           #. Second item is the Y coordinate of the current region vertex.

        This array is equivalent to the :attr:`voronoi_verts` array flattened
        over the first dimension of that array.
    voronoi_verts : ndarray
        Three-dimensional Numpy array of the vertex coordinates of all
        polygonal regions in the Voronoi diagram producing this cell cluster,
        whose dimensions are structured as those of the :attr:`cell_verts`
        attribute (replacing "cell" with "Voronoi region").
    '''

    # ..................{ INITIALIZORS                      }..................
    #FIXME: Refactor this method as follows:
    #
    #* Remove the optional "p" parameter, which is no longer (and should
    #  ideally *NEVER*) be required here.
    #* Refactor all calls to this method to avoid passing this parameter.

    # Avoid circular import dependencies.
    @type_check
    def __init__(self, p: 'betse.science.parameters.Parameters') -> None:
        '''
        Initialize this cell cluster with the passed simulation configuration.

        Parameters
        ----------
        p : betse.science.parameters.Parameters
            Current simulation configuration.
        '''

        pass

    # ..................{ MAKERS                            }..................
    MAKE_WORLD_PROGRESS_TOTAL = 5
    '''
    Cumulative number of times that each call of the :meth:`make_world` method
    calls either the :meth:`SimCallbacksBC.progressed` callback or
    higher-level callbacks calling that callback (e.g.,
    :meth:`SimCallbacksBC.progressed_next`).

    This magic number *must* be manually synchronized with the
    implementation of both the :meth:`make_world` method and methods
    transitively called by that method. Failure to do so *will* result in fatal
    exceptions. There exists no reasonable means of enforcing this constraint.
    '''


    @type_check
    def make_world(self, phase: SimPhase) -> None:
        '''
        Construct the entire world (including both the cell cluster and
        extracellular environment) to be simulated for the passed seed
        simulation phase.

        Parameters
        --------
        phase : SimPhase
            Current simulation phase.
        '''

        # Notify the sink callback of the current state of progress.
        phase.callbacks.progress_stated(
            'Creating pseudo-random cell cluster...')

        # If this is *NOT* the seed phase, raise an exception.
        phase.die_unless_kind_seed()

        # Create the cell lattice, which serves as the seed grid underlying all
        # subsequent data structures.
        if phase.p.svg_override is True: # If user requests cell centres and clip from svg files
            self.parse_svg(phase.p)

        else: # Otherwise make cell centre points algorithmically
            self._make_cell_lattice(phase.p)


        # Define the initial seed point collection:
        # seed_points_o = np.vstack((self.clust_xy, self.bbox))
        # seed_points_o = np.vstack((self.clust_xy))

        # Create the initial Voronoi diagram after making the initial seed
        # points of the cell lattice.
        self.make_voronoi(phase)

        # Clean the Voronoi diagram of empty data structures.
        # self._clean_voronoi(phase.p)

        # Calculate the centroids of all Voronoi cells.
        self.cell_index(phase.p)

        # If optionally refining the Voronoi mesh, do so.
        # if phase.p.refine_mesh:
        #     self._refine_voronoi(phase)

        # Notify the sink callback of the current state of progress.
        phase.callbacks.progressed_next(
            status='Calculating geometric cell properties...')

        # Create individual cell polygon vertices.
        self.cellVerts(phase.p)

        # Create a variety of matrices used in routine cells calculations.
        self.cellMatrices(phase.p)

        # Notify the sink callback of the current state of progress.
        phase.callbacks.progressed_next(
            status='Calculating geometric cell volumes...')

        # Calculate the volume of each cell and its internal regions.
        self.cell_vols(phase.p)

        # Calculate membrane nearest neighbours, ECM interaction, boundary
        # tags, and so forth.
        self.mem_processing(phase.p)

        # Calculate the nearest neighbour array for each cell.
        self.near_neigh(phase.p)
        # self.voronoiGrid(phase.p)

        # Notify the sink callback of the current state of progress.
        phase.callbacks.progressed_next(
            status='Creating extracellular matrix (ECM) grid...')

        # Create the ECM grid.
        self.makeECM(phase.p)

        # Fefine features of the ECM grid.
        self.environment(phase.p)
        self.grid_len = len(self.xypts)
        self.make_maskM(phase.p)
        # self.env_weighting(p)

        # Notify the sink callback of the current state of progress.
        phase.callbacks.progressed_next(
            status='Creating environmental voltage Poisson solver...')

        # Avoid storing the non-inverse matrix, which only consumes memory.
        bdic = {'N': 'value', 'S': 'value', 'E': 'value', 'W': 'value'}
        _, self.lapENVinv = self.grid_obj.makeLaplacian(bound=bdic)

        # Notify the sink callback of the current state of progress.
        phase.callbacks.progressed_next(
            status='Creating environmental current Poisson solver...')

        # Avoid storing the non-inverse matrix, which only consumes memory.
        bdic = {'N': 'flux', 'S': 'flux', 'E': 'flux', 'W': 'flux'}
        _, self.lapENV_P_inv = self.grid_obj.makeLaplacian(bound=bdic)

        # Set all Laplacian matrices to "None" to allow for flexible creation
        # of Laplacians and inverses on the cell grid (i.e., two boundary
        # conditions sets).
        self.lapGJinv = None
        self.lapGJ = None

        # Other matrices.
        self.M_sum_mem_to_ecm = None   # used for deformation
        self.gradMem = None  # used for electroosmosis

        # Coerce this to an array for advanced indexing and assignments.
        self.cell_i = np.asarray(self.cell_i)
        self.gj_default_weights = np.ones(len(self.mem_i))

    # ..................{ DEFORMERS                         }..................
    def deformWorld(self, p, ecm_verts) -> None:
        '''
        Recalculates the current world to accomodate a mechanical deformation.

        Note that each cell is assumed to be incompressible by default.
        Therefore, cell volumes and total surface area are not updated in a
        deformation routine.
        '''

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

        # Finish up by creating indices vectors and converting to Numpy arrays
        # where needed:
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

    @type_check
    def _make_cell_lattice(
        # Avoid circular import dependencies.
        self, p: 'betse.science.parameters.Parameters') -> None:
        '''
        Create the randomly scattered lattice of **cell seed points** (i.e.,
        cell centres) defined on a 2D world space, with dimensions supplied by
        the :attr:`wsx` variable in meters.

        Specifically, this method defines the following object attributes:

        * ``xmin``, ``xmax``, ``ymin``, ``ymax``, the minimum and maximum points
          of the global world space.
        * ``centre``, the centre of the global world space.
        * ``clust_xy``, the list of `(x, y)` Cartesian coordinates of seed
          points.

        The amount of deviation from a square grid is specified by the
        :attr:`cell_lattice_disorder` variable, defined from 0 (perfect square
        grid) to 1 (full noise).

        Parameters
        -----------
        p : Parameters
            Current simulation configuration.
        '''

        # If a standard square cell lattice is requested, create this lattice.
        if p.cell_lattice_type is CellLatticeType.SQUARE:
            # Log this creation.
            logs.log_info('Creating square cell lattice...')

            # Linear vectors which are the "ticks" of the X and Y dimensions.
            x_v = np.linspace(0, (p.nx - 1) * (p.d_cell + p.ac), p.nx)  # create lattice vector x
            y_v = np.linspace(0, (p.ny - 1) * (p.d_cell + p.ac), p.ny)  # create lattice vector y

            # next define a 2d array of lattice points using the x- and y- vectors
            x_2d, y_2d = np.meshgrid(x_v, y_v)  # create 2D array of lattice points

            # now create a matrix of points that will add a +/- deviation to each point centre
            x_rnd = p.cell_lattice_disorder * p.d_cell * (np.random.rand(p.ny, p.nx) - 0.5)  # create a mix of random deltas x dir
            y_rnd = p.cell_lattice_disorder * p.d_cell * (np.random.rand(p.ny, p.nx) - 0.5)  # create a mix of random deltas x dir

            # Add the noise effect to the world point matrices.
            x_2d = x_2d + x_rnd
            y_2d = y_2d + y_rnd

            # Redefine the results.
            xypts = np.vstack((x_2d.ravel(), y_2d.ravel())).T

        # If a hexagonal square cell lattice is requested, create this lattice.
        elif p.cell_lattice_type is CellLatticeType.HEX:
            # Log this creation.
            logs.log_info('Creating hexagonal cell lattice...')

            # Recalculate the number of lattice sites.
            # r = (p.d_cell/2)*np.sqrt(3)
            r = p.d_cell/2

            # "Blue" and "red" simply refer to the different seed point arrays
            # used to create the hexagonal lattice.
            a_hex = np.sqrt(3) * r  # distance flat radius (y-axis hex radius)
            b_hex = 3 * r  # x-spacing between midpoint of blue and red grids

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
            x_rnd = p.cell_lattice_disorder * p.d_cell * (np.random.rand(len(x_2d)) - 0.5)  # create a mix of random deltas x dir
            y_rnd = p.cell_lattice_disorder * p.d_cell * (np.random.rand(len(y_2d)) - 0.5)  # create a mix of random deltas x dir

            # # add the noise effect to the world point matrices and redefine the results
            x_2d = x_2d + x_rnd
            y_2d = y_2d + y_rnd

            xypts = np.vstack((x_2d, y_2d)).T

        # Else, this lattice type is unrecognized. Raise us up an exception!
        else:
            raise BetseSimConfException(
                'Cell lattice type "{}" unrecognized.'.format(
                    p.cell_lattice_type))

        # Define a data structure that holds [x,y] coordinate points of each 2d
        # grid-matrix entry.

        # define geometric limits and centre for the cluster of points
        self.xmin = np.min(xypts[:,0])
        self.xmax = np.max(xypts[:,0])
        self.ymin = np.min(xypts[:,1])
        self.ymax = np.max(xypts[:,1])

        bbox = np.asarray(
            [[self.xmin, self.ymin],
             [self.xmax, self.ymin],
             [self.xmax, self.ymax],
             [self.xmin, self.ymax]])

        self.centre = xypts.mean(axis=0)
        # self.clust_xy = xypts

        self.clust_xy = np.vstack((xypts, bbox))

        # alter the actual bounding box to create a nice clipping and plotting polygon
        self.bbox = np.asarray(
            [[self.xmin - p.cell_radius, self.ymin - p.cell_radius],
             [self.xmax + p.cell_radius, self.ymin - p.cell_radius],
             [self.xmax + p.cell_radius, self.ymax + p.cell_radius],
             [self.xmin - p.cell_radius, self.ymax + p.cell_radius]])

        self.xmin = self.xmin - p.cell_radius
        self.ymin = self.ymin - p.cell_radius
        self.xmax = self.xmax + p.cell_radius
        self.ymax = self.ymax + p.cell_radius

        # define axis bounds for easy plotting:
        self.xyaxis = [self.xmin, self.xmax, self.ymin, self.ymax]


    @type_check
    def make_voronoi(self, phase: SimPhase) -> None:
        '''
        Calculate, close, and clip the Voronoi diagram from the cell lattice
        previously created by the :func:`_make_cell_lattice` method.

        Specifically, this method (in order):

        #. Calculates this diagram.
        #. Closes this diagram at the square boundaries of the global
           environment.
        #. Removes all cells from this diagram to define a cluster shape.

        Parameters
        --------
        phase : SimPhase
            Current simulation phase.

        Creates
        ---------
        self.ecm_verts              x,y points of Voronoi cell vertices (nulled at world creation endpoint)
        self.ecm_verts_unique       x,y points of unique Voronoi cell vertices (nulled at world creation endpoint)
        self.ecm_polyinds           indices into self.ecm_verts_unique list defining each Voronoi polygon
                                    (nulled at world creation endpoint)
        self.cluster_mask           Matrix of booleans defining masked shape of cell cluster
        self.msize                  Size of bitmap (side pixel number)
        '''

        # Log this action.
        logs.log_info('Creating Voronoi geometry... ')

        # Cell cluster image picker, producing the cell cluster image mask.
        from betse.science.tissue.picker.tispickimage import TissuePickerImage
        image_picker = TissuePickerImage(
            filename=phase.p.tissue_default.picker_image_filename,
            dirname=phase.p.conf_dirname)

        # Cell cluster image mask, clipping the cell cluster against a
        # user-defined image file.
        image_mask = image_picker.get_image_mask(cells=self)

        # Dictionary mapping from cell lattice enumeration members to
        # corresponding low-level DEC-specific strings.
        lattice_to_mesh_type = {
            CellLatticeType.HEX:    'tri',
            CellLatticeType.SQUARE: 'rect',
        }

        if phase.p.single_cell: # if simulating only a single cell:

            self.mesh = DECMesh(cell_radius=phase.p.cell_radius,
                                seed_points=None,
                                use_alpha_shape=False,
                                allow_merging=False,
                                single_cell_noise=0.1,
                                single_cell_sides=6,
                                image_mask=None,
                                make_all_operators=False, # FIXME: later fix this; for now set to False
                                mesh_type=lattice_to_mesh_type[phase.p.cell_lattice_type],
                                center = self.centre) # Optional center point for cluster (used to center a single cell)

            self.mesh.init_mesh()


        else: # otherwise, create a whole mesh of cells:

            if phase.p.svg_override: # If using an svg file, make the mesh directly from the points:

                self.mesh = DECMesh(cell_radius = None,
                         seed_points=self.clust_xy,
                         alpha_shape = phase.p.alpha_shape,
                         use_centroids = phase.p.use_centroids,
                         use_alpha_shape = True,
                         single_cell_noise = 0.5,
                         single_cell_sides = 6,
                         image_mask = None,
                         allow_merging = True,
                         merge_thresh = 0.2,
                         make_all_operators=False,  # FIXME: later fix this; for now set to False
                         mesh_type = lattice_to_mesh_type[phase.p.cell_lattice_type])

                self.mesh.init_mesh()

            else: # Otherwise, make and clip a point cloud to the base shape:
                # Initialize the Discrete Exterior Calculus (DEC) mesh object:
                self.mesh = DECMesh(
                    cell_radius=phase.p.cell_radius,
                    seed_points=self.clust_xy,
                    alpha_shape=phase.p.alpha_shape,
                    use_centroids=phase.p.use_centroids,
                    use_alpha_shape=True,
                    single_cell_noise=0.5,
                    single_cell_sides=6,
                    image_mask=None,
                    allow_merging=True,
                    merge_thresh=0.1,
                    close_thresh=0.1,
                    make_all_operators=False,  # FIXME: later fix this; for now set to False
                    mesh_type=lattice_to_mesh_type[phase.p.cell_lattice_type],
                )

                # mesh refinement:
                self.mesh.clip_and_refine(image_mask, smoothing=None, refinement=phase.p.refine_mesh,
                                     max_steps=phase.p.maximum_voronoi_steps,
                                     convergence=phase.p.voronoi_convergence,
                                     fix_bounds=True)


        # Clip the Voronoi cluster to the shape of the clipping bitmap -------------------------------------------------
        self.ecm_verts = self.mesh.vcell_verts # voronoi verts of clipped cluster, nested as polygons defining each cell

        # self.voronoi_verts = []  # track all voronoi cells, even those not in cluster (used as grid for masking)

        self.cluster_mask = image_mask.clipping_matrix  # keep track of cluster mask and its size
        self.msize = image_mask.msize

        self.ecm_verts_unique = self.mesh.vor_verts  # convert to numpy array

    def cell_index(self, p) -> None:
        '''
        (Re)define the array of all cell centres (i.e., :attr:`cell_centres`)
        from every corresponding Voronoi polygon.

        This method is typically called *after* creating, closing, and clipping
        the Voronoi diagram for this cell cluster. Doing so creates an ordering
        of cells consistent with all Voronoi polygons, membrane domains, and
        extracellular matrix (ECM)-driven polygons and segments.
        '''

        self.cell_centres = self.mesh.vor_cents

        self.cell_vol = self.mesh.vor_sa*p.cell_height

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

        #FIXME: This logic has been duplicated in both deformWorld() and
        #quickVerts() methods, which... isn't the best. Consider centralizing
        #these duplicated variable assignments into a common private method.

        # Finish up by creating indices vectors and converting to Numpy arrays where needed:
        self.cell_i = [x for x in range(0,len(self.cell_centres))]
        self.mem_i  = [x for x in range(0,len(self.mem_mids_flat))]

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

        #Indexing by "[self.mem_i]" is superfluous here. Lively jumping snakes!
        self.mem_to_cells = indmap_mem[:,0]   # gives cell index for each mem_i index placeholder

        # construct a mapping giving membrane index for each cell_i------------------------------------------------
        self.cell_to_mems = []

        for cell_index in self.cell_i:
            # One-dimensional Numpy array of the indices of all membranes of
            # the current cell. Note that:
            #
            # * "self.mem_to_cells == cell_index" is a Numpy boolean array
            # * ndarray.nonzero() returns an n-tuple of the indices of all items
            #   of this boolean array evaluating to True, where "n" is the
            #   number of dimensions of the "mem_to_cells" array. Since this
            #   array is one-dimensional, this is guaranteed to be a 1-tuple
            #   whose first item providing this array is indexed at "[0]".
            #
            #It is confusing. It is Numpy.
            mems_index = (self.mem_to_cells == cell_index).nonzero()[0]
            self.cell_to_mems.append(mems_index)

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
        cellVertTree = cKDTree(self.mem_verts)

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
        self.mem_i =  [x for x in range(0,len(self.mem_mids_flat))]

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

    def cellMatrices(self, p) -> None:
        '''
        Create the main matrices used in routine calculations on the Cell Grid.
        '''

        # Log this action.
        logs.log_info(
            'Creating computational matrices for cell-cell transfers...')

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

        self.M_sum_mems_inv = np.linalg.pinv(self.M_sum_mems)  # matrix inverse of M_sum_mems for div-free cell calcs
        self.num_mems = np.asarray(self.num_mems)  # number of membranes per cell
        self.mem_distance = p.cell_space + 2*p.tm # distance between two adjacent intracellluar spaces
        self.cell_number = self.cell_centres.shape[0]

        # matrix for calculating gradients around the cell circumference:
        self.gradTheta = np.zeros((len(self.mem_i), len(self.mem_i)))

        # matrix for averaging points from mem-mids to companion polygon midpoints:
        # self.aveTheta = np.zeros((len(self.mem_i), len(self.mem_i)))

        # matrix storing the radial length:
        self.radial_len = np.zeros(len(self.mem_i))

        for cell_i, mem_i in enumerate(self.cell_to_mems):

            mem_io = np.roll(mem_i, 1)

            # distance between points:
            li = self.mem_mids_flat[mem_i] - self.mem_mids_flat[mem_io]
            lm = np.sqrt(li[:, 0] ** 2 + li[:, 1] ** 2)

            self.radial_len[mem_i] = np.abs(lm)

            self.gradTheta[mem_i, mem_i] = 1 / lm
            self.gradTheta[mem_i, mem_io] = -1 / lm
            #
            # self.aveTheta[mem_i, mem_i] = 1/2
            # self.aveTheta[mem_i, mem_io] = 1/2

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

    def cell_vols(self, p) -> None:
        '''
        Divide each cell into a spider-web like volume with an outer "pie-box"
        region volume associated with each membrane, and a region around the
        central point (centroid region).

        This function calculates the volumes of cells, membrane "pie-boxes" and
        centroids, and associated parameters of **cell chords** (i.e., the
        distance between the cell center and each membrane midpoint).
        '''

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

    def mem_processing(self, p) -> None:
        '''
        Calculate the membrane-specific nearest neighbour pairs, interaction
        with the extracellular data points of the Voronoi grid lattice, and
        membranes on the boundary.
        '''

        # Log this action.
        logs.log_info('Creating gap junctions... ')

        #-- find nearest neighbour cell-cell junctions via adjacent membranes-------------------------------------------
        sc = 3.0*(1-p.scale_cell)*p.cell_radius
        memTree = cKDTree(self.mem_mids_flat)

        mem_nn_o = memTree.query_ball_point(self.mem_mids_flat, sc)
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

        # Perfect bflags mems:

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

        # linear k index:
        self.index_k = [x for x in range(0,len(self.xypts))]

        # FIXME temporary addition of block of sphaghetti code for later use/improvement:
        #--------------------------------------------------------------
        # Create singular ecm points and data-structures mapping from cell to ecm nn
        # where 'ecm' is different than the environmental grid points defined in self.xypts and corresponds
        # to the shared membrane midpoint between two cells (e.g. a true 'extracellular' matrix point).

        bflags_ecm_o = [] # flags for ecm points at boundaries
        ecm_points_o = [] # the points, prior to processing

        for i, (mem_i, mem_j) in enumerate(self.mem_nn): # for each membrane mid pair in mem near neighbours...

            xmem = self.mem_mids_flat[:, 0]
            ymem = self.mem_mids_flat[:, 1]

            xi = (xmem[mem_i] + xmem[mem_j])/2 # average the x and y points
            yi = (ymem[mem_i] + ymem[mem_j])/2

            ecm_points_o.append([xi, yi]) # append the pair to the points list

            if mem_i == mem_j:         # if membrane indices are equal, they must be a boundary:
                bflags_ecm_o.append([xi, yi])

        ecm_points_unique = np.unique(ecm_points_o, axis = 0)   # get only the unique points from the above
        bflags_ecm_unique = np.unique(bflags_ecm_o, axis = 0)

        ecm_list = ecm_points_unique.tolist()     # convert them back to a list
        bflags_list = bflags_ecm_unique.tolist()

        bflags_ecm = []

        for bf in bflags_list:

            bflags_ecm.append(ecm_list.index(bf))   # get the proper index for the boundary flag from points list

        self.ecm_points = np.asarray(ecm_points_unique)  # assign final data structures
        self.bflags_ecm = np.asarray(bflags_ecm)

        self.ecm_i = np.asarray([ii for ii in range(len(self.ecm_points))])

        xmem_mids = (self.mem_mids_flat[self.mem_nn[:,0]]
                     + self.mem_mids_flat[self.mem_nn[:,1]])/2

        points_tree_ecm = cKDTree(self.ecm_points)


        mem_to_ecm = []
        for xm in xmem_mids:

            _,ei = points_tree_ecm.query(xm)
            mem_to_ecm.append(ei)

        self.mem_to_ecm = np.asarray(mem_to_ecm)


        cell_to_ecm = []
        for ci, mem_is in enumerate(self.cell_to_mems):

            mem_pts = xmem_mids[mem_is] # get membrane mids surrounding this cell
            _, eci = points_tree_ecm.query(mem_pts)
            cell_to_ecm.append(eci.tolist())

        self.cell_to_ecm = np.asarray(cell_to_ecm)

        self.all_points = np.vstack((self.cell_centres, self.ecm_points))
        self.all_points_imap = np.hstack((self.cell_i, self.ecm_i))
        self.all_i = np.asarray([ii for ii in range(len(self.all_points))])
        self.all_points_ecm_i = len(self.cell_i) + self.ecm_i
        self.all_points_cell_i = self.cell_i
        self.adl = len(self.all_points)

        points_tree_mems = cKDTree(self.mem_mids_flat)
        ecm_to_mems_o = points_tree_mems.query_ball_point(self.ecm_points, r = 10*p.cell_space)

        ecm_to_mems = []
        ecm_to_cells = []
        for i, ecm_is in enumerate(ecm_to_mems_o):

            ecm_to_mems.append(ecm_is)
            ecm_to_cells.append(self.mem_to_cells[ecm_is])

        self.ecm_to_mems = np.asarray(ecm_to_mems)
        self.ecm_to_cells = np.asarray(ecm_to_cells)

        # Surface area and volume of ecm spaces:
        self.ecm_vol = []
        self.ecm_sa = []

        for mis in self.ecm_to_mems:
            esa = self.mem_sa[mis[0]]
            self.ecm_vol.append(esa*p.cell_space)
            self.ecm_sa.append(esa*1)

        self.ecm_vol = np.asarray(self.ecm_vol)
        self.ecm_sa = np.asarray(self.ecm_sa)

        self.ecmdl = len(self.ecm_points)

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
        self.points_tree = cKDTree(self.xypts)

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

        # create the matrix that allows individual membrane normal fluxes to be mapped to each ecm square:
        # If Fmem is the normal component of a vector field wrt individual membranes,
        # the result of M_divmap_mem2ecm *dot* Fmem  is the divergence of the flux wrt the environment.
        self.M_divmap_mem2ecm = np.zeros((len(self.xypts), len(self.mem_i)))

        for mem_i, ecm_i in enumerate(self.map_mem2ecm):
            mem_sa = self.mem_sa[mem_i]
            self.M_divmap_mem2ecm[ecm_i, mem_i] += (mem_sa)
            # self.M_divmap_mem2ecm[ecm_i, mem_i] += (mem_sa) / (p.cell_height*(self.delta**2))

    def graphLaplacian(self, p) -> None:
        '''
        Define an abstract inverse Laplacian used to solve Poisson's equation
        on the irregular Voronoi grid of this cell cluster.

        Parameters
        ----------
        p : Parameters
            Current simulation configuration.

        Creates
        ----------
        self.lapGJinv          Solver for Poisson equation with Neumann (zero gradient) boundary
        '''

        # Log this action.
        logs.log_debug('Creating cell network Poisson solver...')


        #----DEC matrix creation
        # Hodge star for edge length ratios:
        star_eij = np.diag(self.mesh.vor_edge_len/self.mesh.tri_edge_len)

        # First term in the laplacian matrix:
        L1 = np.dot(star_eij, self.mesh.delta_tri_0)

        star_a = np.diag(1/self.mesh.vor_sa)

        L2 = np.dot(star_a, -self.mesh.delta_tri_0.T)

        # Inverse terms:
        # Hodge star for edge length ratios:
        star_eij_inv = np.diag(self.mesh.tri_edge_len/self.mesh.vor_edge_len)

        # First term in the laplacian matrix:
        L1_inv = np.dot(self.mesh.delta_tri_0_inv, star_eij_inv)

        star_a_inv = np.diag(self.mesh.vor_sa)

        L2_inv = np.dot(-self.mesh.delta_tri_0_inv.T, star_a_inv)


        self.lapGJinv = np.dot(L1_inv, L2_inv)

        # if p.td_deform is True:
        #     # if time0dependent deformation is selected, also save the direct Laplacian operator:
        self.lapGJ = np.dot(L2, L1)

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

    @type_check
    def redo_gj(self, phase: SimPhase) -> None:
        '''
        (Re)create the gap junction connection network after assessing tissue
        profile requests for the passed simulation phase.

        Parameters
        ----------
        phase : SimPhase
            Current simulation phase.
        '''

        # Localize frequently referenced phase variables for convenience.
        dyna = phase.dyna
        p    = phase.p

        flag_cell_nn = [ [] for x in range(0,len(self.cell_i))]

        for tissue_name, tissue_profile in dyna.tissue_name_to_profile.items():
            # Step through gj's and find cases where connection is split between
            # cells in different tissues.
            if tissue_profile.is_gj_insular:
                # Get the cell target inds for this tissue.
                cell_targets = dyna.cell_target_inds[tissue_name]

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

        self.cell_nn_connected = np.asarray(new_cell_nn)

        # Redo the number and average nearest neighbours per cell:
        self.num_nn = []  # initialize a list that will hold number of nns to a cell

        for indices in self.cell_nn:
            self.num_nn.append(len(indices))

        self.average_nn = (sum(self.num_nn)/len(self.num_nn))
        self.num_nn = np.asarray(self.num_nn)
        self.gj_default_weights = np.ones(len(self.mem_i))

        for cell_i, nn_cell_i_set in enumerate(self.cell_nn_connected):
            mem_i_set = self.cell_to_mems[cell_i]  # get all the membranes for this cell

            for mem_i in mem_i_set:
                mem_j = self.nn_i[mem_i]  # get the current neighbour mem and cell...

                if mem_j != mem_i:  # if we're not on a boundary
                    cell_j = self.mem_to_cells[mem_j]

                    if cell_j not in nn_cell_i_set:  # if the partner cell is no longer listed as a nn...
                        # ...then the gap junction weight to zero to fully inhibit coupling:
                        self.gj_default_weights[mem_i] = 0.0

                        #...then set both the membrane and cell neighbour spot to "self":
                        # self.nn_i[mem_i] = mem_i
                        # self.cell_nn_i[mem_i] = [cell_i,cell_i]

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

                self.nn_len.append(-1) # FIXME -- this seems like a horrific idea...

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

    @type_check
    def save_cluster(self, phase: SimPhase) -> None:
        '''
        Pickle (i.e., save) this cell cluster to the file configured by the
        passed seed simulation phase.

        Parameters
        ----------
        phase : SimPhase
            Current simulation phase.
        '''

        # Log this pickling.
        logs.log_debug(
            'Pickling cell cluster to: %s', phase.p.seed_pickle_filename)

        # Pickle this cell cluster.
        datadump = [self, phase.p]
        fh.saveSim(phase.p.seed_pickle_filename, datadump)

    def make_maskM(self,p):
        """
        Create structures for plotting interpolated data on cell centres
        and differentiating between the cell cluster and environment.
        """

        voronoiTree = cKDTree(self.xypts)
        _, self.map_voronoi2ecm = voronoiTree.query(self.ecm_verts_unique)

        self.voronoi_mask = np.zeros(len(self.xypts))
        self.voronoi_mask[self.map_voronoi2ecm]=1

        xv = np.linspace(self.xmin,self.xmax,p.plot_grid_size)
        yv = np.linspace(self.xmin,self.xmax,p.plot_grid_size)

        X,Y = np.meshgrid(xv,yv)

        self.Xgrid = X
        self.Ygrid = Y

        self.maskM = interp.griddata(
            (self.xypts[:,0],self.xypts[:,1]),
            self.voronoi_mask,(self.Xgrid,self.Ygrid),
            method='linear',fill_value=0)

        self.maskM = ndimage.filters.gaussian_filter(self.maskM, 1, mode='nearest')
        self.maskM = np.round(self.maskM,0)

        self.maskECM = interp.griddata(
            (X.ravel(),Y.ravel()),self.maskM.ravel(),
            (self.X, self.Y),
            method='linear',fill_value=0)
        self.maskECM = ndimage.filters.gaussian_filter(self.maskECM, 1, mode='nearest')
        self.maskECM = np.round(self.maskECM,0)

        self.inds_env = list(*(self.maskECM.ravel() == 0).nonzero())
        self.inds_clust = list(*(self.maskECM.ravel() == 1).nonzero())

        self.voronoi_centres = self.xypts

        vertTree = cKDTree(self.voronoi_centres)
        _, self.cell_to_grid = vertTree.query(self.cell_centres)

    def intra_updater(self,p):
        """
        Calculates a matrix that takes values on membrane midpoints,
        interpolates them to cell vertices, and adds them all together as half
        of a finite volume integration for the pie-box regions. The other half
        will come from the centroid region value.
        """

        self.M_int_mems = np.zeros((len(self.mem_i), len(self.mem_i)))

        #FIXME: "i" appears to be unused here. In optimistic theory, this
        #implies that this loop should be reducible to:
        #    for inds in self.cell_to_mems:
        #Maybe? May the misty dawn exhale its hot breath upon you!
        for i, inds in enumerate(self.cell_to_mems):

            # get the set of indices for the cell:
            inds = np.asarray(inds)

            inds_p1 = np.roll(inds, 1)
            inds_o = np.roll(inds, 0)
            inds_n1 = np.roll(inds, -1)

            self.M_int_mems[inds_o, inds_o] = (1/3)
            self.M_int_mems[inds_o, inds_p1] = (1/12)
            self.M_int_mems[inds_o, inds_n1] = (1/12)

    def deform_tools(self,p):

        # Data structures specific for deformation option------------------------------

        logs.log_info('Creating computational tools for mechanical deformation... ')

        # create a data structure that will allow us to repackage ecm_verts and re-build the
        # cells world after deforming ecm_verts_unique:

        # first get the search-points tree:
        ecmTree = cKDTree(self.ecm_verts_unique)

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

        #FIXME: "i" appears to be unused here. In optimistic theory, this
        #implies that this loop should be reducible to:
        #    for inds in self.cell_to_mems:
        #Maybe? Unshroud the penultimate technique, Dagalfor!
        for i, inds in enumerate(self.cell_to_mems):

            inds = np.asarray(inds)

            inds_p1 = np.roll(inds,1)
            inds_o = np.roll(inds,0)

            dist = self.mem_mids_flat[inds_p1] - self.mem_mids_flat[inds_o]
            len_mem = np.sqrt(dist[:,0]**2 + dist[:,1]**2)

            self.gradMem[inds_o,inds_p1] = 1/len_mem.mean()
            self.gradMem[inds_o,inds_o] = -1/len_mem.mean()

    # SVG processing functions-----------------------------------------------------------------------------------------
    @type_check
    def parse_svg(self, p: 'betse.science.parameters.Parameters') -> None:
        '''
        Parse the passed SVG-formatted file into a clipping bitmask, cell
        centers, tissue profiles, and related cell cluster metadata.

        Parameters
        ----------
        fname_cells : str
            Absolute or relative filename of the SVG-formatted file defining cell centres to be parsed.

        fname_clip : str
            Absolute or relative filename of the SVG-formatted file defining clipping curve to be parsed.
        '''

        doc = minidom.parse(p.svg_cells_fname)
        circle_objects = doc.getElementsByTagName('circle')

        circle_x = np.asarray([np.float64(co.getAttribute('cx')) for co in circle_objects])
        circle_y = np.asarray([np.float64(co.getAttribute('cy')) for co in circle_objects])

        circle_x = (circle_x) / p.svg_size
        circle_y = (p.svg_size - circle_y) / p.svg_size

        xypts = np.column_stack((circle_x*p.wsx, circle_y*p.wsx))

        fill_regex = regexes.compile_regex(r'^.*?(?:^|;)fill:#([0-9a-fA-F]{6})(?:$|;)')

        self.seed_fills = nparray.from_iterable(
            regexes.get_match_group_first(
                text=co.getAttribute('style'), regex=fill_regex)
            for co in circle_objects
        )

        doc.unlink()

        # define geometric limits and centre for the cluster of points
        self.xmin = 0
        self.xmax = p.wsx
        self.ymin = 0
        self.ymax = p.wsx

        self.bbox = np.asarray(
            [[self.xmin, self.ymin],
             [self.xmax, self.ymin],
             [self.xmax, self.ymax],
             [self.xmin, self.ymax]])

        self.centre = np.asarray([0.5*p.wsx, 0.5*p.wsx])
        self.clust_xy = xypts

        self.xyaxis = [self.xmin, self.xmax, self.ymin, self.ymax]

        #-----Utility functions--------------------------------------------------------------------------------------------

    def gradient(self, SS):
        """
        Calculates the gradient based on differences of a property SS between cell centres.
        The gradient is defined at each membrane.
        """
        gSS = (SS[self.cell_nn_i[:, 1]] - SS[self.cell_nn_i[:, 0]]) / (self.nn_len)

        gSx = gSS * self.mem_vects_flat[:, 2]
        gSy = gSS * self.mem_vects_flat[:, 3]

        return gSS, gSx, gSy

    def meanval(self, SS):
        """
        Calculates the mean value based on average of property SS between nearest-neighbour cell centres.
        """
        if SS.shape[0] == len(self.cell_i): # if parameter defined on cell centres:
            mSS = (SS[self.cell_nn_i[:, 1]] + SS[self.cell_nn_i[:, 0]]) / 2

        elif SS.shape[0] == len(self.mem_i): # if parameter defined on membrane mids:

            mSS = (SS[self.nn_i] + SS[self.mem_i])/2

        return mSS

    def average_vector(self, Smx, Smy):
        """
        Takes vector quantity (Smx, Smy) defined at membranes and calculates the averaged
        single vector at the cell centre.
        """
        Scx = np.dot(self.M_sum_mems, Smx * self.mem_sa) / self.cell_sa
        Scy = np.dot(self.M_sum_mems, Smy * self.mem_sa) / self.cell_sa

        return Scx, Scy

    def mag(self, Sx, Sy):
        """
        Returns the magnitude of vector (Sx, Sy) defined on any system.
        """

        So = np.sqrt(Sx ** 2 + Sy ** 2)

        return So

    def mem_normal_component(self, Sx, Sy):
        """
        Returns the normal component of a vector (Sx, Sy) at each membrane.
        """

        if Sx.shape[0] == len(self.cell_i):

            Sxm = Sx[self.mem_to_cells]
            Sym = Sy[self.mem_to_cells]

        elif Sx.shape[0] == len(self.mem_i):

            Sxm = Sx
            Sym = Sx

        else:
            logs.log_error("Shape of input is wrong!")

        Sn = Sxm * self.mem_vects_flat[:, 2] + Sym * self.mem_vects_flat[:, 3]

        return Sn

    def div(self, gSx, gSy, cbound = False):
        """
        Calculates the divergence of a vector (gSx, gSy) defined on membranes.
        """

        gSn = gSx * self.mem_vects_flat[:, 2] + gSy * self.mem_vects_flat[:, 3]

        if cbound is True: # close the boundary (zero-flux boundary condition)
            gSn[self.bflags_mems] = 0.0

        divS = np.dot(self.M_sum_mems, gSn * self.mem_sa) / self.cell_vol

        return divS

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
            raise BetseSimConfException("Input to cells.curl not defined properly."
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

        Phi = np.dot(self.lapGJinv, div_F + rho)

        gPhi = (Phi[self.cell_nn_i[:, 1]] - Phi[self.cell_nn_i[:, 0]]) / (self.nn_len)

        gPx = gPhi*self.nn_tx
        gPy = gPhi*self.nn_ty

        # make the field divergence-free:
        Fx = fxo - gPx
        Fy = fyo - gPy

        Fn = Fx*self.cell_vects_flat[:,2] + Fy*self.cell_vects_flat[:,3]


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

        # if bounds_closed is True:

        AA = np.dot(self.lapGJinv, -curlF)

        Ax, Ay, _ = self.curl(0, 0, AA)

        # -----------------

        if rot_only is False:

            Fxm = self.interp_to_mem(Fx)
            Fym = self.interp_to_mem(Fy)

            Fn = Fxm * nx + Fym * ny

            divF = np.dot(self.M_sum_mems, Fn * self.mem_sa) / self.cell_vol

            BB = np.dot(self.lapGJinv, divF)

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

    def single_cell_div_free(self, uxo, uyo):
        # now, make the transport field divergence-free wrt individual cells (divergence-free is the way to be!
        divU = self.div(uxo, uyo, cbound=False)  # divergence of the field at each membrane
        Pi = np.dot(self.M_sum_mems_inv, divU) * (
            self.cell_vol[self.mem_to_cells] / self.mem_sa)  # 'pressure" field to create div-free case

        ux = uxo - Pi * self.mem_vects_flat[:, 2]  # corrected vector at the membrane
        uy = uyo - Pi * self.mem_vects_flat[:, 3]  # corrected vector at the membrane

        return ux, uy

    # ..........{ PROPERTIES ~ membrane                  }.....................
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

    # ..........{ PROPERTIES ~ mappers                   }.....................
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
        ``n`` are as defined above.
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
        * ``1/k`` if this cell contains this membrane, where ``k`` is the
          number of membranes this cell contains.

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
        # * "self.M_sum_mems.T", the transpose of the "M_sum_mems" matrix.
        #   Since this matrix is of size m x n for m the number of cells and n
        #   is the number of membranes, this transpose is a matrix of size n x
        #   m. Each element of this transpose is either:
        #   * 0 if this cell does *NOT* contain this membrane.
        #   * 1 if this cell contains this membrane.
        # * "... / self.num_mems", normalizing each cell membrane element of
        #   this matrix by the number of membranes in that cell. Since
        #   "num_mems" is a row vector of length m whose elements are the
        #   number of membranes in that cell, each column of this transpose is
        #   divided by the corresponding element of this row vector.
        return self.M_sum_mems.T / self.num_mems

    # ..........{ MAPPERS                                }.....................
    #FIXME: To reduce code duplication:
    #
    #    # Globally replace all instances of this...
    #    np.dot(cells.M_sum_mems, some_array) / cells.num_mems
    #
    #    # ...with this.
    #    cells.calculate_divergence(some_array)
    #
    #    # ...which should internally do something like this:
    #    np.dot(some_array, cells.M_sum_mems.T) / cells.num_mems
    #
    #After doing so, consider replacing the "M_sum_mems" array with its
    #tranpspose, which is a # the more general-purpose and efficient array.
    #Unlike the former, the latter is applicable to both arrays and matrices.
    @type_check
    def map_membranes_midpoint_to_cells_centre(
        self, membranes_midpoint_data: SequenceTypes) -> ndarray:
        """
        Convert the passed one- or two-dimensional sequence of arbitrary data
        spatially situated at cell membrane midpoints into a Numpy array of the
        same dimensionality and data spatially situated at cell centres.

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
        membranes_midpoint_data = nparray.from_iterable(membranes_midpoint_data)

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

    # ..........{ MAPPERS ~ cells centre                  }.....................
    def map_cells_centre_to_grids_centre(
        self, *args, **kwargs) -> ndarray:
        """
        Convert the passed one- or two-dimensional sequence of arbitrary data
        spatially situated at cell centres into a Numpy array of the same
        dimensionality and data spatially situated at environmental (i.e.,
        extracellular) grid space centres.

        See Also
        -----------
        :meth:`map_cells_centre_to_points`
            Further details.
        """

        return self.map_cells_centre_to_points(
            *args,

            # 2-tuple of the target X and Y coordinates of all grid space
            # centres to interpolate this data onto.
            target_points=(self.X, self.Y),

            # Mask nullifying all output data mapped onto extracellular spaces.
            # While the "fill_value" parameter internally passed to the
            # griddata() function by this call should already do so, assumptions
            # are a dangerous mistress.
            data_factor=self.maskECM,

            # All remaining arguments.
            **kwargs
        )


    def map_cells_centre_to_membranes_midpoint(
        self, *args, **kwargs) -> ndarray:
        """
        Convert the passed one- or two-dimensional sequence of arbitrary data
        spatially situated at cell centres into a Numpy array of the same
        dimensionality and data spatially situated at cell membrane midpoints.

        See Also
        -----------
        :meth:`map_cells_centre_to_points`
            Further details.
        """

        return self.map_cells_centre_to_points(
            *args,

            # 2-tuple of the target X and Y coordinates of all cell membrane
            # midpoints to interpolate this data onto.
            target_points=(
                self.mem_mids_flat[:, 0], self.mem_mids_flat[:, 1]),

            # All remaining arguments.
            **kwargs
        )


    def map_cells_centre_to_points(
        self,

        # Mandatory parameters.
        cells_centre_data: SequenceTypes,
        target_points: SequenceTypes,

        # Optional parameters.
        interp_method: str = 'linear',
        data_factor: NumericOrSequenceTypes = 1,
    ) -> ndarray:
        """
        Convert the passed one- or two-dimensional sequence of arbitrary data
        spatially situated at cell centres into a Numpy array of the same
        dimensionality and data spatially situated at the passed target points.

        Parameters
        -----------
        cells_centre_data : SequenceTypes
            One- or two-dimensional sequence of length the number of cells such
            that each element is arbitrary cell data spatially situated at the
            centre of that cell.
        target_points : SequenceTypes
            Two-dimensional sequence of the coordinates of all target points to
            interpolate this cell data onto, whose:
            #. First dimension indexes first the X and then Y axis.
            #. Second dimension indexes each target point such that each element
               is the coordinate of that point along this axis.
            For each target point, the output data returned for source data
            defined at cells whose membranes overlap this target point is
            defined as follows:
            * If this target point resides *inside* the convex hull of this cell
              cluster, this output data is spatially interpolated from the
              centres of these cells onto this target point.
            * Else, 0. In this case, this target point resides *outside* the
              convex hull of this cell cluster. The source data is spatially
              situated at cell centres and hence perfectly intracellular. No
              extracellular source data exists with which to interpolate onto
              extracellular target points.
        interp_method : optional[str]
            Interpolation type to pass to the
            :func:`scipy.interpolate.gridddata` function (e.g., ``nearest``,
            ``linear``, ``cubic``). Defaults to ``linear``.
        data_factor : NumericOrSequenceTypes
            Integer, float, or one-dimensional sequence of integers or floats
            by which to multiply all elements of the returned array. Defaults to
            1, in which case this array is returned as is. Specifically, if this
            array is:
            * One-dimensional, this array itself is multiplied by this factor.
            * Two-dimensional, each element of the first dimension of this array
              is multiplied by this factor.

        Returns
        -----------
        ndarray
            One- or two-dimensional Numpy array of length the number of target
            points such that each element is arbitrary cell data spatially
            interpolated onto the nearast target point.
        """

        # Numpy arrays converted from the passed sequences.
        cells_centre_data = nparray.from_iterable(cells_centre_data)

        # If this source data is neither one- nor two-dimensional, raise an
        # exception.
        if len(cells_centre_data.shape) not in {1, 2}:
            raise BetseSequenceException(
                'Source data dimensionality {} neither one- nor '
                'two-dimensional.'.format(len(cells_centre_data.shape)))

        # If the last dimension of this source data does *NOT* index all cells
        # in this cluster, raise an exception.
        if cells_centre_data.shape[-1] != self.cell_centres.shape[0]:
            raise BetseSequenceException(
                'Source data not spatially situated at cell centres '
                '(i.e., last dimension length {} not {}).'.format(
                    cells_centre_data.shape[-1], self.cell_centres.shape[0]))

        # If the first dimension of these target points are *NOT* coordinate
        # pairs, raise an exception.
        if len(target_points) != 2:
            raise BetseSequenceException(
                'Target points not coordinate pairs '
                '(i.e., first dimension length {} not 2).'.format(
                    len(target_points)))

        # 2-tuple of the X and Y coordinates of all cell centres.
        cell_centres = (self.cell_centres[:, 0], self.cell_centres[:, 1])

        # Dictionary of all keyword arguments to pass to all calls to the
        # griddata() function performed below.
        griddate_kwargs = {
            # 2-tuple of all source X and Y coordinates to interpolate from.
            'points': cell_centres,

            # 2-tuple of all target X and Y coordinates to interpolate onto.
            'xi': target_points,

            # Machine-readable string specifying the interpolation type.
            'method': interp_method,

            # Default data value to assign all output points residing outside
            # the convex hull of the cell centres being interpolated from, which
            # are thus non-interpolatable. For safety, this data is nullified.
            # The default "fill_value" is NaN, which is absurdly unsafe.
            'fill_value': 0,
        }

        # If this source data is one-dimensional, map this data from cell
        # centres onto target points via a single interpolation call.
        if len(cells_centre_data.shape) == 1:
            return data_factor * interp.griddata(
                values=cells_centre_data, **griddate_kwargs)
        # Else, this data is two-dimensional. Since SciPy requires the second
        # parameter to the interp.griddata() function be one-dimensional, this
        # is mappable from cell centres onto target points *ONLY* via multiple
        # interpolation calls. While inefficient, SciPy offers no alternatives.
        else:
            # Two-dimensional output list aggregating all interpolation results.
            cells_centre_data_interpolated = []

            # For each one-dimensional input array of source data in this
            # two-dimensional input array of such data...
            for cells_centre_data_one in cells_centre_data:
                # One-dimensional output array interpolated from this array.
                cells_centre_data_one_interpolated = (
                    data_factor * interp.griddata(
                        values=cells_centre_data_one, **griddate_kwargs))

                # Append this output array to this output list.
                cells_centre_data_interpolated.append(
                    cells_centre_data_one_interpolated)

            # Return this output list converted back to an output array.
            return nparray.from_iterable(cells_centre_data_interpolated)
