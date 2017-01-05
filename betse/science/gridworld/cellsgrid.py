#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

"""
This module contains the class GridWorld, which is specific to quasi-continuous simulations.
GridWorld holds all data structures relating to the size of the environment,
the extent of the cell cluster, the co-ordinates of cell
centre points, and all kinds of data relating to individual cell properties.

GridWorld is defined on a rectilinear grid with equal dimensions in the x and y directions.

GridWorld assumes each tick in the grid is the potential location of a cell and its extracellular
matrix surrounds.

The initialization method of the Cells class sets-up
and crops the cell cluster to an optional user-defined geometry input
(a set of points arranged in counter-clockwise order and
defining a closed polygon). Other methods define the cell centres and create cell-cell gap junctions (GJs).



"""

# FIXME when running gridworld, make sure to set p.sim_ECM to false so that everything works OK in downstream
# modules...!
#FIXME do an ecm volume

import math
import os
import os.path

import numpy as np
import scipy.spatial as sps

from betse.science import filehandling as fh
from betse.science import toolbox as tb
from betse.science.tissue.bitmapper import BitMapper
from betse.util.io.log import logs


class GridWorld(object):

    def __init__(self,p):

        self.do_once_cavity = True  # boolean to ensure any cavity preparation happens only once...
        self.do_once_cuts = True      #boolean to ensure cut lines are only prepared once...

        # initialize some parameters that may or may not be used...
        self.cavity_inds = []
        self.true_env_inds = []
        self.cavity_volume = []

        self.fileInit(p)
        self.makeSeeds(p)
        self.makeCluster(p)
        # self.near_neigh(p)
        self.environment(p)

        self.generalMask = self.makeMask(mask_type='exterior bound')
        self.cellMask = self.makeMask(mask_type = 'cluster bound')

        logs.log_info('Creating environmental Poisson solver...')

        self.Ainv = self.makeLaplacian(self.generalMask)

        logs.log_info('Creating cluster-specific Poisson solver...')
        self.Ainv_cells = self.makeLaplacian(self.cellMask)

        # loggers.log_info('Cell cluster creation complete!')

        # # save the cell cluster
        # loggers.log_info('Saving the cell cluster... ')
        #
        # # celf = copy.deepcopy(self)
        # datadump = [self,p]
        # fh.saveSim(self.savedWorld,datadump)
        # message = 'Cell cluster saved to' + ' ' + self.savedWorld
        # loggers.log_info(message)

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

    def makeSeeds(self,p):

        """
        makeSeeds returns an regular scatter
        of points defined on a world space
        with dimensions wsx, wsy in [m].


        Creates
        -------
        self.xypts      numpy array listing [x,y] of world seed points

        self.xmin, self.xmax      dimensions of world grid
        self.ymin, self.ymax

        self.centre     [x,y] coordinate of world centre

        Notes
        -------
        Uses Numpy arrays

        """

        # base parameter definitions
        self.delta = p.d_cell + p.cell_space  # spacing between grid points
        self.grid_n = int(p.wsx/self.delta) # number of grid points

        # linear vectors for spatial dimensions
        self.x_v = np.linspace(0, p.wsx, self.grid_n)  # create lattice vector x
        self.y_v = np.linspace(0, p.wsx, self.grid_n)  # create lattice vector y

        # next define a 2d array of lattice points using the x- and y- vectors
        self.X, self.Y = np.meshgrid(self.x_v, self.y_v)  # create 2D array of lattice points

        # define a data structure that holds [x,y] coordinate points of each 2d grid-matrix entry
        self.xypts = []
        self.map_ij2k = []

        for i, row in enumerate(self.X):
            for j, val in enumerate(row):
                xpt = self.X[i,j]
                ypt = self.Y[i,j]
                self.xypts.append([xpt,ypt])
                self.map_ij2k.append([i,j])

        # linear k index:
        self.index_k = [x for x in range(0,len(self.xypts))]

        # map for converting between k linear index and 2d i,j matrix
        self.map_ij2k = np.asarray(self.map_ij2k)
        self.xypts = np.asarray(self.xypts)
        self.grid_len = len(self.xypts)

        # define geometric limits and centre for the cluster of points
        self.xmin = np.min(self.X)
        self.xmax = np.max(self.X)
        self.ymin = np.min(self.Y)
        self.ymax = np.max(self.Y)

        self.centre = self.xypts.mean(axis=0)

    def makeCluster(self,p):
        """
        Uses bitmaps or circle mask to define a cluster of cells on the x,y grid.
        """

        # load the bitmap used to clip the cell cluster and create a clipping function
        logs.log_info('Clipping geometry to cluster shape... ')
        bitmasker = BitMapper(
            p.clipping_bitmap_matcher,
            self.xmin, self.xmax, self.ymin, self.ymax)

        self.cell_bools_k = []
        for point in self.xypts:
            cell_truth = float(*bitmasker.clipping_function(point[0],point[1]))
#           cell_truth = 1.0
            self.cell_bools_k.append(cell_truth)

        self.cluster_mask = bitmasker.clipping_matrix
        self.msize = bitmasker.msize
        self.cell_bools_k = np.asarray(self.cell_bools_k)

        for k, ind in enumerate(self.cell_bools_k):
            if ind > 0.0:
                self.cell_bools_k[k] = 1

        cell_inds = (self.cell_bools_k > 0.0).nonzero()
        self.cell_centres = self.xypts[cell_inds]

        # create an index for only the cell centres
        self.cell_i = [x for x in range(0,len(self.cell_centres))]

        # properties of individual cells
        self.cell_vol = (4/3)*math.pi*p.rc**3   # cell volume (unifrom)
        self.cell_sa = 4*math.pi*p.rc**2    # whole cell surface area (uniform)

        # calculate volume of local ecm space for each cell
        rext = p.rc + p.cell_space
        self.ecm_vol = ((4/3)*math.pi*(rext**3 - p.rc**3))/2

        # calculating centre, min, max of cell cluster
        self.clust_centre = np.mean(self.cell_centres)
        self.clust_x_max = np.max(self.cell_centres[:,0])
        self.clust_x_min = np.min(self.cell_centres[:,0])
        self.clust_y_max = np.max(self.cell_centres[:,1])
        self.clust_y_min = np.min(self.cell_centres[:,1])

        self.cell_number = int(np.sum(self.cell_bools_k))

    def near_neigh(self,p):

        """
        Calculate the nearest neighbours for each cell centre in the cluster and return a numpy
        array of nn indices with an index consistent with all other data lists for the cluster.

        Creates
        -------
        self.cell_nn            A nested list defining the indices of all nearest neighbours to each cell
        self.num_nn             A list enumerating the number of nn each cell connects with
        self.gj_vects           An array of gj midpoints and corresponding tangent vectors to each gj
        self.gap_jun_i          A list of index pairs to self.cell_centres, each pair defining a unique cell-cell GJ
        self.gj_i               The index to unique gap junctions (ordering of self.gj_vects and self.gap_jun_i)
        self.gj_len             The length of indiviual gap junctions
        self.gjMatrix           Allows matrix multiplication to deliver fluxes from each cell to neighbours
        self.cell_to_gj         A mapping giving gj_i indices for each position (ordered to cell_i)

        Notes
        -------
        Uses numpy arrays
        Uses scipy spatial KDTree search algorithm

        """

        logs.log_info('Creating gap junctions... ')

        cell_tree = sps.KDTree(self.cell_centres)
        self.cell_nn=cell_tree.query_ball_point(self.cell_centres,p.search_d*p.d_cell)

        self.num_nn = []  # initialize a list that will hold number of nns to a cell

        for indices in self.cell_nn:
            self.num_nn.append(len(indices) -1)  # minus one because query cell is included in each nn list

        self.average_nn = (sum(self.num_nn)/len(self.num_nn))

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

        self.gap_jun_i = np.asarray(self.gap_jun_i)


        # define an index to gap junction values
        self.gj_i = [x for x in range(0,len(self.gap_jun_i))]

        # calculate lengths of gap junctions
        seg1=self.cell_centres[self.gap_jun_i[:,0]]
        seg2=self.cell_centres[self.gap_jun_i[:,1]]
        nn_diff_gj = (seg2 - seg1)**2
        self.gj_len = np.sqrt(nn_diff_gj[:,0] + nn_diff_gj[:,1])

        # compute mapping between cell and gj:
        self.cell_to_gj =[[] for x in range(0,len(self.cell_i))]

        for i, inds in enumerate(self.gap_jun_i):
            ind1 = inds[0]
            ind2 = inds[1]
            self.cell_to_gj[ind1].append(i)
            self.cell_to_gj[ind2].append(i)

        self.cell_to_gj = np.asarray(self.cell_to_gj)

        # calculating matrix for gap junction flux calculation between cells
        # divide by the number of nearest neighbours to properly distribute each flux
        # over a fraction of the total cell surface area
        self.gjMatrix = np.zeros((len(self.gj_i),len(self.cell_i)))
        for igj, pair in enumerate(self.gap_jun_i):
            ci = pair[0]
            cj = pair[1]
            nn_i = self.num_nn[0]
            nn_j = self.num_nn[1]
            self.gjMatrix[igj,ci] = -1
            self.gjMatrix[igj,cj] = 1

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

        logs.log_info('Tagging environmental boundary points... ')

        con_hull = tb.alpha_shape(points, alpha/p.d_cell)  # get the concave hull for the membrane midpoints
        con_hull = np.asarray(con_hull)

        bflags = np.unique(con_hull)    # get the value of unique indices from segments

        bmask = np.array((points[:,0]))
        bmask[:] = 0
        bmask[bflags] = 1

        return bflags, bmask

    def environment(self,p):

        """
        Defines conditions for points in contact with the global environment at the outer boundary
        of the square world.

        Note: this is how to access all elements in a main-grid format from the k-vector:
        Z[cells.map_ij2k[:,0],cells.map_ij2k[:,1]]

        or access elements in a main-grid format to a subset of the k-vector:

        accessing bottom boundary:
        Z[cells.map_ij2k[bBot_k][:,0],cells.map_ij2k[bBot_k][:,1]]

        accessing all cells:
        Z[cells.map_ij2k[map_cell2ecm][:,0],cells.map_ij2k[map_cell2ecm][:,1]] = 10.0

        """

        logs.log_info('Setting global environmental conditions... ')

        # first obtain a structure to map to total xypts vector index:
        points_tree = sps.KDTree(self.xypts)

        # define a mapping between a cell and its ecm space in the full list of xy points for the world:
        self.map_cell2ecm = list(points_tree.query(self.cell_centres))[1]

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

        self.bBot_k = list(points_tree.query(bBot_pts))[1]
        self.bTop_k = list(points_tree.query(bTop_pts))[1]
        self.bL_k = list(points_tree.query(bL_pts))[1]
        self.bR_k = list(points_tree.query(bR_pts))[1]


        # get a list of k indices to the four 1st interior (global) boundaries of the rectangular world:
        bBot_xp = self.X[1,1:-1]
        bTop_xp = self.X[-2,1:-1]
        bL_xp = self.X[1:-1,1]
        bR_xp = self.X[1:-1,-2]

        bBot_yp = self.Y[1,1:-1]
        bTop_yp = self.Y[-2,1:-1]
        bL_yp = self.Y[1:-1,1]
        bR_yp = self.Y[1:-1,-2]

        bBot_pts_p = np.column_stack((bBot_xp, bBot_yp))
        bTop_pts_p = np.column_stack((bTop_xp, bTop_yp))
        bL_pts_p = np.column_stack((bL_xp, bL_yp))
        bR_pts_p = np.column_stack((bR_xp, bR_yp))

        self.bBot_kp = list(points_tree.query(bBot_pts_p))[1]
        self.bTop_kp = list(points_tree.query(bTop_pts_p))[1]
        self.bL_kp = list(points_tree.query(bL_pts_p))[1]
        self.bR_kp = list(points_tree.query(bR_pts_p))[1]

        # # for purposes of Laplacian solver, need to remove boundary elements and define a "core"
        # # of interior points:
        # core_x = self.X[1:-1,1:-1]
        # core_y = self.Y[1:-1,1:-1]
        # core_pts = np.column_stack((core_x.ravel(), core_y.ravel()))
        #
        # self.core_k = list(points_tree.query(core_pts))[1]
        #
        # self.core_len = len(self.core_k)
        #
        # self.core_n = self.grid_n - 2
        #
        # # create a mapping from the ip, jp indexing to the kp index of the unraveled core matrix:
        # M = []
        # for i in range(self.core_n):
        #     row = []
        #     for j in range(self.core_n):
        #         row.append([i,j])
        #
        #     M.append(row)
        #
        # self.map_ij2k_core = []
        # for row in M:
        #     for inds in row:
        #         self.map_ij2k_core.append(inds)

        # create a list of indices to the environmental space

        self.map_envSpace = []

        for ind in self.index_k:

            try:
                self.map_cell2ecm.tolist().index(ind)

            except:
                self.map_envSpace.append(ind)

        # finally, get a list of indices to points on the cell cluster boundary:

        bpoints, _ = self.boundTag(self.xypts[self.map_cell2ecm],p)

        self.bound_pts_xy = self.xypts[self.map_cell2ecm[bpoints]]

        self.bound_pts_k = list(points_tree.query(self.bound_pts_xy))[1]

    def makeMask(self, mask_type = 'exterior bound'):

        if mask_type == 'cluster bound':

            maskM = np.zeros(self.X.shape)

            maskM[self.map_ij2k[self.map_cell2ecm][:,0], self.map_ij2k[self.map_cell2ecm][:,1]] =1
            maskM[self.map_ij2k[self.bound_pts_k][:,0], self.map_ij2k[self.bound_pts_k][:,1]] =-1

        elif mask_type == 'exterior bound':
            maskM = np.ones(self.X.shape)
            maskM[1,1:-1]= -1
            maskM[-2,1:-1] =-1
            maskM[1:-1,1]=-1
            maskM[1:-1,-2] =-1

            maskM[0,:]= 0
            maskM[-1,:] = 0
            maskM[:,0]= 0
            maskM[:,-1] =0

        return maskM

    def makeLaplacian(self,maskM):
        """
        Generate the discrete finite central difference 2D Laplacian operator based on:

        d2 Uij / dx2 + d2 Uij / dy2 = (Ui+1,j + Ui-1,j + Ui,j+1, Ui,j-1 - 4Uij)*1/(delta_x)**2

        To solve the Poisson equation:

        A*U = fxy

        by:

        U = Ainv*fxy

        For a domain with an irregular boundary which is *embedded* within a 2D Cartesian grid.

        Parameters
        -----------
        maskM               A numpy nd array of the world grid, where 0 stands for external point where U = 0,
                            1 for an internal point where U is calculated, and -1 for a point on a boundary.

        Returns
        --------
        Ainv                Creates a laplacian matrix solver for the specified geometry (including boundary conditions)
        bound_pts_k         A list of k-indices for the boundary points (specific to the maskM) for modifying source

        Notes
        -------
        ip and jp are related to the original matrix index scheme by ip = i +1 and j = jp + 1.

        """

        # initialize the laplacian operator matrix:
        A = np.zeros((self.grid_len,self.grid_len))
        bound_pts_k = []

        for k, (i, j) in enumerate(self.map_ij2k):

            maskVal = maskM[i, j]

            if maskVal == 1:  # interior point, do full algorithm...

                k_ip1_j = self.map_ij2k.tolist().index([i + 1,j])
                k_in1_j = self.map_ij2k.tolist().index([i-1,j])
                k_i_jp1 = self.map_ij2k.tolist().index([i,j+1])
                k_i_jn1 = self.map_ij2k.tolist().index([i,j-1])

                A[k, k_ip1_j] = 1
                A[k, k_in1_j] = 1
                A[k, k_i_jp1] = 1
                A[k, k_i_jn1] = 1
                A[k,k] = -4

            elif maskVal == -1:  # boundary point, treat as such
                A[k,k] = 1
                bound_pts_k.append(k)

            elif maskVal == 0:  # exterior point, treat as such
                A[k,k] = 1

        # complete the laplacian operator by diving through by grid spacing squared:
        Ai = A/self.delta**2

        # calculate the inverse, which is stored for solution calculation of Laplace and Poisson equations
        Ainv = np.linalg.pinv(Ai)

        return Ainv

    def redo_gj(self,dyna,p,savecells =True):

        # profile_names = list(p.profiles.keys())  # names of each tissue profile...
        profile_names = dyna.tissue_profile_names
        new_gj_nn_i = []
        new_gj_vects = []
        removal_flags = np.zeros(len(self.gj_i))

        for name in profile_names:

            cell_targets = dyna.cell_target_inds[name]   # get the cell target inds for this tissue
            insular_flag = p.profiles[name]['insular gj']

            # step through gj's and find cases where connection is split between cells in different tissues:
            if insular_flag is True:

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

        if savecells is True:

            # save the cell cluster
            logs.log_info('Saving the cell cluster... ')

            # celf = copy.deepcopy(self)
            datadump = [self,p]
            fh.saveSim(self.savedWorld,datadump)
            message = 'Cell cluster saved to' + ' ' + self.savedWorld
            logs.log_info(message)

    # def makeLaplacian(self):
    #     """
    #     Generate the discrete finite central difference 2D Laplacian operator based on:
    #
    #     d2 Uij / dx2 + d2 Uij / dy2 = (Ui+1,j + Ui-1,j + Ui,j+1, Ui,j-1 - 4Uij)*1/(delta_x)**2
    #
    #     To solve the Poisson equation:
    #
    #     A*U = fxy
    #
    #     by:
    #
    #     U = Ainv*fxy
    #
    #     Creates
    #     --------
    #     self.Ainv
    #
    #     Notes
    #     -------
    #     ip and jp are related to the original matrix index scheme by ip = i +1 and j = jp + 1.
    #
    #
    #     """
    #     # initialize the laplacian operator matrix:
    #     A = np.zeros((self.core_len,self.core_len))
    #
    #     # # create a temporary Python list version of the full index map
    #     # map_ij2k = self.map_ij2k.tolist()
    #
    #     # lists to store the indices of
    #     self.bBot_kp = []
    #     self.bTop_kp = []
    #     self.bL_kp = []
    #     self.bR_kp = []
    #
    #
    #     for kp, (ip, jp) in enumerate(self.map_ij2k_core):
    #
    #         if ip + 1 < self.core_n:
    #
    #             k_ip1_j = self.map_ij2k_core.index([ip+1,jp])
    #             A[kp, k_ip1_j] = 1
    #
    #         else:
    #             # store the kp index of the bottom boundary for off the cuff modification of sol vector:
    #             self.bTop_kp.append(kp)
    #
    #             # deal with the bounary value by moving to the RHS of the equation:
    #             # bval = self.bBot[jp+1]
    #             # self.FF[k] = self.FF[k] - (bval/self.delta**2)
    #
    #         if ip - 1 >= 0:
    #
    #             k_in1_j = self.map_ij2k_core.index([ip-1,jp])
    #             A[kp, k_in1_j] = 1
    #
    #         else:
    #             # deal with the bounary value by moving to the RHS of the equation:
    #             self.bBot_kp.append(kp)
    #             # bval = self.bTop[jp+1]
    #             # self.FF[k] = self.FF[k] - (bval/self.delta**2)
    #
    #         if jp + 1 < self.core_n:
    #
    #             k_i_jp1 = self.map_ij2k_core.index([ip,jp+1])
    #             A[kp, k_i_jp1] = 1
    #
    #         else:
    #             # deal with the bounary value by moving to the RHS of the equation:
    #             self.bL_kp.append(kp)
    #
    #             # bval = self.bR[ip + 1]
    #             # self.FF[k] = self.FF[k] - (bval/self.delta**2)
    #
    #         if jp -1 >= 0:
    #             k_i_jn1 = self.map_ij2k_core.index([ip,jp-1])
    #             A[kp, k_i_jn1] = 1
    #
    #         else:
    #             # deal with the bounary value by moving to the RHS of the equation:
    #             self.bR_kp.append(kp)
    #
    #             # bval = self.bL[ip+1]
    #             # self.FF[k] = self.FF[k] - (bval/self.delta**2)
    #
    #         A[kp,kp] = -4
    #
    #     # complete the laplacian operator by diving through by grid spacing squared:
    #
    #     A = A/self.delta**2
    #
    #     loggers.log_info('Creating Poisson solver... ')
    #
    #     # calculate the inverse, which is stored for solution calculation of Laplace and Poisson equations
    #     self.Ainv = np.linalg.pinv(A)
