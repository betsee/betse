#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.


import numpy as np
import math
import scipy.ndimage
import scipy.spatial as sps


class FiniteDiffSolver(object):
    """
    Provides methods to establish a rectangular Marker and Cell (MACs) grid, compute discrete Laplacian operators
    and solve the Poisson, Laplacian, Heat, or Unsteady Stokes Flow (linearized Navier-Stokes) equations in 2D using
    the Central Finite Difference method.

    """

    def __init__(self):

        pass

    def makeGrid(self,grid_n, xmin = -1, xmax = 1, ymin = -1, ymax = 1):

        """
        Create linear vectors, rectangular grids and expanded grids of spatial coordinates
        and solution vectors.This method produces a standard grid.

        delta:          spacing between grid points (same in the x and y dimensions)
        xmin:           minimum x-value  [m]
        xmax:           max x-value [m]
        ymin:           minimum y-value  [m]
        ymax:           max y value [m]

        """

        # Create linear vectors and rectangular grids of spatial coordinates (indexed by i and j)
        self.x = np.linspace(xmin,xmax,grid_n)
        self.y = np.linspace(ymin,ymax,grid_n)

        self.X,self.Y = np.meshgrid(self.x,self.y)

        self.delta_x = self.x[1] - self.x[0]
        self.delta_y = self.y[1] - self.y[0]

        # define a data structure that holds [x,y] coordinate points of each 2d grid-matrix entry
        self.xypts = []
        self.map_ij2k = []

        for i, row in enumerate(self.X):
            for j, val in enumerate(row):
                xpt = self.X[i,j]
                ypt = self.Y[i,j]
                self.xypts.append([xpt,ypt])
                self.map_ij2k.append([i,j])

        self.map_ij2k = np.asarray(self.map_ij2k)

    def cell_grid(self,grid_delta, xmin, xmax, ymin, ymax):

        """
        Make a staggered n x n grid for use with finite volume or MACs method.
        One n x n grid contains the cell centres.
        One n+1 x n+1 grid contains the cell corners (vertices)
        One n+1 x n grid contains the u-coordinate points for velocity and x-coordinate for fluxes
        One n x n+1 grid contains the v-coordinate points for velocity and y-coordinate for fluxes

        """

        # max-min dimensions of the 2D world space:
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

        # spacing of the grid points (uniform in x and y dimensions):
        self.delta = grid_delta

        # number of points in the x and y directions
        self.grid_nx = int((self.xmax - self.xmin)/grid_delta)
        self.grid_ny = int((self.ymax - self.ymin)/grid_delta)

        # creation of the first mesh -- these are the corners of each MACs grid square:
        xv = np.linspace(xmin,xmax,self.grid_nx+1)
        yv = np.linspace(ymin,ymax,self.grid_ny+1)
        self.verts_X, self.verts_Y = np.meshgrid(xv,yv)

        x_verts = self.verts_X.ravel()
        y_verts = self.verts_Y.ravel()
        self.xy_verts = np.column_stack((x_verts,y_verts))
        self.verts_shape = self.verts_X.shape

        # grid and coordinate vectors defining the centre points of each cell:
        self.cents_X = (self.verts_X[0:-1,0:-1] + self.verts_X[0:-1,1:])/2
        self.cents_Y = (self.verts_Y[0:-1,0:-1] + self.verts_Y[1:,0:-1])/2

        x_cent = self.cents_X.ravel()
        y_cent = self.cents_Y.ravel()

        # unravelled x,y points for cell centres:
        self.xy_cents = np.column_stack((x_cent,y_cent))

        self.cents_shape = self.cents_X.shape

        # define a mapping between the (i,j) centers mesh and the linear unravelled k index:
        self.map_ij2k_cents = []

        for i, row in enumerate(self.cents_X):
            for j, val in enumerate(row):
                self.map_ij2k_cents.append([i,j])

        self.map_ij2k_cents = np.asarray(self.map_ij2k_cents)

        # Next define grids for the x and y cordinates as midpoints of each MACs cell:
        # cell side midpoints in the x (u) and x (v) directions:
        mx = self.cents_X[0,:]
        my = self.verts_Y[:,0]
        self.v_X, self.v_Y = np.meshgrid(mx,my)

        x_m = self.v_X.ravel()
        y_m = self.v_Y.ravel()

        # points on which y co-ordinate of velocity is defined:
        self.v_pts = np.column_stack((x_m,y_m))
        self.v_shape = self.v_X.shape

        self.map_ij2k_v = []

        for i, row in enumerate(self.v_X):
            for j, val in enumerate(row):
                self.map_ij2k_v.append([i,j])

        self.map_ij2k_v = np.asarray(self.map_ij2k_v)

        nx = self.verts_X[0,:]
        ny = self.cents_Y[:,0]
        self.u_X, self.u_Y = np.meshgrid(nx,ny)

        x_n = self.u_X.ravel()
        y_n = self.u_Y.ravel()

        # points on which x co-ordinate of velocity is defined:
        self.u_pts = np.column_stack((x_n,y_n))
        self.u_shape = self.u_X.shape

        self.map_ij2k_u = []

        for i, row in enumerate(self.u_X):
            for j, val in enumerate(row):
                self.map_ij2k_u.append([i,j])

        self.map_ij2k_u = np.asarray(self.map_ij2k_u)

    def makeLaplacian(self, bound = {'N':'value','S':'value','E':'value','W':'value'}):
        """
        Calculate a Laplacian operator matrix suitable for solving a 2D Poisson equation
        on a regular Cartesian grid with square boundaries. Note: the graph must have
        equal spacing in the x and y directions (the same delta in x and y directions).

        """

        size_rows = self.cents_shape[0]
        size_cols = self.cents_shape[1]

        sze = size_rows*size_cols

        A = np.zeros((sze,sze))

        for k, (i,j) in enumerate(self.map_ij2k_cents):

            # if we're not on a main boundary:
            if i != 0 and j != 0 and i != size_rows-1 and j != size_cols-1:

                k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1,j])
                k_in1_j = self.map_ij2k_cents.tolist().index([i-1,j])
                k_i_jp1 = self.map_ij2k_cents.tolist().index([i,j+1])
                k_i_jn1 = self.map_ij2k_cents.tolist().index([i,j-1])

                A[k, k_ip1_j] = 1
                A[k, k_in1_j] = 1
                A[k, k_i_jp1] = 1
                A[k, k_i_jn1] = 1
                A[k,k] = -4


            elif i == 0 and j != 0 and j != size_cols -1: # if on the bottom (South) boundary:

                if bound['S'] == 'flux':

                    k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1,j])
                    k_i_jp1 = self.map_ij2k_cents.tolist().index([i,j+1])
                    k_i_jn1 = self.map_ij2k_cents.tolist().index([i,j-1])

                    A[k, k_ip1_j] = 1
                    A[k, k_i_jp1] = 1
                    A[k, k_i_jn1] = 1

                    A[k,k] = -3

                elif bound['S'] == 'value':

                    k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1,j])

                    A[k,k] = 1
                    # A[k,k_ip1_j] = 1

            elif i == size_rows -1 and j != 0 and j != size_cols -1: # if on the top (North) boundary:

                if bound['N'] == 'flux':

                    k_in1_j = self.map_ij2k_cents.tolist().index([i-1,j])
                    k_i_jp1 = self.map_ij2k_cents.tolist().index([i,j+1])
                    k_i_jn1 = self.map_ij2k_cents.tolist().index([i,j-1])

                    A[k, k_in1_j] = 1
                    A[k, k_i_jp1] = 1
                    A[k, k_i_jn1] = 1

                    A[k,k] = -3

                elif bound['N'] == 'value':

                    k_in1_j = self.map_ij2k_cents.tolist().index([i-1,j])

                    A[k,k] = 1
                    # A[k,k_in1_j] = 1

            elif j == 0 and i != 0 and i != size_rows -1:  # if on the left (West) boundary:

                if bound['W'] == 'flux':

                    k_i_jp1 = self.map_ij2k_cents.tolist().index([i,j+1])
                    k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1,j])
                    k_in1_j = self.map_ij2k_cents.tolist().index([i-1,j])


                    A[k, k_i_jp1] = 1
                    A[k, k_ip1_j] = 1
                    A[k, k_in1_j] = 1

                    A[k,k] = -3

                elif bound['W'] == 'value':

                    k_i_jp1 = self.map_ij2k_cents.tolist().index([i,j+1])

                    A[k,k] = 1
                    # A[k,k_i_jp1] = 1

            elif j == size_cols -1 and i != 0 and i != size_rows -1: # if on the right (East) boundary:

                if bound['E'] == 'flux':

                    k_i_jn1 = self.map_ij2k_cents.tolist().index([i,j-1])
                    k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1,j])
                    k_in1_j = self.map_ij2k_cents.tolist().index([i-1,j])

                    A[k, k_i_jn1] = 1
                    A[k, k_ip1_j] = 1
                    A[k, k_in1_j] = 1

                    A[k,k] = -3

                elif bound['E'] == 'value':

                    k_i_jn1 = self.map_ij2k_cents.tolist().index([i,j-1])

                    A[k,k] = 1
                    # A[k,k_i_jn1] = 1

            # corners:
            elif i == 0 and j == 0: # SW corner

                if bound['S'] == 'flux':

                    k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1,j])
                    k_i_jp1 = self.map_ij2k_cents.tolist().index([i,j+1])

                    A[k, k_i_jp1] = 1
                    A[k, k_ip1_j] = 1
                    A[k,k] = -2

                else:
                    A[k,k] = 1

            elif i == size_rows - 1 and j == 0: # NW corner

                if bound['N'] == 'flux':

                    k_in1_j = self.map_ij2k_cents.tolist().index([i - 1,j])
                    k_i_jp1 = self.map_ij2k_cents.tolist().index([i,j+1])

                    A[k, k_i_jp1] = 1
                    A[k, k_in1_j] = 1
                    A[k,k] = -2

                else:
                    A[k,k] = 1

            elif i == 0 and j == size_cols - 1: # SE corner

                if bound['E'] == 'flux':

                    k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1,j])
                    k_i_jn1 = self.map_ij2k_cents.tolist().index([i,j-1])

                    A[k, k_i_jn1] = 1
                    A[k, k_ip1_j] = 1
                    A[k,k] = -2

                else:
                    A[k,k] = 1

            elif i == size_rows - 1 and j == size_cols - 1: # NE corner

                if bound['E'] == 'flux':

                    k_in1_j = self.map_ij2k_cents.tolist().index([i - 1,j])
                    k_i_jn1 = self.map_ij2k_cents.tolist().index([i,j-1])

                    A[k, k_i_jn1] = 1
                    A[k, k_in1_j] = 1
                    A[k,k] = -2

                else:
                    A[k,k] = 1


        A = A/(self.delta**2)

        # calculate the inverse, which is stored for solution calculation of Laplace and Poisson equations
        Ainv = np.linalg.pinv(A)

        return A, Ainv

    def makeScreenedLaplacian(self, ko=1.0e9, bound = {'N':'value','S':'value','E':'value','W':'value'}):
        """
        Calculate a Laplacian operator matrix suitable for solving a 2D Poisson equation
        on a regular Cartesian grid with square boundaries. Note: the graph must have
        equal spacing in the x and y directions (the same delta in x and y directions).

        """

        size_rows = self.cents_shape[0]
        size_cols = self.cents_shape[1]

        sze = size_rows*size_cols

        A = np.zeros((sze,sze))

        for k, (i,j) in enumerate(self.map_ij2k_cents):

            # if we're not on a main boundary:
            if i != 0 and j != 0 and i != size_rows-1 and j != size_cols-1:

                k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1,j])
                k_in1_j = self.map_ij2k_cents.tolist().index([i-1,j])
                k_i_jp1 = self.map_ij2k_cents.tolist().index([i,j+1])
                k_i_jn1 = self.map_ij2k_cents.tolist().index([i,j-1])

                A[k, k_ip1_j] = 1
                A[k, k_in1_j] = 1
                A[k, k_i_jp1] = 1
                A[k, k_i_jn1] = 1
                A[k,k] = -4 - (self.delta**2)*(ko**2)


            elif i == 0 and j != 0 and j != size_cols -1: # if on the bottom (South) boundary:

                if bound['S'] == 'flux':

                    k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1,j])
                    k_i_jp1 = self.map_ij2k_cents.tolist().index([i,j+1])
                    k_i_jn1 = self.map_ij2k_cents.tolist().index([i,j-1])

                    A[k, k_ip1_j] = 1
                    A[k, k_i_jp1] = 1
                    A[k, k_i_jn1] = 1

                    A[k,k] = -3 - (self.delta**2)*(ko**2)

                elif bound['S'] == 'value':

                    k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1,j])

                    A[k,k] = 1
                    # A[k,k_ip1_j] = 1

            elif i == size_rows -1 and j != 0 and j != size_cols -1: # if on the top (North) boundary:

                if bound['N'] == 'flux':

                    k_in1_j = self.map_ij2k_cents.tolist().index([i-1,j])
                    k_i_jp1 = self.map_ij2k_cents.tolist().index([i,j+1])
                    k_i_jn1 = self.map_ij2k_cents.tolist().index([i,j-1])

                    A[k, k_in1_j] = 1
                    A[k, k_i_jp1] = 1
                    A[k, k_i_jn1] = 1

                    A[k,k] = -3 - (self.delta**2)*(ko**2)

                elif bound['N'] == 'value':

                    k_in1_j = self.map_ij2k_cents.tolist().index([i-1,j])

                    A[k,k] = 1
                    # A[k,k_in1_j] = 1

            elif j == 0 and i != 0 and i != size_rows -1:  # if on the left (West) boundary:

                if bound['W'] == 'flux':

                    k_i_jp1 = self.map_ij2k_cents.tolist().index([i,j+1])
                    k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1,j])
                    k_in1_j = self.map_ij2k_cents.tolist().index([i-1,j])


                    A[k, k_i_jp1] = 1
                    A[k, k_ip1_j] = 1
                    A[k, k_in1_j] = 1

                    A[k,k] = -3 - (self.delta**2)*(ko**2)

                elif bound['W'] == 'value':

                    k_i_jp1 = self.map_ij2k_cents.tolist().index([i,j+1])

                    A[k,k] = 1
                    # A[k,k_i_jp1] = 1

            elif j == size_cols -1 and i != 0 and i != size_rows -1: # if on the right (East) boundary:

                if bound['E'] == 'flux':

                    k_i_jn1 = self.map_ij2k_cents.tolist().index([i,j-1])
                    k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1,j])
                    k_in1_j = self.map_ij2k_cents.tolist().index([i-1,j])

                    A[k, k_i_jn1] = 1
                    A[k, k_ip1_j] = 1
                    A[k, k_in1_j] = 1

                    A[k,k] = -3 - (self.delta**2)*(ko**2)

                elif bound['E'] == 'value':

                    k_i_jn1 = self.map_ij2k_cents.tolist().index([i,j-1])

                    A[k,k] = 1
                    # A[k,k_i_jn1] = 1

            # corners:
            elif i == 0 and j == 0: # SW corner

                if bound['S'] == 'flux':

                    k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1,j])
                    k_i_jp1 = self.map_ij2k_cents.tolist().index([i,j+1])

                    A[k, k_i_jp1] = 1
                    A[k, k_ip1_j] = 1
                    A[k,k] = -2 - (self.delta**2)*(ko**2)

                else:
                    A[k,k] = 1

            elif i == size_rows - 1 and j == 0: # NW corner

                if bound['N'] == 'flux':

                    k_in1_j = self.map_ij2k_cents.tolist().index([i - 1,j])
                    k_i_jp1 = self.map_ij2k_cents.tolist().index([i,j+1])

                    A[k, k_i_jp1] = 1
                    A[k, k_in1_j] = 1
                    A[k,k] = -2 - (self.delta**2)*(ko**2)

                else:
                    A[k,k] = 1

            elif i == 0 and j == size_cols - 1: # SE corner

                if bound['E'] == 'flux':

                    k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1,j])
                    k_i_jn1 = self.map_ij2k_cents.tolist().index([i,j-1])

                    A[k, k_i_jn1] = 1
                    A[k, k_ip1_j] = 1
                    A[k,k] = -2 - (self.delta**2)*(ko**2)

                else:
                    A[k,k] = 1

            elif i == size_rows - 1 and j == size_cols - 1: # NE corner

                if bound['E'] == 'flux':

                    k_in1_j = self.map_ij2k_cents.tolist().index([i - 1,j])
                    k_i_jn1 = self.map_ij2k_cents.tolist().index([i,j-1])

                    A[k, k_i_jn1] = 1
                    A[k, k_in1_j] = 1
                    A[k,k] = -2 - (self.delta**2)*(ko**2)

                else:
                    A[k,k] = 1


        A = A/(self.delta**2)

        # calculate the inverse, which is stored for solution calculation of Laplace and Poisson equations
        Ainv = np.linalg.pinv(A)

        return A, Ainv

    def makeLaplacian_u(self, bound = {'N':'value','S':'value','E':'value','W':'value'}):
        """
        Calculate a Laplacian operator matrix suitable for solving a 2D Poisson equation
        on a regular Cartesian grid with square boundaries. Note: the graph must have
        equal spacing in the x and y directions (the same delta in x and y directions).

        """

        size_rows = self.u_shape[0]
        size_cols = self.u_shape[1]

        sze = size_rows*size_cols

        A = np.zeros((sze,sze))

        for k, (i,j) in enumerate(self.map_ij2k_u):

            # if we're not on a main boundary:
            if i != 0 and j != 0 and i != size_rows-1 and j != size_cols-1:

                k_ip1_j = self.map_ij2k_u.tolist().index([i + 1,j])
                k_in1_j = self.map_ij2k_u.tolist().index([i-1,j])
                k_i_jp1 = self.map_ij2k_u.tolist().index([i,j+1])
                k_i_jn1 = self.map_ij2k_u.tolist().index([i,j-1])

                A[k, k_ip1_j] = 1
                A[k, k_in1_j] = 1
                A[k, k_i_jp1] = 1
                A[k, k_i_jn1] = 1
                A[k,k] = -4


            if i == 0: # if on the bottom (South) boundary:

                if bound['S'] == 'flux':

                    k_ip1_j = self.map_ij2k_u.tolist().index([i + 1,j])
                    A[k, k_ip1_j] = 1

                    A[k,k] = -1

                elif bound['S'] == 'value':

                    k_ip1_j = self.map_ij2k_u.tolist().index([i + 1,j])

                    A[k,k] = 1
                    A[k,k_ip1_j] = 1

            if i == size_rows -1: # if on the top (North) boundary:

                if bound['N'] == 'flux':

                    k_in1_j = self.map_ij2k_u.tolist().index([i-1,j])
                    A[k, k_in1_j] = 1

                    A[k,k] = -1

                elif bound['N'] == 'value':

                    k_in1_j = self.map_ij2k_u.tolist().index([i-1,j])

                    A[k,k] = 1
                    A[k,k_in1_j] = 1

            if j == 0: # if on the left (West) boundary:

                if bound['W'] == 'flux':

                    k_i_jp1 = self.map_ij2k_u.tolist().index([i,j+1])
                    A[k, k_i_jp1] = 1

                    A[k,k] = -1

                elif bound['W'] == 'value':

                    k_i_jp1 = self.map_ij2k_u.tolist().index([i,j+1])

                    A[k,k] = 1
                    A[k,k_i_jp1] = 1

            if j == size_cols -1: # if on the right (East) boundary:

                if bound['E'] == 'flux':

                    k_i_jn1 = self.map_ij2k_u.tolist().index([i,j-1])
                    A[k, k_i_jn1] = 1

                    A[k,k] = -1

                elif bound['E'] == 'value':

                    k_i_jn1 = self.map_ij2k_u.tolist().index([i,j-1])

                    A[k,k] = 1
                    A[k,k_i_jn1] = 1

        A = A/(self.delta**2)

        # calculate the inverse, which is stored for solution calculation of Laplace and Poisson equations
        Ainv = np.linalg.pinv(A)

        return A, Ainv

    def makeLaplacian_v(self, bound = {'N':'value','S':'value','E':'value','W':'value'}):
        """
        Calculate a Laplacian operator matrix suitable for solving a 2D Poisson equation
        on a regular Cartesian grid with square boundaries. Note: the graph must have
        equal spacing in the x and y directions (the same delta in x and y directions).

        """

        size_rows = self.v_shape[0]
        size_cols = self.v_shape[1]

        sze = size_rows*size_cols

        A = np.zeros((sze,sze))

        for k, (i,j) in enumerate(self.map_ij2k_v):

            # if we're not on a main boundary:
            if i != 0 and j != 0 and i != size_rows-1 and j != size_cols-1:

                k_ip1_j = self.map_ij2k_v.tolist().index([i + 1,j])
                k_in1_j = self.map_ij2k_v.tolist().index([i-1,j])
                k_i_jp1 = self.map_ij2k_v.tolist().index([i,j+1])
                k_i_jn1 = self.map_ij2k_v.tolist().index([i,j-1])

                A[k, k_ip1_j] = 1
                A[k, k_in1_j] = 1
                A[k, k_i_jp1] = 1
                A[k, k_i_jn1] = 1
                A[k,k] = -4


            if i == 0: # if on the bottom (South) boundary:

                if bound['S'] == 'flux':

                    k_ip1_j = self.map_ij2k_v.tolist().index([i + 1,j])
                    A[k, k_ip1_j] = 1

                    A[k,k] = -1

                elif bound['S'] == 'value':

                    k_ip1_j = self.map_ij2k_v.tolist().index([i + 1,j])

                    A[k,k] = 1
                    A[k,k_ip1_j] = 1

            if i == size_rows -1: # if on the top (North) boundary:

                if bound['N'] == 'flux':

                    k_in1_j = self.map_ij2k_v.tolist().index([i-1,j])
                    A[k, k_in1_j] = 1

                    A[k,k] = -1

                elif bound['N'] == 'value':

                    k_in1_j = self.map_ij2k_v.tolist().index([i-1,j])

                    A[k,k] = 1
                    A[k,k_in1_j] = 1

            if j == 0: # if on the left (West) boundary:

                if bound['W'] == 'flux':

                    k_i_jp1 = self.map_ij2k_v.tolist().index([i,j+1])
                    A[k, k_i_jp1] = 1

                    A[k,k] = -1

                elif bound['W'] == 'value':

                    k_i_jp1 = self.map_ij2k_v.tolist().index([i,j+1])

                    A[k,k] = 1
                    A[k,k_i_jp1] = 1

            if j == size_cols -1: # if on the right (East) boundary:

                if bound['E'] == 'flux':

                    k_i_jn1 = self.map_ij2k_v.tolist().index([i,j-1])
                    A[k, k_i_jn1] = 1

                    A[k,k] = -1

                elif bound['E'] == 'value':

                    k_i_jn1 = self.map_ij2k_v.tolist().index([i,j-1])

                    A[k,k] = 1
                    A[k,k_i_jn1] = 1

        A = A/(self.delta**2)

        # calculate the inverse, which is stored for solution calculation of Laplace and Poisson equations
        Ainv = np.linalg.pinv(A)

        return A, Ainv

    def makeIntegrator(self):
        """
        Calculate a Finite Difference integration operator matrix for the
        2D grid. Stencil takes 1/2 of central point and 1/8 of points interpolated
        to the side of each imaginary grid square.

        """

        size_rows = self.cents_shape[0]
        size_cols = self.cents_shape[1]

        sze = size_rows*size_cols

        A = np.zeros((sze,sze))

        for k, (i, j) in enumerate(self.map_ij2k_cents):

            # if we're not on a main boundary:
            if i != 0 and j != 0 and i != size_rows - 1 and j != size_cols - 1:

                k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1, j])
                k_in1_j = self.map_ij2k_cents.tolist().index([i - 1, j])
                k_i_jp1 = self.map_ij2k_cents.tolist().index([i, j + 1])
                k_i_jn1 = self.map_ij2k_cents.tolist().index([i, j - 1])

                A[k, k_ip1_j] = 1 / 16
                A[k, k_in1_j] = 1 / 16
                A[k, k_i_jp1] = 1 / 16
                A[k, k_i_jn1] = 1 / 16
                A[k, k] = (1 / 2) + (1 / 4)


            elif i == 0 and j != 0 and j != size_cols - 1:  # if on the bottom (South) boundary:

                k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1, j])
                k_i_jp1 = self.map_ij2k_cents.tolist().index([i, j + 1])
                k_i_jn1 = self.map_ij2k_cents.tolist().index([i, j - 1])

                A[k, k_ip1_j] = 1 / 16
                A[k, k_i_jp1] = 1 / 16
                A[k, k_i_jn1] = 1 / 16

                A[k, k] = (1 / 2) + (3 / 16)


            elif i == size_rows - 1 and j != 0 and j != size_cols - 1:  # if on the top (North) boundary:

                k_in1_j = self.map_ij2k_cents.tolist().index([i - 1, j])
                k_i_jp1 = self.map_ij2k_cents.tolist().index([i, j + 1])
                k_i_jn1 = self.map_ij2k_cents.tolist().index([i, j - 1])

                A[k, k_in1_j] = 1 / 16
                A[k, k_i_jp1] = 1 / 16
                A[k, k_i_jn1] = 1 / 16

                A[k, k] = (1 / 2) + (3 / 16)

            elif j == 0 and i != 0 and i != size_rows - 1:  # if on the left (West) boundary:

                k_i_jp1 = self.map_ij2k_cents.tolist().index([i, j + 1])
                k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1, j])
                k_in1_j = self.map_ij2k_cents.tolist().index([i - 1, j])

                A[k, k_i_jp1] = 1 / 16
                A[k, k_ip1_j] = 1 / 16
                A[k, k_in1_j] = 1 / 16

                A[k, k] = (1 / 2) + (3 / 16)

            elif j == size_cols - 1 and i != 0 and i != size_rows - 1:  # if on the right (East) boundary:

                k_i_jn1 = self.map_ij2k_cents.tolist().index([i, j - 1])
                k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1, j])
                k_in1_j = self.map_ij2k_cents.tolist().index([i - 1, j])

                A[k, k_i_jn1] = 1 / 16
                A[k, k_ip1_j] = 1 / 16
                A[k, k_in1_j] = 1 / 16

                A[k, k] = (1 / 2) + (3 / 16)

            # corners:
            elif i == 0 and j == 0:  # SW corner

                k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1, j])
                k_i_jp1 = self.map_ij2k_cents.tolist().index([i, j + 1])

                A[k, k_i_jp1] = 1 / 16
                A[k, k_ip1_j] = 1 / 16
                A[k, k] = (1 / 2) + (1 / 8)

            elif i == size_rows - 1 and j == 0:  # NW corner

                k_in1_j = self.map_ij2k_cents.tolist().index([i - 1, j])
                k_i_jp1 = self.map_ij2k_cents.tolist().index([i, j + 1])

                A[k, k_i_jp1] = 1 / 16
                A[k, k_in1_j] = 1 / 16
                A[k, k] = A[k, k] = (1 / 2) + (1 / 8)

            elif i == 0 and j == size_cols - 1:  # SE corner

                k_ip1_j = self.map_ij2k_cents.tolist().index([i + 1, j])
                k_i_jn1 = self.map_ij2k_cents.tolist().index([i, j - 1])

                A[k, k_i_jn1] = 1 / 16
                A[k, k_ip1_j] = 1 / 16
                A[k, k] = (1 / 2) + (1 / 8)

            elif i == size_rows - 1 and j == size_cols - 1:  # NE corner

                k_in1_j = self.map_ij2k_cents.tolist().index([i - 1, j])
                k_i_jn1 = self.map_ij2k_cents.tolist().index([i, j - 1])

                A[k, k_i_jn1] = 1 / 16
                A[k, k_in1_j] = 1 / 16
                A[k, k] = (1 / 2) + (1 / 8)

        return A

    def stokes_kernel(self):
        """
        Calculate the linearized Navier-Stokes equations using
        the Finite Difference method on a Marker and Cell (MACs) grid.

        """
        # parameters of the liquid
        rho = 1e3   # density
        visc = 0.1  # visocity

        # body force on the liquid, with components at the x and y coordinates of velocity:
        Fy = np.zeros(self.v_shape)
        Fy[5:-5,5:-5]=-9.81*rho

        Fx = np.zeros(self.u_shape)

        Ainv = self.makeLaplacian()

        time_step = 1e-5
        end_time = 1e-3
        time_points = end_time/time_step
        time = np.linspace(0,end_time,time_points)

        # Unsteady Stokes Flow Solver on MACs grids

        # Initial conditions:

        # Velocity u= ui + vj
        u = np.zeros(self.u_shape)
        v = np.zeros(self.v_shape)

        P_time = []
        u_time = []
        v_time = []

        for t in time:

            # reinforce boundary conditions
            #left
            u[:,0] = 0
            # right
            u[:,-1] = 0
            # top
            u[-1,:] = 0
            # bottom
            u[0,:] = 0

            # left
            v[:,0] = 0
            # right
            v[:,-1] = 0
            # top
            v[-1,:] = 0
            # bottom
            v[0,:] = 0

            # calculate the flow, omitting the pressure term:
            lap_u = laplacian(u,self.delta_u_x,self.delta_u_y)
            lap_v = laplacian(v,self.delta_v_x,self.delta_v_y)

            u = u + (time_step/rho)*(visc*lap_u + Fx)

            v = v + (time_step/rho)*(visc*lap_v + Fy)

            # take the divergence of the flow field using a forward difference
            # that creates a matrix the same size as the pressure matrix:

            u_dx = (u[:,1:] - u[:,0:-1])/self.delta_u_x
            v_dy = (v[1:,:] - v[0:-1,:])/self.delta_v_y

            div_u = u_dx + v_dy

            #...solve for the pressure in terms of existing divergence:

            source = (rho/time_step)*div_u.ravel()
    #         source[bvals] = 0

            P = np.dot(Ainv, source)
            P = P.reshape(self.cents_shape)

            # enforce zero gradient boundary conditions on P:
            P[:,0] = P[:,1]
            P[:,-1] = P[:,-2]
            P[0,:] = P[1,:]
            P[-1,:] = P[-2,:]

            # Take the gradient of the pressue:
            gPxo, gPyo = gradient(P,self.delta_cents_x,self.delta_cents_y)

            gPx = np.zeros(self.u_shape)
            gPx[:,0:-1] = gPxo
            gPx[:,-1] = gPxo[:,-1]

            gPy = np.zeros(self.v_shape)
            gPy[0:-1,:] = gPyo
            gPy[-1,:] = gPyo[-1,:]

            # subtract the pressure from the solution to yeild a divergence-free flow field
            u = u - gPx*(time_step/rho)
            v = v - gPy*(time_step/rho)

            # interpolate u and v values at the centre for easy plotting:
            u_at_c = u[:,0:-1]

            v_at_c = v[0:-1,:]

            # reinforce boundary conditions
            #left
            u_at_c[:,0] = 0
            # right
            u_at_c[:,-1] = 0
            # top
            u_at_c[-1,:] = 0
            # bottom
            u_at_c[0,:] = 0

            # left
            v_at_c[:,0] = 0
            # right
            v_at_c[:,-1] = 0
            # top
            v_at_c[-1,:] = 0
            # bottom
            v_at_c[0,:] = 0

            P_time.append(P)
            u_time.append(u_at_c)
            v_time.append(v_at_c)

    def grid_gradient(self,P,bounds='closed'):
        """
        Calculates a gradient for the MACs grid for a property P
        defined on MACs cell centres.

        Allows boundary conditions to be set: 'closed' means no flux out of boundary
        'open' means an exchange of flux can occur.

        """

        # gradient for ecm:
        gPx = np.zeros(self.u_shape)
        gPy = np.zeros(self.v_shape)

        gPxo = (P[:,1:] - P[:,0:-1])/(2*self.delta)
        gPyo = (P[1:,:] - P[0:-1,:])/(2*self.delta)

        gPx[:,1:-1] = gPxo
        gPy[1:-1,:] = gPyo

        # decide on boundary conditions:
        if bounds == 'open':
            # no "acceleration" on any boundary:
            gPx[:,0] = gPx[:,1]
            gPx[:,-1] = gPx[:,-2]
            gPx[0,:] = gPx[1,:]
            gPx[-1,:] = gPx[-2,:]

            gPy[0,:] = gPy[1,:]
            gPy[-1,:] = gPy[-2,:]
            gPy[:,0] = gPy[:,1]
            gPy[:,-1] = gPy[:,-2]

        elif bounds == 'closed':
            # no flux:
            gPx[:,0] = 0
            gPx[:,-1] = 0
            gPy[0,:] = 0
            gPy[-1,:] = 0

        elif bounds == 'none':
            pass

        return gPx, gPy

    def grid_int(self,F,bounds='closed'):

        # interpolate F to x and y midpoints:
        Fx = np.zeros(self.u_shape)
        Fy = np.zeros(self.v_shape)

        Fxo = (F[:,1:] + F[:,0:-1])/2
        Fyo = (F[1:,:] + F[0:-1,:])/2

        Fx[:,1:-1] = Fxo
        Fy[1:-1,:] = Fyo

        # set boundary conditions:
        if bounds == 'open':
            # zero "gradient" on any boundary:
            Fx[:,0] = Fx[:,1]
            Fx[:,-1] = Fx[:,-2]
            # Fx[0,:] = Fx[1,:]
            # Fx[-1,:] = Fx[-2,:]

            Fy[0,:] = Fy[1,:]
            Fy[-1,:] = Fy[-2,:]
            # Fy[:,0] = Fy[:,1]
            # Fy[:,-1] = Fy[:,-2]

        elif bounds == 'closed':
            # no flux:
            Fx[:,0] = 0
            Fx[:,-1] = 0
            Fy[0,:] = 0
            Fy[-1,:] = 0

        elif bounds == 'none':
            pass

        F_int = np.zeros(F.shape)

        eP = Fx[:,1:] # east midpoints
        wP = Fx[:,0:-1] # west midpoints
        nP = Fy[1:,:] # north midpoints
        sP = Fy[0:-1,:] # south midpoints

        # do the numerical integration:
        F_int[:,:] = (1/2)*F
        F_int[:,:] = F_int[:,:] + (1/8)*nP
        F_int[:,:] = F_int[:,:] + (1/8)*sP
        F_int[:,:] = F_int[:,:] + (1/8)*eP
        F_int[:,:] = F_int[:,:] + (1/8)*wP

        return F_int

def jacobi(A,b,N=50,x=None):
    """
    Solves the equation Ax=b via the Jacobi iterative method.
    Courtesy of Michael Halls-Moore
    https://www.quantstart.com/articles/Jacobi-Method-in-Python-and-NumPy
    """
    # Create an initial guess if needed
    if x is None:
        x = np.zeros(len(A[0]))

    # Create a vector of the diagonal elements of A
    # and subtract them from A
    D = np.diag(A)
    R = A - np.diagflat(D)

    # Iterate for N times
    for i in range(N):
        xo = x
        x = (b - np.dot(R,x)) / D
        er = np.sum((x-xo)**2)/len(x)

    return x, er

def laplacian(F,delx,dely=None):
    """
    Uses discrete central finite difference to calculate the action of a Laplacian operator on a
    dataset F.

    Parameters
    ------------
    F:          A Numpy ndarray defined as F(x,y) on a grid
    dx:         Spacing of the grid (must be uniform in x and y directions)

    Returns
    -----------
    del2F        The Laplacian of F.

    """

    if dely is None:
        dely = delx

    ddFx_interior = (F[:,2:] - 2*F[:,1:-1] + F[:,0:-2])/(delx*dely)
    ddFy_interior = (F[2:,:] - 2*F[1:-1,:] + F[0:-2,:])/(delx*dely)

    ddF_B = (F[-3,:] - 2*F[-2,:] + F[-1,:])/(dely*dely)
    ddF_T = (F[0,:] - 2*F[-1,:] + F[-2,:])/(dely*dely)
    ddF_L = (F[:,2] - 2*F[:,1] + F[:,0])/(delx*delx)
    ddF_R = (F[:,-3] - 2*F[:,-2] + F[:,-1])/(delx*delx)

    lapF_interior = ddFx_interior[1:-1,:] + ddFy_interior[:,1:-1]

    # initialize the dFx and dFy arrays:
    ddF = np.zeros(F.shape)

    ddF[1:-1,1:-1] = lapF_interior

    ddF[:,0] = ddF_L
    ddF[:,-1] = ddF_R

    ddF[0,:] = ddF_T
    ddF[-1,:] = ddF_B


    return ddF

def gradient(F,delx,dely=None):
    # gradient using numpy slicing:

    if dely is None:
        dely = delx

    # calculate the discrete central first derivatives on the internal mesh points:
    dF_interior_y = -(F[:-2,:] - F[2:,:])/(2*dely)
    dF_interior_x = -(F[:,:-2] - F[:,2:])/(2*delx)

    # calculate the discrete forward or backward first derivatives on the boundary points:
    dF_B = (F[1,:] - F[0,:])/dely
    dF_T = (F[-1,:] - F[-2,:])/dely
    dF_L = (F[:,1] - F[:,0])/delx
    dF_R = (F[:,-1] - F[:,-2])/delx

    # initialize the dFx and dFy arrays:
    dFx = np.zeros(F.shape)
    dFy = np.zeros(F.shape)

    # build the final dFx and dFy arrays by splicing together internal and boundary derivatives:
    dFx[:,1:-1] = dF_interior_x
    dFy[1:-1,:] = dF_interior_y

    dFx[:,0] = dF_L
    dFx[:,-1] = dF_R

    dFy[0,:] = dF_B
    dFy[-1,:] = dF_T

    return dFx, dFy

def diff(F,delx,axis=0):
    # dertivative using numpy slicing:

    if axis == 1:
        # calculate the discrete central first derivatives on the internal mesh points:
        dF_interior = -(F[:-2,:] - F[2:,:])/(2*delx)

        # calculate the discrete forward or backward first derivatives on the boundary points:
        dF_B = -(F[1,:] - F[0,:])/delx
        dF_T = -(F[-1,:] - F[-2,:])/delx

        dF = np.zeros(F.shape)

        dF[1:-1,:] = dF_interior

        dF[0,:] = dF_B
        dF[-1,:] = dF_T


    elif axis == 0:
        # calculate the discrete central first derivatives on the internal mesh points:
        dF_interior = -(F[:,:-2] - F[:,2:])/(2*delx)

        # calculate the discrete forward or backward first derivatives on the boundary points:
        dF_L = (F[:,0] - F[:,1])/delx
        dF_R = (F[:,-2] - F[:,-1])/delx

        dF = np.zeros(F.shape)

        dF[:,1:-1] = dF_interior

        dF[:,0] = dF_L
        dF[:,-1] = dF_R


    return dF

def divergence(Fx,Fy,delx,dely):

    gx = diff(Fx, delx, axis=0)
    gy = diff(Fy, dely, axis=1)
    div = gx + gy

    return div

def curl(Fx, Fy, delx,dely):

    """
    Calculates the curl of a 2D vector field.
    The returned parameter is the z-axis component (in or out of the page).

    """

    gyx = diff(Fy,delx,axis=0)
    gxy = diff(Fx,dely,axis=1)

    curlF = gyx - gxy

    return curlF

def makeMask(M,xy_pts,X,Y,delta,sensitivity=1.0):
    """
    Creates an nd-array matrix mask with in-boundary
    elements tagged as 1, on-boundary elements tagged
    as -1 and environmental points tagged as 0.

    Parameters
    -----------
    M               A matrix with non-zero data and zero environment
    xy_pts          unravelled xy points for the Cells of M
    X               2D X matrix for the Cells of M
    Y               2D Y matrix for the Cells of M
    delta           Spacing (average x and y) for world of M grid
    sensitivity     1.0 is average, lower selects less boundary points, higher more

    Returns
    -----------
    maskM           The mask matrix
    """

    maskM = np.zeros(M.shape)
    inds = (M > 0).nonzero()
    maskM[inds] = 1

    x_bpts = X[inds].ravel()
    y_bpts =Y[inds].ravel()
    bpts = np.column_stack((x_bpts,y_bpts))

    bound_inds = boundTag(bpts,delta,alpha=sensitivity)
    points_tree = sps.KDTree(xy_pts)
    b_inds = points_tree.query(bpts[bound_inds])[1]

    mm = maskM.ravel()
    mm[b_inds] = -1
    maskM = mm.reshape(M.shape)

    return maskM

def boundTag(points,delta,alpha=1.0):

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
    alpha = alpha/delta

    tri = sps.Delaunay(points)
    tri_edges = []
    circum_r_list = []

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
        area = math.sqrt(abs(s*(s-a)*(s-b)*(s-c)))

        if area > 0:
            circum_r = a*b*c/(4.0*area)

        if area == 0:
            circum_r = a*b*c/(4.0*1e-25)

        # Here's the radius filter:
        if circum_r < 1.0/alpha:
            #self.tri_edges_append([ia, ib])
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

    i = 0
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

#     con_hull = tb.alpha_shape(points, alpha/delta)  # get the concave hull for the membrane midpoints
    concave_hull = np.asarray(concave_hull)

    bflags = np.unique(concave_hull)    # get the value of unique indices from segments

    return bflags

def integrator(P):

    F = np.zeros(P.shape)

    eP = P[:,1:] # east midpoints
    wP = P[:,0:-1] # west midpoints
    nP = P[1:,:] # north midpoints
    sP = P[0:-1,:] # south midpoints

    F[:,:] = (1/2)*P
    F[0:-1,:] = (1/16)*F[0:-1,:] + (1/16)*nP
    F[1:,:] = (1/16)*F[1:,:] + (1/16)*sP
    F[:,0:-1] = (1/16)*F[:,0:-1] + (1/16)*eP
    F[:,1:] = (1/16)*F[:,1:] + (1/16)*wP

    return F

# flux summer
def flux_summer(U,V,P):

    F = np.zeros(P.shape)

    eU = U[:,1::1] # east midpoints
    wU = U[:,0:-1:1] # west midpoints
    nV = V[1::1,:] # north midpoints
    sV = V[0:-1:1,:] # south midpoints

    # do the sum:
    F[:,:] = nV - sV + eU - wU

    return F

def flux_summer_sym(U,V):

    F = np.zeros(U.shape)

    eP = U[:,1:] # east midpoints
    wP = U[:,0:-1] # west midpoints
    nP = V[1:,:] # north midpoints
    sP = V[0:-1,:] # south midpoints

    F[0:-1,:] = nP
    F[1:,:] = F[1:,:] - sP
    F[:,0:-1] = F[:,0:-1] + eP
    F[:,1:] = F[:,1:] - wP

    return F


#-----------------------------------------------------------------------------------------------------------------------
#                                                        WASTELANDS
#-----------------------------------------------------------------------------------------------------------------------

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
    #     Modifies
    #     --------
    #     self.FF
    #
    #
    #     """
    #
    #     A = np.zeros((self.unraveled_length,self.unraveled_length))
    #
    #     for k, (ip, jp) in enumerate(self.map_ij2k):
    #
    #         if ip + 1 < self.core_ny:
    #
    #             k_ip1_j = self.map_ij2k.index([ip+1,jp])
    #             A[k, k_ip1_j] = 1
    #
    #         else:
    #             # deal with the bounary value by moving to the RHS of the equation:
    #             bval = self.bBot[jp+1]
    #             self.FF[k] = self.FF[k] - (bval/self.delta**2)
    #
    #         if ip - 1 >= 0:
    #
    #             k_in1_j = self.map_ij2k.index([ip-1,jp])
    #             A[k, k_in1_j] = 1
    #
    #         else:
    #             # deal with the bounary value by moving to the RHS of the equation:
    #             bval = self.bTop[jp+1]
    #             self.FF[k] = self.FF[k] - (bval/self.delta**2)
    #
    #         if jp + 1 < self.core_nx:
    #
    #             k_i_jp1 = self.map_ij2k.index([ip,jp+1])
    #             A[k, k_i_jp1] = 1
    #
    #         else:
    #             # deal with the bounary value by moving to the RHS of the equation:
    #             bval = self.bR[ip + 1]
    #             self.FF[k] = self.FF[k] - (bval/self.delta**2)
    #
    #         if jp -1 >= 0:
    #             k_i_jn1 = self.map_ij2k.index([ip,jp-1])
    #             A[k, k_i_jn1] = 1
    #
    #         else:
    #             # deal with the bounary value by moving to the RHS of the equation:
    #             bval = self.bL[ip+1]
    #             self.FF[k] = self.FF[k] - (bval/self.delta**2)
    #
    #         A[k,k] = -4
    #
    #     # complete the laplacian operator by diving through by grid spacing squared:
    #
    #     A = A/self.delta**2
    #
    #     # calculate the inverse, which is stored for solution calculation of Laplace and Poisson equations
    #     self.Ainv = np.linalg.pinv(A)

    # def makeLaplacian(grid_len,shape,map_ij2k,delx,dely,maskM=None):
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
#     For a domain with an irregular boundary which is *embedded* within a 2D Cartesian grid.
#
#     Parameters
#     -----------
#     grid_len            The length of the unravelled 2D array
#     map_ij2k            An array where each index k corresponding to unravelled array
#                         maps to the i,jth element of the 2D array
#     delx                Grid spacing in the x-dimension
#     dely                Grid spacing in the y-dimension
#     maskM               A numpy nd array of the world grid, where 0 stands for external point where U = 0,
#                         1 for an internal point where U is calculated, and -1 for a point on a boundary.
#
#     Returns
#     --------
#     Ainv                Creates a laplacian matrix solver for the specified geometry (including boundary conditions)
#     bound_pts_k         A list of k-indices for the boundary points (specific to the maskM) for modifying source
#
#     Notes
#     -------
#     ip and jp are related to the original matrix index scheme by ip = i +1 and j = jp + 1.
#
#     """
#     if maskM is None:
#
#         maskM = np.ones(shape)
#
#         maskM[1,1:-1]= -1
#         maskM[-2,1:-1] =-1
#         maskM[1:-1,1]=-1
#         maskM[1:-1,-2] =-1
#
#         maskM[0,:]= 0
#         maskM[-1,:] = 0
#         maskM[:,0]= 0
#         maskM[:,-1] = 0
#
#
#     # initialize the laplacian operator matrix:
#     A = np.zeros((grid_len,grid_len))
#     bound_pts_k = []
#
#     for k, (i, j) in enumerate(map_ij2k):
#
#         maskVal = maskM[i, j]
#
#         if maskVal == 1:  # interior point, do full algorithm...
#
#             k_ip1_j = map_ij2k.tolist().index([i + 1,j])
#             k_in1_j = map_ij2k.tolist().index([i-1,j])
#             k_i_jp1 = map_ij2k.tolist().index([i,j+1])
#             k_i_jn1 = map_ij2k.tolist().index([i,j-1])
#
#             A[k, k_ip1_j] = 1
#             A[k, k_in1_j] = 1
#             A[k, k_i_jp1] = 1
#             A[k, k_i_jn1] = 1
#             A[k,k] = -4
#
#         elif maskVal == -1:  # boundary point, treat as such
#             A[k,k] = 1
#             bound_pts_k.append(k)
#
#         elif maskVal == 0:  # exterior point, treat as such
#             A[k,k] = 1
#
#     # complete the laplacian operator by diving through by grid spacing squared:
#     Ai = A/(delx*dely)
#
#     # calculate the inverse, which is stored for solution calculation of Laplace and Poisson equations
#     Ainv = np.linalg.pinv(Ai)
#
#     return Ainv, bound_pts_k







