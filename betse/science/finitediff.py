#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

import numpy as np
import math
import scipy.ndimage


class FiniteDiffSolver(object):
    """
    Provides methods to establish a square recti-liniear grid, compute discrete Laplacian operators
       and solve the Poisson, Laplacian, or Heat equation in 2D using
       the Central Finite Difference method.

       Laplace Equation:
       d2U/dx2 + d2U/dy2 = 0

       Poisson Equation:
       d2U/dx2 + d2U/dy2 = fxy

       Heat Equation (spatially varying diffusion constant D):
       U = dD/dx*dC/dx + dD/dy*dC/dy + D*(d2C/dx2 + d2C/dy2)

    """

    def __init__(self):

        pass


    def makeGrid(self,delta = 0.05, xmin = -1, xmax = 1, ymin = -1, ymax = 1):

        """
        Create linear vectors, rectangular grids and expanded grids of spatial coordinates
        and solution vectors.

        delta:          spacing between grid points (same in the x and y dimensions)
        xmin:           minimum x-value  [m]
        xmax:           max x-value [m]
        ymin:           minimum y-value  [m]
        ymax:           max y value [m]

        """
        # check
        if xmax < xmin or ymax < ymin:
            print('error!')   # fixme do this properly

        # calculate the number of points in the x and y directions
        self.grid_nx = int((xmax - xmin)/delta)
        self.grid_ny = int((ymax - ymin)/delta)

        # Create linear vectors and rectangular grids of spatial coordinates (indexed by i and j)
        self.x = np.linspace(xmin,xmax,self.grid_nx)
        self.y = np.linspace(ymin,ymax,self.grid_ny)

        self.X,self.Y = np.meshgrid(self.x,self.y)

        self.delta = delta


    def modGrid(self, X, Y, F, bTop = 0.0, bBot = 0.0, bL=0.0, bR = 0.0, fTop = 0.0, fBot = 0.0, fL = 0.0,
        fR = 0.0, closedBound = True):

        """
        Already have X, Y and F(x,y) grids for the computation? Use this
        function to modify the existing grids for computation.

        Parameters
        -----------
        X, Y            2D numpy ndarrays of gridded spatial coordinates
        F               2D numpy ndarray of solution data values
        bTop       Neumann boundary value for solution at top of rectangular shape
        bBot       Neumann boundary value for solution at bottom of rectangular shape
        bL         Neumann boundary value for solution at left boundary
        bR         Neumann boundary value for solution at right boundary

        """

        self.grid_nx = X.shape[1]
        self.grid_ny = X.shape[0]

        self.delta = X[0,1] - X[0,0]

        self.xmin = np.min(X)
        self.xmax = np.max(X)
        self.ymin = np.min(Y)
        self.ymax = np.max(Y)

        self.F = F

        # Reduced x-y grid removing the boundary points:
        self.core_nx = self.grid_nx - 2
        self.core_ny = self.grid_ny - 2

        self.X_core = self.X[1:-1,1:-1]

        self.Y_core = self.Y[1:-1,1:-1]

        self.F_core = self.F[1:-1,1:-1]

        self.F_pad_top = self.F[0,:]
        self.F_pad_bot = self.F[-1,:]
        self.F_pad_L = self.F[:,0]
        self.F_pad_R = self.F[:,-1]

        self.bTop = self.F_pad_top[:]
        self.bTop[:] = bTop
        self.bBot = self.F_pad_bot[:]
        self.bBot[:] = bBot
        self.bL = self.F_pad_L[:]
        self.bL[:] = bL
        self.bR = self.F_pad_R[:]
        self.bR[:] = bR

        self.fTop = fTop
        self.fBot = fBot
        self.fL = fL
        self.fR = fR

        self.closed_bound = closedBound

        # Unravel the x-y grid into a long vector of dimensions (grid_nx - 2)(grid_ny -2) indexed by k:
        self.FF = self.F_core.ravel()
        self.unraveled_length = len(self.FF)

        # create a mapping from the ip, jp indexing to the k index of the unraveled core matrix:
        M = []
        for i in range(self.core_ny):
            row = []
            for j in range(self.core_nx):
                row.append([i,j])

            M.append(row)

        self.map_ij2k = []
        for row in M:
            for inds in row:
                self.map_ij2k.append(inds)

    def makeLaplacian(self):
        """
        Generate the discrete finite central difference 2D Laplacian operator based on:

        d2 Uij / dx2 + d2 Uij / dy2 = (Ui+1,j + Ui-1,j + Ui,j+1, Ui,j-1 - 4Uij)*1/(delta_x)**2

        To solve the Poisson equation:

        A*U = fxy

        by:

        U = Ainv*fxy

        Creates
        --------
        self.Ainv

        Modifies
        --------
        self.FF


        """

        A = np.zeros((self.unraveled_length,self.unraveled_length))

        for k, (ip, jp) in enumerate(self.map_ij2k):

            if ip + 1 < self.core_ny:

                k_ip1_j = self.map_ij2k.index([ip+1,jp])
                A[k, k_ip1_j] = 1

            else:
                # deal with the bounary value by moving to the RHS of the equation:
                bval = self.bBot[jp+1]
                self.FF[k] = self.FF[k] - (bval/self.delta**2)

            if ip - 1 >= 0:

                k_in1_j = self.map_ij2k.index([ip-1,jp])
                A[k, k_in1_j] = 1

            else:
                # deal with the bounary value by moving to the RHS of the equation:
                bval = self.bTop[jp+1]
                self.FF[k] = self.FF[k] - (bval/self.delta**2)

            if jp + 1 < self.core_nx:

                k_i_jp1 = self.map_ij2k.index([ip,jp+1])
                A[k, k_i_jp1] = 1

            else:
                # deal with the bounary value by moving to the RHS of the equation:
                bval = self.bR[ip + 1]
                self.FF[k] = self.FF[k] - (bval/self.delta**2)

            if jp -1 >= 0:
                k_i_jn1 = self.map_ij2k.index([ip,jp-1])
                A[k, k_i_jn1] = 1

            else:
                # deal with the bounary value by moving to the RHS of the equation:
                bval = self.bL[ip+1]
                self.FF[k] = self.FF[k] - (bval/self.delta**2)

            A[k,k] = -4

        # complete the laplacian operator by diving through by grid spacing squared:

        A = A/self.delta**2

        # calculate the inverse, which is stored for solution calculation of Laplace and Poisson equations
        self.Ainv = np.linalg.inv(A)



    def completeData(self,U):
        """
        Put back boundary values that were removed.

        """

        # U = np.vstack((U, self.F_pad_bot[0:-2]))
        # U = np.vstack((self.F_pad_top[0:-2], U))
        # U = np.column_stack((U, self.F_pad_R))
        # U = np.column_stack((self.F_pad_L, U))

        U = np.vstack((U, self.bBot[0:-2]))
        U = np.vstack((self.bTop[0:-2], U))
        U = np.column_stack((U, self.bR))
        U = np.column_stack((self.bL, U))

        U_padded = U

        return U_padded


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

    def makeHeatEq(self,dt, D=None):

        if D == None:
            self.D = np.zeros(self.X.shape)
            self.D[:,:] = 1e-9
        else:
            self.D = D[:]

        D = scipy.ndimage.filters.gaussian_filter(D, 1.0)

        self.grad_Dx, self.grad_Dy = gradient(D,self.delta)
        self.dt = dt
        self.Z_total = []
        self.Z_time = []

    def heatKernel(self,Z):
        """
         Computational kernel for a heat/diffusion equation with spatially varying diffusion
         constant D:

          dZ       dD   dZ     dD   dZ          d2Z    d2Z
         ----- =  ----*---- + ----*----  +  D*( ---- + ----)
          dt       dx   dx     dy   dy          dx2    dy2

         Parameters:
         ------------
         Z                A data matrix of Z(x,y) values defined as a Numpy ndarray


         Z1               The next iteration of Z(x,y)

        """

        grad_Zx, grad_Zy = gradient(Z,self.delta)
        del2_Z = laplacian(Z,self.delta)

        delta_Z = (self.grad_Dx*grad_Zx + self.grad_Dy*grad_Zy) + self.D*del2_Z

        Z = Z + self.dt*(delta_Z)

        #Dirichlet boundary condition (concentration at boundary):
    #     Z[0,:] = 1.0

        # Neumann boundary condition (flux at boundary)
        # zero flux boundaries
    #     Z[:,-1] = Z[:,-2]
    #     Z[:,0] = Z[:,1]
    #     Z[0,:] = Z[1,:]
    #     Z[-1,:] = Z[-2,:]

        sub_inds = (Z < 0).nonzero()
        Z[sub_inds] = 0

        self.Z_total.append(np.sum(Z))

        self.Z_time.append(Z[:])

def laplacian(F,delx):
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
    ddFx_interior = (F[:,2:] - 2*F[:,1:-1] + F[:,0:-2])/delx**2
    ddFy_interior = (F[2:,:] - 2*F[1:-1,:] + F[0:-2,:])/delx**2

    ddF_T = -(F[2,:] - 2*F[1,:] + F[0,:])/delx**2
    ddF_B = -(F[-3,:] - 2*F[-2,:] + F[-1,:])/delx**2
    ddF_L = -(F[:,2] - 2*F[:,1] + F[:,0])/delx**2
    ddF_R = -(F[:,-3] - 2*F[:,-2] + F[:,-1])/delx**2

    lapF_interior = ddFx_interior[1:-1,:] + ddFy_interior[:,1:-1]

    # initialize the dFx and dFy arrays:
    ddF = np.zeros(F.shape)

    ddF[1:-1,1:-1] = lapF_interior

    ddF[:,0] = ddF_L
    ddF[:,-1] = ddF_R

    ddF[0,:] = ddF_T
    ddF[-1,:] = ddF_B

    return ddF

def gradient(F,delx):
    # gradient using numpy slicing:

    # calculate the discrete central first derivatives on the internal mesh points:
    dF_interior_y = (F[:-2,:] - F[2:,:])/(2*delx)
    dF_interior_x = (F[:,:-2] - F[:,2:])/(2*delx)

    # calculate the discrete forward or backward first derivatives on the boundary points:
    dF_T = (F[0,:] - F[1,:])/delx
    dF_B = (F[-2,:] - F[-1,:])/delx
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




