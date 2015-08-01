#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

import numpy as np
import math
import scipy.ndimage
import scipy.spatial as sps


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

    def makeGrid(self,grid_n, xmin = -1, xmax = 1, ymin = -1, ymax = 1):

        """
        Create linear vectors, rectangular grids and expanded grids of spatial coordinates
        and solution vectors.

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

    def cell_grid(self,grid_n, xmin, xmax, ymin, ymax):

        """
        Make a staggered n x n grid for use with finite volume or MACs method.
        One n x n grid contains the cell centres.
        One n+1 x n+1 grid contains the cell corners (vertices)
        One n+1 x n grid contains the u-coordinate points for veloctiy
        One n x n+1 grid contains the v-coordinate points for velocity

        """
        self.grid_n = grid_n
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

        xv = np.linspace(xmin,xmax,grid_n+1)
        yv = np.linspace(ymin,ymax,grid_n+1)
        X_verts, Y_verts = np.meshgrid(xv,yv)

        # grid and coordinate vectors defining the vertices of each cell:
        x_verts = X_verts.ravel()
        y_verts = Y_verts.ravel()
        self.xy_verts = np.column_stack((x_verts,y_verts))
        self.verts_shape = X_verts.shape

        # grid and coordinate vectors defining the centre points of each cell:
        self.cents_X = (X_verts[0:-1,0:-1] + X_verts[0:-1,1:])/2
        self.cents_Y = (Y_verts[0:-1,0:-1] + Y_verts[1:,0:-1])/2

        x_cent = self.cents_X.ravel()
        y_cent = self.cents_Y.ravel()

        self.xy_cents = np.column_stack((x_cent,y_cent))

        self.cents_shape = self.cents_X.shape

        self.map_ij2k_cents = []

        for i, row in enumerate(self.cents_X):
            for j, val in enumerate(row):
                self.map_ij2k_cents.append([i,j])

        self.map_ij2k_cents = np.asarray(self.map_ij2k_cents)

        # cell side midpoints in the x (u) and x (v) directions:
        mx = self.cents_X[0,:]
        my = Y_verts[:,0]
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

        nx = X_verts[0,:]
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

        # calculate the delta values for different grids:
        ccx = np.unique(self.xy_cents[:,0])
        ccy = np.unique(self.xy_cents[:,1])

        self.delta_cents_x = ccx[1] - ccx[0]
        self.delta_cents_y = ccy[1] - ccy[0]

        mmx = np.unique(self.u_pts[:,0])
        mmy = np.unique(self.u_pts[:,1])
        self.delta_mids_x = mmx[1] - mmx[0]
        self.delta_mids_y = mmy[1] - mmy[0]

def makeLaplacian(grid_len,shape,map_ij2k,delx,dely,maskM=None):
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
    grid_len            The length of the unravelled 2D array
    map_ij2k            An array where each index k corresponding to unravelled array
                        maps to the i,jth element of the 2D array
    delx                Grid spacing in the x-dimension
    dely                Grid spacing in the y-dimension
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
    if maskM is None:

        maskM = np.ones(shape)

        maskM[1,1:-1]= -1
        maskM[-2,1:-1] =-1
        maskM[1:-1,1]=-1
        maskM[1:-1,-2] =-1

        maskM[0,:]= 0
        maskM[-1,:] = 0
        maskM[:,0]= 0
        maskM[:,-1] = 0


    # initialize the laplacian operator matrix:
    A = np.zeros((grid_len,grid_len))
    bound_pts_k = []

    for k, (i, j) in enumerate(map_ij2k):

        maskVal = maskM[i, j]

        if maskVal == 1:  # interior point, do full algorithm...

            k_ip1_j = map_ij2k.tolist().index([i + 1,j])
            k_in1_j = map_ij2k.tolist().index([i-1,j])
            k_i_jp1 = map_ij2k.tolist().index([i,j+1])
            k_i_jn1 = map_ij2k.tolist().index([i,j-1])

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
    Ai = A/(delx*dely)

    # calculate the inverse, which is stored for solution calculation of Laplace and Poisson equations
    Ainv = np.linalg.inv(Ai)

    return Ainv, bound_pts_k

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

    # if dely is None:
    #     dely = delx
    #
    # ddFx_interior = (F[:,2:] - 2*F[:,1:-1] + F[:,0:-2])/(delx*dely)
    # ddFy_interior = (F[2:,:] - 2*F[1:-1,:] + F[0:-2,:])/(delx*dely)
    #
    # ddF_T = (F[-3,:] - 2*F[-2,:] + F[-1,:])/(dely*dely)
    # ddF_B = (F[0,:] - 2*F[-1,:] + F[-2,:])/(dely*dely)
    # ddF_L = (F[:,2] - 2*F[:,1] + F[:,0])/(delx*delx)
    # ddF_R = (F[:,-3] - 2*F[:,-2] + F[:,-1])/(delx*delx)
    #
    # lapF_interior = ddFx_interior[1:-1,:] + ddFy_interior[:,1:-1]
    #
    # # initialize the dFx and dFy arrays:
    # ddF = np.zeros(F.shape)
    #
    # ddF[1:-1,1:-1] = lapF_interior
    #
    # ddF[:,0] = ddF_L
    # ddF[:,-1] = ddF_R
    #
    # ddF[-1,:] = ddF_T
    # ddF[0,:] = ddF_B

    #------------------------------

    gx, gy = gradient(F, delx,dely)
    gxx = diff(gx, delx,axis=0)
    gyy = diff(gy,dely,axis=1)
    ddF = gxx + gyy

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
    dF_L = (F[:,0] - F[:,1])/delx
    dF_R = (F[:,-2] - F[:,-1])/delx

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
        dF_B = (F[1,:] - F[0,:])/delx
        dF_T = (F[-1,:] - F[-2,:])/delx

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
    xy_pts          unravelled xy points for the World of M
    X               2D X matrix for the World of M
    Y               2D Y matrix for the World of M
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
    #     self.Ainv = np.linalg.inv(A)







