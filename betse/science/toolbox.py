#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

"""
The toolbox module contains a number of functions that are used throughout the
BETSE project.
"""


import numpy as np
import scipy.spatial as sps
from scipy import interpolate as interp
# from scipy.special import expit
import math
import copy

#FIXME: For disambiguity, rename to flatten_list_of_lists(). Whereas the
#np.flatten() function flattens an array of arbitrary dimensionality, this
#function only flattens the first and second dimensions of a list of lists.
def flatten(ls_of_ls: list) -> tuple:
    """
    Flattens (i.e., un-nests) a doubly-nested Python "list of lists."

    Parameters
    ----------
    ls_of_ls        a nested list of lists, as in: [[a,b,c],[d,e,f],[g,h,i]]

    Returns
    -------
    ls_flat        a flattened version of the input list, as in:
                   [a,b,c,d,e,f,g,h,i]

    ind_map        returns the indices of the original nested list-of-lists at
                   the index of the new list, as in:
                    ind_map[5] = [0,5]   which would yield the same value for ls_flat[5] and ls_of_ls[0][5]

    Examples
    ----------
    >>> from betse.science import toolbox as tb
    >>> tb.flatten([['a', 'b', 'c'], ['d', 'e', 'f'], ['g', 'h', 'i']])
    (['a',   'b',   'c',   'd',   'e',   'f',   'g',   'h',   'i'],
     [[0,0], [0,1], [0,2], [1,0], [1,1], [1,2], [2,0], [2,1], [2,2]])

    Notes
    -------
    Requires python nested lists of lists. Numpy arrays have their own tools for
    this.
    """
    ls_flat = []
    ind_map =[]
    rind_map = copy.deepcopy(ls_of_ls)   # make a deepcopy of the nested list to get the right shape

    for i, sublist in enumerate(ls_of_ls):    # flatten the array and make a mapping: inds nest --> inds flat
        for j, val in enumerate(sublist):
            ls_flat.append(val)
            ind_map.append([i,j])

    for j, vals in enumerate(ind_map):   # go through the mapping and create a reverse mapping: ind flat --> inds nest
        rind_map[vals[0]][vals[1]] = j

    return ls_flat, ind_map, rind_map   # return the flattened list and the map and reverse map

def area(p):
    """
    Calculates the area of an arbitrarily shaped polygon defined by a set of
    counter-clockwise oriented points in 2D.

    Parameters
    ----------
    p               xy list of polygon points

    Returns
    -------
    area            area of a polygon in square meters

    Notes
    -------
    The algorithm is an application of Green's theorem for the functions -y and
    x, exactly in the way a planimeter works.
    """
    return 0.5 * abs(sum(x0*y1 - x1*y0 for ((x0, y0), (x1, y1)) in zip(p, p[1:] + [p[0]])))

def clip(subjectPolygon, clipPolygon):
    '''
    Clips the subject polygon to the boundary defined by the clip polygon.

    Notes
    -------
    This function implements the Sutherland-Hodgman polygon clipping algorithm,
    and largely comes courtesy Rosetta Code:
    http://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#Python
    '''
    def inside(p):
        return(cp2[0]-cp1[0])*(p[1]-cp1[1]) > (cp2[1]-cp1[1])*(p[0]-cp1[0])

    def computeIntersection():
        dc = [ cp1[0] - cp2[0], cp1[1] - cp2[1] ]
        dp = [ s[0] - e[0], s[1] - e[1] ]
        n1 = cp1[0] * cp2[1] - cp1[1] * cp2[0]
        n2 = s[0] * e[1] - s[1] * e[0]
        n3 = 1.0 / (dc[0] * dp[1] - dc[1] * dp[0])
        return [(n1*dp[0] - n2*dc[0]) * n3, (n1*dp[1] - n2*dc[1]) * n3]

    assert isinstance(subjectPolygon, list)
    assert isinstance(clipPolygon, list)
    assert len(subjectPolygon)
    assert len(clipPolygon)

    outputList = subjectPolygon
    cp1 = clipPolygon[-1]

    for clipVertex in clipPolygon:
        cp2 = clipVertex
        inputList = outputList
        outputList = []
        s = inputList[-1]

        for subjectVertex in inputList:
            e = subjectVertex
            if inside(e):
                if not inside(s):
                    outputList.append(computeIntersection())
                outputList.append(e)
            elif inside(s):
                outputList.append(computeIntersection())
            s = e
        cp1 = cp2

    return(outputList)

def alpha_shape(points, alpha):
    """
    Calculate the alpha_shape of a cluster of points in 2D.

    Parameters
    ----------
    points              A numpy array listing [[x1,y1],[x2,y2]...] for a collection of 2D points

    alpha               The filtering parameter (gauges which triangles to remove from Delaunay triangulation)
                        Note, for this application alpha = 1/d_cell is ideal.

    Returns
    --------
    concave_hull    A list of the indices to vertices in the points structure which define the concave hull
                    (these are all of the points on the boundary).

    Notes
    --------
    Unlike the convex hull, the alpha shape method and concave hull work for complex, concave geometries.
    The result depends on the alpha parameter. A value of alpha = 1/d_cell gives suitable results.

    """

    tri = sps.Delaunay(points)
    tri_edges = []
    circum_r_list = []

    # loop over triangles:
    # ia, ib, ic = indices of corner points of the
    # triangle
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

    return concave_hull

def sigmoid(x,g,y_sat):
    """
    A sigmoidal function (logistic curve) allowing user
    to specify a saturation level (y_sat) and growth rate (g).

    Parameters
    ----------
    x            Input values, may be numpy array or float
    g            Growth rate
    y_sat        Level at which growth saturates

    Returns
    --------
    y            Numpy array or float of values

    """
    y = (y_sat*np.exp(g*x))/(y_sat + (np.exp(g*x)-1))

    if np.isnan(y):
        y= 0

    return y

def hill(x,K,n):

    """
    The Hill equation (log-transformed sigmoid). Function ranges
    from y = 0 to +1.

    Parameters
    ----------
    x            Input values, may be numpy array or float. Note all x>0 !
    K            Value of x at which curve is 1/2 maximum (y=0.5)
    n            Hill co-efficient n<1 negative cooperativity, n>1 positive.

    Returns
    --------
    y            Numpy array or float of values

    """
    # assert x.all() > 0

    y = x**n/((K**n)+(x**n))

    return y

def step(t,t_on,t_change):
    """
    A step function (bounded by 0 and 1) based on a logistic curve
    and allowing user to specify time for step to come on (t_on) and time for
    change from zero to one to happen.

    Parameters
    ----------
    t            Input values, may be numpy array or float
    t_on         Time step turns on
    t_change     Time for change from 0 to 1 (off to on)

    Returns
    --------
    y            Numpy array or float of values

    """
    g = (1/t_change)*10
    y = 1/(1 + (np.exp(-g*(t-t_on))))
    # x = g*(t-t_on)
    # y = expit(x)

    # if np.isnan(y):
    #     y= 0

    return y

def pulse(t,t_on,t_off,t_change):
    """
    A pulse function (bounded by 0 and 1) based on logistic curves
    and allowing user to specify time for step to come on (t_on) and time for
    change from zero to one to happen, and time for step to come off (t_change).

    Parameters
    ----------
    t            Input values, may be numpy array or float
    t_on         Time step turns on
    t_off        Time step turns off
    t_change     Time for change from 0 to 1 (off to on)

    Returns
    --------
    y            Numpy array or float of values

    """
    g = (1/t_change)*10
    x1 = g*(t-t_on)
    x2 = g*(t-t_off)
    y1 = 1/(1 + (np.exp(-g*(t-t_on))))
    y2 = 1/(1 + (np.exp(-g*(t-t_off))))
    # y1 = expit(x1)
    # y2 = expit(x2)

    y = y1 - y2

    # if np.isnan(y):
    #     y=0

    return y

def H(x):

    y = 0.5*(np.sign(x) +1)

    return y

def emptyDict(dic):
    """
    Tests a dictionary to see if all keys are zero
    and returns True if so.

    Parameters:
    --------------
    dic            a dictionary

    :return:
    """
    # search_success = False
    #
    # while search_success is False:
    for entry in dic.values():

        if entry == 0:
            zero_dic = True

        elif type(entry) == list or type(entry)==str or entry != 0:
            zero_dic = False
            break


    return zero_dic

def grid_vector_data(xpts,ypts,zdata_x,zdata_y,cells,p):

    """
    Takes irregularly spaced vector data in the form of linear arrays of x,y,ux,uy and
    returns grids of X, Y, and Z_x and Z_y. Suitable for streamline plots (streamplot)

    Parameters
    ------------
    xpts                        Linear array of x-coordinates
    ypts                        Linear array of y-coordinates
    zdata_x                     Linear array of vector x-components
    zdata_y                     Linear array of vector y-components
    cells                       Instance of a Cells object

    Returns
    ---------
    X, Y, zi_x, zi_y            Arrays corresponding to xpts, ypts and the vector data points
    """
    # x_full = np.linspace(cells.xmin,cells.xmax,cells.msize)
    # y_full = np.linspace(cells.ymin,cells.ymax,cells.msize)

    xgrid = np.linspace(cells.xmin,cells.xmax,p.isamples)
    ygrid = np.linspace(cells.ymin,cells.ymax,p.isamples)
    X, Y = np.meshgrid(xgrid,ygrid)
    # xgrid = cells.x_v
    # ygrid = cells.y_v
    #
    # X = cells.x_2d
    # Y = cells.y_2d

    # create an interpolation function to resample the cluster mask matrix:
    # mask_funk = interp.interp2d(x_full,y_full,cells.cluster_mask)
    # mask_funk = interp.RectBivariateSpline(x_full,y_full,cells.cluster_mask)
    #
    # new_mask = mask_funk.ev(xgrid,ygrid)

    zi_x = interp.griddata((xpts,ypts),zdata_x,(X,Y))
    zi_x = np.nan_to_num(zi_x)
    # zi_x = np.multiply(zi_x,new_mask)

    zi_y = interp.griddata((xpts,ypts),zdata_y,(X,Y))
    zi_y = np.nan_to_num(zi_y)
    # zi_y = np.multiply(zi_y,new_mask)

    return X,Y,zi_x,zi_y

def griddata(xpts,ypts,zdata,gridsize,cells):

    """
    Takes irregularly spaced data in the form of linear arrays of x,y,z and
    returns grids of X, Y, and Z.

    Parameters
    ------------
    xpts                        Linear array of x-coordinates
    ypts                        Linear array of y-coordinates
    zdata                       Linear array of data values at x-y coordinate points
    gridsize                    Resolution of the interpolation grid (recommended 100x100)

    Returns
    ---------
    X, Y, zi_m                  Arrays corresponding to xpts, ypts and the zdata values. Note
                                zi_m is a masked array.

    """

    # x_full = np.linspace(cells.xmin,cells.xmax,cells.msize)
    # y_full = np.linspace(cells.ymin,cells.ymax,cells.msize)

    xlin = np.linspace(cells.xmin,cells.xmax,gridsize)
    ylin = np.linspace(cells.ymin,cells.ymax,gridsize)

    # create an interpolation function to resample the cluster mask matrix:
    # mask_funk =  interp.RectBivariateSpline(x_full,y_full,cells.cluster_mask)
    # new_mask = mask_funk.ev(xlin,ylin)

    # xmin = np.min(xpts)
    # xmax = np.max(xpts)
    # ymin = np.min(ypts)
    # ymax = np.max(ypts)
    #
    # xlin = np.linspace(xmin,xmax,gridsize)
    # ylin = np.linspace(ymin,ymax,gridsize)
    #
    X,Y = np.meshgrid(xlin,ylin)
#
    zi = interp.griddata((xpts,ypts),zdata,(X,Y))
    zi = np.nan_to_num(zi)
    # zi = np.multiply(new_mask,zi)

    return X, Y, zi

def makegrid(xpts,ypts,gridsize,cells):

    """
    Takes irregularly spaced data in the form of linear arrays of x,y and
    returns grids of X, Y.

    Parameters
    ------------
    xpts                        Linear array of x-coordinates
    ypts                        Linear array of y-coordinates
    gridsize                    Resolution of the interpolation grid (recommended 100x100)

    Returns
    ---------
    X, Y                   Arrays corresponding to xpts, ypts.
    dx, dY                 Array spacing vectors

    """

    # xmin = np.min(xpts)
    # xmax = np.max(xpts)
    # ymin = np.min(ypts)
    # ymax = np.max(ypts)

    # xlin = np.linspace(xmin,xmax,gridsize)
    # ylin = np.linspace(ymin,ymax,gridsize)

    xlin = np.linspace(cells.xmin,cells.xmax,gridsize)
    ylin = np.linspace(cells.ymin,cells.ymax,gridsize)

    X,Y = np.meshgrid(xlin,ylin)

    dx = np.gradient(xlin)
    dy = np.gradient(ylin)

    return X, Y, dx, dy

def RK4(f):

    """
    This rather confusing looking function was obtained from
    Rossetta Code's Python RK4 implementation.

    It is modified for a differential equation that does not
    depend on time in its rate equation.

    As an example, suppose there is a differential equation
    dc/dt = sqrt(c)

    To solve using RK4, we would define:

    dc = RK4(lambda c: sqrt(c))

    And update it as:
    c1 = co + dc(co,dt)

    Where co is the original value, dt is the time-step

    """
    return lambda y, dt: (
        lambda dy1: (
            lambda dy2: (
                lambda dy3: (
                    lambda dy4: (dy1 + 2 * dy2 + 2 * dy3 + dy4) / 6
                )(dt * f(y + dy3))
            )(dt * f(y + dy2 / 2))
        )(dt * f(y + dy1 / 2))
    )(dt * f(y))







