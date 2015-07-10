#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

"""
The toolbox module contains a number of functions that are used throughout the
BETSE project.
"""

# FIXME SESS: To my scruffy chagrin, line 315 of this module just threw up the
# following overflow error during a recent "beste try" run:
#
#    /home/leycec/py/betse/betse/science/toolbox.py:315: RuntimeWarning: overflow encountered in exp
#      y2 = 1/(1 + (np.exp(-g*(t-t_off))))
#
# Overflows are *REALLY* bad, so I'm a wee surprised that Python didn't just
# carpet-bomb the whole thing with a fatal exception. "What is overflow?", you
# are probably wondering to your beautiful self. Overflow is an arithmetic error
# specific to computer science. Computers can only hold so much data, right? In
# most older programming languages (e.g., C, C++, Java), integer and floating
# point variables can only hold a finite number of bytes. This means that the
# range of values representable by such variables is bounded. 32-bit integer
# variables, for example, are only able to store values in the range
# [-2147483647, 2147483647]. That probably seems like a really large number of
# possible values. Who would need to represent anything larger than 2 billion!?
# In the real-world, though, it's depressingly common for this range to be
# exceeded. As you may have astutely surmised, this error condition is referred
# to as "overflow."
#
# And it's *REALLY*, *REALLY* bad. The reason why is that overflow causes values
# to wrap around to the other side of the range of supported values. For
# example, consider the value 2147483648. As the range above implies,
# 2147483648 is not representable as a 32-bit integer variable. What happens,
# exactly, when you try to assign 2147483648 to a 32-bit integer variable in
# older programming languages -- like, say, Python 2? Let's find out:
#
#     >>> uh_oh  = 2147483648
#     >>> oh_god = 2147483649
#     >>> print "Oh by God! It's crazy overflow: " + uh_oh + oh_god
#     Oh by God! It's crazy overflow: -2147483648 -2147483647
#
# Yes, you read that correctly. Thanks to the horrifying magic of overflow,
# the positive values 2147483648 and 2147483649 wrap around to the other side of
# the range of supported values... and magically become *VERY*, *VERY* negative.
# We're probably beginning to perceive why overflow may not be a good thing.
#
# In Python 2, integers were bounded just like they were in older programming
# languages. Longs, however, were *NOT*. Longs are an old-school arithmetic
# type. (Just roll with me here, O.K.?) In older languages, longs were typically
# bounded 64-bit integers. In Python 2, longs were what-we-call "arbitrary
# precision" or "unbounded" integers. As the crude jargon implies, longs in
# Python 2 could represent any possible integer value! Sweet, right? Except
# that integers were still bounded.
#
# Everyone thought this was pretty stupid. I mean, why support two different
# integer types -- one of which is bounded? Why not just have a single unbounded
# integer type? Exactly, bro! And that's exactly why we use Python 3 rather than
# Python 2. I probably never mentioned this, but Python 3 made everyone's (and
# especially every scientist and mathematician's) life much easier by just
# getting rid of bounded integers entirely. In Python 3, there are only
# unbounded integers -- which is why you probably haven't run into too many
# oveflow issues.
#
# Unfortunately, what applies to integers does *NOT* apply to floats. In both
# Python 2 and 3, all floats are bounded and hence susceptible to overflow. Yes,
# this sucks. It turns out that it's much, much harder to implement unbounded
# floats than it is unbounded integers. In fact, there's a whole body of
# research devoted to the problem of implementing arbitrary precision floats.
# This means that you have to be *REALLY*, *REALLY* careful when performing
# arithmetic operations on floats -- especially operations that are likely to
# substantially increase or decrease the values of floats. Exponentiation is the
# usual offender here.
#
# Unsurprisingly, this is exactly the operation that's causing overflow above:
# np.exp(). This appears to be a fairly common issue, as shown by googling the
# following search terms:
#
#     numpy "overflow encountered in exp"
#
# Let's revisit the problematic subexpression:
#
#     np.exp(-g*(t-t_off))
#
# Note the negative sign. "Wait a second here, fat bub!" Shouldn't the
# exponentation of a negative number produce a really small number? Why would
# *THAT* overflow? To overflow, we need to produce a really big number, right?
# O.K., so here's where we introduce a new term: "underflow." Underflow is
# basically the inverse of overflow. With overflow, the issue is that the
# resulting number is too large to be represented. With underflow, the issue is
# that the resulting number is too *SMALL* to be represented. Yup! But doesn't
# that mean that the above error message should say "RuntimeWarning: underflow
# encountered in exp" rather than "RuntimeWarning: overflow encountered in exp".
#
# Yes. Yes it does. But even numpy isn't perfect. It said "overflow," but it
# really meant "underflow." But it doesn't really matter, does it? The end
# result is the same: a *REALLY* bjorked-up number possibly bjorking-up the
# entire simulation run.
#
# Enough useless theory already! "Just tell me how to fix it, you."
#
# Hey, hey: don't shoot the ugly-faced messenger. And... well, I don't know how
# to fix it. Not exactly, anyway. Overflow and underflow are both squirmy issues
# that have long vexed software engineers since Time Itself Beshat Itself.
#
# Here are a few helpful suggestions that we should probably consider pursuing:
#
# * Force numpy to raise fatal exceptions rather than printing non-fatal
#   warnings on both overflows and underflows. Happily, numpy already provides a
#   function for doing just that: numpy.seterr(). Overflows and underflows are
#   sufficiently terrible that we really want to halt the entire simulation if
#   even a single one occurs. I'd be a giddy schoolboy to make this happen. Just
#   let me know! For my own rumination, the numpy.seterr() docs live at:
#   http://docs.scipy.org/doc/numpy/reference/generated/numpy.seterr.html
# * Use the external "bigfloat" library. This is an external Python library
#   supporting unbounded floats. Unfortunately, there are two problems with this
#   approach:
#   * Numpy and bigfloat don't know about each other. While bigfloat *DOES*
#     provide an exp() function, that function only accepts single float values
#     rather than the numpy arrays accepted by the numpy.exp() function. While
#     the problematic line of code that's underflowing *COULD* be converted to
#     call bigfloat.exp() rather than numpy.exp(), we'd need to manually loop
#     over all array values to do so. In other words: super-bad slowdown.
#   * Bigfloat has some pretty hefty requirements. Like scipy, bigfloat is
#     simply a wrapper for an external library implemented in a completely
#     different language. This means that installing bigfloat is non-trivial --
#     especially under Windows and OS X, our primary target platforms.
# * Rewrite the problematic expression so as to prevent underflow. This is
#   usually the best approach, but it's also the hardest. In some cases, it's
#   infeasible. This may be one of them. The usual way that this is done is by
#   taking the log of your very small or very large input values (thus producing
#   much more reasonable values that can be worked with without worrying about
#   underflow or overflow) and then exponentiating the resulting output value
#   (thus "undoing" the prior logging). I hab no ideaz!
# * Ignore underflow. This isn't necessarily the worst option. Overflow is
#   really terrible, as documented above, because it fundamentally alters the
#   resulting value in unexpected (and usually catastrophic) ways. Underflow, by
#   comparison, isn't quite as bad. In most cases, the resulting value just gets
#   truncated to 0. This *CAN* be bad, particularly in expressions like
#   "1/np.exp(-blah_blah)" which, given an underflow, would reduce to a killer
#   "1/0". Uhps! This implies that underflow shouldn't be ignored for general-
#   purpose calculations. In the above case, however, the addition by 1
#   guarantees that we won't be dividing by 0 even if np.exp() causes underflow:
#
#      y2 = 1/(1 + (np.exp(-g*(t-t_off))))
#
#   This means that we can probably do something like this:
#
#      old_settings = np.seterror(under = 'ignore')
#      y2 = 1/(1 + (np.exp(-g*(t-t_off))))
#      np.seterr(**old_settings)
#
#   What this does is temporarily ignore underflow for that specific operation
#   and that specific operation only. We'll still want to ensure that overflows
#   and underflows raise exceptions everywhere else -- probably by adding logic
#   to "betse.ignition" resembling:
#
#      np.seterror(all = 'raise')
#
# But wait! There's more. I see with my perscipacious eye that the problematic
# line is actually calculating the logistic function. Now, wouldn't you know it,
# but the "scipy.special" module already provides the logistic function in a
# manner supporting both single floats and numpy arrays. For some inane reason,
# it's called expit(). (Rolls off the tongue. "Just expit, baby!") See:
#
#     http://docs.scipy.org/doc/scipy-0.15.0/reference/generated/scipy.special.expit.html
#
# Given that, you almost definitely want to refactor each statement of the form
# "y = 1/(1 + (np.exp(-...)))" in the BETSE codebase to a statement of the form
# "y = spsp.expit(...)", assuming an earlier "import scipy.special as spsp"
# statement. Using PyCharm, it shouldn't be *TOO* meddlesome to find all places
# where the exp() function is being called. Not all of them will be used for
# calculating the logistic function, of course -- but many of them probably
# will. For example, the problematic line above can be refactored:
#
#     # From this...
#     y2 = 1/(1 + (np.exp(-g*(t-t_off))))
#
#     # ...to this. Suck it, verbosity. Note the conspicuous lack of a negative
#     # sign here. We'll want to be supremely careful when refactoring
#     # statements in the form above to statements in the form below.
#     import scipy.special as spsp
#     y2 = spsp.expit(g*(t-t_off))
#
# The latter is almost guaranteed to avoid both overflow and underflow issues
# *AND* to be substantially faster than the former. Why? Because scipy rocks my
# skinny Casbah. Most scipy developers are aware of this sort of thing, so
# they've (hopefully) gone to great lengths to avoid overflow and underflow with
# awesome mathematical cleverness. We'll want to test that, but I'm intuitively
# confident that it's true. It has to be, or I'll eat my Tilley Hat.
#
# Yay! Yaay!!

import numpy as np
import scipy.spatial as sps
from scipy import interpolate as interp
import math
import copy

def flatten(ls_of_ls):
    """
    Flattens (i.e. un-nests) a nested python "list of lists".

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
        area = math.sqrt(s*(s-a)*(s-b)*(s-c))

        if area > 0:
            circum_r = a*b*c/(4.0*area)

        if area == 0:
            circum_r = a*b*c/(4.0*1e-25)

        # circum_r = a*b*c/(4.0*area)
        # circum_r_list.append(circum_r)

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


    # Sess' improved concave hull finding algorithm.
    # self.exterior_edges = set()
    # self.interior_edges = set()
    #
    # def tri_edges_append(self, edge_list):
    #     edge_tuple = tuple(edge_list)
    #     if edge_tuple in self.exterior_edges:
    #         self.exterior_edges.remove(edge_tuple)
    #         self.interior_edges.add(edge_tuple)
    #     elif edge_tuple not in self.interior_edges:
    #         self.exterior_edges.add(edge_tuple)
    #
    # concave_hull = []
    # for edge_tuple in self.exterior_edges:
    #     concave_hull.append(list(edge_tuple))
    #
    # self.exterior_edges = None
    # self.interior_edges = None

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
    y1 = 1/(1 + (np.exp(-g*(t-t_on))))
    y2 = 1/(1 + (np.exp(-g*(t-t_off))))
    y = y1 - y2
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
    # while search_success == False:
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
    cells                       Instance of a World object

    Returns
    ---------
    X, Y, zi_x, zi_y            Arrays corresponding to xpts, ypts and the vector data points
    """
    x_full = np.linspace(cells.xmin,cells.xmax,cells.msize)
    y_full = np.linspace(cells.ymin,cells.ymax,cells.msize)

    xgrid = np.linspace(cells.xmin,cells.xmax,p.isamples)
    ygrid = np.linspace(cells.ymin,cells.ymax,p.isamples)
    X, Y = np.meshgrid(xgrid,ygrid)
    # xgrid = cells.x_v
    # ygrid = cells.y_v
    #
    # X = cells.x_2d
    # Y = cells.y_2d

    # create an interpolation function to resample the cluster mask matrix:
    mask_funk = interp.interp2d(x_full,y_full,cells.cluster_mask)
    new_mask = mask_funk(xgrid,ygrid)

    zi_x = interp.griddata((xpts,ypts),zdata_x,(X,Y))
    zi_x = np.nan_to_num(zi_x)
    zi_x = np.multiply(zi_x,new_mask)

    zi_y = interp.griddata((xpts,ypts),zdata_y,(X,Y))
    zi_y = np.nan_to_num(zi_y)
    zi_y = np.multiply(zi_y,new_mask)

    return X,Y,zi_x,zi_y

def griddata(xpts,ypts,zdata,gridsize):

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

    xmin = np.min(xpts)
    xmax = np.max(xpts)
    ymin = np.min(ypts)
    ymax = np.max(ypts)

    xlin = np.linspace(xmin,xmax,gridsize)
    ylin = np.linspace(ymin,ymax,gridsize)

    X,Y = np.meshgrid(xlin,ylin)

    zi = interp.griddata((xpts,ypts),zdata,(X,Y))

    return X, Y, zi

def makegrid(xpts,ypts,gridsize):

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

    xmin = np.min(xpts)
    xmax = np.max(xpts)
    ymin = np.min(ypts)
    ymax = np.max(ypts)

    xlin = np.linspace(xmin,xmax,gridsize)
    ylin = np.linspace(ymin,ymax,gridsize)

    X,Y = np.meshgrid(xlin,ylin)

    dx = np.gradient(xlin)
    dy = np.gradient(ylin)

    return X, Y, dx, dy

def periodic(cells,sim,p):

    if p.run_sim == True:
            tt = np.linspace(0,p.sim_tsteps*p.dt,p.sim_tsteps)
    else:
        tt = np.linspace(0,p.init_tsteps*p.dt,p.init_tsteps)   # timestep vector

    yy = np.sin(tt*2*np.pi*p.periodic_properties['frequency'] + p.periodic_properties['phase'])

    y = interp.interp1d(tt,yy)

    return y

def gradient_x(cells,sim,p):


    x = np.linspace(cells.clust_x_min,cells.clust_x_max,100)

    yy = x*p.gradient_x_properties['slope'] + p.gradient_x_properties['offset']

    yy = yy/np.max(yy)

    y = interp.interp1d(x,yy)

    return y

def gradient_y(cells,sim,p):

    z = np.linspace(cells.clust_y_min,cells.clust_y_max,100)

    yy = z*p.gradient_y_properties['slope'] + p.gradient_y_properties['offset']

    yy = yy/np.max(yy)

    y = interp.interp1d(z,yy)

    return y

