#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME figure out why scipy can't read pngs!

import numpy as np
from betse.exceptions import BetseExceptionMethod, BetseExceptionSimulation
from betse.util.path import files, paths
from scipy import interpolate as interp
from scipy import misc


class GeometryBitmap(object):
    """
    Finds a designated bitmap, loads it, makes it into an interpolation
    function, and allows the user to screen a set of points in the space defined
    by `p.wsx` to see if they fall within the colored area of this bitmap.

    All bitmaps loaded must be square, with equal dimensions (pixels). It is
    recommended that the bitmaps used be 500x500 pixels. The bitmap should be a
    completely threshholded image, with black defining the area to be used as
    the clipping mask for the cell cluster, or defining the area for a tissue or
    boundary profile.

    Attributes
    ----------------------------
    filename : str
        Absolute path of this bitmap.
    clipping_matrix : ndarray
        Numpy matrix defining this bitmap's threshholded image.

    Attributes (makeClippingFunctions)
    ----------------------------
    The following attributes are available _only_ after calling the
    `makeClippingFunctions()` method.

    clipping_function : func
        SciPy interpolation function accepting an `(x, y)` point and
        returning `1.0` if that point resides outside this bitmap's colored area
        _or_ `0.0` otherwise.
    clipping_function_fast : func
        Fast variant of `clipping_function` otherwise sharing the same API.

    Attributes (clipPoints)
    ----------------------------
    The following attributes are available _only_ after calling the
    `clipPoints()` method.

    good_inds : ndarray
        #FIXME: Document us up the `GeometryBitmap` bomb.
    good_points : ndarray
        Numpy matrix listing all points `(x, y)` residing inside this bitmap's
        colored area.
    """

    #FIXME: We currently ignore "threshhold_val". Should we just remove it, for
    #limplicity's sake? ("Think of the simplistic children!")
    def __init__(self, filename, dirname, threshhold_val = 0.0):
        """
        Loads, initializes, and creates a threshholded interpolation matrix from
        the passed bitmap file.

        Parameters
        ----------------------------
        filename : str
            Absolute or relative path of the bitmap to be loaded. If relative
            (i.e., _not_ prefixed by a directory separator), this path will be
            canonicalized into an absolute path relative to the directory
            containing our source configuration file.
        dirname : str
            Absolute path of the directory containing the path of the bitmap to
            be loaded (i.e., `filename`). If that path is relative, that path
            will be prefixed by this path to convert that path into an absolute
            path; otherwise, this path will be ignored.
        threshhold_val : int
            The value of the pixel to threshhold to. Greyscale runs from 0
            (black) to 255 (white). The default is to get black pixels (0).
        """
        assert isinstance(filename, str), '{} not a string.'.format(filename)
        assert isinstance( dirname, str), '{} not a string.'.format(dirname)

        # If this is a relative path, convert this into an absolute path
        # relative to the directory containing the source configuration file.
        if paths.is_relative(filename):
            filename = paths.join(dirname, filename)

        # If this absolute path is *NOT* an existing file, raise an exception.
        files.die_unless_file(filename)

        # Store this absolute path.
        self.filename = filename

        # Load this bitmap as a flattened (i.e., grayscale) Numpy array.
        bitmap = misc.imread(filename, flatten=1)

        if bitmap.shape[0] != bitmap.shape[1]:
            raise BetseExceptionSimulation(
                'Bitmap "{}" dimensions not square '
                '(i.e., of the same width and height).'.format(filename))

        self.msize = bitmap.shape[0]

        # find the black pixels (a really basic threshholding!)
        point_inds = (bitmap == 0).nonzero()

        # define a new matrix the same shape as the image and set values to 0 or 1:
        self.clipping_matrix = np.zeros((self.msize, self.msize))
        self.clipping_matrix[point_inds] = 1.0
        self.clipping_matrix = np.flipud(self.clipping_matrix)

    def makeClippingFunctions(self, xmin, xmax, ymin, ymax):
        """
        Initialize the `clipping_function` and `clipping_function_fast`
        attributes of this object to interpolation functions generated for the
        passed range of `(x, y)` points.

        Parameters
        -----------
        xmin : float
            Minimum `x` coordinate accepted by these interpolation functions.
        xmax : float
            Maximum `x` coordinate accepted by these interpolation functions.
        ymin : float
            Minimum `y` coordinate accepted by these interpolation functions.
        ymax : float
            Maximum `y` coordinate accepted by these interpolation functions.
        """
        # Create spatial data vectors that span the extent of the cell seeds and
        # match bitmap pixel number.
        xpts = np.linspace(xmin, xmax, self.msize)
        ypts = np.linspace(ymin, ymax, self.msize)

        # Create an interpolation function that returns zero if the query point
        # is outside the mask and 1 if the query point is in the mask.
        self.clipping_function = interp.interp2d(
            xpts, ypts, self.clipping_matrix)
        self.clipping_function_fast = interp.RectBivariateSpline(
            xpts, ypts, self.clipping_matrix)

    def clipPoints(self, point_list_x, point_list_y):
        """
        Initialize the `good_points` and `good_inds` attributes of this object
        to the subset of the passed list or vector of points residing in this
        bitmap's colored area by calling the clipping function previously
        initialized for this bitmap.

        This method must be called _after_ the `makeClippingFunctions()` method
        has been called.

        Parameters
        -----------
        point_list_x : {list, ndarray}
            List or Numpy vector of x coordinates of points.
        point_list_y : {list, ndarray}
            List or Numpy vector of y coordinates of points.
        """

        # If the makeClippingFunctions() method has not yet been called, raise
        # an exception.
        if not hasattr(self, 'clipping_function'):
            raise BetseExceptionMethod(
                'GeometryBitmap.makeClippingFunctions() not called before calling '
                'GeometryBitmap.clipPoints().')

        self.good_points = []
        self.good_inds = []

        for i, (x, y) in enumerate(zip(point_list_x,point_list_y)):
            #FIXME: The following two lines appear to reduce to a single line:
            #    if self.clipping_function(x,y) != 0.0:
            z = self.clipping_function(x,y)
            if z != 0.0:
                pt = [x,y]
                self.good_points.append(pt)
                self.good_inds.append(i)

        self.good_points = np.asarray(self.good_points)
        self.good_inds = np.asarray(self.good_inds)
