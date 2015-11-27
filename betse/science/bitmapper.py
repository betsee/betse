#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME figure out why scipy can't read pngs!

import numpy as np
from betse.exceptions import BetseExceptionSimulation
from betse.util.path import files, paths
from scipy import interpolate as interp
from scipy import misc

class BitMapper(object):
    """
    Finds a designated bitmap, loads it, makes it into an interpolation
    function, and allows the user to screen a set of points in the space defined
    by `p.wsx`, to see if they fall within the colored area of the bitmap.

    All bitmaps loaded must be square, with equal dimensions (pixels). It is
    recommended that the bitmaps used be 500x500 pixels. The bitmap should be a
    completely threshholded image, with black defining the area to be used as
    the clipping mask for the cell cluster, or defining the area for a tissue or
    boundary profile.

    Attributes
    ----------------------------
    bitmapFile : str
        Absolute path of this bitmap.
    clippingMatrix : ndarray
        Numpy matrix defining this bitmap's threshholded image.
    clipping_function : func
        SciPy interpolation function object accepting a point `(x, y)` and
        returning `1.0` if that point resides inside this bitmap's colored area.
    good_points : ndarray
        Numpy matrix listing all points `(x, y)` residing inside this bitmap's
        colored area.
    """

    #FIXME: We currently ignore "threshhold_val". Should we just remove it, for
    #limplicity's sake? ("Think of the simplistic children!")
    def __init__(self,
        p, bitmap_filename, xmin, xmax, ymin, ymax, threshhold_val = 0.0):
        """
        Loads, initializes, and creates a threshholded interpolation matrix from
        the passed bitmap file.

        Parameters
        ----------------------------
        p : Parameters
            Instance of the `Parameters` object.
        bitmap_filename : str
            Absolute or relative path of the bitmap to be loaded. If relative
            (i.e., _not_ prefixed by a directory separator), this path will be
            canonicalized into an absolute path relative to the directory
            containing our source configuration file.
        threshhold_val : int
            The value of the pixel to threshhold to. Greyscale runs from 0
            (black) to 255 (white). The default is to get black pixels (0).
        """

        # If this is a relative path, convert this into an absolute path
        # relative to the directory containing the source configuration file.
        if paths.is_relative(bitmap_filename):
            bitmap_filename = paths.join(p.config_dirname, bitmap_filename)

        # If this bitmap does *NOT* exist, raise an exception.
        self.bitmapFile = bitmap_filename
        files.die_unless_file(self.bitmapFile)

        # Load this bitmap as a flattened (i.e., grayscale) array.
        bitmap = misc.imread(self.bitmapFile, flatten=1)

        if bitmap.shape[0] != bitmap.shape[1]:
            raise BetseExceptionSimulation(
                'Bitmap "{}" dimensions not square '
                '(i.e., of the same width and height).'.format(self.bitmapFile))

        self.msize = bitmap.shape[0]

        # find the black pixels (a really basic threshholding!)
        point_inds = (bitmap == 0).nonzero()

        # define a new matrix the same shape as the image and set values to 0 or 1:
        self.clippingMatrix = np.zeros((self.msize,self.msize))
        self.clippingMatrix[point_inds]= 1.0
        self.clippingMatrix = np.flipud(self.clippingMatrix)

        # create spatial data vectors that span the extent of the cell seeds and match bitmap pixel number:
        xpts = np.linspace(xmin,xmax,self.msize)
        ypts = np.linspace(ymin,ymax,self.msize)

        # create an interpolation function that returns zero if the query point is outside the mask and
        # 1 if the query point is in the mask:
        self.clipping_function = interp.interp2d(xpts,ypts,self.clippingMatrix)
        self.clipping_function_fast = interp.RectBivariateSpline(xpts,ypts,self.clippingMatrix)

    def clipPoints(self, point_list_x, point_list_y):
        """
        Initialize the `good_points` attribute to the subset of the passed list
        or vector of points residing in this bitmap's colored area by calling
        the clipping function previously initialized for this bitmap.

        Parameters
        -----------
        point_list_x : {list, ndarray}
            List or Numpy vector of x coordinates of points.
        point_list_y : {list, ndarray}
            List or Numpy vector of y coordinates of points.
        """

        self.good_points = []
        self.good_inds = []

        for i, (x, y) in enumerate(zip(point_list_x,point_list_y)):
            z = self.clipping_function(x,y)
            if z != 0.0:
                pt = [x,y]
                self.good_points.append(pt)
                self.good_inds.append(i)

        self.good_points = np.asarray(self.good_points)
        self.good_inds = np.asarray(self.good_inds)
