#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.exceptions import BetseSimException
from betse.util.path import files, paths
from scipy import interpolate as interp
from scipy import misc


class BitMapper(object):
    '''
    Object finding, loading, and converting a passed bitmap into a SciPy-based
    interpolation function.

    All bitmaps loaded must be square, with equal dimensions (pixels). It is
    recommended that the bitmaps used be 500x500 pixels. The bitmap should be a
    completely threshholded image, with black defining the area to be used as
    the clipping mask for the cell cluster, or defining the area for a tissue or
    boundary profile.

    Attributes
    ----------------------------
    clipping_matrix : ndarray
        Numpy matrix defining this bitmap's threshholded image.
    clipping_function : func
        SciPy-based interpolation function accepting an `(x, y)` point and
        returning `1.0` if that point resides outside this bitmap's colored
        pixel area or `0.0` otherwise. This function permits callers to filter
        a passed set of points in the space defined by `p.wsx` for the subset
        residing within this area.
    clipping_function_fast : func
        Fast variant of `clipping_function` otherwise sharing the same API.

    Attributes (clipPoints)
    ----------------------------
    The following attributes are available _only_ after calling the
    `clipPoints()` method.

    good_inds : ndarray
        #FIXME: Document us up the `BitMapper` bomb.
    good_points : ndarray
        Numpy matrix listing all points `(x, y)` residing inside this bitmap's
        colored area.
    '''

    def __init__(self,
        bitmap_matcher, xmin, xmax, ymin, ymax):
        '''
        Load, initialize, and create a threshholded interpolation matrix from
        the passed bitmap file.

        Parameters
        ----------------------------
        bitmap_matcher : TissuePickerABC
            Low-level BETSE-specific object describing this bitmap.
        xmin : float
            Minimum `x` coordinate accepted by these interpolation functions.
        xmax : float
            Maximum `x` coordinate accepted by these interpolation functions.
        ymin : float
            Minimum `y` coordinate accepted by these interpolation functions.
        ymax : float
            Maximum `y` coordinate accepted by these interpolation functions.

        '''
        # Avoid circular import dependencies.
        from betse.science.tissue.tissuepick import TissuePickerBitmap
        assert isinstance(bitmap_matcher, TissuePickerBitmap),\
            '{} not a BETSE-formatted bitmap.'.format(bitmap_matcher)

        # Load this bitmap as a flattened (i.e., grayscale) Numpy array.
        bitmap = misc.imread(bitmap_matcher.filename, flatten=1)

        # bitmap = np.asarray(bitmap, dtype=np.int)

        if bitmap.shape[0] != bitmap.shape[1]:
            raise BetseSimException(
                'Bitmap "{}" dimensions non-square '
                '(i.e., not of the same width and height).'.format(
                    bitmap_matcher.filename))

        # find the black pixels (a really basic threshholding!)
        point_inds = (bitmap != 255).nonzero()

        # define a new matrix the same shape as the image and set values to 0 or 1:
        self.msize = bitmap.shape[0]
        self.clipping_matrix = np.zeros((self.msize, self.msize))
        self.clipping_matrix[point_inds] = 1.0
        self.clipping_matrix = np.flipud(self.clipping_matrix)

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
        '''
        Initialize the `good_points` and `good_inds` attributes of this object
        to the subset of the passed list or vector of points residing in this
        bitmap's colored area by calling the clipping function previously
        initialized for this bitmap.

        Parameters
        -----------
        point_list_x : {list, ndarray}
            List or Numpy vector of x coordinates of points.
        point_list_y : {list, ndarray}
            List or Numpy vector of y coordinates of points.
        '''

        self.good_points = []
        self.good_inds = []

        for i, (x, y) in enumerate(zip(point_list_x, point_list_y)):
            if self.clipping_function(x,y) != 0.0:
                pt = [x,y]
                self.good_points.append(pt)
                self.good_inds.append(i)

        self.good_points = np.asarray(self.good_points)
        self.good_inds = np.asarray(self.good_inds)
