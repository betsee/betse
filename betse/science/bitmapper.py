#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME figure out why scipy can't read pngs!

import numpy as np
from betse.exceptions import BetseExceptionSimulation
from betse.util.path import files, paths
from scipy import interpolate as interp
from scipy import misc

class Bitmapper(object):
    """
    Finds a designated bitmap, loads it, makes it into an interpolation
    function, and allows the user to screen a set of points in the space defined
    by `p.wsx`, to see if they fall within the colored area of the bitmap.

    All bitmaps loaded must be square, with equal dimensions (pixels). It is
    recommended that the bitmaps used be 500x500 pixels. The bitmap should be a
    completely threshholded image, with black defining the area to be used as
    the clipping mask for the cell cluster, or defining the area for a tissue or
    boundary profile.

    Parameters
    ----------------------------
    p                      The instance of the Parameters object used in the simulation.
    desired_bitmap         The designation of the desired bitmap as a string. The string
                           may be:
                           'clipping'   to specify the cluster clipping bitmap is desired
                           'tissue profile 1' to specify the bitmap defining tissue profile 1 is desired
                           'tissue profile 2'     "                           "             2 is desired
                           'boudnary profile 1' to specify the bitmap defining the boundary profile 1

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

    Methods
    ----------------------------
    self.bitmapInit        Loads, initializes, and creates an interpolation object from the bitmap
    self.clipPoints        Takes a set of points as parameters and returns the list of points in the bitmap mask
    """

    def __init__(self, p, desired_bitmap,xmin,xmax,ymin,ymax):
        self.bitmapInit(p, desired_bitmap,xmin,xmax,ymin,ymax)

    def bitmapInit(self,p,desired_bitmap,xmin,xmax,ymin,ymax, threshhold_val = 0.0):
        """
        Finds the appropriate bitmap from the library,
        initializes file loading directory, and loads the
        file to a threshholded matrix.

        Parameters
        ----------------------------
        p                   An instance of the parameters object

        desired_bitmap      This is a string.
                            BETSE has a code for bitmap designation in the
                            model. 'clipping' means the bitmap will be used
                            to clip the cell cluster to a final shape.
                            'tissue profile 1', or 'boundary profile 1'
                            can be used to indicate the bitmap will be used to
                            define the appropriate tissue and boundary profiles, respectively.

        threshhold_val      The value of the pixel to threshhold to. Greyscale runs from 0.0 (black)
                            to 255 (white). The default is to get black pixels (0.0).
        """

        # Absolute or relative path of the bitmap to be loaded.
        self.bitmapFile = p.bitmap_profiles[desired_bitmap]

        # If this is a relative path, convert this into an absolute path
        # relative to the directory containing the source configuration file.
        if paths.is_relative(self.bitmapFile):
            self.bitmapFile = paths.join(p.config_dirname, self.bitmapFile)

        # If this bitmap does *NOT* exist, raise an exception.
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

    def clipPoints(self,point_list_x,point_list_y):

        """
        Uses the self.clipping_function defined for the bitmap
        to screed a list/vector of points, returning those
        that are within the bitmap's colored area.

        Parameters
        -----------
        point_list_x        list or numpy vector of x coordinates of points
        point_list_y        list or numpy vector of y coordinates of points

        Creates
        ------------
        self.good_points    the points falling within the colored area of the bitmap

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








