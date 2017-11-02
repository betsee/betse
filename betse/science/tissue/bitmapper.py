#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.exceptions import BetseSimException
from betse.util.type.types import type_check, NumericTypes, SequenceTypes
from scipy import interpolate
from scipy import misc
from scipy.spatial import ConvexHull

# ....................{ CLASSES                            }....................
#FIXME: Generalize this class to support images of arbitrary (possibly
#non-square) dimensions.
#FIXME: Rename this class to "TissueImageMask".
#FIXME: Document all undocumented attributes.
#FIXME: Document why this and the intrinsically related "TissuePickerBitmap"
#class are separate -- notably, to reduce space consumption by instantiating
#instances of this class as local variables of methods in the latter class.
class BitMapper(object):
    '''
    Object finding, loading, and converting a passed bitmap into a SciPy-based
    interpolation function.

    Attributes
    ----------
    clipping_matrix : ndarray
        Numpy matrix defining this bitmap's threshholded image.
    clipping_function : func
        SciPy-based interpolation function accepting an ``(x, y)`` point and
        returning ``1.0`` if that point resides outside this bitmap's colored
        pixel area or ``0.0`` otherwise. This function permits callers to filter
        a passed set of points in the space defined by :attr:`p.wsx` for the
        subset residing within this area.
    clipping_function_fast : func
        Fast variant of :attr:`clipping_function`, but otherwise sharing the
        same API.

    Attributes (clip_points)
    ----------
    The following attributes are guaranteed to be ``None`` until the
    :meth:`clipPoints` method is externally called.

    good_inds : ndarray
        #FIXME: Document us up the `BitMapper` bomb.
    good_points : ndarray
        Numpy matrix listing all points ``(x, y)`` residing inside this bitmap's
        colored area.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,

        #FIXME: Nonsense. Reduce this parameter to simply:
        #    filename: str,
        #This eliminates the circular dependency and simplifies life for all.
        # Avoid circular import dependencies.
        bitmap_matcher: 'betse.science.tissue.tissuepick.TissuePickerBitmap',

        #FIXME: Rename these parameters to "x_min", "x_max", etc.
        xmin: NumericTypes,
        xmax: NumericTypes,
        ymin: NumericTypes,
        ymax: NumericTypes,
    ) -> None:
        '''
        Load, initialize, and create a threshholded interpolation matrix from
        the passed image mask file.

        Constraints
        ----------
        This image should ideally be completely threshholded, with black pixels
        defining either:

        * The area to be used as the clipping mask for the cell cluster.
        * The area for a tissue or boundary profile.

        This image *must* additionally:

        * Be square (i.e., have equal pixel dimensions). For consistency, image
          dimensions of 500x500 pixels are recommended.
        * Not contain an alpha transparency layer, a SciPy-based constraint and
          hence *not* amenable to change.

        Parameters
        ----------
        bitmap_matcher : TissuePickerBitmap
            Low-level BETSE-specific object describing this bitmap.
        xmin : NumericTypes
            Minimum X coordinate accepted by these interpolation functions.
        xmax : NumericTypes
            Maximum X coordinate accepted by these interpolation functions.
        ymin : NumericTypes
            Minimum Y coordinate accepted by these interpolation functions.
        ymax : NumericTypes
            Maximum Y coordinate accepted by these interpolation functions.
        '''

        # Nullify all uninitialized instance variables for safety.
        self.good_points = None
        self.good_inds = None

        #FIXME: Shift this and the following bitmap validation functionality
        #into a new top-level utility function of this submodule: e.g.,
        #
        #    @type_check
        #    def load_bitmap_mask(filename: str) -> ndarray:
        #FIXME: Additionally validate this image to *NOT* contain an alpha
        #transparency layer, a SciPy-based constraint.

        # Load this bitmap as a flattened (i.e., grayscale) Numpy array.
        bitmap = misc.imread(bitmap_matcher.filename, flatten=1)
        # bitmap = np.asarray(bitmap, dtype=np.int)

        # If this bitmap has non-square dimensions, raise an exception.
        if bitmap.shape[0] != bitmap.shape[1]:
            #FIXME: Consider raising a less ambiguous "BetseBitmapException".
            raise BetseSimException(
                'Bitmap "{}" dimensions non-square '
                '(i.e., width {} != height {}).'.format(
                    bitmap_matcher.filename,
                    bitmap.shape[0],
                    bitmap.shape[1],
                ))

        # Find the black pixels. (This is a really basic threshholding!)
        point_inds = (bitmap != 255).nonzero()

        # New matrix the same shape as the image and set values to 0 or 1.
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
        self.clipping_function = interpolate.interp2d(
            xpts, ypts, self.clipping_matrix)
        self.clipping_function_fast = interpolate.RectBivariateSpline(
            xpts, ypts, self.clipping_matrix)

        # Store some additional information relating to bounding polygon of the
        # clipping image.
        Xclip, Yclip = np.meshgrid(xpts, ypts)

        indsN = (self.clipping_matrix == 1.0).nonzero()

        ptsx = Xclip[indsN].ravel()
        ptsy = Yclip[indsN].ravel()

        # points stack
        ppts = np.column_stack((ptsx, ptsy))

        # Calculate convex hull of shape:
        bclip = ConvexHull(ppts)

        # calculate
        bverts = bclip.vertices

        # store boundary points of the clipping poly representing the image shape:
        bx = ptsx[bverts]
        by = ptsy[bverts]

        # store the points of the clipping poly curve:
        self.clipcurve = np.column_stack((bx, by))

    # ..................{ CLIPPERS                           }..................
    #FIXME: Rename this method to clip_points().
    #FIXME: Rename the "point_list_x" and "point_list_y" parameters to
    #"points_x" and "points_y" (respectively).
    @type_check
    def clipPoints(
        self,
        point_list_x: SequenceTypes,
        point_list_y: SequenceTypes,
    ) -> None:
        '''
        Initialize the :attr:`good_points` and :attr:`good_inds` attributes of
        this object to the subset of the passed list or vector of points
        residing in this bitmap's colored area by calling the clipping function
        previously initialized for this bitmap.

        Parameters
        -----------
        point_list_x : SequenceTypes
            Sequence of X coordinates of points.
        point_list_y : SequenceTypes
            Sequence of Y coordinates of points.
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
