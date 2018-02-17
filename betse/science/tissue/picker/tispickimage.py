#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Class hierarchy collectively implementing various methods for assigning a subset
of the total cell population to the corresponding tissue profile.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.exceptions import BetseImageException
from betse.lib.pil import pilnumpy
from betse.lib.pil.pilnumpy import ImageModeType
from betse.science.tissue.picker.tispickcls import TissuePickerABC
from betse.util.path import files, pathnames
from betse.util.type.types import type_check, NumericSimpleTypes, SequenceTypes
from numpy import ndarray
from scipy import interpolate
from scipy.spatial import ConvexHull

# ....................{ CLASSES ~ utility                  }....................
#FIXME: Generalize this class to support images of arbitrary (possibly
#non-square) dimensions.
#FIXME: Document all undocumented attributes.
class TissuePickerImageMask(object):
    '''
    Cell profile-specific **image mask** (i.e., low-level object converting an
    well-formatted image into multiple SciPy-based interpolation functions
    clipping passed points to that image's colored pixel area).

    Design
    ----------
    This utility class is principally intended to be instantiated by the
    :class:`TissuePickerImage` class, but can (in theory) be instantiated
    directly by sufficiently careful external callers.

    This utility class is typically large in terms of memory consumption. Hence,
    the :class:`TissuePickerImage` class only locally instantiates instances of
    this class for the duration of method calls rather than permanently
    persisting such instances as instance variables of that class.

    Attributes
    ----------
    clipping_matrix : ndarray
        Numpy matrix defining this threshholded image.
    clipping_function : func
        SciPy-based interpolation function accepting an ``(x, y)`` point and
        returning ``1.0`` if that point resides outside this bitmap's colored
        pixel area or ``0.0`` otherwise. This function permits callers to filter
        a passed set of points in the space defined by :attr:`p.wsx` for the
        subset residing within this area.
    clipping_function_fast : func
        Fast variant of :attr:`clipping_function`, but otherwise sharing the
        same API.
    msize : int
        Size in pixels of each square dimension of this image, equivalent to
        both the width and height of this image.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        filename: str,
        x_min: NumericSimpleTypes,
        x_max: NumericSimpleTypes,
        y_min: NumericSimpleTypes,
        y_max: NumericSimpleTypes,
    ) -> None:
        '''
        Load, initialize, and create a threshholded interpolation matrix from
        the passed image mask file.

        Format
        ----------
        This image should ideally be completely threshholded. Each pure black
        pixel (i.e., with red, green, and blue color components all 0) of this
        image collectively defines either:

        * The area to be used as the clipping mask for the cell cluster.
        * The area for a tissue, cut, or boundary profile.

        This image *must* additionally:

        * Be square (i.e., have equal pixel dimensions). For consistency, image
          dimensions of 500x500 pixels are recommended.
        * Contain *no* alpha transparency layer, a SciPy-based constraint and
          hence *not* amenable to change.

        Parameters
        ----------
        filename : str
            Absolute or relative filename of this image.
        x_min : NumericSimpleTypes
            Minimum X coordinate accepted by these interpolation functions.
        x_max : NumericSimpleTypes
            Maximum X coordinate accepted by these interpolation functions.
        y_min : NumericSimpleTypes
            Minimum Y coordinate accepted by these interpolation functions.
        y_max : NumericSimpleTypes
            Maximum Y coordinate accepted by these interpolation functions.

        Raises
        ----------
        BetseFileException
            If this image does *not* exist.
        BetseImageeException
            If this image has either:
            * Non-square dimensions (i.e., differing width and height).
            * An alpha transparency layer.
        '''

        # Load this bitmap as a grayscale Numpy array, preserving floating point
        # precision across this necessarily lossy reduction.
        bitmap = pilnumpy.load_image(
            filename=filename, mode=ImageModeType.GRAYSCALE_FLOAT)

        # If this bitmap has non-square dimensions, raise an exception.
        if bitmap.shape[0] != bitmap.shape[1]:
            raise BetseImageException(
                'Image "{}" dimensions non-square '
                '(i.e., width {} != height {}).'.format(
                    filename, bitmap.shape[0], bitmap.shape[1]))

        # Find the black pixels. (This is a really basic threshholding!)
        point_inds = (bitmap != 255).nonzero()

        # New matrix the same shape as the image and set values to 0 or 1.
        self.msize = bitmap.shape[0]
        self.clipping_matrix = np.zeros((self.msize, self.msize))
        self.clipping_matrix[point_inds] = 1.0
        self.clipping_matrix = np.flipud(self.clipping_matrix)

        # Create spatial data vectors that span the extent of the cell seeds and
        # match bitmap pixel number.
        xpts = np.linspace(x_min, x_max, self.msize)
        ypts = np.linspace(y_min, y_max, self.msize)

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

    # ..................{ GETTERS                            }..................
    @type_check
    def get_clipped_points_index(
        self,
        points_x: SequenceTypes,
        points_y: SequenceTypes,
    ) -> ndarray:
        '''
        One-dimensional Numpy array of the indices of all passed points that
        reside inside and are thus clipped by this image mask's colored area.

        Parameters
        -----------
        points_x : SequenceTypes
            Sequence of the X coordinates of all points to be clipped.
        points_y : SequenceTypes
            Sequence of the Y coordinates of all points to be clipped.

        Returns
        -----------
        ndarray
            One-dimensional Numpy array of the indices of all clipped points.
        '''

        # List of the indices of all passed points clipped by this method.
        clipped_points_index = []

        # For the index and X and Y coordinate of each passed point...
        for points_index, (x, y) in enumerate(zip(points_x, points_y)):
            # If this point resides inside this image mask, clip this point.
            if self.clipping_function(x,y) != 0.0:
                clipped_points_index.append(points_index)

        # Coerce this possibly empty list into a guaranteed integer Numpy array.
        # Since Numpy has no means of deciding whether an empty list should be
        # coerced into an integer or floating point Numpy array, Numpy sensibly
        # defaults to the latter functionality. While that's typically what one
        # wants, that is *NOT* what callers expect in this case. For caller
        # sanity, Numpy *MUST* be instructed to produce an integer array.
        return np.asarray(clipped_points_index, dtype=int)

# ....................{ CLASSES ~ picker                   }....................
class TissuePickerImage(TissuePickerABC):
    '''
    Image-based tissue picker, matching all cells residing inside the colored
    pixel area defined by an associated on-disk image mask file.

    Attributes
    ----------
    filename : str
        Absolute path of this image mask.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, filename: str, dirname: str) -> None:
        '''
        Initialize this tissue picker.

        Parameters
        ----------
        filename : str
            Absolute or relative filename of this image mask. If relative (i.e.,
            *not* prefixed by a directory separator), this filename will be
            canonicalized into an absolute filename relative to the passed
            absolute dirname.
        dirname : str
            Absolute dirname of the directory containing this image mask. If the
            ``filename`` parameter is relative, that filename will be prefixed
            by this path to convert that filename into an absolute path; else,
            this dirname is ignored.
        '''

        # If this is a relative path, convert this into an absolute path
        # relative to the directory containing the source configuration file.
        if pathnames.is_relative(filename):
            filename = pathnames.join(dirname, filename)

        # If this absolute path is *NOT* an existing file, raise an exception.
        files.die_unless_file(filename)

        # Classify this parameter.
        self.filename = filename

    # ..................{ PICKERS                            }..................
    @type_check
    def pick_cells(
        self,
        cells: 'betse.science.cells.Cells',
        p:     'betse.science.parameters.Parameters',
    ) -> SequenceTypes:

        # Cell profile-specific image mask.
        image_mask = self.get_image_mask(cells)

        # One-dimensional Numpy array of the indices of all cells whose centres
        # reside inside and are thus clipped by this image mask's colored area.
        clipped_points_index = image_mask.get_clipped_points_index(
            points_x=cells.cell_centres[:,0],
            points_y=cells.cell_centres[:,1])

        # Return this array.
        return clipped_points_index

    # ..................{ GETTERS                            }..................
    @type_check
    def get_image_mask(
        self, cells: 'betse.science.cells.Cells') -> TissuePickerImageMask:
        '''
        Cell profile-specific image mask, abstracting the low-level pixels of
        this image file into higher-level Python objects.

        Parameters
        ----------
        cells : Cells
            Current cell cluster.

        Returns
        ----------
        TissuePickerImageMask
            Cell profile-specific image mask implementing this picker.
        '''

        return TissuePickerImageMask(
            filename=self.filename,
            x_min=cells.xmin,
            x_max=cells.xmax,
            y_min=cells.ymin,
            y_max=cells.ymax,
        )
