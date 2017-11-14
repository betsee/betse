#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Class hierarchy collectively implementing various methods for assigning a subset
of the total cell population to the corresponding tissue profile.
'''

#FIXME: *UGH.* The scipy.misc.imread() has been officially deprecated and will
#be removed as of SciPy 1.2.0. SciPy now recommends the "imageio" package as a
#replacement dependency providing an equivalent function. We therefore need to
#perform the following:
#
#* Select a minimum version of "imageio" to be required. That would appear to be
#  "imageio >= 2.0.1", which is approximately a year old and a clear demarcation
#  point. Note that imageio 2.0.0 was released a day prior to 2.0.1 and that no
#  distributions (including Gentoo) appear to provide it, suggesting that
#  imageio 2.0.0 was a failed release best forgot. "imageio >= 2.0.1" it is!
#* Add "imageio" as a mandatory dependency to "betse.metadeps".
#* Remove "Pillow" as a mandatory dependency from "betse.metadeps". This
#  dependency was only required for the scipy.misc.imread(). Since this
#  dependency is extremely heavyweight (by compare to "imageio"), this is good.
#* Replace all existing calls to scipy.misc.imread() with calls to
#  imageio.imread(). Although the function name remains the same, we'll need to
#  manually validate that the API remains the same as well.
#* Revise all of the following to install "imageio" in lieue of Pyllow:
#  * "README.md" documentation. Since Anaconda does *NOT* provide an official
#    "conda" package for "imageio", end users will need to manually install
#    "imageio" from a third-party channel. Trivial, but tedious. We should
#    probably consider the following:
#    * Decide which channel to leverage. Cursor research strongly suggests the
#      "conda-forge" channel, whose "imageio" package boasts an overwhelming
#      143,969 total downloads (!). To do so, simply run:
#      $ conda install -c conda-forge imageio
#    * Shift a portion of the "bin/install" directory from BETSEE into BETSE
#      itself. Perhaps the portion that is specific only to BETSE? This will
#      then require that the BETSEE installation script internally:
#      * Perform a "wget" command to download the remote BETSE installation
#        script into a local directory... somewhere. Say, "tmp"?
#      * Run this now local BETSE installation script with the same directories
#        as were passed to the BETSEE installation script. *sigh*
#    * Define a new BETSE-specific "bin/install/conda.bash" shell script
#      performing an Anaconda-based installation -- complete with Git repository
#      cloning stripped straight out of our existing Ubuntu installer. This
#      script should probably *NOT* attempt to install Anaconda but simply fail
#      with a fatal error in the absence of "conda" in the current ${PATH}.
#    * Document both our Ubuntu >= 16.04 and Anaconda installers (in that order,
#      as Ubuntu users would certainly prefer a process integrating with their
#      existing package manager) at the head of "README.md".
#    * Document more explicitly that Windows users should *ABSOLUTELY* install
#      Bash-on-Ubuntu-on-Windows. This should literally be the first item in
#      this Markdown list of installation instructions, as it then permits
#      Windows users to install BETSE via either:
#      * The Ubuntu shell script (strongly recommended).
#      * The Anaconda shell script (feasible but less recommended).
#  * "doc/md/INSTALL.md" documentation. Note that, although recent versions of
#    Ubuntu provide the requisite "python3-imageio" package, that this package
#    has *NOT* been backported to Ubuntu 16.04 (Xenial Xerus). Ergo, this
#    package will need to be manually installed via "pip3" here.
#  * "bin/install" scripts.
#  * Gentoo packages. Helpfully, Gentoo provides an "imageio" package. Yes!

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.exceptions import BetseImageException
from betse.science.tissue.picker.tispickcls import TissuePickerABC
from betse.util.path import files, pathnames
from betse.util.type.types import type_check, NumericTypes, SequenceTypes
from numpy import ndarray
from scipy import interpolate, misc
from scipy.spatial import ConvexHull

# ....................{ CLASSES ~ utility                  }....................
#FIXME: Generalize this class to support images of arbitrary (possibly
#non-square) dimensions.
#FIXME: Document all undocumented attributes.
#FIXME: Document why this and the intrinsically related "TissuePickerImage"
#class are separate -- notably, to reduce space consumption by instantiating
#instances of this class as local variables of methods in the latter class.
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
        x_min: NumericTypes,
        x_max: NumericTypes,
        y_min: NumericTypes,
        y_max: NumericTypes,
    ) -> None:
        '''
        Load, initialize, and create a threshholded interpolation matrix from
        the passed image mask file.

        Caveats
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
        filename : str
            Absolute or relative filename of this image.
        x_min : NumericTypes
            Minimum X coordinate accepted by these interpolation functions.
        x_max : NumericTypes
            Maximum X coordinate accepted by these interpolation functions.
        y_min : NumericTypes
            Minimum Y coordinate accepted by these interpolation functions.
        y_max : NumericTypes
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

        # If this image does *NOT* exist, raise an exception.
        files.die_unless_file(filename)

        #FIXME: Additionally validate this image to *NOT* contain an alpha
        #transparency layer, a SciPy-based constraint.

        # Load this bitmap as a flattened (i.e., grayscale) Numpy array.
        bitmap = misc.imread(filename, flatten=1)
        # bitmap = np.asarray(bitmap, dtype=np.int)

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
