#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level Numpy-based image facilities.

Design
----------
Most functions defined by this submodule are thin wrappers around similar
functionality provided by an arbitrary third-party dependency. As example, the
:func:`load_image` function:

* Currently wraps the :func:`Pillow.Image.load` method.
* Previously wrapped the :func:`scipy.misc.imread` function, which SciPy 1.0.0
  formally deprecated and SciPy 1.2.0 permanently removed. Indeed, SciPy 1.2.0
  broke backward compatibility with downstream consumers.

This submodule principally exists to mitigate the costs associated with similar
catastrophic upstream API changes in the future. We are ready this time.

Two principal choices exist:

* Pillow, the canonical image I/O library for Python. Pillow is a long-standing
  fork of PIL, whose stable API Pillow has inherited and carefully curated.
* imageio, a higher-level image I/O library for Python wrapping lower-level
  libraries (e.g., Pillow, FreeImage), whose API is required to generalize to
  multiple underlying libraries and is hence slightly less stable and
  feature-filled than Pillow. Moreover, imageio currently fails to support all
  fine-grained configurability required by this application (e.g., the capacity
  to convert images to grayscale).

Since the stable Pillow API provides all requisite features and the unstable
imageio API does *not*, this submodule is currently implemented in terms of
Pillow. Callers should *not* assume this to be the case, however.
'''

#FIXME: Revisit imageio when the following feature request is resolved in full:
#    https://github.com/imageio/imageio/issues/289

# ....................{ IMPORTS                           }....................
from PIL import Image
from betse.lib.numpy import nparray
from betse.util.io.log import logs
from betse.util.path import files, pathnames
from betse.util.type.types import type_check, NoneType, NumpyArrayType
from enum import Enum
from numpy import array

# ....................{ ENUMERATIONS                      }....................
#FIXME: Define the remainder of all Pillow modes. See the URL given below.
class ImageModeType(Enum):
    '''
    Enumeration of all types of **image modes** (i.e., :mod:`PIL`-specific
    identifier defining both the type and depth of all pixels in an image)
    supported by functions defined by this submodule.

    This enumeration wraps the non-typesafe, non-human-readable ``mode`` string
    parameter commonly accepted by :mod:`PIL` functions (e.g.,
    :func:`PIL.Image.new`) with typesafe, human-readable variants. The value of
    each member of this enumeration is the corresponding ``mode`` string.

    Attributes
    ----------
    COLOR_RGB : enum
        Integer-based true color mode *without* alpha transparency, providing
        8-bit red, green, and blue pixels.
    COLOR_RGBA : enum
        Integer-based true color mode *with* alpha transparency, providing
        8-bit red, green, blue, and alpha transparency pixels.
    GRAYSCALE_INT : enum
        Integer-based grayscale mode, providing 8-bit black and white pixels.
        Converting color images *not* already in this mode to this mode reduces
        the resulting images to grayscale according to the ITU-R 601-2 luma
        transform in a lossy integer manner (i.e., ``L = R * 299/1000 + G *
        587/1000 + B * 114/1000`` for each pixel such that the fractional
        portion of each such ``L`` is discarded).
    GRAYSCALE_FLOAT : enum
        Floating point-based grayscale mode, providing 32-bit floating point
        pixels. Converting color images *not* already in this mode to this mode
        reduces the resulting images to grayscale according to the ITU-R 601-2
        luma transform in a lossy floating point manner (i.e., ``L = R *
        299/1000 + G * 587/1000 + B * 114/1000`` for each pixel such that the
        fractional portion of each such ``L`` is preserved). This conversion
        generalizes the integer-based conversion to grayscale performed by
        converting to the :attr:`GRAYSCALE_INT` mode by preserving rather than
        discarding the fractional portion of each luma value. Note that the
        only file format currently supporting images of this mode is TIFF.

    See Also
    ----------
    http://pillow.readthedocs.io/en/3.4.x/handbook/concepts.html#modes
        Most recent official documentation for Pillow modes.
    '''

    COLOR_RGB       = 'RGB'
    COLOR_RGBA      = 'RGBA'
    GRAYSCALE_INT   = 'I'
    GRAYSCALE_FLOAT = 'F'

# ....................{ TYPES                             }....................
ImageModeOrNoneTypes = (ImageModeType, NoneType)
'''
Tuple of the type of all image mode enumeration members *and* of the singleton
``None`` object.
'''

# ....................{ CONVERTERS                        }....................
@type_check
def load_image(
    # Mandatory parameters.
    filename: str,

    # Optional parameters.
    is_signed: bool = True,
    mode: ImageModeOrNoneTypes = None,
) -> NumpyArrayType:
    '''
    Load the raw pixel data from the image with the passed filename into a
    multi-dimensional Numpy array and return this array.

    This array is guaranteed to be at least three-dimensional. Specifically, if
    this image is:

    * Greyscale, this array is three-dimensional such that:
      * The first dimension indexes each row of this image.
      * The second dimension indexes each column of the current row.
      * The third dimension indexes the grayscale value of the current pixel.
    * RGB or RGBA, this array is four-dimensional such that:
      * The first dimension indexes each row of this image.
      * The second dimension indexes each column of the current row.
      * The third dimension is either a 3-tuple ``(R, G, B)`` or 4-tuple
        ``(R, G, B, A)`` indexing each color component of the current pixel.
      * The fourth dimension indexes the value of the current color component.

    Caveats
    ----------
    When attempting to load user-defined images of arbitrary filetype, callers
    should pass the following parameters:

    * ``mode``.
    * ``is_signed`` to ``True``. Since this is (and *always* will be) the
      default, *not* passing this parameter satisfies this suggestion.

    Failure to do so invites subtle issues in computations falsely assuming the
    data type and shape of a returned array to be sane, which is *not* the case
    in common edge cases. Thanks to the heterogeneity of image file formats,
    returned arrays may exhibit anomalous features if any of these paremeters
    are *not* passed as suggested. In particular, the ``is_signed`` parameter
    should typically either *not* be passed or be passed as ``True``.

    Failure to do so instructs Pillow to produce an array with data type
    automatically corresponds to that of the input image. Since most (but *not*
    necessarily all) images reside in the :attr:`ImageModeType.COLOR_RGB` and
    :attr:`ImageModeType.COLOR_RGBA` colour spaces whose three- and
    four-channel pixel data is homogenously constrained onto unsigned bytes,
    most arrays returned by this function when explicitly passed an
    ``is_signed`` parameter of ``False`` will be **unsigned byte arrays**
    (i.e., arrays whose data types are :attr:`np.uint8`).

    Is that a subtle problem? **It is.**

    Python silently coerces scalar types as needed to preserve precision across
    operations that change precision. The canonical example is integer
    division.  In Python, dividing two integers that are *not* simple integer
    multiples of one another implicitly expands precision by producing a real
    number rather than integer (e.g., ``1 / 2 == 0.5`` rather than
    ``1 / 2 == 0``).

    On the other hand:

    * For all **signed Numpy arrays** (i.e., arrays whose data types are
      implicitly signed rather than explicitly unsigned), Numpy silently
      coerces the data types of these arrays as needed to preserve precision
      across precision-modifying operations.
    * For all **unsigned Numpy arrays** (i.e., arrays whose data types are
      explicitly unsigned rather than implicitly signed), Numpy silently
      preserves the unsigned facet of these arrays as needed by wrapping all
      numerical results to the integer range of these unsigned data types, thus
      discarding precision across precision-modifying operations.

    The canonical example is integer addition and substraction applied to
    unsigned byte arrays. Since unsigned bytes are confined to the integer
    range ``[0, 255]``, attempting to perform even seemingly trivial
    computation with unsigned byte arrays silently wraps results exceeding this
    range onto this range. The resulting arrays typically contain so-called
    "garbage data." As the following example shows, applying integer
    subtraction to signed but *not* unsigned Numpy arrays produces expected
    results:

        >>> import numpy as np
        >>> unsigned_garbage = np.array(((1,2), (3,4)), dtype=np.uint8)
        >>> unsigned_garbage[:,0] - unsigned_garbage[:,1]
        ... array([255, 255], dtype=uint8)
        >>> signed_nongarbage = np.array(((1,2), (3,4)))
        >>> signed_nongarbage[:,0] - signed_nongarbage[:,1]
        ... array([-1, -1])

    Design
    ----------
    This utility function is a thin wrapper around a similar function provided
    by some unspecified third-party dependency. This function currently wraps
    the :meth:`PIL.Image.open` method but previously wrapped the:

    * :func:`imageio.imread` function, which failed to expose support for
      colourspace conversion provided by Pillow.
    * :func:`scipy.misc.imread` function, which SciPy 1.0.0 formally deprecated
      and SciPy 1.2.0 permanently killed. Thus, SciPy 1.2.0 broke backward
      compatibility with downstream applications (notably, *this* application)
      requiring that API.

    This utility function principally exists to mitigate the costs associated
    with similar upstream API changes in the future. (We are ready this time.)

    Parameters
    ----------
    filename : str
        Absolute or relative filename of this image.
    is_signed : optional[bool]
        ``True`` only if converting the possibly unsigned array loaded from
        this image into a signed array. Defaults to ``True`` for the reasons
        detailed above. Since explicitly setting this to ``False`` invites
        errors in computations employing the returned array, callers should do
        so *ONLY* where these issues are acknowledged and handled
        appropriately.
    mode : ImageModeOrNoneTypes
        Type and depth of all pixels in the array loaded from this image,
        converted from this image's pixel data according to industry-standard
        image processing transforms implemented by :mod:`PIL`. Note that this
        is *not* the type and depth of all pixels in the input image, which
        :mod:`PIL` implicitly detects and hence requires no explicit
        designation. Defaults to ``None``, in which case no such conversion is
        performed (i.e., this image's pixel data is returned as is).

    Returns
    ----------
    ndarray
        Numpy array loaded from this image.
    '''

    # Log this load attempt.
    logs.log_debug('Loading image "%s"...', pathnames.get_basename(filename))

    # If this image does *NOT* exist, raise an exception.
    files.die_unless_file(filename)

    # In-memory image object loaded from this on-disk image file.
    image = Image.open(filename)

    # If the caller requests a mode conversion *AND* this image is not already
    # of the required mode, do this conversion.
    if mode is not None and mode != image.mode:
        image = image.convert(mode.value)

    # Numpy array converted from this image to be returned.
    image_array = array(image)

    # If converting unsigned to signed arrays, do so.
    if is_signed:
        image_array = nparray.to_signed(image_array)

    # Return this array.
    return image_array
