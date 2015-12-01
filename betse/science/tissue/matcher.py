#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# ....................{ IMPORTS                            }....................
from abc import ABCMeta
from betse.exceptions import BetseExceptionParameters
from betse.util.path import files, paths
from betse.util.type import types

# ....................{ BASE                               }....................
class TissueMatcher(object, metaclass = ABCMeta):
    '''
    Abstract base class of all tissue matching classes.

    Instances of this class assign a subset of all cells and/or cell membranes
    matching subclass-specific criteria (e.g., explicit indexing, randomized
    selection, spatial location) to the corresponding tissue profile.
    '''
    pass

class TissueMatcherEverything(TissueMatcher):
    '''
    Ubiquitous (i.e., all-inclusive) tissue matcher.

    This matcher unconditionally assigns _all_ cells and/or cell membranes to
    the corresponding tissue profile.
    '''
    pass

# ....................{ BITMAP                             }....................
class TissueMatcherBitmap(TissueMatcher):
    '''
    Bitmap-specific tissue matcher.

    This matcher assigns all cells residing inside a bitmap's colored pixel area
    to the corresponding tissue profile.

    Attributes
    ----------------------------
    filename : str
        Absolute path of this bitmap.
    '''

    def __init__(self, filename, dirname):
        '''
        Initialize this matcher.

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
        '''
        assert types.is_str(filename), types.assert_is_nonstr(filename)
        assert types.is_str( dirname), types.assert_is_nonstr( dirname)

        # If this is a relative path, convert this into an absolute path
        # relative to the directory containing the source configuration file.
        if paths.is_relative(filename):
            filename = paths.join(dirname, filename)

        # If this absolute path is *NOT* an existing file, raise an exception.
        files.die_unless_file(filename)

        # Persist this path.
        self.filename = filename

# ....................{ INDICES                            }....................
class TissueMatcherIndices(TissueMatcher):
    '''
    Indices-specific tissue matcher.

    This matcher assigns all cells with the listed indices to the corresponding
    tissue profile.

    Attributes
    ----------------------------
    indices : collections.Sequence
        Sequence (e.g., list, tuple) of the indices of all cells to be assigned
        to this tissue profile.
    '''

    def __init__(self, indices):
        '''
        Initialize this matcher.

        Parameters
        ----------------------------
        indices : collections.Sequence
            See the class docstring.
        '''
        assert types.is_sequence_nonstr(indices),\
            types.assert_is_not_sequence_nonstr(indices)
        self.indices = indices

# ....................{ RANDOM                             }....................
class TissueMatcherRandom(TissueMatcher):
    '''
    Randomized cell matcher.

    This matcher randomly assigns a predetermined percentage of all cells and/or
    cell membranes to the corresponding tissue profile.

    Attributes
    ----------------------------
    percentage : float
        Percentage of the cell population to be assigned to this profile. For
        sanity, this should be in the range `[0,0, 100.0]`.
    '''

    def __init__(self, percentage):
        '''
        Initialize this matcher.

        Parameters
        ----------------------------
        percentage : float
            See the class docstring.
        '''
        assert types.is_numeric(percentage),\
            types.assert_is_not_numeric(percentage)

        # If this is not a valid percentage, raise an exception. This is
        # important enough to always test rather than defer to assertions.
        if not 0.0 <= percentage <= 100.0:
            raise BetseExceptionParameters(
                '{} not in the range [0.0, 100.0].'.format(percentage))

        self.percentage = percentage
