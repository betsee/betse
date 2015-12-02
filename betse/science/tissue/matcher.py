#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# ....................{ IMPORTS                            }....................
from abc import ABCMeta
from betse.exceptions import (
    BetseExceptionMethodUnimplemented, BetseExceptionParameters)
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

class TissueMatcherAll(TissueMatcher):
    '''
    All-inclusive tissue matcher.

    This matcher unconditionally matches _all_ cells and/or cell membranes.
    '''
    pass

# ....................{ BITMAP                             }....................
class TissueMatcherBitmap(TissueMatcher):
    '''
    Bitmap-specific tissue matcher.

    This matcher matches all cells residing inside the colored pixel area
    defined by an associated bitmap file.

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
            Absolute or relative path of the desired bitmap. If relative (i.e.,
            _not_ prefixed by a directory separator), this path will be
            canonicalized into an absolute path relative to the directory
            containing the current simulation's configuration file.
        dirname : str
            Absolute path of the directory containing the path of the bitmap to
            be loaded (i.e., `filename`). If that path is relative, that path
            will be prefixed by this path to convert that path into an absolute
            path; otherwise, this path will be ignored.
        '''
        assert types.is_str(filename), types.assert_nonstr(filename)
        assert types.is_str( dirname), types.assert_nonstr( dirname)

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

    This matcher matches all cells with the listed indices.

    Attributes
    ----------------------------
    indices : collections.Sequence
        Sequence (e.g., list, tuple) of the indices of all cells to be matched.
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
            types.assert_not_sequence_nonstr(indices)
        self.indices = indices

# ....................{ RANDOM                             }....................
class TissueMatcherRandom(TissueMatcher):
    '''
    Randomized cell matcher.

    This matcher randomly matches a percentage of cells.

    Attributes
    ----------------------------
    percentage : {int, float}
        Percentage of the total cell population to be randomly matched as an
        integer or float in the range `[0,0, 100.0]`.
    '''

    def __init__(self, percentage):
        '''
        Initialize this matcher.

        Parameters
        ----------------------------
        percentage : {int, float}
            See the class docstring.
        '''
        assert types.is_numeric(percentage),\
            types.assert_not_numeric(percentage)

        # If this is not a valid percentage, raise an exception. This is
        # important enough to always test rather than defer to assertions.
        if not 0.0 <= percentage <= 100.0:
            raise BetseExceptionParameters(
                '{} not in the range [0.0, 100.0].'.format(percentage))

        self.percentage = percentage
