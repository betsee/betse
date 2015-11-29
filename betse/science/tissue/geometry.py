#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# ....................{ IMPORTS                            }....................
from abc import ABCMeta
from betse.util.path import files, paths

# ....................{ BASE                               }....................
class Geometry(object, metaclass = ABCMeta):
    '''
    Abstract base class of all geometry descriptor classes.

    Geometry descriptors are objects describing which subset of the cell
    population will be assigned to a particular tissue profile.
    '''
    pass

class GeometryUbiquity(Geometry):
    '''
    Ubiquitous (i.e., all-inclusive) geometry descriptor.

    This descriptor unconditionally assigns _all_ cells to a particular tissue
    profile.
    '''
    pass

# ....................{ BITMAP                             }....................
class GeometryBitmap(Geometry):
    '''
    Bitmap-specific geometry descriptor.

    This descriptor assigns all cells residing inside a bitmap's colored pixel
    area to a particular tissue profile.

    Attributes
    ----------------------------
    filename : str
        Absolute path of this bitmap.
    '''

    def __init__(self, filename, dirname):
        '''
        Initialize this descriptor.

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
        assert isinstance(filename, str), '{} not a string.'.format(filename)
        assert isinstance( dirname, str), '{} not a string.'.format(dirname)

        # If this is a relative path, convert this into an absolute path
        # relative to the directory containing the source configuration file.
        if paths.is_relative(filename):
            filename = paths.join(dirname, filename)

        # If this absolute path is *NOT* an existing file, raise an exception.
        files.die_unless_file(filename)

        # Persist this path.
        self.filename = filename

# ....................{ RANDOM                             }....................
class GeometryRandom(Geometry):
    '''
    Randomized geometry descriptor.

    This descriptor randomly assigns a predetermined percentage of the cell
    population to a particular tissue profile.

    Attributes
    ----------------------------
    percentage : float
        Percentage of the cell population to be assigned to this profile. For
        sanity, this should be in the range `[0,0, 100.0]`.
    '''

    def __init__(self, percentage):
        '''
        Initialize this descriptor.

        Parameters
        ----------------------------
        percentage : float
            See the class docstring.
        '''
        assert isinstance(percentage, (int, float)),\
            '{} not an integer or float.'.format(percentage)
        self.percentage = percentage
