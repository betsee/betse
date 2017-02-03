#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Class hierarchy collectively implementing various methods for assigning a subset
of the total cell population to the corresponding tissue profile.
'''

# ....................{ IMPORTS                            }....................
import random
from abc import ABCMeta, abstractmethod
from betse.exceptions import BetseSimConfigException
from betse.science.math import toolbox
from betse.util.path import files, paths
from betse.util.type.types import (
    type_check, MappingType, NumericTypes, SequenceTypes)

# ....................{ SUPERCLASS                         }....................
class TissuePickerABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all tissue matching classes.

    Instances of this class assign a subset of all cells matching
    subclass-specific criteria (e.g., explicit indexing, randomized selection,
    spatial location) to the corresponding tissue profile.
    '''

    # ..................{ ABSTRACT                           }..................
    @abstractmethod
    def get_cell_indices(
        self,
        cells: 'betse.science.cells.Cells',
        p:     'betse.science.parameters.Parameters',
        ignoreECM: bool = False,
    ) -> SequenceTypes:
        '''
        One-dimensional Numpy array of the indices of all cells in the passed
        cell cluster selected by this tissue picker.

        Parameters
        ---------------------------------
        cells : Cells
            Current cell cluster.
        p : Parameters
            Current simulation configuration.
        ignoreECM : bool
            ``True`` if extracellular spaces are to be ignored; ``False`` if
            extracellular spaces are to be simulated. Defaults to ``False``.

        Returns
        ---------------------------------
        ndarray
            See method synopsis above.
        '''

        pass

# ....................{ SUBCLASSES                         }....................
class TissuePickerAll(TissuePickerABC):
    '''
    All-inclusive tissue picker.

    This picker unconditionally matches *all* cells.
    '''

    # ..................{ GETTERS                            }..................
    @type_check
    def get_cell_indices(
        self,
        cells: 'betse.science.cells.Cells',
        p:     'betse.science.parameters.Parameters',
        ignoreECM: bool = False,
    ) -> SequenceTypes:

        #FIXME: The following logic is reapeated throughout this submodule. DRY!

        # If either not simulating *OR* ignoring extracellular spaces, get cell
        # indices only.
        if ignoreECM:
            target_inds = cells.cell_i
        # Else, simulate extracellular spaces.
        else:
            target_inds = cells.cell_to_mems[cells.cell_i]
            target_inds, _, _ = toolbox.flatten(target_inds)

        return target_inds

# ....................{ SUBCLASSES ~ bitmap                }....................
class TissuePickerBitmap(TissuePickerABC):
    '''
    Bitmap-specific tissue picker.

    This matcher matches all cells residing inside the colored pixel area
    defined by an associated bitmap file.

    Attributes
    ----------------------------
    filename : str
        Absolute path of this bitmap.
    '''

    # ..................{ PUBLIC                             }..................
    @type_check
    def __init__(self, filename: str, dirname: str) -> None:
        '''
        Initialize this tissue picker.

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
            path; else, this path is ignored.
        '''

        # If this is a relative path, convert this into an absolute path
        # relative to the directory containing the source configuration file.
        if paths.is_relative(filename):
            filename = paths.join(dirname, filename)

        # If this absolute path is *NOT* an existing file, raise an exception.
        files.die_unless_file(filename)

        # Persist this path.
        self.filename = filename

    # ..................{ GETTERS                            }..................
    @type_check
    def get_cell_indices(
        self,
        cells: 'betse.science.cells.Cells',
        p:     'betse.science.parameters.Parameters',
        ignoreECM: bool = False,
    ) -> SequenceTypes:

        # Calculate the indices of all cells residing inside this bitmap.
        bitmask = self.get_bitmapper(cells)
        target_inds = bitmask.good_inds

        # If simulating electromagnetism and at least one cell matches...
        if ignoreECM is False and len(target_inds):
            target_inds = cells.cell_to_mems[target_inds]
            target_inds,_,_ = toolbox.flatten(target_inds)

        return target_inds


    @type_check
    def get_bitmapper(self, cells: 'betse.science.cells.Cells'):
        '''
        :class:`BitMapper` object providing the indices of all cells residing
        inside this bitmap.

        Parameters
        ---------------------------------
        cells : Cells
            Current cell cluster.
        '''

        # Avoid circular import dependencies.
        from betse.science.tissue.bitmapper import BitMapper

        # Create and return the desired bitmap. (Note this object is typically
        # large and hence intentionally *NOT* cached as an object attribute.)
        bitmapper = BitMapper(
            self, cells.xmin, cells.xmax, cells.ymin, cells.ymax)
        bitmapper.clipPoints(cells.cell_centres[:,0], cells.cell_centres[:,1])

        return bitmapper

# ....................{ SUBCLASSES ~ indices               }....................
class TissuePickerIndices(TissuePickerABC):
    '''
    Indices-specific tissue picker.

    This matcher matches all cells with the listed indices.

    Attributes
    ----------------------------
    indices : SequenceTypes
        SequenceTypes (e.g., list, tuple) of the indices of all cells to be
        matched.
    '''

    @type_check
    def __init__(self, indices: SequenceTypes) -> None:
        '''
        Initialize this tissue picker.

        Parameters
        ----------------------------
        indices : SequenceTypes
            See the class docstring.
        '''

        self.indices = indices


    @type_check
    def get_cell_indices(
        self,
        cells: 'betse.science.cells.Cells',
        p:     'betse.science.parameters.Parameters',
        ignoreECM: bool = False,
    ) -> SequenceTypes:

        # If either not simulating *OR* ignoring extracellular spaces, get cell
        # indices only.
        if ignoreECM:
            target_inds = self.indices
        # Else if simulating extracellular spaces, get membrane indices as well.
        else:
            target_inds = cells.cell_to_mems[self.indices]
            target_inds,_,_ = toolbox.flatten(target_inds)

        return target_inds

# ....................{ SUBCLASSES ~ random                }....................
class TissuePickerRandom(TissuePickerABC):
    '''
    Randomized cell picker.

    This picker randomly matches a percentage of cells.

    Attributes
    ----------------------------
    percentage : NumericTypes
        Percentage of the total cell population to be randomly matched as an
        integer or float in the range ``[0,0, 100.0]``.
    '''

    @type_check
    def __init__(self, percentage: NumericTypes) -> None:
        '''
        Initialize this tissue picker.

        Parameters
        ----------------------------
        percentage : NumericTypes
            See the class docstring.
        '''

        # If this is not a valid percentage, raise an exception. This is
        # important enough to always test rather than defer to assertions.
        if not 0.0 <= percentage <= 100.0:
            raise BetseSimConfigException(
                '{} not in the range [0.0, 100.0].'.format(percentage))

        # Classify this parameter.
        self.percentage = percentage


    @type_check
    def get_cell_indices(
        self,
        cells: 'betse.science.cells.Cells',
        p:     'betse.science.parameters.Parameters',
        ignoreECM: bool = False,
    ) -> SequenceTypes:

        data_length = len(cells.cell_i)
        data_fraction = int((self.percentage/100)*data_length)
        cell_i_copy = cells.cell_i[:]
        random.shuffle(cell_i_copy)
        target_inds_cell = [cells.cell_i[x] for x in range(0,data_fraction)]

        # If either not simulating *OR* ignoring extracellular spaces, get cell
        # indices only.
        if ignoreECM:
            target_inds = target_inds_cell
        # Else if simulating extracellular spaces, get membrane indices.
        else:
            target_inds = cells.cell_to_mems[target_inds_cell]
            target_inds,_,_ = toolbox.flatten(target_inds)

        return target_inds

# ....................{ MAKERS                             }....................
@type_check
def make(
    p: 'betse.science.parameters.Parameters',
    conf: MappingType,
) -> TissuePickerABC:
    '''
    Create and return a concrete instance of this abstract base class from the
    passed picker-specific dictionary and simulation parameters.

    Parameters
    ----------------------------
    p : Parameters
        Current simulation configuration.
    conf : MappingType
        Dictionary describing the type and contents of the tissue picker to be
        created via the following key-value pairs:
        * ``type``, a string enumeration.

    Returns
    ----------------------------
    TissuePickerABC
        Concrete instance of this abstract base class.
    '''

    picker = None
    picker_type = conf['type']

    if picker_type == 'all':
        picker = TissuePickerAll()
    elif picker_type == 'bitmap':
        picker = make_bitmap(p=p, conf=conf['bitmap'])
    elif picker_type == 'indices':
        picker = TissuePickerIndices(conf['indices'])
    elif picker_type == 'random':
        picker = TissuePickerRandom(conf['random'])
    else:
        raise BetseSimConfigException(
            'Tissue picker type "{}"' 'unrecognized.'.format(picker_type))

    return picker


@type_check
def make_bitmap(
    p: 'betse.science.parameters.Parameters',
    conf: MappingType,
) -> TissuePickerBitmap:
    '''
    Create and return an instance of this class from the passed bitmap-specific
    dictionary and simulation parameters.

    Parameters
    ----------------------------
    p : Parameters
        Current simulation configuration.
    conf : MappingType
        Dictionary describing the type and contents of the bitmap tissue picker
        to be created via the following key-value pairs:
        * ``file``, the absolute or relative path of this bitmap's file.

    Returns
    ----------------------------
    TissuePickerBitmap
        Instance of this class.
    '''

    return TissuePickerBitmap(
        filename=conf['file'], dirname=p.config_dirname)
