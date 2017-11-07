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
from betse.util.type.types import type_check, NumericTypes, SequenceTypes

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

        #FIXME: Remove this parameter for brevity. It's *NEVER* used anywhere.
        p:     'betse.science.parameters.Parameters',

        #FIXME: Refactor as follows:
        #
        #* Make this parameter mandatory rather than optional. Fortunately, all
        #  calls to this method in the codebase *ALWAYS* pass this parameter.
        #* Rename to "is_ecm_handled".
        #* Invert all boolean logic referencing this boolean. Double negatives
        #  make my weary head blisters ache.
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
