#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Cell cluster profile** (i.e., collection of parameters associated with a
subset of the cells in a cell cluster) functionality.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta
from betse.science.tissue.picker.tispickcls import TissuePickerABC
from betse.util.type.types import type_check, MappingType

# ....................{ SUPERCLASS                         }....................
class CellsProfileABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all **cell cluster profile** (i.e., collection of
    parameters associated with a subset of the cells in a cell cluster)
    subclasses.

    Instances of this class assign a subset of all cells matching
    subclass-specific criteria (e.g., explicit indexing, randomized selection,
    spatial location) to the corresponding tissue profile.

    Attributes
    ----------
    name : str
        Unique name of the current profile.
    z_order : int
        1-based order in which this profile will be plotted with respect to all
        other plotted profiles. Specifically, profiles with larger z order will
        be plotted over profiles with smaller z order.
    picker : TissuePickerABC
        Object assigning a subset of the cell population to this profile.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, name: str, z_order: int, picker: TissuePickerABC) -> None:
        '''
        Initialize this cell cluster region.

        Parameters
        ----------
        name : str
            Unique name of this profile.
        z_order : int
            1-based order in which this profile will be plotted with respect to
            all other plotted profiles. Specifically, profiles with larger z
            order will be plotted over profiles with smaller z order.
        picker : TissuePickerABC
            Object assigning a subset of the cell population to this profile.
        '''

        # Classify all passed parameters.
        self.name = name
        self.z_order = z_order
        self.picker = picker

# ....................{ SUBCLASSES                         }....................
class TissueProfile(CellsProfileABC):
    '''
    **Tissue profile** (i.e., cell cluster region unconditionally assigning all
    cells in this region the same initial parameters).

    Since these parameters are typically derived from real-world observation of
    biological tissue, the spacialization of these parameters collectively
    defines a "tissue" and is hence referred to as a **tissue profile**.

    Attributes
    ----------
    is_gj_insular : bool
        ``True`` only if the gap junctions of all cells originating in this
        tissue are **insular** (i.e., prevented from connecting to cells in
        other tissues), implying these gap junctions to be strictly
        intra-tissue.
    mem_diffusion_name_to_const : MappingType
        Dictionary mapping from the name of each membrane diffusion constant
        applied to all cells in this tissue to a floating point number
        specifying that constant.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        is_gj_insular: bool,
        mem_diffusion_name_to_const : MappingType,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this tissue profile.

        Parameters
        ----------
        is_gj_insular : bool
            ``True`` only if the gap junctions of all cells originating in this
            tissue are **insular** (i.e., prevented from connecting to cells in
            other tissues), implying these gap junctions to be strictly
            intra-tissue.
        mem_diffusion_name_to_const : MappingType
            Dictionary mapping from the name of each membrane diffusion constant
            applied to all cells in this tissue to a floating point number
            specifying that constant.

        All remaining parameters are passed as is to the superclass constructor.
        '''

        # Initialize our superclass with all remaining parameters.
        super().__init__(*args, **kwargs)

        # Classify all passed parameters.
        self.is_gj_insular = is_gj_insular
        self.mem_diffusion_name_to_const = mem_diffusion_name_to_const


#FIXME: Actually do something here.
class CutProfile(CellsProfileABC):
    '''
    Profile identifying all cells to be permanently removed by a cutting
    event subsequently triggered during the current tissue simulation.

    There exists a many-to-one relation between cut profiles and cutting events.
    That is to say, each cutting event references zero or more cut profiles.
    While unwieldy, this disconnection permits cutting event (but _not_ cut
    profile) parameters to be modified without requiring simulation reseeding or
    reinitialization. Welcome to the code jungle.
    '''

    pass
