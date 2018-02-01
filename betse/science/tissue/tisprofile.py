#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
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
        ``True`` only if gap junctions originating at cells in this tissue are
        **insular** (i.e., prevented from connecting to cells in other tissues),
        implying these gap junctions to be strictly intra-tissue.

    Attributes (Membrane Diffusion)
    ----------
    Dm_Na : float
        Sodium (Na+) membrane diffusion constant in m2/s.
    Dm_K : float
        Potassium (K+) membrane diffusion constant in m2/s.
    Dm_Cl : float
        Chloride (Cl-) membrane diffusion constant in m2/s.
    Dm_Ca : float
        Calcium (Ca2+) membrane diffusion constant in m2/s.
    Dm_H : float
        Hydrogen (H+) membrane diffusion constant in m2/s.
    Dm_M : float
        Charge balance anion (M-) membrane diffusion constant in m2/s.
    Dm_P : float
        Protein (P-) membrane diffusion constant in m2/s.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        is_gj_insular: bool,
        Dm_Na: float,
        Dm_K: float,
        Dm_Cl: float,
        Dm_Ca: float,
        Dm_M: float,
        Dm_P: float,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this tissue profile.

        Parameters
        ----------
        is_gj_insular : bool
            ``True`` only if gap junctions originating at cells in this tissue
            are **insular** (i.e., prevented from connecting to cells in other
            tissues), implying these gap junctions to be strictly intra-tissue.
        Dm_Na : float
            Sodium (Na+) membrane diffusion constant in m2/s.
        Dm_K : float
            Potassium (K+) membrane diffusion constant in m2/s.
        Dm_Cl : float
            Chloride (Cl-) membrane diffusion constant in m2/s.
        Dm_Ca : float
            Calcium (Ca2+) membrane diffusion constant in m2/s.
        Dm_H : float
            Hydrogen (H+) membrane diffusion constant in m2/s.
        Dm_M : float
            Charge balance anion (M-) membrane diffusion constant in m2/s.
        Dm_P : float
            Protein (P-) membrane diffusion constant in m2/s.

        All remaining parameters are passed as is to the superclass constructor.
        '''

        # Initialize our superclass with all remaining parameters.
        super().__init__(*args, **kwargs)

        # Classify all passed parameters.
        self.is_gj_insular = is_gj_insular
        self.Dm_Na = Dm_Na
        self.Dm_K  = Dm_K
        self.Dm_Cl = Dm_Cl
        self.Dm_Ca = Dm_Ca
        self.Dm_M  = Dm_M
        self.Dm_P  = Dm_P


#FIXME: Actually do something here.
class CutProfile(CellsProfileABC):
    '''
    **Cut profile** (i.e., cell cluster region to be permanently removed by a
    cutting event triggered during the simulation phase).

    There exists a many-to-one relation between cut profiles and cutting events.
    That is to say, each cutting event references zero or more cut profiles.
    While unwieldy, this disconnection permits cutting event (but *not* cut
    profile) parameters to be modified without requiring simulation reseeding or
    reinitialization. Welcome to the code jungle.
    '''

    pass
