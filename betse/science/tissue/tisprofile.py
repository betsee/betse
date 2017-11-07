#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level classes aggregating all parameters pertaining to tissue and cut
profiles.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta
from betse.exceptions import BetseSimConfigException
from betse.science.tissue.picker import tispickcls
from betse.science.tissue.picker.tispickcls import TissuePickerABC
from betse.util.type.types import type_check, MappingType

# ....................{ SUPERCLASS                         }....................
class TissueABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all tissue-centric profile classes.

    Instances of this class assign a subset of all cells matching
    subclass-specific criteria (e.g., explicit indexing, randomized selection,
    spatial location) to the corresponding tissue profile.

    Attributes
    ----------------------------
    name : str
        Unique name of the current profile.
    z_order : int
        1-based order in which this profile will be plotted with respect to all
        other plotted profiles. Specifically, profiles with larger z order will
        be plotted over profiles with smaller z order.
    picker : TissuePickerABC
        Object assigning a subset of the cell population to this profile.
    '''

    # ..................{ CONCRETE                           }..................
    @type_check
    def __init__(self, name: str, z_order: int, picker: TissuePickerABC) -> None:
        '''
        Initialize this tissue.

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

        # Classify the passed parameters.
        self.name = name
        self.z_order = z_order
        self.picker = picker

# ....................{ SUBCLASSES                         }....................
#FIXME: Actually do something here.
class TissueProfile(TissueABC):
    '''
    Parameters associated with a subset of the cell population.

    Since these parameters are typically derived from real-world observation of
    biological tissue, the spacialization of these parameters collectively
    defines a "tissue" and is hence referred to as a **tissue profile**.
    '''

    pass


#FIXME: Actually do something here.
class TissueCut(TissueABC):
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
