#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level classes aggregating all parameters pertaining to simulation events.
'''

# ....................{ IMPORTS                            }....................
from betse.science.tissue.event.tiseveabc import SimEventSpikeABC
# from betse.util.io.log import logs
from betse.util.type.types import type_check, SequenceTypes

# ....................{ SUBCLASSES                         }....................
#FIXME: Refactor the TissueHandler.removeCells() method into a new
#fire() method of this subclass.
class SimEventCut(SimEventSpikeABC):
    '''
    **Cutting event** (i.e., event removing a region of the current cluster at
    some time step during the simulation phase).

    Attributes
    ----------
    profile_names : list
        List of the names of all applicable cut profiles, each describing a
        subset of the cell population to be removed.
    '''

    # ..................{ PUBLIC                             }..................
    @type_check
    def __init__(
        self,
        profile_names: SequenceTypes,
        *args, **kwargs
    ) -> None:

        # Initialize our superclass with all remaining parameters..
        super().__init__(*args, **kwargs)

        # Classify all passed parameters.
        self.profile_names = profile_names
