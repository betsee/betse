#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level classes aggregating all parameters pertaining to simulation events.
'''

# ....................{ IMPORTS                           }....................
from betse.science.tissue.event.tiseveabc import SimEventSpikeABC
# from betse.util.io.log import logs
from betse.util.type.types import type_check, SequenceTypes

# ....................{ SUBCLASSES                        }....................
#FIXME: Refactor the TissueHandler.removeCells() method into a new
#fire() method of this subclass.
class SimEventCut(SimEventSpikeABC):
    '''
    **Cutting event** (i.e., event removing one or more cells of the current
    cell cluster at some time step during the simulation phase).

    Attributes
    ----------
    '''

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def __init__(
        self,
        *args, **kwargs
    ) -> None:

        # Initialize our superclass with all remaining parameters..
        super().__init__(*args, **kwargs)
