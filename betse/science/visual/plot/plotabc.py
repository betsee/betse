#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all post-simulation plot subclasses.
'''

# ....................{ IMPORTS                            }....................
from betse.science.phase.phasecls import SimPhase
from betse.science.visual.visabc import VisualCellsABC
from betse.util.type.types import type_check

# ....................{ CLASSES                            }....................
class PlotCellsAfterSolving(VisualCellsABC):
    '''
    Abstract base class of all post-simulation cell plot subclasses, plotting
    simulation data over the cell cluster _after_ rather than _during_
    simulation modelling.
    '''

    @type_check
    def __init__(
        self,
        phase: SimPhase,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this post-simulation plot.

        Parameters
        ----------
        phase: SimPhase
            Current simulation phase.

        See the superclass `__init__()` method for all remaining parameters.
        '''

        # Initialize our superclass.
        super().__init__(
            *args,

            # Pass this simulation phase as is to our superclass.
            phase=phase,

            # Save and show this post-simulation plot only if this configuration
            # enables doing so.
            is_save=phase.p.plot.is_after_sim_save,
            is_show=phase.p.plot.is_after_sim_show,

            # Save all post-simulation plots to the same parent directory.
            save_dir_parent_basename='plot',

            # Pass all remaining arguments as is to our superclass.
            **kwargs
        )
