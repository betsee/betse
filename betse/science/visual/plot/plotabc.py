#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all post-simulation plot subclasses.
'''

# ....................{ IMPORTS                            }....................
from betse.science.parameters import Parameters
from betse.science.visual.visualabc import VisualCellsABC
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
        p: Parameters, #'betse.science.parameters.Parameters',
        *args, **kwargs
    ) -> None:
        '''
        Initialize this post-simulation plot.

        Parameters
        ----------
        p : Parameters
            Current simulation configuration.

        See the superclass `__init__()` method for all remaining parameters.
        '''

        # Initialize our superclass.
        super().__init__(
            # Pass this simulation configuration as is to our superclass.
            p=p,

            # Save and show this post-simulation plot only if this configuration
            # enables doing so.
            is_save=p.plot.is_after_sim_save,
            is_show=p.plot.is_after_sim_show,

            # Save all post-simulation plots to the same parent directory.
            save_dir_parent_basename='plot',

            # Pass all remaining arguments as is to our superclass.
            *args, **kwargs
        )
