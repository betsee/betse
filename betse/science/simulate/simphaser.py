#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level functionality pertaining to **simulation phases** (i.e., requisite
simulation steps ordered in a typical computational pipeline).
'''

# ....................{ IMPORTS                            }....................
from betse.science.cells import Cells
from betse.science.parameters import Parameters
from betse.science.sim import Simulator

# ....................{ CLASSES ~ phaser                   }....................
class SimPhaser(object):
    '''
    High-level simulation phase class encapsulating the three principal objects
    needed for running a single simulation phase.

    Attributes
    ----------
    _cells : Cells
        Current cell cluster.
    _p : Parameters
        Current simulation configuration.
    _sim : Simulator
        Current simulation.
    '''

    # ..................{ INITIALIZORS                       }..................
    def __init__(
        self, cells: Cells, p: Parameters, sim: Simulator) -> None:
        '''
        Initialize this simulation phase instance.

        Parameters
        ----------
        cells : Cells
            Current cell cluster.
        p : Parameters
            Current simulation configuration.
        sim : Simulator
            Current simulation.
        '''

        # Classify the passed parameters.
        self.cells = cells
        self.p = p
        self.sim = sim
