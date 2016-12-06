#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level functionality pertaining to **simulation phases** (i.e., requisite
simulation steps ordered in a typical computational pipeline).
'''

# ....................{ IMPORTS                            }....................
from betse.science.cells import Cells
from betse.science.parameters import Parameters
from betse.science.sim import Simulator
from betse.util.type.obj import objsize
from betse.util.type.obj.objsize import SizeProfilableABC

# ....................{ CLASSES ~ phaser                   }....................
class SimPhaser(SizeProfilableABC):
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

    # ..................{ SUPERCLASS                         }..................
    #FIXME: Consider generalizing the current approach to support arbitrary
    #objects to which an arbitrary number of non-builtin variables are bound.
    #After doing so:
    #
    #* The "SizeProfilableABC" class could be removed entirely.
    #* The objsize.get_size_profile() function could conceivably be improved to
    #  support a new optional "depth" keyword argument defaulting to 0. If 1,
    #  the logic below would be performed instead (in generalized form). If any
    #  other number, an exception would be raised.
    #* The profilers._profile_callable_size() function could conceivably be
    #  improved to unconditionally pass "depth=1" to its unconditional call to
    #  objsize.get_size_profile().
    #
    #Doing so is hardly essential, but would simplify and streamline the
    #implementation burden of memory profiling throughout the codebase.

    def get_size_profile(self, *args, **kwargs) -> str:

        # Number of non-builtin variables bound to this object.
        VAR_COUNT = 3

        # Maximum number of the largest non-builtin variables bound to this
        # object to synopsize the sizes of if any or None otherwise.
        vars_max = kwargs.get('vars_max', None)

        # If such a maximum is requested, average this maximum over this number
        # of non-builtin bound variables (rounding down to the nearest integer).
        if vars_max is not None:
            vars_max = round(vars_max / VAR_COUNT)

        # Human-readable strings synopsizing these objects' memory consumption.
        size_profile_cells = objsize.get_size_profile(
            obj=self.cells, vars_max=vars_max)
        size_profile_p = objsize.get_size_profile(
            obj=self.p, vars_max=vars_max)
        size_profile_sim = objsize.get_size_profile(
            obj=self.sim, vars_max=vars_max)

        # Return the concatenation of these synopses.
        return 'Cells {}\n\nSimulator {}\n\nParameters {}'.format(
            size_profile_cells, size_profile_p, size_profile_sim)
