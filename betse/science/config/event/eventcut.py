#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level classes aggregating all parameters pertaining to simulation events.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseMethodUnimplementedException
from betse.science.config.event.eventabc import SimEventSpikeABC
from betse.util.io.log import logs
from betse.util.type.types import (
    type_check, NoneType, NumericTypes, SequenceTypes)

# ....................{ SUBCLASSES                         }....................
class SimEventCut(SimEventSpikeABC):
    '''
    Event removing a subset of the cell population at some simulation time step.

    Attributes
    ----------------------------
    profile_names : list
        List of the names of all applicable cut profiles, each describing a
        subset of the cell population to be removed.
    '''

    # ..................{ PUBLIC                             }..................
    @type_check
    def __init__(
        self,
        profile_names: SequenceTypes,
        time: NumericTypes,
    ) -> None:

        super().__init__(time)

        self.profile_names = profile_names


    #FIXME: Refactor the handler.removeCells() function into this method. Before
    #we do so, note that this will require refactoring this method's signature
    #everywhere to resemble:
    #    def fire(
    #        self,
    #        sim: 'Simulation',
    #        cells: 'Cells',
    #        p: 'Parameters',
    #        t: float) -> None:
    @type_check
    def fire(self, sim: 'betse.science.sim.Simulator', t: NumericTypes) -> None:
        raise BetseMethodUnimplementedException()

# ....................{ MAKERS                             }....................
@type_check
def make(p: 'betse.science.parameters.Parameters') -> (SimEventCut, NoneType):
    '''
    Create and return a new :class:`SimEventCut` instance if enabled by the
    passed simulation configuration *or* ``None`` otherwise.

    Parameters
    ----------------------------
    p : Parameters
        Current simulation configuration.

    Returns
    ----------------------------
    SimEventCut, NoneType
        Either:
        * If enabled by the passed simulation configuration, a new
          :class:`SimEventCut` instance.
        * Else, ``None``.
    '''

    # Object to be returned, defaulting to nothing.
    action = None

    ce = p._conf['cutting event']

    # If this event is enabled, create an instance of this class.
    if bool(ce['event happens']):
        # If profiles are enabled, parse this event.
        if len(p.profiles):
            action = SimEventCut(
                # Time step at which to cut. For simplicity, this is coerced
                # to be the start of the simulation.
                time=0.0,
                profile_names=ce['apply to'],
            )
        # Else, log a non-fatal warning.
        else:
            logs.log_warning(
                'Ignoring cutting event, as cut profiles are disabled.')

    return action
