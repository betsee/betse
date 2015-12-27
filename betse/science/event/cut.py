#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level classes aggregating all parameters pertaining to simulation events.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseExceptionParameters
from betse.science.event.abc import Event
from betse.science.tissue.picker import TissuePicker
from betse.util.io import loggers
from betse.util.type import types

# ....................{ EVENT                              }....................
#FIXME: Actually use us, please!
class EventCut(Event):
    '''
    Event permanently removing a subset of the cell population at the start
    (i.e., first time step) of the current tissue simulation.

    Attributes
    ----------------------------
    time : float
        Time step at which to cut.
    picker : TissuePicker
        Object selecting the subset of the cell population to be removed.
    '''

    # ..................{ PUBLIC ~ static                    }..................
    @staticmethod
    def make(params: 'Parameters') -> 'EventCut':
        assert types.is_parameters(params), types.assert_not_parameters(params)

        # Object to be returned, defaulting to nothing.
        event = None

        # To permit a subset of event parameters to be modified without
        # requiring simulation reinitialization, these parameters are subdivided
        # into multiple top-level configuration keys. Such is the code ghetto.
        ce = params.config['cutting event']
        tpd = params.config['tissue profile definition']

        # If this event is enabled, create an instance of this class.
        if bool(ce['event happens']):
            # If tissue profiles are enabled, parse this event.
            if tpd['profiles enabled']:
                event = EventCut(
                    picker=TissuePicker.make(
                        tpd['cut profile']['cell targets'], params),
                )
            # Else, print a non-fatal warning.
            else:
                loggers.log_warning(
                    'Ignoring cutting event, as tissue profiles are disabled.')

            # # If extracellular spaces are enabled, parse this event.
            # if params.sim_ECM:

        return event

    # ..................{ PUBLIC                             }..................
    def __init__(
        self,
        time: float,
        picker: 'TissuePicker',
    ) -> None:
        assert types.is_numeric(time), types.assert_not_numeric(time)
        assert types.is_tissue_picker(picker), \
            types.assert_not_tissue_picker(picker)

        self.time = time
        self.picker = picker

    # def fire(self, sim: 'Simulation', t: float) -> None:
    #     pass
