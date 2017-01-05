#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level classes aggregating all parameters pertaining to simulation events.
'''

# ....................{ IMPORTS                            }....................
from betse.science.event.abc import Action
from betse.util.io.log import logs
from betse.util.type import types

# ....................{ EVENT                              }....................
class ActionCut(Action):
    '''
    Event permanently removing a subset of the cell population at some time step
    of the simulation.

    Attributes
    ----------------------------
    profile_names : list
        List of the names of all applicable cut profiles, each describing a
        subset of the cell population to be removed.

    '''

    # ..................{ PUBLIC ~ static                    }..................
    @staticmethod
    def make(params: 'Parameters') -> 'ActionCut':
        assert types.is_parameters(params), types.assert_not_parameters(params)

        # Object to be returned, defaulting to nothing.
        action = None

        ce = params.config['cutting event']

        # If this event is enabled, create an instance of this class.
        if bool(ce['event happens']):
            # If profiles are enabled, parse this event.
            if len(params.profiles):
                action = ActionCut(
                    # Time step at which to cut. For simplicity, this is coerced
                    # to be the start of the simulation.
                    time=0.0,
                    profile_names=ce['apply to'],
                )

            # Else, print a non-fatal warning.
            else:
                logs.log_warning(
                    'Ignoring cutting event, as cut profiles are disabled.')

        return action

    # ..................{ PUBLIC                             }..................
    def __init__(
        self,
        profile_names: list,
        time: float
    ) -> None:
        assert types.is_sequence_nonstr(profile_names), (
            types.assert_not_sequence_nonstr(profile_names))


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
    def fire(self, sim: 'Simulation', t: float) -> None:
        raise TypeError('Not implemented yet!')
