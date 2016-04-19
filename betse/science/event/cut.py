#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
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
    is_hurt_cells_leaky : bool
        `True` if the membranes of unremoved cells adjacent to removed cells are
        damaged by such removal and hence open (i.e., leak) to the environment.
    hurt_cell_leakage : float
        Extent to which damaged membranes leak to the environment if
        `is_hurt_cells_leaky` is `True` or ignored otherwise. In the former case, all
        free diffusion constants will be multiplied by this value and hence:
        * Increased if this value is greater than 1.0.
        * Decreased if this value is less than 1.0.
        * Unmodified if this value is 1.0.
        For example, the resulting Na+ membrane diffusion constant for all
        damaged membranes will be `hurt_cell_leakage * tissue_profile.Dm_Na`.
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
                    is_hurt_cells_leaky=bool(ce['membranes damaged']),
                    hurt_cell_leakage=float(ce['membrane leakage']),
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
        time: float,
        is_hurt_cells_leaky: bool,
        hurt_cell_leakage: float,
    ) -> None:
        assert types.is_sequence_nonstr(profile_names), (
            types.assert_not_sequence_nonstr(profile_names))
        assert types.is_bool(is_hurt_cells_leaky), (
            types.assert_not_bool(is_hurt_cells_leaky))
        assert types.is_numeric(hurt_cell_leakage), (
            types.assert_not_numeric(hurt_cell_leakage))

        super().__init__(time)

        self.profile_names = profile_names
        self.is_hurt_cells_leaky = is_hurt_cells_leaky
        self.hurt_cell_leakage = hurt_cell_leakage


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
