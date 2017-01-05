#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level classes aggregating all parameters pertaining to simulation events.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseSimConfigException
from betse.science import toolbox
from betse.science.event.abc import Pulse
from betse.util.io.log import logs
from betse.util.type import types

# ....................{ EVENT                              }....................
class PulseVoltage(Pulse):
    '''
    Event applying a directed voltage to the environmental boundary during some
    time period of the simulation.

    A positive voltage will be applied to one boundary edge during this period;
    a negative voltage will be applied to another.

    Attributes
    ----------------------------
    peak_voltage : float
        Maximum voltage (V) to be applied.
    positive_voltage_boundary : str
        String identifying the boundary edge to apply this positive voltage to.
        Valid values include:
        * `T`, the top boundary.
        * `B`, the bottom boundary.
        * `L`, the left boundary.
        * `R`, the right boundary.
    negative_voltage_boundary : str
        String identifying the boundary edge to apply this negative voltage to.
        Valid values are as for `side_positive_voltage` above.
    '''

    # ..................{ PUBLIC ~ static                    }..................
    @staticmethod
    def make(params: 'Parameters') -> 'PulseVoltage':
        assert types.is_parameters(params), types.assert_not_parameters(params)

        # Object to be returned, defaulting to nothing.
        event = None

        # If this event is enabled, create an instance of this class.
        aev = params.config['apply external voltage']
        if bool(aev['event happens']):
            # If extracellular spaces are enabled, parse this event.
            if params.sim_ECM:
                event = PulseVoltage(
                    start_time=float(aev['change start']),
                    stop_time=float(aev['change finish']),
                    step_rate=float(aev['change rate']),
                    peak_voltage=float(aev['peak voltage']),
                    positive_voltage_boundary=\
                        _convert_boundary_str_to_char(
                            aev['positive voltage boundary']),
                    negative_voltage_boundary=\
                        _convert_boundary_str_to_char(
                            aev['negative voltage boundary']),
                )
            # Else, print a non-fatal warning.
            else:
                logs.log_warning(
                    'Ignoring voltage event, '
                    'as extracellular spaces are disabled.')

        return event

    # ..................{ PUBLIC                             }..................
    def __init__(
        self,
        start_time: float,
        stop_time: float,
        step_rate: float,
        peak_voltage: float,
        positive_voltage_boundary: float,
        negative_voltage_boundary: float,
    ) -> None:
        assert types.is_numeric(peak_voltage), (
            types.assert_not_numeric(peak_voltage))
        assert types.is_char(positive_voltage_boundary), (
            types.assert_not_char(positive_voltage_boundary))
        assert types.is_char(negative_voltage_boundary), (
            types.assert_not_char(negative_voltage_boundary))

        super().__init__(
            start_time=start_time,
            stop_time=stop_time,
            step_rate=step_rate,)

        self.peak_voltage = peak_voltage
        self.positive_voltage_boundary = positive_voltage_boundary
        self.negative_voltage_boundary = negative_voltage_boundary


    def fire(self, sim: 'Simulation', t: float) -> None:
        effector = toolbox.pulse(
            t, self.start_time, self.stop_time, self.step_rate)

        sim.bound_V[self.positive_voltage_boundary] = \
             self.peak_voltage * effector
        sim.bound_V[self.negative_voltage_boundary] = \
            -self.peak_voltage * effector

# ....................{ CONVERTERS                         }....................
def _convert_boundary_str_to_char(side: str) -> str:
        '''
        Convert the passed human-readable string constant identifying an
        environmental boundary edge (e.g., `top`) into the corresponding
        machine-readable character constant (e.g., `T`).
        '''
        assert types.is_str(side), types.assert_not_str(side)

        if   side == 'top':    return 'T'
        elif side == 'bottom': return 'B'
        elif side == 'left':   return 'L'
        elif side == 'right':  return 'R'
        else:
            raise BetseSimConfigException(
                'Boundary edge "{}" unrecognized.'.format(side))
