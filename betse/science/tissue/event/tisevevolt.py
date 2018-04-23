#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level classes aggregating all parameters pertaining to simulation events.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseSimConfException
from betse.science.tissue.event.tiseveabc import SimEventPulseABC
from betse.science.math import toolbox
from betse.util.io.log import logs
from betse.util.type.types import type_check, NoneType, NumericSimpleTypes

# ....................{ SUBCLASSES                         }....................
class SimEventPulseVoltage(SimEventPulseABC):
    '''
    Event applying a directed voltage to the environmental boundary for some
    range of simulation time steps.

    A positive voltage will be applied to one boundary edge during this period;
    a negative voltage will be applied to another.

    Attributes
    ----------------------------
    peak_voltage : NumericSimpleTypes
        Maximum voltage (V) to be applied.
    positive_voltage_boundary : str
        Character identifying the boundary edge to apply this positive voltage
        to. Valid values include:
        * ``T``, the top boundary.
        * ``B``, the bottom boundary.
        * ``L``, the left boundary.
        * ``R``, the right boundary.
    negative_voltage_boundary : str
        Character identifying the boundary edge to apply this negative voltage
        to. Valid values are as for ``side_positive_voltage`` above.
    '''

    # ..................{ PUBLIC                             }..................
    @type_check
    def __init__(
        self,
        start_time: NumericSimpleTypes,
        stop_time: NumericSimpleTypes,
        step_rate: NumericSimpleTypes,
        peak_voltage: NumericSimpleTypes,
        positive_voltage_boundary: str,
        negative_voltage_boundary: str,
    ) -> None:

        # Initialize our superclass with some passed parameters.
        super().__init__(
            start_time=start_time,
            stop_time=stop_time,
            step_rate=step_rate,
        )

        # Classify all remaining parameters.
        self.peak_voltage = peak_voltage
        self.positive_voltage_boundary = positive_voltage_boundary
        self.negative_voltage_boundary = negative_voltage_boundary


    @type_check
    def fire(self, sim: 'betse.science.sim.Simulator', t: NumericSimpleTypes) -> None:

        effector = toolbox.pulse(
            t, self.start_time, self.stop_time, self.step_rate)

        sim.bound_V[self.positive_voltage_boundary] = (
             self.peak_voltage * effector)
        sim.bound_V[self.negative_voltage_boundary] = (
            -self.peak_voltage * effector)

# ....................{ MAKERS                             }....................
@type_check
def make(p: 'betse.science.parameters.Parameters') -> (SimEventPulseVoltage, NoneType):
    '''
    Create and return a new :class:`SimEventPulseVoltage` instance if enabled by the
    passed simulation configuration *or* ``None`` otherwise.

    Parameters
    ----------------------------
    p : Parameters
        Current simulation configuration.

    Returns
    ----------------------------
    SimEventPulseVoltage, NoneType
        Either:
        * If enabled by the passed simulation configuration, a new
          :class:`SimEventPulseVoltage` instance.
        * Else, ``None``.
    '''

    # Object to be returned, defaulting to nothing.
    event = None

    # If this event is enabled, create an instance of this class.
    aev = p._conf['apply external voltage']

    if bool(aev['event happens']):
        # If extracellular spaces are enabled, parse this event.
        if p.is_ecm:
            event = SimEventPulseVoltage(
                start_time=float(aev['change start']),
                stop_time=float(aev['change finish']),
                step_rate=float(aev['change rate']),
                peak_voltage=float(aev['peak voltage']),
                positive_voltage_boundary=(
                    _convert_boundary_str_to_char(
                        aev['positive voltage boundary'])),
                negative_voltage_boundary=(
                    _convert_boundary_str_to_char(
                        aev['negative voltage boundary'])),
            )
        # Else, log a non-fatal warning.
        else:
            logs.log_warning(
                'Ignoring voltage event, '
                'as extracellular spaces are disabled.')

    return event

# ....................{ CONVERTERS                         }....................
#FIXME: Utter rubbish. Refactor this function to leverage a private global
#dictionary constant instead: e.g.,
#
#    BOUNDARY_STR_TO_CHAR = {
#        'top':    'T',
#        'bottom': 'B',
#        'left':   'L',
#        'right':  'R',
#    }
#    boundary_char = BOUNDARY_STR_TO_CHAR.get(side, None)
#    if boundary_char is None:
#        raise BetseSimConfException(
#            'Boundary edge "{}" unrecognized.'.format(side))
#    return boundary_char
#
#The sane means of doing so would be to refactor the implementations of the
#"positive_voltage_boundary" and "positive_voltage_boundary" attributes above
#into @property_cached-style properties. Possibly? Or not. That does rather
#seem like overkill. The current function-based approach is probably superior.
@type_check
def _convert_boundary_str_to_char(side: str) -> str:
    '''
    Convert the passed human-readable string constant identifying an
    environmental boundary edge (e.g., `top`) into the corresponding
    machine-readable character constant (e.g., `T`).
    '''

    if   side == 'top':    return 'T'
    elif side == 'bottom': return 'B'
    elif side == 'left':   return 'L'
    elif side == 'right':  return 'R'
    else:
        raise BetseSimConfException(
            'Boundary edge "{}" unrecognized.'.format(side))
