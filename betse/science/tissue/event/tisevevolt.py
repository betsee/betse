#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level classes aggregating all parameters pertaining to simulation events.
'''

# ....................{ IMPORTS                           }....................
from betse.exceptions import BetseSimEventException
from betse.science.tissue.event.tiseveabc import SimEventPulseABC
from betse.science.math import toolbox
from betse.science.phase.phasecls import SimPhase
# from betse.util.io.log import logs
from betse.util.type.types import type_check

# ....................{ SUBCLASSES                        }....................
class SimEventPulseVoltage(SimEventPulseABC):
    '''
    **Voltage event** (i.e., event applying a directed voltage to the
    environmental boundary for a range of simulation time steps).

    A positive voltage will be applied to one boundary edge during this period;
    a negative voltage will be applied to another.

    Attributes
    ----------------------------
    peak_voltage : float
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

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def __init__(self, p: 'betse.science.parameters.Parameters') -> None:
        '''
        Initialize this voltage event for the passed simulation configuration.

        Attributes
        ----------
        p : betse.science.parameters.Parameters
            Current simulation configuration.
        '''

        aev = p._conf['apply external voltage']

        # Initialize our superclass with some passed parameters.
        super().__init__(
            p=p,
            start_time_step=float(aev['change start']),
            stop_time_step=float(aev['change finish']),
            time_step_rate=float(aev['change rate']),
        )

        # Classify all remaining parameters.
        self.peak_voltage = float(aev['peak voltage'])
        self.positive_voltage_boundary = _convert_boundary_str_to_char(
            aev['positive voltage boundary'])
        self.negative_voltage_boundary = _convert_boundary_str_to_char(
            aev['negative voltage boundary'])


    #FIXME: Refactor to resemble the superclass method signature.
    @type_check
    def fire(self, phase: SimPhase, time_step: float) -> None:

        effector = toolbox.pulse(
            time_step,
            self.start_time_step,
            self.stop_time_step,
            self.time_step_rate,
        )

        phase.sim.bound_V[self.positive_voltage_boundary] = (
             self.peak_voltage * effector)
        phase.sim.bound_V[self.negative_voltage_boundary] = (
            -self.peak_voltage * effector)

# ....................{ CONVERTERS                         }....................
#FIXME: Refactor this function to leverage a private global dictionary: e.g.,
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
#FIXME: Actually, why do we even require character constants like "T" here?
#These should all be leveraging enumeration member constants, in which the
#conversion from human-readable names to enumeration member constants will
#already be implicitly defined by the "Enum" class -- which, in turn, will
#enable us to entirely remove this function.
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
        raise BetseSimEventException(
            'Boundary edge "{}" unrecognized.'.format(side))
