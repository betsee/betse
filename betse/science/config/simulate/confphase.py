#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
YAML-backed simulation phase subconfigurations.
'''

# ....................{ IMPORTS                            }....................
from betse.science.config.confabc import SimConfABC, conf_alias
# from betse.util.type.types import type_check, NumericTypes

# ....................{ SUBCLASSES                         }....................
#FIXME: Implement us up.
#FIXME: Refactor the Parameters.set_time_profile() method to leverage this.
class SimConfPhase(SimConfABC):
    '''
    YAML-backed simulation phase subconfiguration, encapsulating the
    configuration of a single phase (e.g., seed, initialization, simulation)
    parsed from the current YAML-formatted simulation configuration file.

    Attributes (Path)
    ----------
    . : .
        .

    Attributes (Time)
    ----------
    . : .
        .
    '''

    # ..................{ ALIASES ~ wut                      }..................
    wut = conf_alias("['wut']['wut']", int)
