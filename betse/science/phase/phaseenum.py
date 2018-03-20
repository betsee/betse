#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **simulation phase enumeration** (e.g., :class:`enum.Enum` subclass
describing different types of simulation phases) functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.enums import EnumOrdered

# ....................{ ENUMS                              }....................
SimPhaseKind = EnumOrdered('SimPhaseKind', ('SEED', 'INIT', 'SIM',))
'''
Ordered enumeration of all possible simulation phases.

Each member of this enumeration is arbitrarily comparable to each other member.
Each member's value is less than that of another member's value if and only if
the former simulation phase is performed _before_ the latter. Specifically,
this enumeration is a total ordering such that:

    >>> SimPhaseKind.SEED < SimPhaseKind.INIT < SimPhaseKind.SIM
    True

Attributes
----------
SEED : enum
    Seed simulation phase, as implemented by the
    :meth:`betse.science.simrunner.SimRunner.seed` method. This phase creates
    the cell cluster and caches this cluster to an output file.
INIT : enum
    Initialization simulation phase, as implemented by the
    :meth:`betse.science.simrunner.SimRunner.init` method. This phase
    initializes the previously created cell cluster from a cached input file
    and caches this initialization to an output file.
SIM : enum
    Proper simulation phase, as implemented by the
    :meth:`betse.science.simrunner.SimRunner.sim` method. This phase simulates
    the previously initialized cell cluster from a cached input file and caches
    this simulation to an output file.
'''
