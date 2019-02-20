#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Simulation modelling enumerations** (i.e., :class:`enum.Enum` subclasses
internally required by simulation modelling phases).
'''

# ....................{ IMPORTS                           }....................
from betse.util.type import enums

# ....................{ ENUMS                             }....................
SimPhaseKind = enums.make_enum(
    class_name='SimPhaseKind',
    member_names=('SEED', 'INIT', 'SIM',),
    is_ordered=True,
    doc='''
Ordered enumeration of all supported types of **simulation phases** (i.e.,
consecutive steps for the process of simulating simulation configurations).

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
''')
