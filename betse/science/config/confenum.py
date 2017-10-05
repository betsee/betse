#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Enumeration types** (i.e., :class:`enum.Enum` subclasses) required to express
typed data in simulation configuration files.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.enums import make_enum

# ....................{ ENUMS                              }....................
CellLatticeType = make_enum(
    class_name='CellLatticeType', member_names=('HEXAGONAL', 'SQUARE',))
'''
Enumeration of all supported types of **base cell lattices** (i.e., uniform grid
to which cells are situated *before* random lattice disorder is applied).

Attributes
----------
HEXAGONAL : enum
    Hexagonal base cell lattice, situating cells along a hexagonal grid.
SQUARE : enum
    Rectilinear base cell lattice, situating cells along a square grid.
'''


IonProfileType = make_enum(
    class_name='IonProfileType', member_names=(
        'BASIC', 'BASIC_CA', 'MAMMAL', 'AMPHIBIAN', 'CUSTOM',))
'''
Enumeration of all supported types of **ion profiles** (i.e., predefined sets of
all extracellular and cytosolic ions enabled by this simulation).

Note that *all* ion profiles unconditionally enable an unidentified
charge-balance anion denoted M-, as required for both simulation stability and
correctness.

Attributes
----------
BASIC : enum
    Ion profile enabling M-, Na+, K+, and proteins- ions. This profile is the
    proper subset of all other predefined ion profiles.
BASIC_CA : enum
    Ion profile enabling M-, Na+, K+, proteins-, and Ca2+ ions. This profile is
    the superset of the :attr:`BASIC` profile enabling Ca2+ ions.
MAMMAL : enum
    Ion profile enabling M-, Na+, K+, proteins-, Ca2+, Cl-, and H+ ions
    expressed in amniotic environmental concentrations, principally intended for
    mammalian cell clusters. This profile is the superset of the
    :attr:`BASIC_CA` profile enabling Cl- and H+ ions.
AMPHIBIAN : enum
    Ion profile enabling M-, Na+, K+, proteins-, Ca2+, Cl-, and H+ ions
    expressed in aquatic environmental concentrations, principally intended for
    amphibian cell clusters. This profile enables the same ions as the
    :attr:`MAMMAL` profile -- albeit in differing concentrations.
CUSTOM : enum
    User-defined ion profile. See the
    :attr:`betse.science.parameters.Parameters.ions_custom` variable.
'''
