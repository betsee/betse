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
