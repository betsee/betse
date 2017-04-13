#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level vector subcache functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.science.math.vector.veccls import VectorCellsCache
from betse.science.simulate.cache.cacheabc import SimPhaseCacheABC
from betse.util.type.call.memoizers import property_cached
# from betse.util.type.types import type_check

# ....................{ SUBCLASSES                         }....................
#FIXME: Refactor all "betse.science.math.vector.vectormake" functions into
#cached properties of this class.
class SimPhaseCacheVector(SimPhaseCacheABC):
    '''
    Simulation phase-specific vector subcache, persisting all previously
    constructed vectors for a single simulation phase.
    '''

    pass
