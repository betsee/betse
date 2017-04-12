#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level **simulation phase cache** (e.g., container persisting previously
constructed large-scale objects for a simulation phase) functionality.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta  #, abstractmethod
from betse.science.simulate.simphase import SimPhase
from betse.util.py import references
from betse.util.type.types import type_check

# ....................{ CLASSES                            }....................
class SimPhaseCache(object):
    '''
    High-level simulation phase cache, persisting previously constructed
    large-scale objects for subsequent reuse by a single simulation phase.

    Attributes
    ----------
    phase : SimPhase
        Parent simulation phase.

    Attributes (Subcache)
    ----------
    cells : SimPhaseSubcacheCells
        Subcache of all cell cluster objects constructed for this phase.
    vector : SimPhaseSubcacheVector
        Subcache of all vectors and vector fields constructed for this phase.
    '''

    # ..................{ INITIALIZORS                       }..................
    @type_check
    def __init__(self, phase: SimPhase) -> None:
        '''
        Initialize this simulation phase cache.

        Parameters
        ----------
        phase : SimPhase
            Parent simulation phase.
        '''

        # Avoid circular import dependencies.
        from betse.science.simulate.cache.cachevec import SimPhaseSubcacheVector

        # Classify all passed parameters as weak rather than strong reference,
        # circumventing circular references and complications thereof.
        self.phase = references.proxy_weak(phase)

        # Classify all subcaches.
        self.cells = SimPhaseSubcacheCells(self)
        self.vector = SimPhaseSubcacheVector(self)

# ....................{ SUPERCLASSES                       }....................
class SimPhaseSubcacheABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all **simulation phase subcache** (i.e., container
    persisting previously constructed large-scale objects for some facet of a
    simulation phase) subclasses.

    Design
    ----------
    Subcaches provide namespace isolation but are otherwise purely superficial.
    Technically, all objects cached by a subcache could simply be cached by the
    parent :class:`SimPhaseCache` object instead. Doing so would quickly become
    cumbersome, however, both for maintenance and reuse.

    Attributes
    ----------
    _cache : SimPhaseCache
        Parent simulation phase cache, permitting this subcache to access
        objects cached by other subcaches.
    _phase : SimPhase
        Parent simulation phase. This attribute is technically accessible as
        :attr:`_cache.phase` but is provided by this class for convenience.
    '''

    # ..................{ INITIALIZORS                       }..................
    @type_check
    def __init__(self, cache: SimPhaseCache) -> None:
        '''
        Initialize this simulation phase subcache.

        Parameters
        ----------
        cache : SimPhaseCache
            Parent simulation phase cache.
        '''

        # Classify all passed parameters as weak rather than strong reference,
        # circumventing circular references and complications thereof.
        self._cache = references.proxy_weak(cache)

        # For convenience, reclassify the current phase as well.
        self._phase = self._cache.phase

# ....................{ SUBCLASSES                         }....................
#FIXME: Shift into a new "cachecells" submodule of this subpackage.
class SimPhaseSubcacheCells(SimPhaseSubcacheABC):
    '''
    Lower-level simulation phase cell cluster subcache, persisting all
    previously constructed cell cluster objects for a single simulation phase.

    Design
    ----------
    All objects persisted by this subcache are required only sporadically by
    one or more isolated features (e.g., plots, animations). Hence, these
    objects are *not* suitable for ownership by the general-purpose
    :class:`Cells` class.
    '''

    pass
