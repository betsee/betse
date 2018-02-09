#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Class hierarchy designing **simulation phase requirements** (i.e., objects
encapsulating the current state of a single simulation feature for a given
simulation phase).
'''

# ....................{ IMPORTS                            }....................
from betse.science.phase.phasecls import SimPhase
from betse.science.phase.require.abc.phasereqabc import SimPhaseRequirementABC
from betse.util.io.log import logs
from betse.util.type import iterables
from betse.util.type.call.memoizers import property_cached
from betse.util.type.text import strs
from betse.util.type.types import type_check, IterableOrNoneTypes, NoneType

# ....................{ SUBCLASSES ~ requirements          }....................
class SimPhaseRequirements(SimPhaseRequirementABC, frozenset):
    '''
    Immutable set of simulation phase requirements, requiring zero or more
    arbitrary requirements to be satisfied.

    Operators
    ----------
    This set strictly conforms to the :class:`frozenset` API and hence supports
    the following set operators for combining multiple instances of this class:

    * ``|``, performing set union.
    * ``&``, performing set intersection.
    * ``-``, performing asymmetric set difference.
    * ``^``, performing symmetric set difference.

    Likewise, this set supports the following set operators for comparing any
    two instances of this class in terms of subset and superset relations:

    * ``a < b``, ``True`` only if set ``a`` is a proper subset of set ``b``.
    * ``a <= b``, ``True`` only if set ``a`` is a subset of set ``b``.
    * ``a > b``, ``True`` only if set ``a`` is a proper superset of set ``b``.
    * ``a >= b``, ``True`` only if set ``a`` is a superset of set ``b``.

    Design
    ----------
    This set satisfies the :class:`SimPhaseRequirementABC` API and is thus
    itself a valid requirement. Indeed, this set permits instances of that
    lower-level API to be efficiently composed together via an immutable
    set-like API and is thus recommended over that API for external usage.

    This set inherits the abstract :class:`collections.abc.Set` mixin rather
    than the concrete :class:`frozenset` type, as the latter is *not* a mixin.
    In all respects, however, this set is usable as a :class:`frozenset`.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __new__(cls, iterable: IterableOrNoneTypes = None) -> (
        'betse.science.phase.require.abc.phasereqset.SimPhaseRequirements'):
        '''
        Create, initialize, and return this requirement.

        Design
        ----------
        To satisfy the :class:`frozenset` API, this method *must* accept only
        a single iterable. To quote `official documentation`_:

            Since some set operations create new sets, the default mixin methods
            need a way to create new instances from an iterable. The class
            constructor is assumed to have a signature in the form
            ``ClassName(iterable)``. That assumption is factored-out to an
            internal classmethod called :meth:`_from_iterable` which calls
            ``cls(iterable)`` to produce a new set. If the :class:`Set` mixin is
            being used in a class with a different constructor signature, you
            will need to override :meth:`_from_iterable` with a classmethod that
            can construct new instances from an iterable argument.

        .. _official documentation:
            https://docs.python.org/3/library/collections.abc.html#collections.abc.AsyncGenerator

        Parameters
        ----------
        iterable : IterableOrNoneTypes
            Iterable of zero or more arbitrary requirements, collectively
            defining this requirement. Since this requirement is an immutable
            set, this iterable is the *only* means of defining this requirement.
            Defaults to ``None``, in which case this set is permanently empty.
        '''

        # If no iterable was passed, default this iterable to the empty tuple.
        if iterable is None:
            iterable = ()

        # If any item of the passed iterable is *NOT* a requirement, raise an
        # exception.
        iterables.die_unless_items_instance_of(
            iterable=iterable, cls=SimPhaseRequirementABC)

        # Create, initialize, and return a new instance of this class containing
        # all unique items of this iterable.
        return super().__new__(cls, iterable)

        # self = super().__new__(cls, iterable)
        # logs.log_info('requirements for "%s": %r', self.name, self)
        # return self

    # ..................{ SUPERCLASS ~ requirement           }..................
    # Abstract properties required to be implemented by the
    # "SimPhaseRequirementABC" superclass.

    @property_cached
    def name(self) -> str:
        '''
        Human-readable name of this requirement, defined as the conjunction of
        the names of all requirements comprising this set (e.g., "full solver,
        fluid flow, and calcium ions (Ca2+)").

        For usability, the first character of this name should typically be
        lower- rather than uppercase.
        '''

        # If this set is empty, return a sane empty name.
        if not self:
            return 'no simulation features'
        # Else, this set contains at least one requirement.

        # Generator yielding the name of each requirement in this set.
        requirement_names = (requirement.name for requirement in self)

        # Create, return, and cache this name to be a human-readable conjunction
        # of these names.
        return strs.join_as_conjunction(*requirement_names)


    @type_check
    def is_satisfied(self, phase: SimPhase) -> bool:
        '''
        ``True`` only if all child requirements composed by this parent
        requirement are satisfied by the passed simulation phase.
        '''

        # This requirement is satisfied if and only if...
        return all(
            # This child requirement is satisfied by this phase.
            requirement.is_satisfied(phase)
            # For each child requirement of this parent requirement...
            for requirement in self)


    @type_check
    def set_satisfied(self, phase: SimPhase) -> None:
        '''
        Modify the passed simulation phase as needed to ensure that all
        child requirements composed by this parent requirement are satisfied
        by this phase.
        '''

        # For each child requirement of this parent requirement...
        for requirement in self:
            # Ensure this child requirement is satisfied by this phase.
            requirement.set_satisfied(phase)

    # ..................{ TESTERS                            }..................
    # For orthogonality with the existing frozenset.isdisjoint() method, this
    # method is intentionally *NOT* named is_intersection().
    def isintersection(self, other: object) -> bool:
        '''
        ``True`` only if the passed object is also an immutable set of
        simulation phase requirements *and* these two sets are **intersecting**
        (i.e., these sets have a non-null intersection), in which case these
        sets share one or more requirements in common.
        '''

        return not self.isdisjoint(other)

    # ..................{ OPERATORS                          }..................
    #FIXME: Something appears to be critically broken here. The
    #frozenset.__or__() method refuses to invoke our
    #SimPhaseRequirements.__new__() method returning an instance of
    #"SimPhaseRequirements" and instead invokes its own frozenset.__new__()
    #method returning an instance of "frozenset". This is a well-known issue
    #with the Python 2.x line purportedly resolved in Python 3.0 over a decade
    #ago: issue #1721812.
    #
    #Unfortunately, this issue still seems to be unresolved -- at least, for us.
    #The only alternative would be to completely refactor this subclass to
    #inherit from the "Hashable" and "Set" APIs of the "collections.abc" module
    #instead. We actually attempted that initially. Unfortunately, that approach
    #suffered several gotchas:
    #
    #* It was considerably more heavyweight than this approach, requiring in
    #  upwards of 400 more lines of code.
    #* It didn't actually work, either. The issue was disturbingly similar to
    #  this issue, albeit different. Under that approach, attempting to take the
    #  union of two instances of this class correctly returned a new instance of
    #  the expected type. For unknown reasons, that union was *ALWAYS* empty.
    #
    #Ultimately, the current approach is the lesser of two evils. Since this is
    #demonstrably awful, however, we should consider submitting a minimal-length
    #example (MLE) as a StackOverflow question. In theory, that should resolve
    #exactly what is going on here.
    #FIXME: *ALL OTHER SET OPERATORS RETURNING NEW SETS MUST ALSO BE OVERLOADED
    #IN THE SAME EXACT MANNER.* Until then, they remain similarly broken. I'm
    #sufficiently discouraged by this, at the moment, that I can't be bothered.
    #Nonetheless, for the record, these methods are:
    #
    #    __and__, __rand__, __ror__, __rsub__, __rxor__, __sub__, __xor__,
    #    add, copy, difference, intersection, symmetric_difference, union

    def __or__(self, other):

        # If "other" is an instance of this class, defer to the corresponding
        # comparison operator of the "frozenset" class; else, fail gracefully.
        return (
            type(self)(frozenset(self) | frozenset(other))
            if isinstance(other, SimPhaseRequirements) else
            NotImplemented)

# ....................{ TYPES                              }....................
SimPhaseRequirementsOrNoneTypes = (SimPhaseRequirements, NoneType)
'''
Tuple of both the immutable requirements set tupe *and* that of the singleton
``None`` object.
'''
