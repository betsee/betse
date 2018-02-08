#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Class hierarchy designing **simulation phase requirements** (i.e., objects
encapsulating the current state of a single simulation feature for a given
simulation phase).
'''

#FIXME: Generalize the requirement globals defined below to be more globally
#usable. Currently, they require the current "phase" object be passed to their
#is_satisfied() methods, which is cumbersome at best. To rectify this, first
#note that whether or not a given phase satisfies a given requirement is a
#constant that should *NEVER* change during that phase. Ergo, we can efficiently
#associate a given phase with all of the following requirements as follows:
#
#* Define a new "SimPhaseRequirements" class in this submodule.
#* Define a SimPhaseRequirements.__init__() method:
#  * Accepting only the current "SimPhase" object.
#  * For each requirement global defined below (e.g., "DEFORM"):
#    * Defining a boolean instance variable of the same name, lowercased and
#      prefixed by "is_", whose value is the value of the corresponding
#      requirement global passed this phase.
#
#For example, as a first-draft implementation:
#
# ....................{ CLASSES                            }....................
# #FIXME: Document us up.
# class SimPhaseRequirements(object):
#
#     def __init__(self, phase: SimPhase) -> None:
#
#         self.is_deform = DEFORM(phase)
#
#Pretty sleek, eh? While these instance variables could also be defined as
#read-only properties, doing so is significantly more cumbersome and slightly
#less efficient. So, why bother? There are too many requirements to make doing
#so worthwhile or sane. (Yay!)

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractproperty
from betse.science.phase.phasecls import SimPhase
# from betse.util.io.log import logs
from betse.util.type import iterables
from betse.util.type.call.memoizers import property_cached
from betse.util.type.text import strs
from betse.util.type.types import (
    type_check,
    CallableTypes,
    EnumMemberType,
    IterableTypes,
    IterableOrNoneTypes,
    MappingOrNoneTypes,
    NoneType,
)

# Despite the ambiguous class name "Set," this abstract mixin is indeed specific
# to immutable rather than mutable sets, as evidenced by the existence of
# "collections.abc.MutableSet".
from collections.abc import Hashable
from collections.abc import Set as ImmutableSet

# ....................{ SUPERCLASSES                       }....................
class SimPhaseRequirementABC(metaclass=ABCMeta):
    '''
    Abstract base class of all **simulation phase requirement** (i.e., object
    encapsulating the current state of a single simulation feature for a given
    simulation phase) subclasses.
    '''

    # ..................{ INITIALIZERS                       }..................
    # Empty initializer, defined merely to simplify subclass initialization.
    @type_check
    def __init__(self, *args, **kwargs) -> None:
        '''
        Initialize this requirement.
        '''

        pass

    # ..................{ PROPERTIES                         }..................
    # Subclasses are required to implement the following abstract properties.

    @abstractproperty
    def name(self) -> str:
        '''
        Human-readable name of this requirement (e.g., "extracellular spaces").

        For usability, the first character of this name should typically be
        lower- rather than uppercase.
        '''

        pass


    @abstractproperty
    def is_satisfied(self) -> CallableTypes:
        '''
        Callable (e.g., function, lambda) passed only the current simulation
        phase, returning ``True`` only if this requirement is enabled in this
        phase.
        '''

        pass


    @abstractproperty
    def set_satisfied(self) -> CallableTypes:
        '''
        Callable (e.g., function, lambda) passed only the current simulation
        phase, enabling this requirement in this phase.
        '''

        pass

# ....................{ SUBCLASSES ~ requirement           }....................
class SimPhaseRequirement(SimPhaseRequirementABC):
    '''
    **Simulation phase requirement** (i.e., high-level object encapsulating the
    current state of a single simulation feature, such as extracellular spaces).

    Caveats
    ----------
    The higher-level :class:`SimPhaseRequirements` class, which permits multiple
    instances of this lower-level class to be efficiently composed together via
    an immutable set-like API, is strongly preferable for external usage.
    Instances of this lower-level class should typically *only* be instantiated
    by the companion :mod:`betse.science.phase.require.phasereqs` submodule.

    Attributes
    ----------
    is_satisfied : CallableTypes
        Callable (e.g., function, lambda) passed the current simulation phase,
        returning ``True`` only if this requirement is satisfied by this phase.
    set_satisfied : CallableTypes
        Callable (e.g., function, lambda) passed the current simulation phase,
        modifying this phase as needed to ensure this requirement is satisfied
        by this phase.
    name : str
        Human-readable name of this requirement (e.g., ``extracellular
        spaces``). For simplicity, the first character of this name should
        typically be lower- rather than uppercase.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        name: str,
        is_satisfied: CallableTypes,
        set_satisfied: CallableTypes,
    ) -> None:
        '''
        Initialize this requirement.

        Parameters
        ----------
        _name : str
            Human-readable name of this requirement (e.g., ``extracellular
            spaces``). For simplicity, the first character of this name should
            typically be lower- rather than uppercase.
        _is_satisfied : CallableTypes
            Callable (e.g., function, lambda) passed only the current
            simulation phase, returning ``True`` only if this requirement is
            enabled in this phase.
        _set_satisfied : CallableTypes
            Callable (e.g., function, lambda) passed only the current
            simulation phase, enabling this requirement in this phase.
        '''

        # Classify all passed parameters.
        self._is_satisfied = is_satisfied
        self._set_satisfied = set_satisfied
        self._name = name

    # ..................{ PROPERTIES                         }..................
    # Subclasses are required to implement the following abstract properties.

    @property
    def name(self) -> str:
        return self._name

    @property
    def is_satisfied(self) -> CallableTypes:
        return self._is_satisfied

    @property
    def set_satisfied(self) -> CallableTypes:
        return self._set_satisfied

# ....................{ SUBCLASSES ~ requirement : body    }....................
class SimPhaseRequirementEmbodied(SimPhaseRequirement):
    '''
    String-based simulation phase requirement, initialized by strings defining
    this requirement's functions rather than actual functions.

    This requirement is a caller convenience simplifying initialization in the
    common case of a requirement whose functions are definable as strings.

    Attributes
    ----------
    _func_bodies : str
        String of Python code with which the :attr:`is_satisfied` and
        :attr:`set_satisfied` functions are defined. Technically, this string
        need *not* be classified as an instance variable. Pragmatically, doing
        so assists debugging elsewhere in the codebase.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,

        # Mandatory parameters.
        is_satisfied_body: str,
        set_satisfied_body: str,

        # Optional parameters.
        body_attrs: MappingOrNoneTypes = None,

        # Parental parameters.
        *args, **kwargs
    ) -> None:
        '''
        Initialize this requirement.

        Parameters
        ----------
        is_satisfied_body: str
            String of zero or more arbitrary complex Python statements comprising
            the body of a dynamically defined function:
            * Passed only the current simulation phase.
            * Returning ``True`` only if this requirement is enabled in this
              phase.
        set_satisfied_body: str
            String of zero or more arbitrary complex Python statements comprising
            the body of a dynamically defined function:
            * Passed only the current simulation phase.
            * Enabling this requirement in this phase.
        body_attrs : MappingOrNoneTypes
            Dictionary mapping from the name to value of each attribute globally
            exposed to the bodies of these dynamically defined functions.
            Defaults to ``None``, in which case these bodies may only reference:
            * The ``phase`` attribute defining the current simulation phase.
            * The ``SimPhase`` class.
            * the ``type_check`` decorator.

        All remaining parameters are passed as is to the superclass
        :meth:`SimPhaseRequirement.__init__` method.
        '''

        # In unpassed, default this parameter to the empty dictionary.
        if body_attrs is None:
            body_attrs = {}

        # Dictionary mapping from the name to value of each attribute globally
        # exposed to the declaration of these dynamically defined functions.
        #
        # Since the exec() statement called below adds key-value pairs providing
        # these functions to this dictionary, this dictionary is copied from the
        # passed dictionary to avoid mutating caller data.
        func_globals = {
            'type_check': type_check,
            'SimPhase': SimPhase,
        }
        func_globals.update(body_attrs)

        # Dictionary mapping from the name to value of each such function.
        funcs = {}

        # Raw string defining the functions to be passed to our superclass,
        # classified to assist debugging elsewhere in the codebase.
        self._func_bodies = '''
@type_check
def is_satisfied(phase: SimPhase) -> bool:
    {is_satisfied_body}

@type_check
def set_satisfied(phase: SimPhase) -> None:
    {set_satisfied_body}
'''.format(
    is_satisfied_body=is_satisfied_body,
    set_satisfied_body=set_satisfied_body,
)

        # Dynamically define these functions.
        exec(self._func_bodies, func_globals, funcs)
        # logs.log_debug('requirement: %s; funcs: %r; functions: %s',
        #     kwargs['name'], funcs, func_bodies,)

        # Initialize our superclass with these functions.
        super().__init__(
            *args,
            is_satisfied=funcs['is_satisfied'],
            set_satisfied=funcs['set_satisfied'],
            **kwargs)


class SimPhaseRequirementBoolExpr(SimPhaseRequirementEmbodied):
    '''
    String-based boolean simulation phase requirement, initialized by a single
    **boolean expression string** (i.e., string dynamically evaluating to a
    single gettable and settable boolean value).

    This requirement is a caller convenience simplifying initialization in the
    common case of a requirement reducing to a single boolean.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, bool_expr: str, *args, **kwargs) -> None:
        '''
        Initialize this requirement.

        Parameters
        ----------
        bool_expr : str
            Arbitrarily complex Python expression relative to the current
            simulation phase evaluating to a gettable and settable attribute
            whose value is this requirement's boolean state (e.g.,
            ``phase.p.is_ecm``). This expression may (and typically should)
            refer to the ``phase`` variable, bound to the current simulation
            phase.

        All remaining parameters are passed as is to the superclass
        :meth:`SimPhaseRequirementEmbodied.__init__` method.
        '''

        # Initialize our superclass with all passed parameters.
        super().__init__(
            *args,
            is_satisfied_body='return ' + bool_expr,
            set_satisfied_body=bool_expr + ' = True',
            **kwargs)


class SimPhaseRequirementEnumExpr(SimPhaseRequirementEmbodied):
    '''
    String-based enumeration member simulation phase requirement, initialized by
    a single **enumeration member expression string** (i.e., string dynamically
    evaluating to a single gettable and settable enumeration member value).

    This requirement is a caller convenience simplifying initialization in the
    common case of a requirement reducing to a single enumeration member.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        enum_expr: str,
        enum_member: EnumMemberType,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this requirement.

        Parameters
        ----------
        enum_expr : str
            Arbitrarily complex Python expression relative to the current
            simulation phase evaluating to a gettable and settable attribute
            whose value is this requirement's enumeration state (e.g.,
            ``phase.p.solver_type``). This expression may (and typically should)
            refer to the ``phase`` variable, bound to the current simulation
            phase.
        enum_member : EnumMemberType
            Member of this enumeration satisfying this requirement (e.g.,
            ``SolverType.FULL``).

        All remaining parameters are passed as is to the superclass
        :meth:`SimPhaseRequirementEmbodied.__init__` method.
        '''

        # Bodies of the functions dynamically defined by our superclass.
        is_satisfied_body = 'return {} is {}'.format(
            enum_expr, enum_member.name)
        set_satisfied_body = '{} = {}'.format(
            enum_expr, enum_member.name)

        # Dictionary mapping from the name to value of this enumeration
        # member, locally exposing this member to the above bodies.
        body_attrs = {enum_member.name: enum_member}

        # Initialize our superclass with all passed parameters.
        super().__init__(
            *args,
            is_satisfied_body=is_satisfied_body,
            set_satisfied_body=set_satisfied_body,
            body_attrs=body_attrs,
            **kwargs)

# ....................{ SUBCLASSES ~ requirements          }....................
class SimPhaseRequirements(Hashable, ImmutableSet, SimPhaseRequirementABC):
    '''
    Immutable set of simulation phase requirements, requiring zero or more
    arbitrary requirements to be satisfied.

    This set strictly conforms to the :class:`frozenset` API and hence supports
    the following set operators for combining multiple instances of this class:

    * `|`, performing set union.
    * `&`, performing set intersection.
    * `-`, performing asymmetric set difference.
    * `^`, performing symmetric set difference.

    Design
    ----------
    This set satisfies the :class:`SimPhaseRequirementABC` API and is thus
    itself a valid requirement. Indeed, this set permits instances of that
    lower-level API to be efficiently composed together via an immutable
    set-like API and is thus recommended over that API for external usage.

    This set inherits the abstract :class:`collections.abc.Set` mixin rather
    than the concrete :class:`frozenset` type, as the latter is *not* a mixin.
    In all respects, however, this set is usable as a :class:`frozenset`.

    Attributes
    ----------
    _requirements : frozenset
        Immutable set of zero or more arbitrary requirements, collectively
        defining this requirement.
    '''

    # ..................{ INSANITY                           }..................
    # Satisfy the "Hashable" mixin via the "ImmutableSet" mixin. While
    # ostensibly crazy, this is the officially documented means of doing so.
    # Failing to do so commonly results in the following exception on attempting
    # to utilize instances of this class elsewhere:
    #     TypeError: unhashable type: 'SimPhaseRequirements'
    __hash__ = ImmutableSet._hash

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, iterable: IterableOrNoneTypes = None) -> None:
        '''
        Initialize this requirement.

        Design
        ----------
        To satisfy the :class:`ImmutableSet` API, this method *must* accept only
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

        # Initialize our superclasses.
        super().__init__()

        # If no iterable was passed, default this iterable to the empty tuple.
        if iterable is None:
            iterable = ()

        # If any item of the passed iterable is *NOT* a requirement, raise an
        # exception.
        iterables.die_unless_items_instance_of(
            iterable=iterable, cls=SimPhaseRequirementABC)

        # Classify this iterable as a set for subsequent efficiency.
        self._requirements = frozenset(iterable)

    # ..................{ SUPERCLASS ~ set                   }..................
    # Abstract methods required to be implemented by the "ImmutableSet" mixin.

    def __iter__(self):
        return iter(self._requirements)

    def __contains__(self, value) -> bool:
        return value in self._requirements

    def __len__(self) -> int:
        return len(self._requirements)

    # ..................{ SUPERCLASS ~ set : compare         }..................
    # Abstract methods recommended to be implemented by the "ImmutableSet"
    # mixin. To quote official documentation: "To override the comparisons
    # (presumably for speed, as the semantics are fixed), redefine __le__() and
    # __ge__(), then the other operations will automatically follow suit."

    def __le__(self, other) -> bool:

        # If "other" is an instance of this class, defer to the corresponding
        # comparison operator of the "frozenset" class; else, fail gracefully.
        return (
            self._requirements <= other._requirements
            if isinstance(other, SimPhaseRequirements) else
            NotImplemented)


    def __ge__(self, other) -> bool:

        # If "other" is an instance of this class, defer to the corresponding
        # comparison operator of the "frozenset" class; else, fail gracefully.
        return (
            self._requirements >= other._requirements
            if isinstance(other, SimPhaseRequirements) else
            NotImplemented)

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
        if not self._requirements:
            return 'no simulation features'
        # Else, this set contains at least one requirement.

        # Generator yielding the name of each requirement in this set.
        requirement_names = (
            requirement.name for requirement in self._requirements)

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
            for requirement in self._requirements)


    @type_check
    def set_satisfied(self, phase: SimPhase) -> None:
        '''
        Modify the passed simulation phase as needed to ensure that all
        child requirements composed by this parent requirement are satisfied
        by this phase.
        '''

        # For each child requirement of this parent requirement...
        for requirement in self._requirements:
            # Ensure this child requirement is satisfied by this phase.
            requirement.set_satisfied(phase)

# ....................{ TYPES                              }....................
SimPhaseRequirementsOrNoneTypes = (SimPhaseRequirements, NoneType)
'''
Tuple of both the immutable requirements set tupe *and* that of the singleton
``None`` object.
'''
