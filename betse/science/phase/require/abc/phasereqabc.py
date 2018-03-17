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
from abc import ABCMeta
from betse.science.phase.phasecls import SimPhase
# from betse.util.io.log import logs
from betse.util.type.decorator.deccls import abstractproperty
from betse.util.type.types import (
    type_check,
    CallableTypes,
    EnumMemberType,
    MappingOrNoneTypes,
)

# ....................{ SUPERCLASSES                       }....................
class SimPhaseRequirementABC(metaclass=ABCMeta):
    '''
    Abstract base class of all **simulation phase requirement** (i.e., object
    encapsulating the current state of a single simulation feature for a given
    simulation phase) subclasses.
    '''

    # ..................{ PROPERTIES ~ abstract              }..................
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

    To preserve hashability and hence usability in dictionaries and sets, this
    requirement is immutable. Callers should avoid breaking encapsulation to
    modify this requirement.

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

        # Initialize our superclass.
        super().__init__()

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
