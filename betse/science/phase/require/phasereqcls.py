#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **simulation pipeline requirement** (i.e., prerequisite
simulation feature required by a runner) functionality.
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
from betse.science.phase.phasecls import SimPhase
# from betse.util.io.log import logs
from betse.util.type.text import strs
from betse.util.type.types import (
    type_check,
    CallableTypes,
    EnumMemberType,
    IterableTypes,
    MappingOrNoneTypes,
)

# ....................{ SUPERCLASSES                       }....................
class SimPhaseRequirement(object):
    '''
    **Simulation pipeline requirement** (i.e., object encapsulating the
    current state of a single simulation feature, such as extracellular spaces,
    required by one or more simulation pipeline runners).

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
        is_satisfied: CallableTypes,
        set_satisfied: CallableTypes,
        name: str,
    ) -> None:
        '''
        Initialize this requirement.

        Parameters
        ----------
        is_satisfied : CallableTypes
            Callable (e.g., function, lambda) passed only the current
            simulation phase, returning ``True`` only if this requirement is
            enabled in this phase.
        set_satisfied : CallableTypes
            Callable (e.g., function, lambda) passed only the current
            simulation phase, enabling this requirement in this phase.
        name : str
            Human-readable name of this requirement (e.g., ``extracellular
            spaces``). For simplicity, the first character of this name should
            typically be lower- rather than uppercase.
        '''

        # Classify all passed parameters.
        self.is_satisfied = is_satisfied
        self.set_satisfied = set_satisfied
        self.name = name


class SimPhaseRequirementEmbodied(SimPhaseRequirement):
    '''
    Simulation pipeline requirement initialized by high-level Python
    function bodies rather than low-level callables.

    This requirement is a caller convenience simplifying initialization in the
    common case of a requirement reducing to two function bodies.

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
            String of one or more arbitrary complex Python statements comprising
            the body of a dynamically defined function:
            * Passed only the current simulation phase.
            * Returning ``True`` only if this requirement is enabled in this
              phase.
        set_satisfied_body: str
            String of one or more arbitrary complex Python statements comprising
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

# ....................{ SUBCLASSES                         }....................
class SimPhaseRequirementBoolExpr(SimPhaseRequirementEmbodied):
    '''
    Simulation pipeline requirement initialized by a high-level **boolean
    Python expression** (i.e., both gettable and settable as a boolean value)
    rather than a pair of low-level callables.

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
    Simulation pipeline requirement initialized by a high-level **enumeration
    Python expression** (i.e., both gettable and settable as an enumeration
    value) rather than a pair of low-level callables.

    This requirement is a caller convenience simplifying initialization in the
    common case of a requirement reducing to a single enumeration.
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

# ....................{ SUBCLASSES ~ all                   }....................
class SimPhaseRequirementAll(SimPhaseRequirement):
    '''
    Parent simulation pipeline requirement requiring two or more child
    requirements to be **conjunctively satisfied** (i.e., to *all* be
    satisfied by a given simulation phase).

    This requirement permits multiple requirements to be composed together into
    an implicit logical conjunction expressing the ``and`` operation.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, requirements: IterableTypes, *args, **kwargs) -> None:
        '''
        Initialize this requirement.

        Parameters
        ----------
        requirements : IterableTypes
            Iterable of two or more child requirements to be conjuctively
            composed into this parent requirement.

        All remaining parameters are passed as is to the superclass
        :meth:`SimPhaseRequirementEmbodied.__init__` method.
        '''

        @type_check
        def is_satisfied(phase: SimPhase) -> bool:
            '''
            ``True`` only if all child requirements composed by this parent
            requirement are satisfied by the passed simulation phase.
            '''

            # This requirement is satisfied if and only if...
            return all(
                # This child requirement is satisfied by this phase.
                requirement.is_satisfied(phase)
                # For each child requirement of this parent requirement...
                for requirement in requirements)


        @type_check
        def set_satisfied(phase: SimPhase) -> None:
            '''
            Modify the passed simulation phase as needed to ensure that all
            child requirements composed by this parent requirement are satisfied
            by this phase.
            '''

            # For each child requirement of this parent requirement...
            for requirement in requirements:
                # Ensure this child requirement is satisfied by this phase.
                requirement.set_satisfied(phase)


        # Initialize our superclass with all passed parameters.
        super().__init__(
            *args,
            is_satisfied=is_satisfied,
            set_satisfied=set_satisfied,
            **kwargs)

# ....................{ SUBCLASSES ~ all : solver          }....................
class SimPhaseRequirementSolverFullAll(SimPhaseRequirementAll):
    '''
    Parent simulation pipeline requirement requiring both the complete BETSE
    solver *and* two or more passed child requirements to be **conjunctively
    satisfied** (i.e., to *all* be satisfied by a given simulation phase).

    This requirement is a caller convenience simplifying initialization in the
    common case of multiple requirements requiring the complete BETSE solver.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, requirements: IterableTypes, *args, **kwargs) -> None:
        '''
        Initialize this requirement.

        If unpassed, the ``name`` parameter defaults to the human-readable
        conjunction of the names of all passed requirements (e.g., "full solver,
        fluid flow, and calcium ions (Ca2+)").

        Parameters
        ----------
        requirements : IterableTypes
            Iterable of two or more child requirements to be conjuctively
            composed into this parent requirement.

        All remaining parameters are passed as is to the superclass
        :meth:`SimPhaseRequirementAll.__init__` method.
        '''

        # Avoid circular import dependencies.
        from betse.science.phase.require.phasereqs import SOLVER_FULL

        # If this requirement's name was *NOT* passed...
        if 'name' not in kwargs:
            # Generator yielding the name of each passed requirement.
            requirement_names = (
                requirement.name for requirement in requirements)

            # Default this requirement's name to a concatenation of these names.
            kwargs['name'] = strs.join_as_conjunction(*requirement_names)

        # Set of all passed requirements extended by a requirement for the
        # complete BETSE solver. Since the passed sequence of requirements need
        # *NOT* be a list or tuple, the general-purpose extend() method is used.
        requirements_full = {SOLVER_FULL,} | set(requirements)

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, requirements=requirements_full, **kwargs)


class SimPhaseRequirementSolverFullAnd(SimPhaseRequirementSolverFullAll):
    '''
    Parent simulation pipeline requirement requiring both the complete BETSE
    solver *and* a single passed child requirement to be **conjunctively
    satisfied** (i.e., to *all* be satisfied by a given simulation phase).

    This requirement is a caller convenience simplifying initialization in the
    common case of a single requirement requiring the complete BETSE solver.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self, requirement: SimPhaseRequirement, *args, **kwargs) -> None:
        '''
        Initialize this requirement.

        For convenience, the name of this requirement defaults to that of the
        passed requirement if the ``name`` parameter is *not* explicitly passed.

        Parameters
        ----------
        requirement : SimPhaseRequirement
            Child requirement to be conjuctively composed into this parent
            requirement.

        All remaining parameters are passed as is to the superclass
        :meth:`SimPhaseRequirementSolverFullAll.__init__` method.
        '''

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, requirements={requirement,}, **kwargs)


class SimPhaseRequirementIon(SimPhaseRequirementSolverFullAnd):
    '''
    Simulation pipeline requirement initialized by a high-level ion name
    rather than low-level callables.

    This requirement is a caller convenience simplifying initialization in the
    common case of a requirement reducing to a single ion. Since simulation of
    ion concentrations requires the complete BETSE solver, this requirement is
    the conjunction of that requirement *and* a custom requirement specific to
    this ion.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, ion_name: str, *args, **kwargs) -> None:
        '''
        Initialize this requirement.

        Parameters
        ----------
        ion_name : str
            Name of the required ion, equivalent to a key of the
            :attr:`Parameters.ions_dict` dictionary (e.g., ``Ca``, requiring
            calcium ions).

        All remaining parameters are passed as is to the superclass
        :meth:`SimPhaseRequirementAll.__init__` method.
        '''

        # Requirement that the passed ion is satisfied by a given phase.
        ion_requirement = SimPhaseRequirementEmbodied(
            name='{} ions'.format(ion_name),
            is_satisfied_body=(
                'return phase.p.ions_dict[{!r}] == 1'.format(ion_name)),
            set_satisfied_body=(
                'phase.p.ions_dict[{!r}] == 1'.format(ion_name)),
        )

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, requirement=ion_requirement, **kwargs)
