#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **simulation pipeline requirement** (i.e., prerequisite
simulation feature required by a runner) functionality.
'''

#FIXME: Generalize this submodule as follows:
#
#* Rename this submodule to "simphasereq".
#* Shift this submodule into the "betse.science.simulate" subpackage.
#* Rename the "SimPipeRunner" prefix of all class names below to simply
#  "SimPhase" (e.g., "SimPhaseRequirement" to "SimPhaseRequirement").
#FIXME: But that's just the beginning. Ultimately, we'd like to permit the
#global requirements defined below to be more readily usable. Currently, they
#require the current "phase" object be passed to their is_satisfied() methods,
#which is cumbersome at best.
#
#To rectify this, first note that whether or not a given phase satisfies a given
#requirement is a constant that should *NEVER* change during that phase. Ergo,
#we can efficiently associate a given phase with all of the following
#requirements as follows:
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
from betse.science.simulate.simphase import SimPhase
from betse.util.type.types import type_check, CallableTypes, SequenceTypes

SimPhase   # ignore IDE warnings

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
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        is_satisfied_body: str,
        set_satisfied_body: str,
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

        All remaining parameters are passed as is to the superclass
        :meth:`SimPhaseRequirement.__init__` method.
        '''

        # Raw string defining the functions to be passed to our superclass.
        func_bodies = '''
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

        # Dictionary mapping from local attribute names to values. Since these
        # functions require no such attributes, the empty dictionary suffices.
        local_attrs = {}

        # Dynamically define these functions.
        exec(func_bodies, globals(), local_attrs)

        # Initialize our superclass with these functions.
        super().__init__(
            *args,
            is_satisfied =local_attrs['is_satisfied'],
            set_satisfied=local_attrs['set_satisfied'],
            **kwargs)

# ....................{ SUBCLASSES                         }....................
class SimPhaseRequirementBoolExpr(SimPhaseRequirementEmbodied):
    '''
    Simulation pipeline requirement initialized by a high-level boolean
    Python expression rather than low-level callables.

    This requirement is a caller convenience simplifying initialization in the
    common case of a requirement reducing to a simple Python expression.
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
    def __init__(self, requirements: SequenceTypes, *args, **kwargs) -> None:
        '''
        Initialize this requirement.

        Parameters
        ----------
        requirements : SequenceTypes
            Sequence of two or more child requirements to be conjuctively
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
    def __init__(self, requirements: SequenceTypes, *args, **kwargs) -> None:
        '''
        Initialize this requirement.

        For convenience, the name of this requirement defaults to that of the
        first passed requirement when the ``name`` parameter is *not* explicitly
        passed.

        Parameters
        ----------
        requirements : SequenceTypes
            Sequence of two or more child requirements to be conjuctively
            composed into this parent requirement.

        All remaining parameters are passed as is to the superclass
        :meth:`SimPhaseRequirementAll.__init__` method.
        '''

        # If no name was passed, default this requirement's name to that of the
        # first passed requirement.
        if 'name' not in kwargs:
            kwargs['name'] = requirements[0].name

        # List of all passed requirements extended by a requirement for the
        # complete BETSE solver. Since the passed sequence of requirements need
        # *NOT* be a list or tuple, the general-purpose extend() method is used.
        requirements_full = [SOLVER_FULL]
        requirements_full.extend(requirements)

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
        super().__init__(*args, requirements=(requirement,), **kwargs)


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

# ....................{ REQUIREMENTS                       }....................
#FIXME: Add a new "FULL" requirement tracking which of the full or fast
#variants of the simulator are enabled. Plots and animations assuming any of the
#following simulation features should be skipped when this requirement is *NOT*
#met: currents, fields, extracellular voltage, and voltage polarity.
#
#Additionally, we need to ensure that the following simulation features are
#disabled when using the fast solver: fluid, deformation, osmosis, ion
#concentrations, pressure.
#FIXME: Exercise this functionality with tests as follows:
#
#* Define a new "test_cli_sim_fast" functional test enabling "p.is_solver_fast".
#* Rename:
#  * "test_cli_sim_noecm" to "test_cli_sim_full_noecm".
#  * "test_cli_sim_ecm" to "test_cli_sim_full_ecm".

SOLVER_FULL = SimPhaseRequirementBoolExpr(
    name='complete BETSE solver', bool_expr='phase.p.is_solver_full')
'''
Requirement that a simulation phase enable the complete BETSE solver.
'''

# ....................{ REQUIREMENTS ~ feature             }....................
DEFORM = SimPhaseRequirementSolverFullAnd(
    requirement=SimPhaseRequirementBoolExpr(
        name='cellular deformation', bool_expr='phase.p.deformation'))
'''
Requirement that a simulation phase enable cellular deformations.
'''


#FIXME: Does this require the full solver? We assume yes, but...
ECM = SimPhaseRequirementBoolExpr(
    name='extracellular spaces', bool_expr='phase.p.is_ecm')
'''
Requirement that a simulation phase enable the extracellular matrix (ECM), also
referred to as "extracellular spaces."
'''


ELECTROOSMOSIS = SimPhaseRequirementSolverFullAnd(
    requirement=SimPhaseRequirementBoolExpr(
        name='electroosmotic flow', bool_expr='phase.p.sim_eosmosis'))
'''
Requirement that a simulation phase enable electroosmotic flow (EOF).
'''


#FIXME: Does this require the full solver? We assume yes, but... Yes!
FLUID = SimPhaseRequirementBoolExpr(
    name='fluid flow', bool_expr='phase.p.fluid_flow')
'''
Requirement that a simulation phase enable fluid flow.
'''

# ....................{ REQUIREMENTS ~ ion                 }....................
ION_CALCIUM = SimPhaseRequirementIon(
    name='calcium ions (Ca2+)', ion_name='Ca')
'''
Requirement that a simulation phase enable calcium ions (Ca2+).
'''


ION_CHLORIDE = SimPhaseRequirementIon(
    name='chloride ions (Cl-)', ion_name='Cl')
'''
Requirement that a simulation phase enable chloride ions (Cl-).
'''


ION_POTASSIUM = SimPhaseRequirementIon(
    name='potassium ions (K+)', ion_name='K')
'''
Requirement that a simulation phase enable potassium ions (K+).
'''


ION_M_ANION = SimPhaseRequirementIon(
    name='M anions (M-)', ion_name='M')
'''
Requirement that a simulation phase enable M anions (M-).
'''


ION_SODIUM = SimPhaseRequirementIon(
    name='sodium ions (Na+)', ion_name='Na')
'''
Requirement that a simulation phase enable sodium ions (Na+).
'''

# ....................{ REQUIREMENTS ~ pressure            }....................
#FIXME: Does this require the full solver? We assume yes, but... Yes!
PRESSURE_OSMOTIC = SimPhaseRequirementBoolExpr(
    name='osmotic pressure', bool_expr='phase.p.deform_osmo',)
'''
Requirement that a simulation phase enable osmotic pressure.
'''


#FIXME: Does this require the full solver? We assume yes, but... Yes!
PRESSURE_TOTAL = SimPhaseRequirement(
    name='total pressure',
    is_satisfied=lambda phase:
        phase.p.deform_osmo or phase.p.scheduled_options['pressure'] != 0,

    # For simplicity, define this requirement to be settable by enabling osmotic
    # pressure. While the mechanical pressure event could also be enabled, doing
    # so is less trivial than the former.
    set_satisfied=lambda phase:
        phase.p.__setattr__('deform_osmo', True),
)
'''
Requirement that a simulation phase enable at least one pressure feature:
namely, osmotic pressure or the mechanical pressure intervention.
'''

# ....................{ REQUIREMENTS ~ voltage             }....................
VOLTAGE_MEMBRANE_GHK = SimPhaseRequirementSolverFullAnd(
    requirement=SimPhaseRequirementBoolExpr(
        name='Goldman-Hodgkin-Katz (GHK) calculation',
        bool_expr='phase.p.GHK_calc'))
'''
Requirement that a simulation phase enable alternative calculation of
transmembrane voltages (Vmem) given the Goldman-Hodgkin-Katz (GHK) equation.
'''


VOLTAGE_POLARITY = SimPhaseRequirementSolverFullAnd(
    requirement=SimPhaseRequirement(
        name='cellular voltage polarizability',
        is_satisfied =lambda phase: phase.p.cell_polarizability > 0,
        set_satisfied=lambda phase: phase.p.__setattr__(
            'cell_polarizability', 1e-4)))
'''
Requirement that a simulation phase enable cellular voltage polarizability.
'''
