#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **simulation pipeline runner requirement** (i.e., prerequisite
simulation feature required by a runner) functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.science.simulate.simphase import SimPhase
from betse.util.type.types import type_check, CallableTypes

if False: SimPhase   # ignore IDE warnings

# ....................{ SUPERCLASSES                       }....................
class SimPipeRunnerRequirement(object):
    '''
    **Simulation pipeline runner requirement** (i.e., object encapsulating the
    current state of a single simulation feature, such as extracellular spaces,
    required by one or more simulation pipeline runners).

    Attributes
    ----------
    is_satisfied : CallableTypes
        Callable (e.g., function, lambda) passed only the current simulation
        phase, returning ``True`` only if this requirement is enabled in this
        phase.
    set_satisfied : CallableTypes
        Callable (e.g., function, lambda) passed only the current simulation
        phase, enabling this requirement in this phase.
    name : str
        Human-readable lowercase name of this requirement (e.g.,
        ``extracellular spaces``).
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        is_satisfied: CallableTypes,
        set_satisfied: CallableTypes,
        name: str,
        **kwargs
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
            Human-readable lowercase name of this requirement (e.g.,
            ``extracellular spaces``).
        '''

        # Classify all passed parameters.
        self.is_satisfied = is_satisfied
        self.set_satisfied = set_satisfied
        self.name = name


class SimPipeRunnerRequirementBodies(SimPipeRunnerRequirement):
    '''
    Simulation pipeline runner requirement initialized by high-level Python
    function bodies rather than low-level callables.

    This subclass is a caller convenience simplifying initialization in the
    common case of a runner requirement reducing to two function bodies.
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
        :meth:`SimPipeRunnerRequirement.__init__` method.
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
class SimPipeRunnerRequirementBoolExpr(SimPipeRunnerRequirementBodies):
    '''
    Simulation pipeline runner requirement initialized by a high-level boolean
    Python expression rather than low-level callables.

    This subclass is a caller convenience simplifying initialization in the
    common case of a runner requirement reducing to a simple Python expression.
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
        :meth:`SimPipeRunnerRequirement.__init__` method.
        '''

        # Initialize our superclass with all passed parameters.
        super().__init__(
            *args,
            is_satisfied_body='return ' + bool_expr,
            set_satisfied_body=bool_expr + ' = True',
            **kwargs)


class SimPipeRunnerRequirementIon(SimPipeRunnerRequirementBodies):
    '''
    Simulation pipeline runner requirement initialized by a high-level ion name
    rather than low-level callables.

    This subclass is a caller convenience simplifying initialization in the
    common case of a runner requirement reducing to a single ion.
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
        :meth:`SimPipeRunnerRequirement.__init__` method.
        '''

        # Initialize our superclass with all passed parameters.
        super().__init__(
            *args,
            is_satisfied_body=(
                'return phase.p.ions_dict[{!r}] == 1'.format(ion_name)),
            set_satisfied_body=(
                'phase.p.ions_dict[{!r}] == 1'.format(ion_name)),
            **kwargs)

# ....................{ REQUIREMENTS                       }....................
DEFORM = SimPipeRunnerRequirementBoolExpr(
    name='cellular deformation', bool_expr='phase.p.deformation')
'''
Requirement that a simulation phase enable cellular deformations.
'''


ECM = SimPipeRunnerRequirementBoolExpr(
    name='extracellular spaces', bool_expr='phase.p.is_ecm')
'''
Requirement that a simulation phase enable the extracellular matrix (ECM), also
referred to as "extracellular spaces."
'''


ELECTROOSMOSIS = SimPipeRunnerRequirementBoolExpr(
    name='electroosmotic flow', bool_expr='phase.p.sim_eosmosis')
'''
Requirement that a simulation phase enable electroosmotic flow (EOF).
'''


FLUID = SimPipeRunnerRequirementBoolExpr(
    name='fluid flow', bool_expr='phase.p.fluid_flow')
'''
Requirement that a simulation phase enable fluid flow.
'''

# ....................{ REQUIREMENTS ~ ion                 }....................
ION_CALCIUM = SimPipeRunnerRequirementIon(
    name='calcium ions (Ca2+)', ion_name='Ca')
'''
Requirement that a simulation phase enable calcium ions (Ca2+).
'''


ION_CHLORIDE = SimPipeRunnerRequirementIon(
    name='chloride ions (Cl-)', ion_name='Cl')
'''
Requirement that a simulation phase enable chloride ions (Cl-).
'''


ION_HYDROGEN = SimPipeRunnerRequirementIon(
    name='hydrogen ions (H+)', ion_name='H')
'''
Requirement that a simulation phase enable hydrogen ions (H+).
'''


ION_POTASSIUM = SimPipeRunnerRequirementIon(
    name='potassium ions (K+)', ion_name='K')
'''
Requirement that a simulation phase enable potassium ions (K+).
'''


ION_M_ANION = SimPipeRunnerRequirementIon(
    name='M anions (M-)', ion_name='M')
'''
Requirement that a simulation phase enable M anions (M-).
'''


ION_SODIUM = SimPipeRunnerRequirementIon(
    name='sodium ions (Na+)', ion_name='Na')
'''
Requirement that a simulation phase enable sodium ions (Na+).
'''

# ....................{ REQUIREMENTS ~ pressure            }....................
PRESSURE_OSMOTIC = SimPipeRunnerRequirementBoolExpr(
    name='osmotic pressure', bool_expr='phase.p.deform_osmo',)
'''
Requirement that a simulation phase enable osmotic pressure.
'''


PRESSURE_TOTAL = SimPipeRunnerRequirement(
    name='total pressure',
    is_satisfied =lambda phase:
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
VOLTAGE_MEMBRANE_GHK = SimPipeRunnerRequirementBoolExpr(
    name='Goldman-Hodgkin-Katz (GHK) calculation', bool_expr='phase.p.GHK_calc')
'''
Requirement that a simulation phase enable alternative calculation of
transmembrane voltages (Vmem) given the Goldman-Hodgkin-Katz (GHK) equation.
'''


VOLTAGE_POLARITY = SimPipeRunnerRequirement(
    name='cellular voltage polarizability',
    is_satisfied =lambda phase: phase.p.cell_polarizability > 0,
    set_satisfied=lambda phase: phase.p.__setattr__('cell_polarizability', 1e-4),
)
'''
Requirement that a simulation phase enable cellular voltage polarizability.
'''
