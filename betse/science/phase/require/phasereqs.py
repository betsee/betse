#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Simulation phase requirement sets** (i.e., immutable sets of simulation phase
requirements, requiring zero or more arbitrary requirements to be satisfied).

Design
----------
This submodule *only* collects submodule-scope constants, all of which are
instances of the :class:`SimPhaseRequirements` class.
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
from betse.science.phase.require import phasereqsmake
from betse.science.phase.require.abc.phasereqset import SimPhaseRequirements
from betse.science.config.confenum import SolverType

# ....................{ REQUIREMENTS                       }....................
NONE = SimPhaseRequirements()
'''
**Null requirements** (i.e., requirements unconditionally satisfied by all
simulation phases regardless of which simulation features those phases enable).
'''


ELECTRIC_CURRENT = NONE
'''
Requirements that a simulation phase enable electrical currents and hence at
least intracellular electrical currents.

Simulation phases satisfying this requirement do *not* necessarily enable
extracellular electrical currents, which additionally require the full solver
to be enabled.
'''


ELECTRIC_FIELD = NONE
'''
Requirements that a simulation phase enable electrical fields and hence at
least intracellular electrical fields.

Simulation phases satisfying this requirement do *not* necessarily enable
extracellular electrical fields, which additionally require the full solver to
be enabled.
'''


MICROTUBULE = NONE
'''
Requirements that a simulation phase enable microtubules.
'''

# ....................{ REQUIREMENTS ~ bool                }....................
# Requirements satisfied by enabling a requirement-specific boolean.

ECM = phasereqsmake.make_requirements_bool_expr(
    name='extracellular spaces', bool_expr='phase.p.is_ecm')
'''
Requirements that a simulation phase enable the extracellular matrix (ECM), also
referred to as "extracellular spaces."

Note that this requirement does *not* itself require the full BETSE solver.
While the equivalent circuit-based solver mostly ignores extracellular spaces,
gene regulatory networks (GRNs) supported by this solver explicitly support
extracellular spaces. Ergo, this requirement is solver-agnostic.
'''

# ....................{ REQUIREMENTS ~ solver              }....................
#FIXME: Add a new "FULL" requirement tracking which of the full or fast
#variants of the simulator are enabled. Plots and animations assuming any of the
#following simulation features should be skipped when this requirement is *NOT*
#met: currents, fields, extracellular voltage, and voltage polarity.
#
#Additionally, we need to ensure that the following simulation features are
#disabled when using the fast solver: fluid, deformation, osmosis, ion
#concentrations, pressure.

# This solver is internally referenced by solver classes instantiated below
# (e.g., "SimPhaseRequirementsIon", "SimPhaseRequirementsSolverFullAnd") and
# *MUST* thus be declared first -- before these classes are instantiated.

SOLVER_FULL = phasereqsmake.make_requirements_enum_expr(
    name='full BETSE solver',
    enum_expr='phase.p.solver_type',
    enum_member=SolverType.FULL,
)
'''
Requirements that a simulation phase enable the full BETSE solver.
'''

# ....................{ REQUIREMENTS ~ solver : equal      }....................
# Requirements satisfied by enabling only the full solver. Technically, these
# requirements are currently identical to the previously defined "SOLVER_FULL"
# requirement. This need *NOT* necessarily be the case, however. Theoretically,
# subsequent versions of this application could expose boolean settings in the
# simulation configuration file permitting these currently unconditional
# simulation features to be selectively disabled.


ELECTRIC_CURRENT_EXTRA = SOLVER_FULL
'''
Requirements that a simulation phase enable extracellular electrical currents.

Simulation phases satisfying this requirement necessarily enable the full solver
but *not* extracellular spaces, which are *not* required to simulate
extracellular electrical currents.
'''


ELECTRIC_CURRENT_MEMBRANE = SOLVER_FULL
'''
Requirements that a simulation phase enable transmembrane electrical currents.

Simulation phases satisfying this requirement necessarily enable the full solver
but *not* extracellular spaces, which are *not* required to simulate
transmembrane electrical currents.
'''


PUMP_NAKATPASE = SOLVER_FULL
'''
Requirements that a simulation phase enable the Na-K-ATPase membrane pump.
'''

# ....................{ REQUIREMENTS ~ solver : equal      }....................
# Requirements satisfied by enabling only both the full solver *AND*
# extracellular spaces. See above for comments.

SOLVER_FULL_ECM = SOLVER_FULL | ECM
'''
Requirements that a simulation phase enable both extracellular spaces *and* a
solver able to fully simulate these spaces (e.g., the full BETSE solver).
'''


ELECTRIC_FIELD_EXTRA = SOLVER_FULL_ECM
'''
Requirements that a simulation phase enable extracellular electrical fields and
hence both extracellular spaces *and* a solver capable of fully simulating these
spaces (e.g., the full BETSE solver).
'''


VOLTAGE_EXTRA = SOLVER_FULL
'''
Requirements that a simulation phase enable extracellular voltages and hence
both extracellular spaces *and* a solver capable of fully simulating these
spaces (e.g., the full BETSE solver).
'''

# ....................{ REQUIREMENTS ~ solver : bool       }....................
# Requirements satisfied by enabling only the full solver and a
# requirement-specific boolean.

DEFORM = SOLVER_FULL | phasereqsmake.make_requirements_bool_expr(
    name='cellular deformation', bool_expr='phase.p.deformation')
'''
Requirements that a simulation phase enable cellular deformations.
'''


ELECTROOSMOSIS = SOLVER_FULL | phasereqsmake.make_requirements_bool_expr(
    name='electroosmotic flow', bool_expr='phase.p.sim_eosmosis')
'''
Requirements that a simulation phase enable electroosmotic flow (EOF).
'''


FLUID = SOLVER_FULL | phasereqsmake.make_requirements_bool_expr(
    name='fluid flow', bool_expr='phase.p.fluid_flow')
'''
Requirements that a simulation phase enable fluid flow and hence at least
intracellular fluid flow.

Simulation phases satisfying this requirement do *not* necessarily enable
extracellular fluid flow, which additionally require extracellular spaces to be
enabled.
'''


PRESSURE_OSMOTIC = SOLVER_FULL | phasereqsmake.make_requirements_bool_expr(
    name='osmotic pressure', bool_expr='phase.p.deform_osmo')
'''
Requirements that a simulation phase enable osmotic pressure.
'''


VOLTAGE_MEMBRANE_GHK = SOLVER_FULL | phasereqsmake.make_requirements_bool_expr(
    name='Goldman-Hodgkin-Katz (GHK) calculation', bool_expr='phase.p.GHK_calc')
'''
Requirements that a simulation phase enable alternative calculation of
transmembrane voltages (Vmem) given the Goldman-Hodgkin-Katz (GHK) equation.
'''

# ....................{ REQUIREMENTS ~ solver : ion        }....................
# Requirements satisfied by enabling only the full solver and a
# requirement-specific ion.

ION_CALCIUM = phasereqsmake.make_requirements_ion(
    name='calcium ions (Ca2+)', ion_name='Ca')
'''
Requirements that a simulation phase enable calcium ions (Ca2+).
'''


ION_CHLORIDE = phasereqsmake.make_requirements_ion(
    name='chloride ions (Cl-)', ion_name='Cl')
'''
Requirements that a simulation phase enable chloride ions (Cl-).
'''


ION_POTASSIUM = phasereqsmake.make_requirements_ion(
    name='potassium ions (K+)', ion_name='K')
'''
Requirements that a simulation phase enable potassium ions (K+).
'''


ION_M_ANION = phasereqsmake.make_requirements_ion(
    name='M anions (M-)', ion_name='M')
'''
Requirements that a simulation phase enable M anions (M-).
'''


ION_SODIUM = phasereqsmake.make_requirements_ion(
    name='sodium ions (Na+)', ion_name='Na')
'''
Requirements that a simulation phase enable sodium ions (Na+).
'''

# ....................{ REQUIREMENTS ~ solver : lambda     }....................
# Requirements satisfied by enabling only the full solver and a pair of
# requirement-specific lambda functions.

PRESSURE_TOTAL = SOLVER_FULL | phasereqsmake.make_requirements_funcs(
    name='total pressure',
    is_satisfied=lambda phase:
        phase.p.deform_osmo or phase.p.scheduled_options['pressure'] != 0,

    # For simplicity, this requirement is defined to be settable by enabling
    # osmotic pressure. While the mechanical pressure event could also be
    # enabled, doing so is less trivial than the former.
    set_satisfied=lambda phase: phase.p.__setattr__('deform_osmo', True))
'''
Requirements that a simulation phase enable at least one pressure feature:
namely, osmotic pressure or the mechanical pressure intervention.
'''


VOLTAGE_POLARITY = SOLVER_FULL | phasereqsmake.make_requirements_funcs(
    name='cellular voltage polarizability',
    is_satisfied =lambda phase: phase.p.cell_polarizability > 0,
    set_satisfied=lambda phase: phase.p.__setattr__(
        'cell_polarizability', 1e-4))
'''
Requirements that a simulation phase enable cellular voltage polarizability.
'''

# ....................{ REQUIREMENTS ~ solver : ecm        }....................
# Requirements satisfied by enabling only the full solver, extracellular spaces,
# and a single arbitrary requirement.

FLUID_EXTRA = FLUID | ECM
'''
Requirements that a simulation phase enable extracellular fluid flow and hence
fluid flow, extracellular spaces, *and* a solver capable of fully simulating
these features (e.g., the full BETSE solver).
'''


ION_CALCIUM_EXTRA = ION_CALCIUM | ECM
'''
Requirements that a simulation phase enable extracellular calcium ions (Ca2+)
and hence calcium ions, extracellular spaces, *and* a solver capable of fully
simulating these features (e.g., the full BETSE solver).
'''
