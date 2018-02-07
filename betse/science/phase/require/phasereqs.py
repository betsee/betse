#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **simulation phase requirement** (i.e., optional simulation feature
required for a given phase by an external caller) functionality.
'''

#FIXME: *WOOPS.* Sadly, our first-pass attempt at an API for composing multiple
#requirements into a single requirement via the "SimPhaseRequirementAll" class
#are an abject failure. While this class does superficially work, it fails at a
#deeper level in several significant ways:
#
#* It's impossible to sanely query. This is the critical failure. Consider a
#  pipeline runner requirements set resembling
#  "requirements={phasereqs.VOLTAGE_EXTRA,}". Ideally, the "VOLTAGE_EXTRA"
#  requirement should itself require the "EXTRA" requirement via the
#  "SimPhaseRequirementAll" class. If we were to do so, however, the above
#  requirements set could no longer be trivially queried as to whether this set
#  enables extracellur spaces. Ergo, we currently write this set as
#  "requirements={phasereqs.VOLTAGE_EXTRA,}" instead. Clearly, this is awful.
#* It's inefficient. Consider a pipeline runner requirements set resembling
#  "requirements={phasereqs.DEFORM, phasereqs.FLUID,}". What's the issue here?
#  The "DEFORM" and "FLUID" requirements both require the "SOLVER_FULL"
#  requirement, thus duplicating that same shared requirement.
#
#Addressing these concerns requires (...get it?) a fundamental rethink of this
#API. In particular, it's increasingly clear that each and every requirement
#should be a *SET* of all "SimPhaseRequirement" objects required by that
#requirement. While some requirements would reduce to sets of only a single
#requirement, most would comprise two or more requirements.
#
#To implement this cleanly, we probably want to pursue the following:
#
#* Rename this submodule to "phasereqs".
#* Define a new "SimPhaseRequirements" subclass in the "phasereqcls" submodule.
#  This subclass should subclass an immutable set superclass. This might be the
#  builtin "set" type *OR* this might be a more appropriate ABC. Investigate.
#* Implement the __init__() method of this subclass to validate that all passed
#  set items are requirements: e.g.,
#
#       from betse.util.type import iterables
#
#       # If any item of this set is *NOT* a requirement, raise an exception.
#       iterables.die_unless_items_instance_of(
#           iterable=requirements, cls=SimPhaseRequirement)
#
#* Refactor all globals defined below to be instances of the
#  "SimPhaseRequirements" subclass whose items are all "SimPhaseRequirements"
#  instances. For example, "DEFORM" would then reduce to:
#
#    DEFORM = SimPhaseRequirements(SOLVER_FULL, SimPhaseRequirementBoolExpr(
#        name='cellular deformation', bool_expr='phase.p.deformation'))
#* Refactor all uses of these globals to acknowledge their set-like nature. For
#  example:
#
#    # Change this...
#    requirements={phasereqs.VOLTAGE_EXTRA, phasereqs.ECM,},
#
#    # ...into this.
#    requirements=phasereqs.VOLTAGE_EXTRA,
#
#Naturally, multiple requirements would then be trivially composable via the set
#union operator "|". This is the future. See to it, please.

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
from betse.science.phase.require.phasereqcls import (
    SimPhaseRequirement,
    SimPhaseRequirementBoolExpr,
    SimPhaseRequirementEnumExpr,
    SimPhaseRequirementIon,
    SimPhaseRequirementSolverFullAnd,
)
from betse.science.config.confenum import SolverType

# ....................{ REQUIREMENTS ~ bool                }....................
# Requirements satisfied by enabling a requirement-specific boolean.

ECM = SimPhaseRequirementBoolExpr(
    name='extracellular spaces', bool_expr='phase.p.is_ecm')
'''
Requirement that a simulation phase enable the extracellular matrix (ECM), also
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
# (e.g., "SimPhaseRequirementIon", "SimPhaseRequirementSolverFullAnd") and
# *MUST* thus be declared first -- before these classes are instantiated.
SOLVER_FULL = SimPhaseRequirementEnumExpr(
    name='full BETSE solver',
    enum_expr='phase.p.solver_type',
    enum_member=SolverType.FULL,
)
'''
Requirement that a simulation phase enable the full BETSE solver.
'''

# ....................{ REQUIREMENTS ~ solver : equal      }....................
# Requirements satisfied by enabling only the full solver. Technically, these
# requirements are currently identical to the previously defined "SOLVER_FULL"
# requirement. This need *NOT* necessarily be the case, however. Theoretically,
# subsequent versions of this application could expose boolean settings in the
# simulation configuration file permitting these currently unconditional
# simulation features to be selectively disabled.

ELECTRIC_CURRENT = SOLVER_FULL
'''
Requirement that a simulation phase enable electrical currents.
'''


ELECTRIC_FIELD = SOLVER_FULL
'''
Requirement that a simulation phase enable electrical fields.
'''


MICROTUBULE = SOLVER_FULL
'''
Requirement that a simulation phase enable microtubules.
'''


PUMP_NAKATPASE = SOLVER_FULL
'''
Requirement that a simulation phase enable the Na-K-ATPase membrane pump.
'''


#FIXME: Improve this parent requirement to also require "ECM" *AFTER* addressing
#the over-arching "FIXME" concern above.
VOLTAGE_EXTRA = SOLVER_FULL
'''
Requirement that a simulation phase enable extracellular voltages and hence both
extracellular spaces *and* a solver capable of fully simulating these spaces.
'''

# ....................{ REQUIREMENTS ~ solver : bool       }....................
# Requirements satisfied by enabling only the full solver and a
# requirement-specific boolean.

DEFORM = SimPhaseRequirementSolverFullAnd(
    requirement=SimPhaseRequirementBoolExpr(
        name='cellular deformation', bool_expr='phase.p.deformation'))
'''
Requirement that a simulation phase enable cellular deformations.
'''


ELECTROOSMOSIS = SimPhaseRequirementSolverFullAnd(
    requirement=SimPhaseRequirementBoolExpr(
        name='electroosmotic flow', bool_expr='phase.p.sim_eosmosis'))
'''
Requirement that a simulation phase enable electroosmotic flow (EOF).
'''


FLUID = SimPhaseRequirementSolverFullAnd(
    SimPhaseRequirementBoolExpr(
        name='fluid flow', bool_expr='phase.p.fluid_flow'))
'''
Requirement that a simulation phase enable fluid flow.
'''


PRESSURE_OSMOTIC = SimPhaseRequirementSolverFullAnd(
    requirement=SimPhaseRequirementBoolExpr(
        name='osmotic pressure', bool_expr='phase.p.deform_osmo'))
'''
Requirement that a simulation phase enable osmotic pressure.
'''


VOLTAGE_MEMBRANE_GHK = SimPhaseRequirementSolverFullAnd(
    requirement=SimPhaseRequirementBoolExpr(
        name='Goldman-Hodgkin-Katz (GHK) calculation',
        bool_expr='phase.p.GHK_calc'))
'''
Requirement that a simulation phase enable alternative calculation of
transmembrane voltages (Vmem) given the Goldman-Hodgkin-Katz (GHK) equation.
'''

# ....................{ REQUIREMENTS ~ solver : lambda     }....................
# Requirements satisfied by enabling only the full solver and a pair of
# requirement-specific lambda functions.

PRESSURE_TOTAL = SimPhaseRequirementSolverFullAnd(
    requirement=SimPhaseRequirement(
        name='total pressure',
        is_satisfied=lambda phase:
            phase.p.deform_osmo or phase.p.scheduled_options['pressure'] != 0,

        # For simplicity, this requirement is defined to be settable by enabling
        # osmotic pressure. While the mechanical pressure event could also be
        # enabled, doing so is less trivial than the former.
        set_satisfied=lambda phase: phase.p.__setattr__('deform_osmo', True)))
'''
Requirement that a simulation phase enable at least one pressure feature:
namely, osmotic pressure or the mechanical pressure intervention.
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

# ....................{ REQUIREMENTS ~ solver : ion        }....................
# Requirements satisfied by enabling only the full solver and a
# requirement-specific ion.

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
