#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Utility functions simplifying creation of **simulation phase requirements**
(i.e., objects encapsulating the current state of a single simulation feature
for a given simulation phase).
'''

# ....................{ IMPORTS                            }....................
from betse.science.phase.require.phasereqcls import (
    SimPhaseRequirementABC,
    SimPhaseRequirementEmbodied,
    SimPhaseRequirementBoolExpr,
    SimPhaseRequirementEnumExpr,
    SimPhaseRequirements,
)
# from betse.util.io.log import logs
from betse.util.type.types import type_check, IterableTypes

# ....................{ MAKERS ~ one                       }....................
@type_check
def make_requirements_one(
    requirement: SimPhaseRequirementABC) -> SimPhaseRequirements:
    '''
    Immutable set of only one arbitrary simulation phase requirement.

    Parameters
    ----------
    requirement : SimPhaseRequirementABC
        Arbitrary requirement defining this set.
    '''

    # Create and return this requirement set from this requirement.
    return SimPhaseRequirements(iterable=(requirement,))


@type_check
def make_requirements_bool_expr(*args, **kwargs) -> SimPhaseRequirements:
    '''
    Immutable set of only one string-based boolean simulation phase requirement.

    All passed parameters are passed to the
    :meth:`SimPhaseRequirementBoolExpr.__init__` method as is.
    '''

    # String-based boolean simulation phase requirement.
    requirement = SimPhaseRequirementBoolExpr(*args, **kwargs)

    # Create and return this requirement set from this requirement.
    return make_requirements_one(requirement)


@type_check
def make_requirements_enum_expr(*args, **kwargs) -> SimPhaseRequirements:
    '''
    Immutable set of only one string-based enumeration member simulation phase
    requirement.

    All passed parameters are passed to the
    :meth:`SimPhaseRequirementEnumExpr.__init__` method as is.
    '''

    # String-based enumeration member simulation phase requirement.
    requirement = SimPhaseRequirementEnumExpr(*args, **kwargs)

    # Create and return this requirement set from this requirement.
    return make_requirements_one(requirement)

# ....................{ MAKERS ~ solver : ecm              }....................
@type_check
def make_requirements_ecm(requirements: IterableTypes) -> SimPhaseRequirements:
    '''
    Immutable set of extracellular spaces-predicate simulation phase
    requirements, requiring both extracellular spaces *and* one or more
    arbitrary requirements to be satisfied.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable of one or more arbitrary requirements defining this set.
    '''

    # Avoid circular import dependencies.
    from betse.science.phase.require.phasereqs import ECM

    # List requiring both extracellular spaces *AND* all passed requirements.
    requirements_full = [ECM]
    requirements_full.extend(requirements)

    # Create and return this requirement set from this list.
    return SimPhaseRequirements(requirements_full)


@type_check
def make_requirements_ecm_and(
    requirement: SimPhaseRequirementABC) -> SimPhaseRequirements:
    '''
    Immutable set of extracellular spaces-based simulation phase requirements,
    requiring both extracellular spaces *and* a single arbitrary requirement to
    be satisfied.

    Parameters
    ----------
    requirement : SimPhaseRequirementABC
        Arbitrary requirement defining this set.
    '''

    return make_requirements_ecm(requirements=(requirement,))

# ....................{ MAKERS ~ solver : full             }....................
@type_check
def make_requirements_solver_full(
    requirements: IterableTypes) -> SimPhaseRequirements:
    '''
    Immutable set of full solver-predicate simulation phase requirements,
    requiring both the full BETSE solver *and* one or more arbitrary
    requirements to be satisfied.

    Parameters
    ----------
    iterable : IterableTypes
        Iterable of one or more arbitrary requirements defining this set.
    '''

    # Avoid circular import dependencies.
    from betse.science.phase.require.phasereqs import SOLVER_FULL

    # List requiring both the full BETSE solver *AND* all passed requirements.
    requirements_full = [SOLVER_FULL]
    requirements_full.extend(requirements)

    # Create and return this requirement set from this list.
    return SimPhaseRequirements(requirements_full)


@type_check
def make_requirements_solver_full_and(
    requirement: SimPhaseRequirementABC) -> SimPhaseRequirements:
    '''
    Immutable set of full solver-based simulation phase requirements, requiring
    both the full BETSE solver *and* a single arbitrary requirement to be
    satisfied.

    Parameters
    ----------
    requirement : SimPhaseRequirementABC
        Arbitrary requirement defining this set.
    '''

    return make_requirements_solver_full(requirements=(requirement,))


@type_check
def make_requirements_solver_full_and_bool_expr(
    *args, **kwargs) -> SimPhaseRequirements:
    '''
    Immutable set of full solver-based simulation phase requirements, requiring
    both the full BETSE solver *and* a single string-based boolean
    requirement to be satisfied.

    All passed parameters are passed to the
    :meth:`SimPhaseRequirementBoolExpr.__init__` method as is.
    '''

    # String-based boolean simulation phase requirement.
    requirement = SimPhaseRequirementBoolExpr(*args, **kwargs)

    return make_requirements_solver_full(requirements=(requirement,))


@type_check
def make_requirements_ion(name: str, ion_name: str) -> SimPhaseRequirements:
    '''
    Immutable set of ionic simulation phase requirements, requiring both the
    full BETSE solver *and* a single ion type to be satisfied.

    Since simulation of ion concentrations requires the full BETSE solver, this
    requirement is actually the conjunction of that requirement *and* a custom
    requirement specific to this ion.

    Parameters
    ----------
    name : str
        Human-readable name of this requirement (e.g., "extracellular spaces").
        For usability, the first character of this name should typically be
        lower- rather than uppercase.
    ion_name : str
        Name of the required ion, equivalent to a key of the
        :attr:`Parameters.ions_dict` dictionary (e.g., ``Ca``, requiring
        calcium ions).
    '''

    # Requirement that the passed ion is satisfied by a given phase.
    ion_requirement = SimPhaseRequirementEmbodied(
        name=name,
        is_satisfied_body=(
            'return phase.p.ions_dict[{!r}] == 1'.format(ion_name)),
        set_satisfied_body=(
            'phase.p.ions_dict[{!r}] = 1'.format(ion_name)),
    )

    # Create and return this requirement set from this requirement.
    return make_requirements_solver_full_and(requirement=ion_requirement)
