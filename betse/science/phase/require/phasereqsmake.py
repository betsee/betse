#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Utility functions simplifying creation of **simulation phase requirement sets**
(i.e., immutable sets of simulation phase requirements, requiring zero or more
arbitrary requirements to be satisfied).
'''

# ....................{ IMPORTS                            }....................
from betse.science.phase.require.abc.phasereqabc import (
    SimPhaseRequirementABC,
    SimPhaseRequirement,
    SimPhaseRequirementEmbodied,
    SimPhaseRequirementBoolExpr,
    SimPhaseRequirementEnumExpr,
)
from betse.science.phase.require.abc.phasereqset import SimPhaseRequirements
# from betse.util.io.log import logs
from betse.util.type.types import type_check

# ....................{ MAKERS                             }....................
@type_check
def make_requirements(
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
def make_requirements_funcs(*args, **kwargs) -> SimPhaseRequirements:
    '''
    Immutable set of only one callable-based simulation phase requirement.

    All passed parameters are passed to the
    :meth:`SimPhaseRequirement.__init__` method as is.
    '''

    # Callable-based simulation phase requirement.
    requirement = SimPhaseRequirement(*args, **kwargs)

    # Create and return this requirement set from this requirement.
    return make_requirements(requirement)


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
    return make_requirements(requirement)


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
    return make_requirements(requirement)


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

    # Avoid circular import dependencies.
    from betse.science.phase.require.phasereqs import SOLVER_FULL

    # Requirement that the passed ion is satisfied by a given phase.
    ion_requirement = SimPhaseRequirementEmbodied(
        name=name,
        is_satisfied_body=(
            'return phase.p.ions_dict[{!r}] == 1'.format(ion_name)),
        set_satisfied_body=(
            'phase.p.ions_dict[{!r}] = 1'.format(ion_name)),
    )

    # Create and return this requirement set from this requirement.
    return SOLVER_FULL | make_requirements(requirement=ion_requirement)
