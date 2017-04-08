#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level **simulation pipeline runner requirement** (i.e., prerequisite
simulation feature required by a runner) functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.science.simulate.simphase import SimPhaseABC
from betse.util.type.cls.expralias import ExprAliasUnbound
from betse.util.type.types import type_check

# ....................{ CLASSES                            }....................
class SimPipeRunnerRequirement(object):
    '''
    Object encapsulating the current boolean state of a simulation feature
    (e.g., extracellular spaces) required by simulation pipeline runners.

    Attributes
    ----------
    name : str
        Human-readable lowercase name of this requirement (e.g.,
        ``extracellular spaces``).
    _expr_alias : ExprAliasUnbound
        Expression alias, dynamically referring to an arbitrarily complex source
        Python expression relative to the current simulation phase evaluating to
        this requirement's boolean state (e.g., ``phase.p.sim_ECM``). Since each
        instance of this class is shared amongst all phases rather than bound to
        a single phase, this alias is unbound rather than bound.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        expr: str,
        name: str,
        **kwargs
    ) -> None:
        '''
        Initialize this requirement.

        Parameters
        ----------
        expr : str
            Arbitrarily complex Python expression relative to the current
            simulation phase evaluating to this requirement's boolean state
            (e.g., ``phase.p.sim_ECM``). This expression may (and typically
            should) refer to the ``phase`` variable, bound to the current
            simulation phase.
        name : str
            Human-readable lowercase name of this requirement (e.g.,
            ``extracellular spaces``).

        All remaining keyword parameters are passed as is to the
        :meth:`ExprAliasUnbound.__init__` method.
        '''

        # Expression alias encapsulating the passed Python expression. Since
        # this alias is only ever accessed with the public methods defined below
        # that are already type-checked, this alias avoids additional validation
        # (e.g., a "cls=bool" parameter).
        self._expr_alias = ExprAliasUnbound(
            expr=expr, obj_name='phase', **kwargs)

        # Classify all remaining parameters.
        self.name = name

    # ..................{ TESTERS                            }..................
    #FIXME: Annotate this method as strictly returning values of type "bool".
    #Sadly, due to the following design deficiencies elsewhere in the codebase,
    #this method occasionally returns non-boolean types:
    #
    #* The ion requirements (e.g., "ION_SODIUM") return integers retrieved from the
    #  "phase.p.ions_dict" dictionary (e.g., 'phase.p.ions_dict["Na"]'), which
    #  are guaranteed to be either 0 or 1 and hence are effectively boolean.

    @type_check
    def is_satisfied(self, phase: SimPhaseABC):
        '''
        ``True`` only if the passed simulation phase satisfies this requirement.

        Specifically, this method returns the current boolean state of the
        simulation feature encapsulated by this requirement.

        Parameters
        ----------
        phase : SimPhaseABC
            Current simulation phase.

        Returns
        ----------
        bool
            ``True`` only if this phase satisfies this requirement.
        '''

        return self._expr_alias.get(phase)

    # ..................{ SETTERS                            }..................
    @type_check
    def set_satisfied(self, phase: SimPhaseABC, is_satisfied: bool) -> bool:
        '''
        Force the passed simulation phase to either satisfy or *not* satisfy
        this requirement.

        Specifically, this method sets the current boolean state of the
        simulation feature encapsulated by this requirement to this boolean.

        Parameters
        ----------
        phase : SimPhaseABC
            Current simulation phase.
        is_satisfied : bool
            ``True`` only if this phase is to satisfy this requirement.
        '''

        return self._expr_alias.set(phase, is_satisfied)

# ....................{ CONSTANTS                          }....................
DEFORM = SimPipeRunnerRequirement(
    expr='phase.p.deformation', name='cellular deformation',)
'''
Requirement that a simulation phase enable cellular deformations.
'''


ECM = SimPipeRunnerRequirement(
    expr='phase.p.sim_ECM', name='extracellular spaces',)
'''
Requirement that a simulation phase enable the extracellular matrix (ECM), also
referred to as "extracellular spaces."
'''


ELECTROOSMOSIS = SimPipeRunnerRequirement(
    expr='phase.p.sim_eosmosis', name='electroosmotic flow',)
'''
Requirement that a simulation phase enable electroosmotic flow (EOF).
'''


FLUID = SimPipeRunnerRequirement(
    expr='phase.p.fluid_flow', name='fluid flow',)
'''
Requirement that a simulation phase enable fluid flow.
'''


GHK = SimPipeRunnerRequirement(
    expr='phase.p.GHK_calc', name='Goldman calculation',)
'''
Requirement that a simulation phase enable alternative calculation of
transmembrane voltages (Vmem) given the Goldman-Hodgkin-Katz (GHK) equation.
'''


ION_CALCIUM = SimPipeRunnerRequirement(
    expr='phase.p.ions_dict["Ca"]', name='calcium (Ca2+) ions',)
'''
Requirement that a simulation phase enable calcium (Ca2+) ions.
'''


ION_HYDROGEN = SimPipeRunnerRequirement(
    expr='phase.p.ions_dict["H"]', name='hydrogen (H+) ions',)
'''
Requirement that a simulation phase enable hydrogen (H+) ions.
'''


ION_POTASSIUM = SimPipeRunnerRequirement(
    expr='phase.p.ions_dict["K"]', name='potassium (K+) ions',)
'''
Requirement that a simulation phase enable potassium (K+) ions.
'''


ION_M_ANION = SimPipeRunnerRequirement(
    expr='phase.p.ions_dict["M"]', name='M anions (M-)',)
'''
Requirement that a simulation phase enable M anions (M-).
'''


ION_SODIUM = SimPipeRunnerRequirement(
    expr='phase.p.ions_dict["Na"]', name='sodium (Na+) ions',)
'''
Requirement that a simulation phase enable sodium (Na+) ions.
'''


PRESSURE_OSMOTIC = SimPipeRunnerRequirement(
    expr='phase.p.deform_osmo', name='osmotic pressure',)
'''
Requirement that a simulation phase enable osmotic pressure.
'''


PRESSURE_TOTAL = SimPipeRunnerRequirement(
    # For simplicity, define this requirement to be settable by enabling osmotic
    # pressure. While the mechanical pressure event could also be enabled, doing
    # so is less trivial than the former.
    expr         ='phase.p.deform_osmo or phase.p.is_event_pressure',
    expr_settable='phase.p.deform_osmo',
    name='total pressure',)
'''
Requirement that a simulation phase enable at least one pressure feature
(namely, osmotic pressure and/or the mechanical pressure intervention).
'''
