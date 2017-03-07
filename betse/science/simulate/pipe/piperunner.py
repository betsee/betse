#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level **simulation pipeline runner** (i.e., simulation activity iteratively
run by its parent pipeline) functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseSimPipelineException
from betse.science.simulate.simphase import SimPhaseABC
from betse.util.type import enums, strs
from betse.util.type.cls.decorators import MethodDecorator
from betse.util.type.cls.expralias import ExprAliasUnbound
from betse.util.type.types import (
    type_check, CallableTypes, SequenceTypes, SetOrNoneTypes,)
from enum import Enum

# ....................{ DECORATORS                         }....................
@type_check
def runner_metadata(
    # Mandatory metadata.
    categories: SequenceTypes,

    # Optional metadata.
    requirements: SetOrNoneTypes = None,
) -> CallableTypes:
    '''
    Decorator annotating simulation pipeline **runners** (i.e., methods of
    :class:`SimPipelinerABC` subclasses with names prefixed by
    :attr:`SimPipelinerABC._RUNNER_METHOD_NAME_PREFIX`) with custom metadata.

    All such runners decorated by this decorator are guaranteed to be instances
    of the :class:`SimPipelineRunner` class, which provides all metadata passed
    to this decorator as instance variables of the same name.

    Caveats
    ----------
    **This decorator is strictly optional.** Runners *not* decorated by this
    decorator are still runnable from simulation pipelines. Since this decorator
    annotates runners with metadata, however, unannotated runners will *not* be
    usable by external interfaces expecting this metadata -- typically, GUIs
    populating interactive widget fields by this metadata.

    **Runner methods decorated by this decorator should not be decorated by
    other decorators.** In particular, decorated methods should *not* also be
    decorated by :func:`@type_check`, which this decorator already internally
    decorates all decorated methods by.

    Parameters
    ----------
    categories : SequenceTypes
        Sequence of one or more human-readable strings iteratively naming all
        arbitrary categories to which this runner belongs (in descending order
        of hierarchical taxonomy). Categories are arbitrary labels accessed
        *only* by external interfaces and are otherwise ignored by the core
        codebase. Specifically:
        * The first string in this sequence names an arbitrary **root category**
          (e.g., root node in a tree view), intended to be shared between
          multiple runners. This string is typically a broadly applicable label
          such as ``Voltage Plots``.
        * The last string in this sequence names an arbitrary **leaf category**
          (e.g., leaf node in a tree view), intended to be unique to a single
          runner. This string is typically a narrowly applicable label such as
          ``Extracellular Voltage Plot``.
        * All other strings in this sequence name arbitrary categories of
          increasingly fine-grained depth, again intended to be shared between
          multiple runners.
    requirements : optional[SetType]
        Set of zero or more :class:`SimPipelineRunnerRequirementType` instances
        specifying all simulation features required by this runner. This
        decorator then decorates this runner by performing the following logic
        immediately *before* calling this runner:
        * For each requirement in this set...
          * If this requirement is unsatisfied by the current simulation phase
            (e.g., as the configuration for this phase disables extracellular
            spaces), this requirement (and hence this runner) is unsatisfied.
            Since this constitutes a fatal error, an
            :class:`BetseSimPipelineUnsatisfiedException` is raised.
          * Else, this runner is run.
    '''

    @type_check
    def _runner_metadata_closure(
        method: CallableTypes) -> SimPipelineRunnerMetadata:
        '''
        Closure both type-checking *and* annotating the passed simulation
        pipeline runner method with the metadata passed to the outer decorator
        defining this closure, returning an instance of the class decorator
        exposing this metadata to external interfaces.

        Parameters
        ----------
        method : CallableTypes
            Unbound method (i.e., function) to be decorated by (in order):
            #. The :func:`@type_check` decorator, type checking this method.
               For efficiency, callers should ensure this method is *not*
               externally decorated by this decorator.
            #. The :class:`SimPipelineRunnerMetadata` class decorator,
               annotating this method with this metadata.

        See Also
        ----------
        :func:`runner_metadata`
            Further details.
        '''

        # Return an instance of the class decorator exposing this metadata.
        return SimPipelineRunnerMetadata(
            # As a caller convenience, ensure this method is type-checked.
            method=type_check(method),
            categories=categories,
            requirements=requirements,
        )

    # Return the closure accepting the method to be decorated.
    return _runner_metadata_closure

# ....................{ DECORATORS ~ alias                 }....................
# For each abstract "SimPipelineABC" subclass (e.g., "SimPipelinerExportABC"),
# alias the @runner_metadata decorator to a name specific to that subclass.

exporter_metadata = runner_metadata
'''
Decorator annotating simulation export pipeline **runners** (i.e., methods of
:class:`SimPipelinerExportABC` subclasses with names prefixed by
:attr:`SimPipelinerExportABC._RUNNER_METHOD_NAME_PREFIX`) with custom metadata.

See Also
----------
:func:`runner_metadata`
    Further details.
'''

# ....................{ CLASSES                            }....................
#FIXME: Rename to simply "SimPipelineRunner".
class SimPipelineRunnerMetadata(MethodDecorator):
    '''
    Class decorator annotating simulation pipeline runners with custom metadata.

    All such runners decorated by the :func:`runner_metadata` decorator are
    guaranteed to be instances of this class, which provides all metadata passed
    to this decorator as instance variables of the same name.

    Attributes
    ----------
    categories : SequenceTypes
        Sequence of one or more human-readable strings iteratively naming all
        arbitrary categories to which this runner belongs (in descending order
        of hierarchical taxonomy).
    method_name : str
        Name of the method implementing this runner.
    requirements : SetType
        Set of zero or more :class:`SimPipelineRunnerRequirement` instances
        specifying all simulation features required by this runner.
    description : str
        Human-readable description of this runner as a **single-line string**
        (i.e., containing no newlines).

    See Also
    ----------
    :func:`runner_metadata`
        Further details.
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        method: CallableTypes,
        categories: SequenceTypes,
        requirements: SetOrNoneTypes,
    ) -> None:
        '''
        Initialize this class decorator.

        Parameters
        ----------
        method: CallableTypes
            Unbound method (i.e., function) to be decorated.
        categories : SequenceTypes
            Sequence of one or more human-readable category names.
        requirements: SetOrNoneTypes
            Set of zero or more :class:`SimPipelineRunnerRequirementType`
            instances.

        Raises
        ----------
        BetseSimPipelineException
            If this method has no docstring.
        '''

        # Initialize our superclass with the passed method.
        super().__init__(method)

        # Set of all requirements, converted from the passed set of enumeration
        # members to the values encapsulated by these members.
        self.requirements = set()

        # For each enumeration member in this set (defaulting to the empty
        # tuple for efficiency)...
        for requirement in requirements or ():
            # If this is *NOT* an enumeration member, raise an exception.
            enums.die_unless_enum_member(
                enum_type=SimPipelineRunnerRequirementType,
                enum_member=requirement)

            # Add the requirement encapsulated by this enumeration member to
            # this set.
            self.requirements.add(requirement.value)

        # Classify all remaining passed parameters.
        self.categories = categories
        self.method_name = method.__name__

        # Default this runner's description to its docstring.
        self.description = method.__doc__

        # If this docstring is empty, raise an exception.
        if not self.description:
            raise BetseSimPipelineException(
                'Runner method {}() has no docstring.'.format(method.__name__))
        # Else, this docstring is non-empty.

        # Transform this docstring into a description by...
        self.description = (
            # Removing all leading and trailing whitespace.
            strs.remove_whitespace_presuffix(
            # Reducing from a (possibly) multi- to single-line string.
            strs.unwrap(self.description)))

    # ..................{ CALLERS                            }..................
    @type_check
    def __call__(
        self,

        # To avoid circular import dependencies, this is type-checked as a
        # fully-qualified class name resolved at runtime.
        pipeline: 'betse.science.simulate.pipe.pipeabc.SimPipelinerABC',
        *args,
        **kwargs
    ) -> object:

        # If this runner is unsatisfied, raise an exception.
        pipeline.die_unless_runner_satisfied(self)
        # Else, this runner is satisfied. The prior call logged the attempt to
        # run this runner; now actually do so.

        # Defer to the superclass implementation to run this runner.
        return super().__call__(pipeline, *args, **kwargs)


class SimPipelineRunnerRequirement(object):
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
    def __init__(self, expr: str, name: str) -> None:
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
        '''

        # Expression alias encapsulating the passed Python expression. Since
        # this alias is only ever accessed with the public methods defined below
        # that are already type-checked, this alias avoids additional validation
        # (e.g., a "cls=bool" parameter).
        self._expr_alias = ExprAliasUnbound(expr=expr, obj_name='phase')

        # Classify all remaining parameters.
        self.name = name

    # ..................{ TESTERS                            }..................
    #FIXME: Annotate this method as strictly returning values of type "bool".
    #Sadly, due to the following design deficiencies elsewhere in the codebase,
    #this method occasionally returns non-boolean types:
    #
    #* The ion requirements (e.g., "ION_NA") return integers retrieved from the
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

# ....................{ ENUMERATIONS                       }....................
#FIXME: Sadly, this has turned out to be extreme overkill. An enumeration isn't
#actually required anywhere; what *IS* required are the constants currently
#embedded within the namespace of this enumeration. Hence, simplify this as
#follows:
#
#* Rename this submodule to "piperun".
#* Define a new "piperunreq" submodule.
#* For each constant declared by this enumeration, shift that constant into that
#  submodule with the same name (e.g., shift
#  "SimPipelineRunnerRequirementType.ECM" to "piperunreq.ECM".
#* Refactor the SimPipelineRunnerMetadata.__init__() method to accept a
#  "requirements" parameter that is a set of "SimPipelineRunnerRequirement"
#  instances rather than a set of "SimPipelineRunnerRequirementType" instances.
#* Remove this enumeration.

class SimPipelineRunnerRequirementType(Enum):
    '''
    Enumeration of all recognized requirements for simulation pipeline runners.

    Each member of this enumeration is an instance of the
    :class:`SimPipelineRunnerRequirement` class, ensuring that each requirement
    is both efficiently testable as an enumeration member *and* associated with
    additional metadata in a type-safe manner.

    Attributes
    ----------
    DEFORM : SimPipelineRunnerRequirement
        Requirement that a simulation phase enable cellular deformations.
    ECM : SimPipelineRunnerRequirement
        Requirement that a simulation phase enable the extracellular matrix
        (ECM), also referred to as "extracellular spaces."
    ELECTROOSMOSIS : SimPipelineRunnerRequirement
        Requirement that a simulation phase enable electroosmotic flow (EOF).
    FLUID : SimPipelineRunnerRequirement
        Requirement that a simulation phase enable fluid flow.
    ION_CA : SimPipelineRunnerRequirement
        Requirement that a simulation phase enable calcium (Ca2+) ions.
    ION_H : SimPipelineRunnerRequirement
        Requirement that a simulation phase enable hydrogen (H+) ions.
    ION_K : SimPipelineRunnerRequirement
        Requirement that a simulation phase enable potassium (K+) ions.
    ION_NA : SimPipelineRunnerRequirement
        Requirement that a simulation phase enable sodium (Na+) ions.
    PRESSURE_MECHANICAL : SimPipelineRunnerRequirement
        Requirement that a simulation phase enable the mechanical pressure
        intervention.
    PRESSURE_OSMOTIC : SimPipelineRunnerRequirement
        Requirement that a simulation phase enable the mechanical pressure
        intervention.
    '''

    DEFORM = SimPipelineRunnerRequirement(
        expr='phase.p.deformation', name='cellular deformation',)
    ECM = SimPipelineRunnerRequirement(
        expr='phase.p.sim_ECM', name='extracellular spaces',)
    ELECTROOSMOSIS = SimPipelineRunnerRequirement(
        expr='phase.p.sim_eosmosis', name='electroosmotic flow',)
    FLUID = SimPipelineRunnerRequirement(
        expr='phase.p.fluid_flow', name='fluid flow',)
    ION_CA = SimPipelineRunnerRequirement(
        expr='phase.p.ions_dict["Ca"]', name='calcium (Ca2+) ions',)
    ION_H = SimPipelineRunnerRequirement(
        expr='phase.p.ions_dict["H"]', name='hydrogen (H+) ions',)
    ION_K = SimPipelineRunnerRequirement(
        expr='phase.p.ions_dict["K"]', name='potassium (K+) ions',)
    ION_NA = SimPipelineRunnerRequirement(
        expr='phase.p.ions_dict["Na"]', name='sodium (Na+) ions',)
    PRESSURE_MECHANICAL = SimPipelineRunnerRequirement(
        expr='phase.p.is_event_pressure', name='mechanical pressure',)
    PRESSURE_OSMOTIC = SimPipelineRunnerRequirement(
        expr='phase.p.deform_osmo', name='osmotic pressure',)

# ....................{ ENUMERATIONS ~ alias               }....................
# For each abstract "SimPipelineABC" subclass (e.g., "SimPipelinerExportABC"),
# alias the @runner_metadata decorator to a name specific to that subclass.

exporter_requirement = SimPipelineRunnerRequirementType
'''
Enumeration of all recognized requirements for simulation export pipeline
runners.

See Also
----------
:class:`SimPipelineRunnerRequirementType`
    Further details.
'''
