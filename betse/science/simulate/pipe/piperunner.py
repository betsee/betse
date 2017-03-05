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
from betse.util.type import strs
from betse.util.type.cls.decorators import MethodDecorator
from betse.util.type.types import (
    type_check, CallableTypes, SequenceTypes, SetOrNoneTypes,)

# ....................{ DECORATORS                         }....................
#FIXME: Revise the "requirements" parameter docstring as follows:
#
#* Define a new "SimPipelineRunnerRequirement" class encapsulating as simple
#  public instance variables all metadata pertaining to a type of simulation
#  feature to be required by pipeline runner, including:
#  * "SimPipelineRunnerRequirement.name", the human-readable lowercase name of
#    that simulation feature (e.g., "extracellular spaces").
#  Likewise, "SimPipelineRunnerRequirement" methods might include:
#  * def is_satisfied(phase: SimPhaseABC) -> bool, returning "True" only if the
#    passed phase enabled this feature. Note that, since a passed parameter is
#    required, this cannot be reduced to a simple property.
#  * def set_satisfied(phase: SimPhaseABC, is_satisfied: bool) -> None, setting
#    the current state of this feature for the passed phase.
#* Define a new "SimPipelineRunnerRequirementType" enumeration type in this
#  submodule. Each member of this enumeration should be an instance of the
#  aforementioned "SimPipelineRunnerRequirement" class, thus associating each
#  requirement with its metadata. Each enumeration member then signifies a type
#  of simulation feature to be required by pipeline runners: e.g.,
#  * "SimPipelineRunnerRequirementType.ECM", signifying extracellular spaces.
#  This type might be defined like so:
#
#    class SimPipelineRunnerRequirementType(Enum):
#        ECM = SimPipelineRunnerRequirement(
#            name='extracellular spaces',
#            is_satisfied=lambda phase:
#               phase.p.sim_ECM,
#            set_satisfied=lambda phase, value:
#               phase.p.__setattr__('sim_ECM', value),
#        )
#
#* Rather than a set of attribute names (e.g., "{'sim_ECM',}"), this parameter
#  should instead accept a set of "SimPipelineRunnerRequirementType" instances.
#
#Awesome, eh? Extremely efficient, typesafe, and sensible. Let's do this.

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
        specifying all simulation features required by this runner.
        #FIXME: Obsoleted, sadly:
        Set of one or more names of arbitrarily nested instance boolean
        variables of the class:`Parameters` instance associated with this
        pipeline, each signifying a simulation feature required by this runner
        (e.g., ``{'sim_ECM',}``, expanding to the ``p.sim_ECM`` boolean
        variable signifying a runner requiring only extracellular spaces). Each
        name may contain one or more ``.`` delimiters, signifying a variable
        contained within one or more parent objects contained within this
        class:`Parameters` instance (e.g., ``{'anim.is_while_sim_save',}``).
        This decorator then decorates this runner by performing the following
        logic immediately *before* calling this runner:
        * For each variable name in this set...
          * If the value of this variable is ``False``, this runner is
            unsatsified by this requirement. Since this constitutes a fatal
            error, an :class:`BetseSimPipelineUnsatisfiedException` is raised.
          * Else, this runner is called.
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
    requirements: SetType
        Set of zero or more :class:`SimPipelineRunnerRequirementType` instances
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

        # Default all unpassed parameters to sane defaults.
        if requirements is None:
            requirements = set()

        # Classify all remaining passed parameters.
        self.categories = categories
        self.requirements = requirements

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
