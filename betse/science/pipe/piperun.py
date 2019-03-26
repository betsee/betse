#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **simulation pipeline runner** (i.e., simulation activity
iteratively run by its parent pipeline) functionality.
'''

#FIXME: Refactor the "SimPipeRunnerMetadata.kind" instance variable from a
#non-type-safe string to a type-safe enumeration. Since the set of all types of
#runner methods supported by a given pipeline is dynamically defined by the
#current implementation of the subclass for that pipeline, we'll need to
#dynamically define one eneumeration type for each pipeline subclass at runtime
#in a presumably clever manner.

# ....................{ IMPORTS                           }....................
import functools
from betse.exceptions import (
    BetseSimPipeException, BetseSimPipeRunnerUnsatisfiedException)
from betse.science.phase.phasecls import SimPhase
from betse.science.phase.require.abc.phasereqset import (
    SimPhaseRequirements, SimPhaseRequirementsOrNoneTypes)
from betse.util.io.log import logs
from betse.util.type.text.string import strs
from betse.util.type.types import type_check, CallableTypes, SequenceTypes

# ....................{ CLASSES                           }....................
class SimPipeRunnerMetadata(object):
    '''
    **Simulation pipeline runner** (i.e., method of a :class:`SimPipeABC`
    subclass decorated by the :func:`piperunner` decorator) metadata.

    This metadata is available via the ``metadata`` instance variable of each
    such runner method.

    Attributes (Pipeline)
    ----------
    noun_singular_lowercase : str
        Human-readable lowercase singular noun synopsizing this type of runner
        (e.g., ``animation``, ``plot``), equivalent to the
        :attr:`SimPipeABC.subclass._noun_singular_lowercase` constant defined
        by the parent pipeline of this runner.
    noun_singular_uppercase : str
        Human-readable uppercase singular noun synopsizing this type of runner
        (e.g., ``animation``, ``plot``), equivalent to the
        :attr:`SimPipeABC.subclass._noun_singular_uppercase` constant defined
        by the parent pipeline of this runner.
    verb_continuous : str
        Human-readable verb in the continuous tense synopsizing the type of
        action performed by this type of runner (e.g., ``Saving``), equivalent
        to the :attr:`SimPipeABC.subclass._VERB_CONTINUOUS` constant defined by
        the parent pipeline of this runner.

    Attributes (Runner)
    ----------
    categories : SequenceTypes
        Sequence of one or more human-readable strings iteratively naming all
        arbitrary categories to which this runner belongs (in descending order
        of hierarchical taxonomy).
    kind : str
        Machine-readable type of this runner as a raw string. Although this
        string is currently equivalent to :attr:`method_name` excluding the
        :attr:`SimPipeABC._RUNNER_METHOD_NAME_PREFIX` substring required by the
        parent pipeline of this runner, callers should *not* assume this
        low-level implementation detail to always be.
    method_name : str
        Machine-readable name of the method implementing this runner.
    requirements : SimPhaseRequirements
        Immutable set of zero or more :class:`SimPhaseRequirement` instances
        specifying all simulation features required by this runner.
    description : str
        Human-readable description of this runner as a **single-line string**
        (i.e., containing no newlines).

    See Also
    ----------
    :func:`piperunner`
        Further details.
    '''

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def __init__(
        self,

        # Mandatory parameters.
        method: CallableTypes,
        categories: SequenceTypes,

        # Optional parameters.
        requirements: SimPhaseRequirementsOrNoneTypes = None,
    ) -> None:
        '''
        Initialize this metadata.

        Parameters
        ----------
        method: CallableTypes
            Unbound method (i.e., function) associated with this metadata.
        categories : SequenceTypes
            Sequence of one or more human-readable category names.
        requirements : SimPhaseRequirementsOrNoneTypes
            Immutable set of zero or more :class:`SimPhaseRequirement`
            instances *or* ``None``, in which case this parameter defaults to
            the empty set of such instances.

        Raises
        ----------
        BetseSimPipeException
            If this method has no docstring.
        '''

        # Initialize our superclass with the passed method.
        super().__init__()

        # Default all unpassed parameters to sane defaults.
        if requirements is None:
            requirements = SimPhaseRequirements()

        # Classify all passed parameters.
        self.requirements = requirements
        self.categories = categories
        self.method_name = method.__name__

        # Nullify all remaining instance variables, subsequently defined by the
        # "SimPipeABCMeta" metaclass at "SimPipeABC" subclass definition time
        # following the definition of this runner method.
        self.kind = None
        self.noun_singular_lowercase = None
        self.noun_singular_uppercase = None
        self.verb_continuous = None

        # Default this runner's description to its docstring.
        self.description = method.__doc__

        # If this docstring is empty, raise an exception.
        if not self.description:
            raise BetseSimPipeException(
                'Pipeline runner method {}() docstring undefined.'.format(
                    self.method_name))
        # Else, this docstring is non-empty.

        # Transform this docstring into a description by...
        self.description = (
            # Removing all leading and trailing whitespace.
            strs.remove_whitespace_presuffix(
            # Reducing from a (possibly) multi- to single-line string.
            strs.unwrap(self.description)))

        # Expose the metadata associated with this runner to callers. Due to
        # Python constraints on class decorators, *ONLY* the bound
        # SimPipeRunner.__call__() method returned by this method below is
        # accessible to callers. Notably, this "SimPipeRunner" instance is not
        # accessible to callers. Ergo, attaching this method to this method is
        # the only means of exposing this metadata to callers. *shrug*
        method.metadata = self

# ....................{ DECORATORS                        }....................
@type_check
def piperunner(
    # Mandatory metadata.
    categories: SequenceTypes,

    # Optional metadata.
    requirements: SimPhaseRequirementsOrNoneTypes = None,
) -> CallableTypes:
    '''
    Decorator annotating **simulation pipeline runners** (i.e.,
    :meth:`SimPipeRunner.__call__` subclasses with names prefixed by
    :attr:`SimPipeABC._RUNNER_METHOD_NAME_PREFIX`) with custom metadata.

    All methods decorated by this decorator are guaranteed to be instances of
    the :class:`SimPipeRunner` class, which provides all metadata passed to
    this decorator as instance variables of the same name.

    Caveats
    ----------
    **This decorator is strictly optional.** Runners *not* decorated by this
    decorator are still runnable from simulation pipelines. Since this
    decorator annotates runners with metadata, however, unannotated runners
    will *not* be usable by external interfaces expecting this metadata --
    typically, GUIs populating interactive widget fields by this metadata.

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

        * The first string in this sequence names an arbitrary **root
          category** (e.g., root node in a tree view), intended to be shared
          between multiple runners. This string is typically a broadly
          applicable label such as ``Voltage Plots``.
        * The last string in this sequence names an arbitrary **leaf category**
          (e.g., leaf node in a tree view), intended to be unique to a single
          runner. This string is typically a narrowly applicable label such as
          ``Extracellular Voltage Plot``.
        * All other strings in this sequence name arbitrary categories of
          increasingly fine-grained depth, again intended to be shared between
          multiple runners.
    requirements : SimPhaseRequirementsOrNoneTypes
        Immutable set of zero or more :class:`SimPhaseRequirement` instances
        specifying all simulation features required by this runner. This
        decorator then decorates this runner by performing the following logic
        immediately *before* calling this runner:

        * For each requirement in this set...

          * If this requirement is unsatisfied by the current simulation phase
            (e.g., as the configuration for this phase disables extracellular
            spaces), this requirement (and hence this runner) is unsatisfied.
            Since this constitutes a fatal error, an
            :class:`BetseSimPipeRunnerUnsatisfiedException` is raised.
          * Else, this runner is run.

        Defaults to ``None``, in which case no such decoration is applied.
    '''

    @type_check
    def _piperunner_method_factory(method: CallableTypes) -> CallableTypes:
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
            #. The :class:`SimPipeRunner` class decorator,
               annotating this method with this metadata.

        See Also
        ----------
        :func:`piperunner`
            Further details.
        '''

        # As a caller convenience, ensure this method is type-checked.
        method_typed = type_check(method)

        @functools.wraps(method_typed)
        @type_check
        def _piperunner_method(
            # To avoid circular import dependencies, this is validated as a
            # fully-qualified class name resolved at runtime.
            #
            # For clarity, this parameter has been renamed from the customary
            # "self" nomenclature for a bound method.
            self_pipeline: 'betse.science.pipe.pipeabc.SimPipeABC',
            phase: SimPhase,
            *args, **kwargs
        ) -> object:
            '''
            Closure validating this simulation pipeline runner to be satisfied
            by the passed simulation pipeline and phase *before* running that
            runner by calling the method passed to the outer decorator defining
            this closure with the metadata passed to the same decorator,
            returning the values returned by that method call.

            Specifically, this closure either:

            * If this runner is **unsatisfied** (i.e., requires one or more
              simulation features disabled by the user-defined simulation
              configuration for this phase), raises an exception.
            * Else, logs an attempt to run this runner.

            Parameters
            ----------
            self_pipeline : SimPipeABC
                Current simulation pipeline.
            phase : SimPhase
                Current simulation phase.

            All remaining parameters are passed as is to the method passed to
            the outer decorator defining this closure.

            Raises
            ----------
            BetseSimVisualUnsatisfiedException
                If this runner is unsatisfied.

            See Also
            ----------
            :func:`piperunner`
                Further details.
            '''

            # Metadata associated with this runner.
            runner_metadata = _piperunner_method.metadata

            # If any runner requirement is unsatisfied, raise an exception.
            for requirement in runner_metadata.requirements:
                if not requirement.is_satisfied(phase):
                    raise BetseSimPipeRunnerUnsatisfiedException(
                        result='{} "{}" requirement unsatisfied'.format(
                            runner_metadata.noun_singular_uppercase,
                            runner_metadata.kind),
                        reason='{} disabled'.format(requirement.name),
                    )
            # Else, all runner requirements are satisfied.

            # Log the subsequent attempt to run this runner.
            logs.log_info(
                '%s %s "%s"...',
                runner_metadata.verb_continuous,
                runner_metadata.noun_singular_lowercase,
                runner_metadata.kind)

            # Else, this runner is satisfied. Since the prior call already
            # logged the attempt to run this runner, avoid redoing so here.
            #
            # Simply call this method to run this runner.
            return method_typed(self_pipeline, phase, *args, **kwargs)

        # Expose this metadata as an instance variable of this method closure.
        _piperunner_method.metadata = SimPipeRunnerMetadata(
            method=method_typed,
            categories=categories,
            requirements=requirements,
        )

        # Return this method closure.
        return _piperunner_method

    # Return this method factory closure.
    return _piperunner_method_factory
