#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **simulation pipeline runner** (i.e., simulation activity iteratively
run by its parent pipeline) functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.exceptions import BetseSimPipeException
from betse.science.phase.require.abc.phasereqset import (
    SimPhaseRequirements, SimPhaseRequirementsOrNoneTypes)
from betse.util.type.decorator.deccls import MethodDecoratorABC
from betse.util.type.text import strs
from betse.util.type.types import (
    type_check, CallableTypes, SequenceTypes,)

# ....................{ CLASSES                            }....................
class SimPipeRunner(MethodDecoratorABC):
    '''
    Class decorator annotating simulation pipeline runners with custom metadata.

    All such runners decorated by the :func:`piperunner` decorator are
    guaranteed to be instances of this class, which provides all metadata passed
    to that decorator as instance variables of the same name.

    Attributes
    ----------
    categories : SequenceTypes
        Sequence of one or more human-readable strings iteratively naming all
        arbitrary categories to which this runner belongs (in descending order
        of hierarchical taxonomy).
    method_name : str
        Name of the method implementing this runner.
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

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,
        method: CallableTypes,
        categories: SequenceTypes,
        requirements: SimPhaseRequirementsOrNoneTypes,
    ) -> None:
        '''
        Initialize this class decorator.

        Parameters
        ----------
        method: CallableTypes
            Unbound method (i.e., function) to be decorated.
        categories : SequenceTypes
            Sequence of one or more human-readable category names.
        requirements : SimPhaseRequirementsOrNoneTypes
            Immutable set of zero or more :class:`SimPhaseRequirement` instances
            *or* ``None``, in which case this parameter defaults to the empty
            immutable set of such instances.

        Raises
        ----------
        BetseSimPipeException
            If this method has no docstring.
        '''

        # Initialize our superclass with the passed method.
        super().__init__(method)

        # Default all unpassed parameters to sane defaults.
        if requirements is None:
            requirements = SimPhaseRequirements()

        # Classify all passed parameters.
        self.requirements = requirements

        # Classify all remaining passed parameters.
        self.categories = categories
        self.method_name = method.__name__

        # Default this runner's description to its docstring.
        self.description = method.__doc__

        # If this docstring is empty, raise an exception.
        if not self.description:
            raise BetseSimPipeException(
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
        pipeline: 'betse.science.phase.pipe.pipeabc.SimPipeABC',
        *args,
        **kwargs
    ) -> object:

        # If this runner is unsatisfied, raise an exception.
        pipeline.die_unless_runner_satisfied(self)
        # Else, this runner is satisfied. Since the prior call already logged
        # the attempt to run this runner, avoid redoing so here.

        # Defer to the superclass implementation to run this runner.
        return super().__call__(pipeline, *args, **kwargs)

# ....................{ DECORATORS                         }....................
@type_check
def piperunner(
    # Mandatory metadata.
    categories: SequenceTypes,

    # Optional metadata.
    requirements: SimPhaseRequirementsOrNoneTypes = None,
) -> CallableTypes:
    '''
    Decorator annotating simulation pipeline **runners** (i.e., methods of
    :class:`SimPipeABC` subclasses with names prefixed by
    :attr:`SimPipeABC._RUNNER_METHOD_NAME_PREFIX`) with custom metadata.

    All methods decorated by this decorator are guaranteed to be instances of
    the :class:`SimPipeRunner` class, which provides all metadata passed to
    this decorator as instance variables of the same name.

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
    def _piperunner_closure(method: CallableTypes) -> SimPipeRunner:
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

        # Return an instance of the class decorator exposing this metadata.
        return SimPipeRunner(
            # As a caller convenience, ensure this method is type-checked.
            method=type_check(method),
            categories=categories,
            requirements=requirements,
        )

    # Return the closure accepting the method to be decorated.
    return _piperunner_closure
