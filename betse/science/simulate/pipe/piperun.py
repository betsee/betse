#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level **simulation pipeline runner** (i.e., simulation activity iteratively
run by its parent pipeline) functionality.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractproperty
from betse.exceptions import BetseSimPipelineException
from betse.science.simulate.pipe.piperunreq import SimPipelineRunnerRequirement
from betse.util.type import strs
from betse.util.type.cls.decorators import MethodDecorator
from betse.util.type.obj import objects
from betse.util.type.types import (
    type_check, CallableTypes, SequenceTypes, SetOrNoneTypes,)

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
        Set of zero or more :class:`SimPipelineRunnerRequirement` instances
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
        Defaults to ``None``, in which case no such decoration is applied.
    '''

    @type_check
    def _runner_metadata_closure(
        method: CallableTypes) -> SimPipelineRunner:
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
            #. The :class:`SimPipelineRunner` class decorator,
               annotating this method with this metadata.

        See Also
        ----------
        :func:`runner_metadata`
            Further details.
        '''

        # Return an instance of the class decorator exposing this metadata.
        return SimPipelineRunner(
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
class SimPipelineRunner(MethodDecorator):
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
            Set of zero or more :class:`SimPipelineRunnerRequirement` instances.

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

        # For each requirement in this set (defaulting to the empty tuple)...
        for requirement in requirements or ():
            # If this is *NOT* a requirement, raise an exception.
            objects.die_unless_instance(
                obj=requirement, cls=SimPipelineRunnerRequirement)

            # Add this requirement to this set.
            self.requirements.add(requirement)

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

# ....................{ INTERFACES                         }....................
class SimPipelineRunnerConf(object, metaclass=ABCMeta):
    '''
    Abstract base class of all subclasses defining a type of **simulation
    pipeline runner arguments** (i.e., simple object encapsulating all input
    parameters to be passed to a method implementing a runner in a
    :class:`SimPipelinerABC` pipeline).

    This class is suitable for use as a multiple-inheritance mixin. To preserve
    the expected method resolution order (MRO) semantics, this class should
    typically be the *last* rather than *first* base class inherited from.

    See Also
    ----------
    :class:`betse.science.config.confabc.SimConfListableABC`
        Class subclassing this base class via multiple inheritance.
    '''

    # ..................{ SUBCLASS                           }..................
    @abstractproperty
    def is_enabled(self) -> bool:
        '''
        ``True`` only if this runner is **enabled** (i.e., present in the
        parent simulation pipeline *and* containing an ``enabled`` boolean set
        to ``True``).
        '''


    @abstractproperty
    def name(self) -> str:
        '''
        Lowercase alphanumeric string uniquely identifying the runner these
        arguments apply to in the parent simulation pipeline (e.g.,
        ``voltage_intra``, signifying an intracellular voltage runner).
        '''
