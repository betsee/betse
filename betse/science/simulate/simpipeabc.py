#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level **simulation pipeline** (i.e., sequence of similar simulation
activities to be iteratively run) functionality.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractproperty
from betse.exceptions import (
    BetseSimConfigException,
    BetseSimPipelineException,
    BetseSimPipelineUnsatisfiedException,
)
from betse.science.simulate.simphase import SimPhaseABC
from betse.util.io.log import logs
from betse.util.type import strs
from betse.util.type.call import callers
from betse.util.type.cls import classes
from betse.util.type.cls.decorators import MethodDecorator
from betse.util.type.obj import objects
from betse.util.type.types import (
    type_check, CallableTypes, GeneratorType, IterableTypes, SequenceTypes)

# ....................{ SUPERCLASSES                       }....................
class SimPipelinerABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all subclasses running a **simulation pipeline**
    (i.e., sequence of similar simulation activities to be iteratively run).

    This class implements the infrastructure for iteratively running all
    currently enabled simulation activities (referred to as "runners")
    comprising the pipeline defined by this subclass, either in parallel *or* in
    series.

    Design
    ----------
    Each subclass is expected to define:

    * One or more public methods with names prefixed by ``run_`` (e.g.,
      ``run_voltage_intra``). For each such method, the name of that method
      excluding the prefix ``run_`` is the name of that method's runner (e.g.,
      ``voltage_intra`` for the method name ``run_voltage_intra``). Each such
      method should:
      * Accept exactly one parameter: ``self``.
      * Return nothing (i.e.,``None``).
    * The abstract :meth:`_runners_conf_enabled` property returning a sequence of the
      names of all runners currently enabled by this pipeline (e.g.,
      ``['voltage_intra', 'ion_hydrogen', 'electric_total']``).

    The :meth:`run` method defined by this base class then dynamically
    implements this pipeline by iterating over the :meth:`_runners_conf_enabled` property
    and, for each enabled runner, calling that runner's method.

    Attributes (Private)
    ----------
    _phase : SimPhaseABC
        Current simulation phase.

    Attributes (Private: Labels)
    ----------
    _label_singular_lowercase : str
        Human-readable lowercase singular noun synopsizing the type of runners
        implemented by this subclass (e.g., ``animation``, ``plot``).
    _label_singular_uppercase : str
        Human-readable singular noun whose first character is uppercase and all
        remaining characters lowercase (e.g., ``Animation``, ``Plot``).
    _label_plural_lowercase : str
        Human-readable lowercase plural noun synopsizing the type of runners
        implemented by this subclass (e.g., ``animations``, ``plots``).
    '''

    # ..................{ CONSTANTS                          }..................
    # Ideally, the following class constants would instead be implemented as
    # class properties. Unfortunately, there exists no @classproperty decorator
    # and no feasible means of implementing that decorator thanks to subtle
    # idiosyncracies in the design of Python's data descriptor protocol. While
    # class properties can technically be manually implemented by a metaclass
    # defining these properties, doing so is cumbersome and quite discouraged.

    _RUNNER_METHOD_NAME_PREFIX = 'run_'
    '''
    Substring prefixing the name of each runner defined by this pipeline.

    This class variable is intended to be overridden by subclasses desiring a
    less ambiguous prefix (e.g., ``export_`` for export pipelines).
    '''

    # ..................{ STATIC ~ iterators                 }..................
    @classmethod
    def iter_runner_methods(cls) -> GeneratorType:
        '''
        Generator yielding each runner method defined by this pipeline subclass.

        For each subclass method whose name is prefixed by
        :attr:`_RUNNER_METHOD_NAME_PREFIX`, this generator yields that method.

        Yields
        ----------
        MethodType
            Each runner method defined by this pipeline subclass.
        '''

        yield from classes.iter_methods_matching(
            cls=cls,
            predicate=lambda method_name: method_name.startswith(
                cls._RUNNER_METHOD_NAME_PREFIX))


    @classmethod
    def iter_runner_names(cls) -> GeneratorType:
        '''
        Generator yielding the name of each runner defined by this pipeline
        subclass.

        For each subclass method whose name is prefixed by
        :attr:`_RUNNER_METHOD_NAME_PREFIX`, this
        generator yields that name excluding the prefixing
        :attr:`_RUNNER_METHOD_NAME_PREFIX`. For example, if this subclass
        defines only two methods ``run_exhra`` and ``run_intra`` whose names are
        prefixed by :attr:`_RUNNER_METHOD_NAME_PREFIX`, this generator first
        yields the string ``extra`` and then the string ``intra`` (in that
        order).

        Yields
        ----------
        str
            Name of each runner defined by this pipeline.
        '''

        return (
            # Yield the name of this method excluding the runner prefix.
            strs.remove_prefix(
                text=runner_method_name, prefix=cls._RUNNER_METHOD_NAME_PREFIX)
            # For the name of each runner method defined by this subclass...
            for runner_method_name, _ in cls.iter_runner_methods()
        )

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,

        # Mandatory parameters.
        phase: SimPhaseABC,
        label_singular: str,

        # Optional parameters.
        label_plural: str = None,
        label_verb: str = 'Running',
    ) -> None:
        '''
        Initialize this pipeline.

        Parameters
        ----------
        phase : SimPhaseABC
            Current simulation phase.
        label_singular : str
            Human-readable singular noun synopsizing the type of runners
            implemented by this subclass (e.g., ``animation``, ``plot``),
            ideally but *not* necessarily lowercase.
        label_plural : optional[str]
            Human-readable plural noun synopsizing the type of runners
            implemented by this subclass (e.g., ``animations``, ``plots``),
            ideally but *not* necessarily lowercase. Defaults to the passed
            ``label_singular`` parameter suffixed by an ``s`` character.
        label_verb: optional[str]
            Human-readable plural noun synopsizing the type of action performed
            by runners implemented by this subclass (e.g., ``Saving``), ideally
            but *not* necessarily capitalized. Defaults to ``Running``.
        '''

        # Default all unpassed parameters.
        if label_plural is None:
            label_plural = label_singular + 's'

        # Classify all passed parameters.
        self._phase = phase
        self._label_singular_lowercase = label_singular
        self._label_plural_lowercase = label_plural
        self._label_verb = label_verb

        # Human-readable capitalized singular noun.
        self._label_singular_uppercase = strs.uppercase_first_char(
            self._label_singular_lowercase)

    # ..................{ RUNNERS                            }..................
    def run(self) -> None:
        '''
        Run all currently enabled pipeline runners if this pipeline itself is
        currently enabled *or* noop otherwise.

        Specifically:

        * If the :meth:`is_enabled` property is ``True`` (implying this pipeline
          to be currently enabled):
          * For each :class:`SimPipelineRunnerConf` instance (corresponding to
            the configuration of a currently enabled pipeline runner) in the
            sequence of these instances provided by the
            :meth:`_runners_conf_enabled` property:
            * Call this method, passed this configuration.
            * If this method reports this runner's requirements to be
              unsatisfied (e.g., due to the current simulation configuration
              disabling extracellular spaces), this runner is ignored with a
              non-fatal warning.
        * Else, log an informative message and return immediately.
        '''

        # If this pipeline is disabled, log this fact and return immediately.
        if not self.is_enabled:
            logs.log_info('Ignoring %s...', self._label_plural_lowercase)
            return
        # Else, this pipeline is enabled.

        # Log this pipeline run.
        logs.log_info(
            '%s %s...', self._label_verb, self._label_plural_lowercase)

        # For the object encapsulating all input arguments to be passed to each
        # currently enabled runner in this pipeline...
        for runner_conf in self._runners_conf_enabled:
            if not isinstance(runner_conf, SimPipelineRunnerConf):
                raise BetseSimPipelineException(
                    '_runners_conf_enabled() item {!r} '
                    'not instance of "SimPipelineRunnerConf".'.format(
                        runner_conf))

            # Name of the method implementing this runner.
            runner_method_name = (
                self._RUNNER_METHOD_NAME_PREFIX + runner_conf.name)

            # Method implementing this runner *OR* None if this runner is
            # unrecognized.
            runner_method = objects.get_method_or_none(
                obj=self, method_name=runner_method_name)

            # If this runner is unrecognized, raise an exception.
            if runner_method is None:
                raise BetseSimPipelineException(
                    '{} "{}" unrecognized.'.format(
                        self._label_singular_uppercase, runner_method_name))
            # Else, this runner is recognized.

            # Attempt to pass this runner these arguments.
            try:
                runner_method(runner_conf)
            # If this runner's requirements are unsatisfied (e.g., due to the
            # current simulation configuration disabling extracellular spaces),
            # ignore this runner with a non-fatal warning and continue.
            except BetseSimPipelineUnsatisfiedException as exception:
                logs.log_warn(
                    'Ignoring %s "%s", as:\n\t%s',
                    self._label_singular_lowercase,
                    runner_conf.name,
                    str(exception))

    # ..................{ PROPERTIES                         }..................
    @property
    def is_enabled(self) -> bool:
        '''
        ``True`` only if the :meth:`run` method should run this pipeline.

        Specifically, if this boolean is:

        * ``False``, the :meth:`run` method reduces to a noop.
        * ``True``, this method behaves as expected (i.e., calls all currently
          enabled runner methods).

        Defaults to ``True``. Pipeline subclasses typically override this
        property to return a boolean derived from the simulation configuration
        file associated with the current phase.
        '''

        return True

    # ..................{ PRIVATE ~ loggers                  }..................
    def _log_run(self) -> None:
        '''
        Log the current attempt to run the calling runner.
        '''

        # Defer to lower-level functionality to do so.
        self._die_unless_intra()

    # ..................{ PRIVATE ~ exceptions               }..................
    @type_check
    def _die_unless(
        self,
        is_satisfied: bool,
        exception_reason: str,
    ) -> None:
        '''
        Raise an exception containing the passed explanation if the passed
        boolean is ``False`` *or* log an attempt to create this animation
        otherwise.

        Parameters
        ----------
        is_satisfied : bool
            ``True`` only if all simulation features required by this runner
            (e.g., extracellular spaces) are currently available.
        exception_reason : str
            Uncapitalized human-readable string to be embedded in the messages
            of exceptions raised by this method, typically explaining all
            features required by this animation.

        Raises
        ----------
        BetseSimVisualUnsatisfiedException
            If this boolean is ``False``.
        '''

        # Name of the runner method either directly or indirectly calling this
        # method (e.g., "run_electric_extra").
        runner_method_name = callers.get_caller_basename_matching(
            predicate=lambda caller_basename:
                caller_basename.startswith(self._RUNNER_METHOD_NAME_PREFIX))

        # Strip the prefixing "run_" from this name, raising a human-readable
        # exception if this name has no such prefix.
        runner_name = strs.remove_prefix(
            text=runner_method_name,
            prefix=self._RUNNER_METHOD_NAME_PREFIX,
            exception_message=(
                'Callable {}() not a "{}"-prefixed runner method.'.format(
                    runner_method_name, self._RUNNER_METHOD_NAME_PREFIX)))

        # If these animation requirements are unsatisfied, raise an exception.
        if not is_satisfied:
            raise BetseSimPipelineUnsatisfiedException(
                '{} "{}" requirements unsatisfied (i.e., {}).'.format(
                    self._label_singular_uppercase,
                    runner_name,
                    exception_reason))

        # Log this attempt to run the calling runner.
        logs.log_info(
            '%s %s "%s"...',
            self._label_verb,
            self._label_singular_lowercase,
            runner_name)

    # ..................{ PRIVATE ~ exceptions : config      }..................
    def _die_unless_intra(self) -> None:
        '''
        Log an attempt to run the calling runner.

        This method is intended to be called by runners requiring that only
        intracellular spaces (i.e., a cell cluster) be enabled. As a cell
        cluster *always* exists, this method reduces to logging this runner.
        '''

        # Log this animation attempt. Since all simulations *ALWAYS* enable
        # support for intracellular spaces, no actual validation is required.
        self._die_unless(
            is_satisfied=True,
            exception_reason='This exception should never be raised.')


    def _die_unless_extra(self) -> None:
        '''
        Raise an exception unless extracellular spaces are enabled *or* log an
        attempt to run the calling runner otherwise.
        '''

        self._die_unless(
            is_satisfied=self._phase.p.sim_ECM,
            exception_reason='extracellular spaces disabled')


    @type_check
    def _die_unless_ion(self, ion_name: str) -> None:
        '''
        Raise an exception unless the ion with the passed name is enabled by the
        current ion profile *or* log an attempt to run the calling runner
        otherwise.

        Parameters
        ----------
        ion_name : str
            Capitalized alphabetic name of the ion required by this animation
            (e.g., ``Ca``, signifying calcium).
        '''

        # If this ion is unrecognized, raise a lower-level exception.
        if ion_name not in self._phase.p.ions_dict:
            raise BetseSimConfigException(
                'Ion "{}" unrecognized.'.format(ion_name))
        # Else, this ion is recognized.

        # Validate whether this ion is enabled or not.
        self._die_unless(
            is_satisfied=self._phase.p.ions_dict[ion_name] != 0,
            exception_reason='ion "{}" disabled'.format(ion_name))

    # ..................{ PRIVATE ~ subclass                 }..................
    @abstractproperty
    def _runners_conf_enabled(self) -> IterableTypes:
        '''
        Iterable of all currently enabled :class:`SimPipelineRunnerConf`
        instances, each encapsulating all input parameters to be passed to the
        method implementing a currently enabled runner in this pipeline.

        Pipeline subclasses typically implement this property to return an
        instance of the :class:``SimConfList`` class listing all runners enabled
        by the simulation configuration file associated with the current phase.
        '''

        pass

# ....................{ SUBCLASSES                         }....................
class SimPipelinerExportABC(SimPipelinerABC):
    '''
    Abstract base class of all **simulation export pipelines** (i.e., subclasses
    iteritavely exporting all variations on a single type of simulation export,
    either in parallel *or* in series).
    '''

    # ..................{ SUPERCLASS ~ constants             }..................
    _RUNNER_METHOD_NAME_PREFIX = 'export_'
    '''
    Substring prefixing the name of each runner defined by this pipeline.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters and an
        # exportation-specific verb.
        super().__init__(*args, label_verb='Exporting', **kwargs)

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
    def name(self) -> None:
        '''
        Lowercase alphanumeric string uniquely identifying the runner these
        arguments apply to in the parent simulation pipeline (e.g.,
        ``voltage_intra``, signifying an intracellular voltage runner).
        '''

# ....................{ DECORATORS                         }....................
@type_check
def runner_metadata(categories: SequenceTypes) -> CallableTypes:
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
        self, method: CallableTypes, categories: SequenceTypes) -> None:
        '''
        Initialize this class decorator.

        Parameters
        ----------
        method: CallableTypes
            Unbound method (i.e., function) to be decorated.
        categories : SequenceTypes
            Sequence of one or more human-readable category names.

        Raises
        ----------
        BetseSimPipelineException
            If this method has no docstring.
        '''

        # Initialize our superclass with the passed method.
        super().__init__(method)

        # Classify all remaining passed parameters.
        self.categories = categories

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
