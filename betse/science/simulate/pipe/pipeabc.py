#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **simulation pipeline** (i.e., sequence of similar simulation
activities to be iteratively run) functionality.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta, abstractproperty
from betse.exceptions import (
    BetseSimPipeException, BetseSimPipeRunnerUnsatisfiedException)
from betse.lib.yaml.abc.yamllistabc import YamlListItemTypedABC
from betse.science.simulate.pipe.piperun import SimPipeRunner
from betse.science.simulate.simphase import SimPhase
from betse.util.io.log import logs
from betse.util.type.obj import objects
from betse.util.type.text import strs
from betse.util.type.types import type_check, GeneratorType, IterableTypes

# ....................{ SUPERCLASSES                       }....................
class SimPipeABC(object, metaclass=ABCMeta):
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
    * The abstract :meth:`_runners_conf` property returning a sequence of the
      names of all runners currently enabled by this pipeline (e.g.,
      ``['voltage_intra', 'ion_hydrogen', 'electric_total']``).

    The :meth:`run` method defined by this base class then dynamically
    implements this pipeline by iterating over the :meth:`_runners_conf` property
    and, for each enabled runner, calling that runner's method.

    Attributes (Private)
    ----------
    _phase : SimPhase
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
    def iter_runners(cls) -> GeneratorType:
        '''
        Generator yielding the 2-tuple ``(runner_name, runner)`` for each
        **runner** (i.e., :class:`SimPipeRunner` instance produced by the
        :func:`piperunner` decorator decorating this runner's method)
        defined by this pipeline subclass.

        This generator excludes all methods defined by this pipeline subclass
        *not* decorated by this decorator.

        Yields
        ----------
        (str, SimPipeRunner)
            2-tuple ``(runner_name, runner)`` where:
            * ``runner_name`` is the name of this runner's underlying method
              excluding the substring :attr:`_RUNNER_METHOD_NAME_PREFIX`
              prefixing this name.
            * ``runner`` is each runner's :class:`SimPipeRunner` instance.
        '''

        # Return a generator comprehension...
        return (
            # Yielding a 2-tuple of:
            #
            # * The name of this runner's method excluding runner prefix.
            # * Each "SimPipeRunner" instance defined on this class.
            (
                strs.remove_prefix(
                    text=runner_method_name,
                    prefix=cls._RUNNER_METHOD_NAME_PREFIX),
                runner
            )
            # For each such runner's method name and instance.
            for runner_method_name, runner in objects.iter_attrs_matching(
                obj=cls,
                predicate=lambda attr_name, attr_value: (
                    isinstance(attr_value, SimPipeRunner))))

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(
        self,

        # Mandatory parameters.
        phase: SimPhase,
        label_singular: str,

        # Optional parameters.
        label_plural: str = None,
        label_verb: str = 'Running',
    ) -> None:
        '''
        Initialize this pipeline.

        Parameters
        ----------
        phase : SimPhase
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

    # ..................{ RUNNERS                            }..................
    def run(self) -> None:
        '''
        Run all currently enabled pipeline runners if this pipeline itself is
        currently enabled *or* noop otherwise.

        Specifically:

        * If the :meth:`is_enabled` property is ``True`` (implying this pipeline
          to be currently enabled):
          * For each :class:`YamlListItemTypedABC` instance (corresponding to
            the subconfiguration of a currently enabled pipeline runner) in the
            sequence of these instances listed by the :meth:`_runners_conf`
            property:
            * Call this method, passed this configuration.
            * If this method reports this runner's requirements to be
              unsatisfied (e.g., due to the current simulation configuration
              disabling extracellular spaces), this runner is ignored with a
              non-fatal warning.
        * Else, log an informative message and return immediately.
        '''

        # If this pipeline is disabled, log this fact and return immediately.
        if not self.is_enabled:
            logs.log_info('Excluding %s...', self._label_plural_lowercase)
            return
        # Else, this pipeline is enabled.

        # Log this pipeline run.
        logs.log_info(
            '%s %s...', self._label_verb, self._label_plural_lowercase)

        # For the object encapsulating all input arguments to be passed to each
        # currently enabled runner in this pipeline...
        for runner_conf in self._runners_conf:
            if not isinstance(runner_conf, YamlListItemTypedABC):
                raise BetseSimPipeException(
                    '_runners_conf() item {!r} '
                    'not instance of "YamlListItemTypedABC".'.format(
                        runner_conf))

            # If this runner is disabled, log this fact and ignore this runner.
            if not runner_conf.is_enabled:
                logs.log_debug(
                    'Ignoring disabled %s "%s"...',
                    self._label_singular_lowercase,
                    runner_conf.name)
                continue
            # Else, this runner is enabled.

            # Name of the method implementing this runner.
            runner_method_name = (
                self._RUNNER_METHOD_NAME_PREFIX + runner_conf.name)

            # Method implementing this runner *OR* None if this runner is
            # unrecognized.
            runner_method = objects.get_method_or_none(
                obj=self, method_name=runner_method_name)

            # If this runner is unrecognized, raise an exception.
            if runner_method is None:
                raise BetseSimPipeException(
                    '{} "{}" unrecognized.'.format(
                        self._label_singular_uppercase, runner_conf.name))
            # Else, this runner is recognized.

            # Attempt to pass this runner these arguments.
            try:
                runner_method(runner_conf)
            # If this runner's requirements are unsatisfied (e.g., due to the
            # current simulation configuration disabling extracellular spaces),
            # ignore this runner with a non-fatal warning and continue.
            except BetseSimPipeRunnerUnsatisfiedException as exception:
                logs.log_warning(
                    'Excluding %s "%s", as %s.',
                    self._label_singular_lowercase,
                    runner_conf.name,
                    exception.reason)

    # ..................{ EXCEPTIONS                         }..................
    @type_check
    def die_unless_runner_satisfied(
        self, runner: SimPipeRunner) -> None:
        '''
        Raise an exception if the passed runner is **unsatisfied** (i.e.,
        requires one or more simulation features disabled for the current
        simulation phase) *or* log an attempt to run this runner otherwise.

        Parameters
        ----------
        runner : SimPipeRunner
            Simulation pipeline runner to be tested.

        Raises
        ----------
        BetseSimVisualUnsatisfiedException
            If this runner is unsatisfied.
        '''

        # Strip the prefixing verb from the name of this runner's method (e.g.,
        # "export_" from "export_voltage_total"), raising a human-readable
        # exception if this name has no such prefix.
        runner_name = strs.remove_prefix(
            text=runner.method_name,
            prefix=self._RUNNER_METHOD_NAME_PREFIX,
            exception_message=(
                'Runner method name "{}" not prefixed by "{}".'.format(
                    runner.method_name, self._RUNNER_METHOD_NAME_PREFIX)))

        # If any runner requirement is unsatisfied, raise an exception.
        for requirement in runner.requirements:
            if not requirement.is_satisfied(phase=self._phase):
                raise BetseSimPipeRunnerUnsatisfiedException(
                    result='{} "{}" requirement unsatisfied'.format(
                        self._label_singular_uppercase, runner_name),
                    reason='{} disabled'.format(requirement.name),
                )
        # Else, all runner requirements are satisfied.

        # Log the subsequent attempt to run this runner.
        logs.log_info(
            '%s %s "%s"...',
            self._label_verb, self._label_singular_lowercase, runner_name)

    # ..................{ PRIVATE ~ subclass                 }..................
    @abstractproperty
    def _runners_conf(self) -> IterableTypes:
        '''
        Iterable of all :class:`YamlListItemTypedABC` instances for the current
        pipeline, each encapsulating all input parameters to be passed to the
        method implementing a runner currently contained in this pipeline.

        Pipeline subclasses typically implement this property to return an
        instance of the :class:``YamlList`` class listing all runners listed
        by the simulation configuration file associated with the current phase.
        '''

        pass

# ....................{ SUBCLASSES                         }....................
class SimPipeExportABC(SimPipeABC):
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
