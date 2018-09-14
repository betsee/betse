#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **simulation pipeline** (i.e., container of related, albeit
isolated, simulation actions to be run iteratively) functionality.
'''

# ....................{ IMPORTS                           }....................
from abc import ABCMeta
from betse.exceptions import (
    BetseSimPipeException, BetseSimPipeRunnerUnsatisfiedException)
from betse.lib.yaml.abc.yamllistabc import YamlListItemTypedABC
from betse.science.phase.phasecls import SimPhase
from betse.science.pipe.piperun import SimPipeRunner
from betse.util.io.log import logs
from betse.util.type.decorator.deccls import abstractmethod, abstractproperty
from betse.util.type.obj import objects
from betse.util.type.text import strs
from betse.util.type.types import type_check, GeneratorType, IterableTypes

# ....................{ SUPERCLASSES                      }....................
#FIXME: Revise docstring to account for the recent large-scale class redesign.
class SimPipeABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all subclasses running a **simulation pipeline**
    (i.e., sequence of similar simulation activities to be iteratively run).

    This class implements the infrastructure for iteratively running all
    currently enabled simulation activities (referred to as "runners")
    comprising the pipeline defined by this subclass, either in parallel *or*
    in series.

    Design
    ----------
    Each subclass is expected to define:

    * One or more public methods with names prefixed by the subclass-specific
      substring :attr:`_RUNNER_METHOD_NAME_PREFIX`, defaulting to ``run_``
      (e.g., ``run_voltage_intra``). For each such method, the name of that
      method excluding that prefix is the name of that method's runner (e.g.,
      ``voltage_intra`` for the method name ``run_voltage_intra``). Each such
      method should:

      * Accept exactly two parameters:

        * ``self``.
        * ``conf``, an object supplying this runner with its configuration.
          This is typically but *not* necessarily an instance of a subclass of
          the :class:`betse.lib.yaml.abc.yamlabc.YamlABC` superclass, enabling
          transparent configuration of runners from YAML-based files.

      * Return nothing (i.e.,``None``).

    * The abstract :meth:`iter_runners_conf` method returning a sequence of the
      names of all runners currently enabled by this pipeline (e.g.,
      ``['voltage_intra', 'ion_hydrogen', 'electric_total']``).

    The :meth:`run` method defined by this base class then dynamically
    implements this pipeline by iterating over the :meth:`iter_runners_conf`
    property and, for each enabled runner, calling that runner's method.

    Attributes (Private)
    ----------
    _phase : SimPhaseOrNoneTypes
        Current simulation phase if this pipeline is currently running (i.e.,
        if the :meth:`run` method is currently being called) *or* ``None``
        otherwise.

    Attributes (Private: Labels)
    ----------
    _noun_singular_lowercase : str
        Human-readable lowercase singular noun synopsizing the type of runners
        implemented by this subclass (e.g., ``animation``, ``plot``).
    _noun_singular_uppercase : str
        Human-readable singular noun whose first character is uppercase and all
        remaining characters lowercase (e.g., ``Animation``, ``Plot``).
    _noun_plural_lowercase : str
        Human-readable lowercase plural noun synopsizing the type of runners
        implemented by this subclass (e.g., ``animations``, ``plots``).
    '''

    # ..................{ CONSTANTS                         }..................
    #FIXME: Actually, the paragraph below no longer strictly applies. Our newly
    #defined @classproperty_readonly decorator now permits read-only class
    #properties to be trivially (and entirely safely) defined with no need for
    #intermediary metaclasses. Ergo, refactor the "_RUNNER_METHOD_NAME_PREFIX"
    #class variable into a comparable _runner_method_name_prefix() class
    #property decorated by the @classproperty_readonly decorator. Thunderstorm!
    #FIXME: Actually, just refactor this into a normal instance property after
    #refactoring the iter_runners() class method into a normal method as well.

    # Ideally, the following class constants would instead be implemented as
    # class properties. Unfortunately, there exists no @classproperty decorator
    # and no feasible means of implementing that decorator due to subtle
    # idiosyncracies in the design of Python's data descriptor protocol. While
    # class properties can technically be manually implemented by a metaclass
    # defining these properties, doing so is cumbersome and quite discouraged.

    _RUNNER_METHOD_NAME_PREFIX = 'run_'
    '''
    Substring prefixing the name of each runner defined by this pipeline.

    This class variable is intended to be overridden by subclasses desiring a
    less ambiguous prefix (e.g., ``export_`` for export pipelines).
    '''

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def __init__(self) -> None:
        '''
        Initialize this pipeline.
        '''

        # Classify all passed parameters.
        self._noun_singular_lowercase = self._noun_singular
        self._noun_singular_uppercase = strs.uppercase_char_first(
            self._noun_singular)
        self._noun_plural_lowercase = self._noun_plural

        # Nullify all remaining instance variables for safety.
        self._phase = None


    def _init_run(self) -> None:
        '''
        Initialize this pipeline for the current call to the :meth:`run`
        method, which calls this method *before* performing subsequent logic.

        Defaults to a noop. Pipeline subclasses may override this method to
        guarantee state assumed by subsequently run pipeline runners (e.g., to
        create external directories required by these runners).
        '''

        pass

    # ..................{ SUBCLASS ~ methods                }..................
    @abstractmethod
    def iter_runners_conf(self, phase: SimPhase) -> IterableTypes:
        '''
        Iterable of all **runner configurations** (i.e.,
        :class:`YamlListItemTypedABC` instances encapsulating all input
        parameters to be passed to the corresponding pipeline runner) for the
        passed simulation phase.

        Pipeline subclasses typically implement this property to return an
        instance of the :class:``YamlList`` class listing all runners listed
        by the simulation configuration file configuring this phase.

        Parameters
        ----------
        phase : SimPhase
            Current simulation phase.
        '''

        pass

    # ..................{ SUBCLASS ~ properties             }..................
    @abstractproperty
    def _noun_singular(self) -> str:
        '''
        Human-readable singular noun synopsizing the type of runners
        implemented by this subclass (e.g., ``animation``, ``plot``), ideally
        but *not* necessarily lowercase.
        '''

        pass

    # ..................{ PROPERTIES                        }..................
    @property
    def name(self) -> str:
        '''
        Human-readable name of this pipeline, intended principally for display
        (e.g., logging) purposes.

        Defaults to suffixing :meth:`_noun_singular` by `` pipeline``.
        '''

        return '{} pipeline'.format(self._noun_singular)

    # ..................{ PROPERTIES ~ private              }..................
    @property
    def _is_enabled(self) -> bool:
        '''
        ``True`` only if the currently called :meth:`run` method should run
        this pipeline.

        Specifically, if this boolean is:

        * ``False``, the :meth:`run` method reduces to a noop.
        * ``True``, this method behaves as expected (i.e., calls all currently
          enabled runner methods).

        Defaults to ``True``. Pipeline subclasses typically override this
        property to return a boolean derived from the simulation configuration
        file associated with the current phase.
        '''

        return True

    # ..................{ PROPERTIES ~ private : str        }..................
    @property
    def _noun_plural(self) -> str:
        '''
        Human-readable plural noun synopsizing the type of runners implemented
        by this subclass (e.g., ``animations``, ``plots``), ideally but *not*
        necessarily lowercase.

        Defaults to suffixing :meth:`_noun_singular` by ``s``.
        '''

        return self._noun_singular + 's'


    @property
    def _verb_continuous(self) -> str:
        '''
        Human-readable verb in the continuous tense synopsizing the type of
        action performed by runners implemented by this subclass (e.g.,
        ``Saving``), ideally but *not* necessarily capitalized.

        Defaults to ``Running``.
        '''

        return 'Running'

    # ..................{ RUNNERS                           }..................
    @type_check
    def run(self, phase: SimPhase) -> None:
        '''
        Run all currently enabled pipeline runners for the passed simulation
        phase if this pipeline is currently enabled in this phase *or* noop
        otherwise.

        Specifically:

        * If the :meth:`_is_enabled` property is ``True`` (implying this
          pipeline to be currently enabled):

          * For each :class:`YamlListItemTypedABC` instance (corresponding to
            the subconfiguration of a currently enabled pipeline runner) in the
            sequence of these instances listed by the :meth:`iter_runners_conf`
            property:

            * Call this method, passed this configuration.
            * If this method reports this runner's requirements to be
              unsatisfied (e.g., due to the current simulation configuration
              disabling extracellular spaces), this runner is ignored with a
              non-fatal warning.

        * Else, log an informative message and return immediately.

        Parameters
        ----------
        phase : SimPhase
            Current simulation phase.
        '''

        # Temporarily classify the passed phase.
        self._phase = phase

        # Initialize this pipeline for the current call to this method *AFTER*
        # classifying this phase, as the former typically requires the latter.
        self._init_run()

        # If this pipeline is disabled, log this fact and return immediately.
        if not self._is_enabled:
            logs.log_info('Excluding %s...', self._noun_plural_lowercase)
            return
        # Else, this pipeline is enabled.

        # Log this pipeline run.
        logs.log_info(
            '%s %s...', self._verb_continuous, self._noun_plural_lowercase)

        # For the object encapsulating all input arguments to be passed to each
        # currently enabled runner in this pipeline...
        for runner_conf in self.iter_runners_conf(self._phase):
            if not isinstance(runner_conf, YamlListItemTypedABC):
                raise BetseSimPipeException(
                    'iter_runners_conf() item {!r} '
                    'not instance of "YamlListItemTypedABC".'.format(
                        runner_conf))

            # If this runner is disabled, log this fact and ignore this runner.
            if not runner_conf.is_enabled:
                logs.log_debug(
                    'Ignoring disabled %s "%s"...',
                    self._noun_singular_lowercase,
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
                        self._noun_singular_uppercase, runner_conf.name))
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
                    self._noun_singular_lowercase,
                    runner_conf.name,
                    exception.reason)

        # Declassify the passed phase for safety.
        self._phase = None

    # ..................{ EXCEPTIONS                        }..................
    @type_check
    def die_unless_runner_satisfied(self, runner: SimPipeRunner) -> None:
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
                        self._noun_singular_uppercase, runner_name),
                    reason='{} disabled'.format(requirement.name),
                )
        # Else, all runner requirements are satisfied.

        # Log the subsequent attempt to run this runner.
        logs.log_info(
            '%s %s "%s"...',
            self._verb_continuous, self._noun_singular_lowercase, runner_name)

    # ..................{ ITERATORS                         }..................
    def iter_runners(self) -> GeneratorType:
        '''
        Generator yielding the 2-tuple ``(runner_name, runner)`` for each
        **runner** (i.e., :class:`SimPipeRunner` instance produced by the
        :func:`piperunner` decorator decorating that runner's method) defined
        by this pipeline subclass.

        This generator excludes all methods defined by this pipeline subclass
        *not* decorated by that decorator.

        Yields
        ----------
        (str, SimPipeRunner)
            2-tuple ``(runner_name, runner)`` where:

            * ``runner_name`` is the name of the method underlying this runner,
              excluding the substring :attr:`_RUNNER_METHOD_NAME_PREFIX`
              prefixing this name.
            * ``runner`` is the :class:`SimPipeRunner` instance defining this
              runner.
        '''

        # Return a generator comprehension...
        return (
            # Iteratively yielding a 2-tuple of:
            #
            # * The name of this runner's method excluding runner prefix.
            # * Each "SimPipeRunner" instance defined on this class.
            (
                strs.remove_prefix(
                    text=runner_method_name,
                    prefix=self._RUNNER_METHOD_NAME_PREFIX,
                ),
                runner
            )

            # For each such runner's method name and value.
            for runner_method_name, runner in objects.iter_attrs_instance_of(
                obj=self, cls=SimPipeRunner)
        )
