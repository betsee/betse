#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **simulation pipeline** (i.e., container of related, albeit
isolated, simulation actions to be run iteratively) functionality.
'''

#FIXME: Metadata duplication between pipeline and pipeline runner classes has
#reached Corinthian heights of inanity. To address both this and the BETSEE
#export substatus issue (i.e., that substatus for exports is currently kludged
#at the BETSE layer and hence untranslateable), consider:
#
#* Define a new "betse.util.type.call.callbacks.CallbacksMetadataABC" abstract
#  base class whose body reduces to simply "pass".
#* Refactor the "betse.util.type.call.callbacks.CallbacksBC" class as follows:
#  * Define a new progress_metadated() method accepting a mandatory
#    "metadata : CallbacksMetadataABC" parameter.
#  * Improve the progressed_next() method to accept an optional parameter of
#    the same type.
#* Define a new "SimPipeMetadataABC" abstract base class in this submodule,
#  subclassing "CallbacksMetadataABC" for improved generality.
#  Shift all of the following "SimPipeABC" class properties into
#  "SimPipeMetadataABC" *INSTANCE* (i.e., standard non-class) properties:
#  * "_NOUN_SINGULAR", renamed to "noun_singular_lower".
#  * "_VERB_CONTINUOUS", renamed to "verb_continuous".
#  * Lastly, define a new "noun_singular_upper" property.
#* Define a new SimPipeABC.metadata() abstract class property. All "SimPipeABC"
#  subclasses should redefine this property by:
#  * Defining a new "SimPipeMetadataABC" concrete subclass.
#  * Overriding SimPipeABC.metadata() to return an instance of that subclass.
#* Remove the SimPipeRunnerMetadata.noun* and verb* instance variables.
#* Define a new "SimPipeRunnerMetadata.pipe" instance variable.
#* Refactor the SimPipeABCMeta.__new__() method to set this variable on each
#  pipeline runner.
#* Refactor the SimPipesExport.export() method to pass the optional "metadata"
#  parameter to the call to phase.callbacks.progressed_next(): e.g.,
#      phase.callbacks.progressed_next(metadata=runner_metadata.pipe)
#
#Voila! Frankly trivial, but fairly involved. A real-world dichotomy in action.

# ....................{ IMPORTS                           }....................
from abc import ABCMeta
from betse.exceptions import BetseSimPipeException
from betse.science.config.export.confexpabc import SimConfExportABC
from betse.science.phase.phasecls import SimPhase
from betse.util.io.log import logs
from betse.util.type.decorator.deccls import abstractmethod #, abstractproperty
from betse.util.type.descriptor.descs import (
    abstractclassproperty_readonly, classproperty_readonly)
from betse.util.type.iterable import itersort
from betse.util.type.obj import objects, objiter, objtest
from betse.util.type.text.string import strs
from betse.util.type.types import (
    type_check,
    ClassType,
    GeneratorType,
    IterableTypes,
    MappingType,
    SequenceTypes,
)

# ....................{ METACLASSES                       }....................
class SimPipeABCMeta(ABCMeta):
    '''
    Metaclass of the abstract :class:`SimPipeABC` base class and all concrete
    subclasses thereof.

    This metaclass dynamically annotates each **simulation pipeline runner**
    (i.e., method decorated by the :func:`piperunner` decorator) on the initial
    definition of each concrete subclass of the :class:`SimPipeABC` superclass
    with metadata pertaining to that subclass, augmenting the core
    runner-specific metadata already annotated on that runner by that
    decorator. (The metadata annotated by this metaclass requires access to the
    class declaring that runner, which the :func:`piperunner` decorator has no
    means of accessing under conventional semantics of method declaration.)

    See Also
    ----------
    https://stackoverflow.com/a/47371003/2809027
        StackOverflow answer mildly inspiring this class.
    '''

    # ..................{ CONSTRUCTORS                      }..................
    def __new__(
        metacls: ClassType,
        class_name: str,
        class_base_classes: IterableTypes,
        class_attrs: MappingType,
        **kwargs
    ) -> ClassType:
        '''
        Annotate each simulation pipeline runner defined by the passed
        simulation pipeline subclass with metadata pertaining to that subclass.
        '''

        # Unsanitized "SimPipeABC" subclass.
        subcls = super().__new__(
            metacls, class_name, class_base_classes, class_attrs)

        # For the method name and method of each runner in this pipeline...
        for runner_method_name, runner_method in subcls.iter_runners_method():
            # Metadata associated with this runner.
            runner_metadata = runner_method.metadata

            # Set this runner's type to the name of this method excluding this
            # pipeline's runner prefix.
            runner_metadata.kind = strs.remove_prefix(
                text=runner_method_name,
                prefix=subcls._RUNNER_METHOD_NAME_PREFIX)

            # Set this runner's human-readable type to that of this pipeline.
            runner_metadata.noun_singular_lowercase = subcls._NOUN_SINGULAR
            runner_metadata.noun_singular_uppercase = (
                subcls._NOUN_SINGULAR.upper())
            runner_metadata.verb_continuous = subcls._VERB_CONTINUOUS

        # Return this sanitized "SimPipeABC" subclass.
        return subcls

# ....................{ SUPERCLASSES                      }....................
#FIXME: Revise docstring to account for the recent large-scale class redesign.
class SimPipeABC(object, metaclass=SimPipeABCMeta):
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
      :attr:`_RUNNER_METHOD_NAME_PREFIX`, defaulting to ``run_`` (e.g.,
      ``run_voltage_intra``). For each such method, the name of that method
      excluding that prefix is the name of that method's runner (e.g.,
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

    # ..................{ SUBCLASS ~ properties             }..................
    # Subclasses are required to implement the following abstract properties.

    @abstractclassproperty_readonly
    def _NOUN_SINGULAR(cls) -> str:
        '''
        Human-readable singular noun synopsizing the type of runners
        implemented by this subclass (e.g., ``animation``, ``plot``), ideally
        but *not* necessarily lowercase.
        '''

        pass

    # ..................{ CLASS ~ properties : concrete     }..................
    @classproperty_readonly
    def _NOUN_PLURAL(cls) -> str:
        '''
        Human-readable plural noun synopsizing the type of runners implemented
        by this subclass (e.g., ``animations``, ``plots``), ideally but *not*
        necessarily lowercase.

        Defaults to suffixing :meth:`_NOUN_SINGULAR` by ``s``.
        '''

        return cls._NOUN_SINGULAR + 's'


    @classproperty_readonly
    def _VERB_CONTINUOUS(cls) -> str:
        '''
        Human-readable verb in the continuous tense synopsizing the type of
        action performed by runners implemented by this subclass (e.g.,
        ``Saving``), ideally but *not* necessarily capitalized.

        Defaults to ``Running``.
        '''

        return 'Running'


    @classproperty_readonly
    def _RUNNER_METHOD_NAME_PREFIX(cls) -> str:
        '''
        Substring prefixing the name of each runner defined by this pipeline.

        Defaults to ``run_``.
        '''

        return  'run_'

    # ..................{ CLASS ~ iterators                 }..................
    @classmethod
    def iter_runners_method(cls) -> GeneratorType:
        '''
        Generator yielding the method name and method of each **simulation
        pipeline runner** (i.e., method bound to this pipeline decorated by the
        :func:`piperunner` decorator).

        This generator excludes all methods defined by this pipeline subclass
        *not* decorated by that decorator.

        Yields
        ----------
        (runner_method_name : str, runner_method : MethodType)
            2-tuple where:

            * ``runner_method_name`` is the name of the method underlying this
              runner, identical to the ``runner_method.__name__`` attribute.
            * ``runner_method`` is this method.
        '''

        # Defer to this generator.
        yield from objiter.iter_methods_prefixed(
            obj=cls, prefix=cls._RUNNER_METHOD_NAME_PREFIX)

    # ..................{ CLASS ~ iterators : metadata      }..................
    @classmethod
    def iter_runners_metadata(cls) -> GeneratorType:
        '''
        Generator yielding the metadata annotating each **simulation pipeline
        runner** (i.e., method bound to this pipeline decorated by the
        :func:`piperunner` decorator).

        This generator excludes all methods defined by this pipeline subclass
        *not* decorated by that decorator.

        Yields
        ----------
        runner_metadata : SimPipeRunnerMetadata
            :class:`SimPipeRunnerMetadata` instance collecting all metadata for
            this runner.
        '''

        # Return a generator comprehension iteratively yielding...
        return (
            # Metadata annotating this runner.
            runner_method.metadata
            # For each runner method defined by this pipeline subclass...
            for _, runner_method in cls.iter_runners_method()
        )


    @classmethod
    @type_check
    def iter_runners_metadata_kind(cls) -> SequenceTypes:
        '''
        Sequence of each **machine-readable type** (i.e., string value of the
        :attr:`SimPipeRunnerMetadata.kind` instance variable) of all runners
        supported by this pipeline (in sorted lexicographic order).
        '''

        # Return a lexicographically sorted sequence of...
        return itersort.sort_ascending(
            # The machine-readable type of this runner.
            runner_metadata.kind
            # For each runner defined by this pipeline subclass...
            for runner_metadata in cls.iter_runners_metadata()
        )

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def __init__(self) -> None:
        '''
        Initialize this pipeline.
        '''

        # Classify all passed parameters.
        self._noun_singular_lowercase = self._NOUN_SINGULAR
        self._noun_singular_uppercase = strs.uppercase_char_first(
            self._NOUN_SINGULAR)
        self._noun_plural_lowercase = self._NOUN_PLURAL


    @type_check
    def init(self, phase: SimPhase) -> None:
        '''
        Initialize this pipeline for the passed simulation phase.

        Defaults to a noop. Pipeline subclasses may override this method to
        guarantee state assumed by subsequently run pipeline runners (e.g., to
        create external directories required by these runners for this phase).

        Parameters
        ----------
        phase : SimPhase
            Current simulation phase.
        '''

        pass

    # ..................{ SUBCLASS ~ methods                }..................
    # Subclasses are required to implement the following abstract methods.

    @abstractmethod
    def iter_runners_conf(self, phase: SimPhase) -> SequenceTypes:
        '''
        Sequence of all **runner configurations** (i.e.,
        :class:`SimConfExportABC` instances encapsulating all input parameters
        passed to the corresponding pipeline runner) for the passed simulation
        phase.

        Pipeline subclasses typically implement this property to return an
        instance of the :class:``YamlList`` class listing all runners listed
        by the simulation configuration file configuring this phase.

        Caveats
        ----------
        The existence of a runner configuration does *not* imply the
        corresponding pipeline runner to be unconditionally enabled. Instead,
        the :attr:`SimConfExportABC.is_enabled` data descriptor defined by all
        runner configurations returned by this method specifies whether or not
        that runner is to be enabled or disabled.

        Parameters
        ----------
        phase : SimPhase
            Current simulation phase.
        '''

        pass

    # ..................{ PROPERTIES                        }..................
    @property
    def name(self) -> str:
        '''
        Human-readable name of this pipeline, intended principally for display
        (e.g., logging) purposes.

        Defaults to suffixing :meth:`_NOUN_SINGULAR` by `` pipeline``.
        '''

        return '{} pipeline'.format(self._NOUN_SINGULAR)

    # ..................{ PROPERTIES ~ private              }..................
    @type_check
    def _is_enabled(self, phase: SimPhase) -> bool:
        '''
        ``True`` only if this pipeline is enabled for the passed simulation
        phase.

        Specifically, if this boolean is:

        * ``False``, the :meth:`run` method reduces to a noop.
        * ``True``, that method behaves as expected (i.e., calls all currently
          enabled runner methods).

        Defaults to ``True``. Pipeline subclasses typically override this
        property to return a boolean derived from the simulation configuration
        file associated with the passed phase.

        Parameters
        ----------
        phase : SimPhase
            Current simulation phase.
        '''

        return True

    # ..................{ ITERATORS                         }..................
    @type_check
    def iter_runners_enabled(self, phase: SimPhase) -> GeneratorType:
        '''
        Generator yielding the method and simulation configuration of each
        **enabled simulation pipeline runner** (i.e., method bound to this
        pipeline decorated by the :func:`piperunner` decorator *and* enabled by
        this simulation configuration) for the passed simulation phase.

        If this simulation phase disables this pipeline, this generator reduces
        to the empty generator (i.e., yields no values); else, this generator
        excludes:

        * All methods defined by this pipeline subclass *not* decorated by the
          :func:`piperunner` decorator.
        * All methods defined by this pipeline subclass decorated by that
          decorator currently disabled by their simulation configuration (e.g.,
          such that the YAML-backed value of the ``enabled`` key in this
          configuration evaluates to ``False``).

        Parameters
        ----------
        phase : SimPhase
            Current simulation phase.

        Yields
        ----------
        (runner_method : MethodType, runner_conf : SimConfExportABC)
            2-tuple where:

            * ``runner_method`` is the method implementing this runner, whose
              method name is guaranteed to be prefixed by the substring
              :attr:`_RUNNER_METHOD_NAME_PREFIX`. Note that this method object
              defines the following custom instance variables:

              * ``metadata``, whose value is an instance of the
                :class:`SimPipeRunnerMetadata` class.

            * ``runner_conf`` is the YAML-backed list item configuring this
              runner, derived from the YAML-formatted simulation configuration
              file associated with the passed simulation phase.

        Raises
        ----------
        BetseSimPipeException
            If any simulation runner configured for this pipeline by the passed
            simulation phase is **unrecognized** (i.e., if this pipeline
            defines no corresponding method).
        '''

        # If this pipeline is disabled, log this fact and return immediately.
        if not self._is_enabled(phase):
            logs.log_debug(
                'Ignoring disabled %s...', self._noun_plural_lowercase)
            return
        # Else, this pipeline is enabled.

        # For each runner configuration specified for this pipeline...
        for runner_conf in self.iter_runners_conf(phase):
            # If this configuration is *NOT* YAML-backed, raise an exception.
            objtest.die_unless_instance(obj=runner_conf, cls=SimConfExportABC)

            # If this runner is disabled, log this fact and ignore this runner.
            if not runner_conf.is_enabled:
                logs.log_debug(
                    'Ignoring disabled %s "%s"...',
                    self._noun_singular_lowercase,
                    runner_conf.kind)
                continue
            # Else, this runner is enabled.

            # Name of the pipeline method implementing this runner.
            runner_method_name = (
                self._RUNNER_METHOD_NAME_PREFIX + runner_conf.kind)

            # Method running this runner if recognized *OR* "None" otherwise.
            runner_method = objects.get_callable_or_none(
                obj=self, callable_name=runner_method_name)

            # If this runner is unrecognized, raise an exception.
            if runner_method is None:
                raise BetseSimPipeException(
                    '{} "{}" unrecognized '
                    '(i.e., method {}.{}() not found).'.format(
                        self._noun_singular_uppercase,
                        runner_conf.kind,
                        objects.get_class_name_unqualified(self),
                        runner_method_name))
            # Else, this runner is recognized.

            # Yield a 2-tuple of this runner's method and configuration.
            yield runner_method, runner_conf
