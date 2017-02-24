#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level **simulation pipeline** (i.e., set of one or more processing ) functionality.
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
    type_check, CallableTypes, GeneratorType, SequenceTypes)

# ....................{ SUPERCLASSES                       }....................
class SimPipelinerABC(object, metaclass=ABCMeta):
    '''
    Abstract base class of all **simulation pipelines** (i.e., subclasses
    iteritavely running all implementations of a single simulation activity,
    either in parallel *or* in series).

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
    * The abstract :meth:`runner_enabled_names` property returning a sequence of the
      names of all runners currently enabled by this pipeline (e.g.,
      ``['voltage_intra', 'ion_hydrogen', 'electric_total']``).

    The :meth:`run` method defined by this base class then dynamically
    implements this pipeline by iterating over the :meth:`runner_enabled_names` property
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
    def iter_runner_method_names(cls) -> GeneratorType:
        '''
        Generator yielding the name of each runner method defined by this
        pipeline subclass.

        For each subclass method whose name is prefixed by
        :attr:`_RUNNER_METHOD_NAME_PREFIX`, this generator yields that method.

        Yields
        ----------
        MethodType
            Each runner method defined by this pipeline subclass.
        '''

        yield from classes.iter_methods_matching(
            cls=cls, predicate=lambda method_name: method_name.startswith(
                cls._RUNNER_METHOD_NAME_PREFIX))


    #FIXME: Is this genuinely useful? We think not, actually. Defining an
    #iterator mapping from each runner method to the human-readable metadata
    #associated with that method would be far more useful.
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

        # For runner method defined by this subclass...
        for anim_method_name in cls.iter_runner_method_names():
            # Yield the name of this method excluding the runner prefix.
            yield strs.remove_prefix(
                text=anim_method_name, prefix=cls._RUNNER_METHOD_NAME_PREFIX)

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
        Run this pipeline.

        Specifically, this method runs all currently enabled runners in this
        pipeline.
        '''

        # Log animation creation.
        logs.log_info(
            '%s %s...', self._label_verb, self._label_plural_lowercase)

        # For the name of each currently enabled runner in this pipeline...
        for runner_name in self.runner_enabled_names:
            # Name of the method implementing this runner.
            runner_method_name = self._RUNNER_METHOD_NAME_PREFIX + runner_name

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

            # Attempt to run this runner.
            try:
                runner_method()
            # If this runner's requirements are unsatisfied (e.g., due to the
            # current simulation configuration disabling extracellular spaces),
            # ignore this runner with a non-fatal warning and continue.
            except BetseSimPipelineUnsatisfiedException as exception:
                logs.log_warn(
                    'Ignoring %s "%s", as:\n\t%s',
                    self._label_singular_lowercase,
                    runner_name,
                    str(exception))

    # ..................{ LOGGERS                            }..................
    def _log_run(self) -> None:
        '''
        Log the current attempt to run the calling runner.
        '''

        # Defer to lower-level functionality to do so.
        self._die_unless_intra()

    # ..................{ EXCEPTIONS                         }..................
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

    # ..................{ EXCEPTIONS ~ config                }..................
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

    # ..................{ SUBCLASS                           }..................
    @abstractproperty
    def runner_enabled_names(self) -> SequenceTypes:
        '''
        Sequence of the names of all currently enabled runners in this pipeline.

        Pipeline subclasses typically implement this property to return the
        user-defined sequence of the names of all runners listed in the
        simulation configuration file associated with the current phase.
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

# ....................{ DECORATORS                         }....................
@type_check
def piperunner(categories: SequenceTypes) -> CallableTypes:
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
    def _piperunner_closure(method: CallableTypes) -> SimPipelineRunner:
        '''
        Closure annotating simulation pipeline runners with custom metadata,
        returning an instance of the class decorator exposing this metadata to
        external interfaces.

        See Also
        ----------
        :func:`piperunner`
            Further details.
        '''

        # Return
        return SimPipelineRunner(method=method, categories=categories)

    # Return the closure accepting the method to be decorated.
    return _piperunner_closure


class SimPipelineRunner(MethodDecorator):
    '''
    Class decorator annotating simulation pipeline runners with custom metadata.

    All such runners decorated by the :func:`piperunner` decorator are
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
    :func:`piperunner`
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
        categories : SequenceTypes
            Sequence of one or more human-readable category names.
        method: CallableTypes
            Unbound method (i.e., function) to be decorated.

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

        # Reduce this docstring from a possibly multiline string to a
        # single-line string containing no newlines.
        self.description = strs.unwrap(self.description)
