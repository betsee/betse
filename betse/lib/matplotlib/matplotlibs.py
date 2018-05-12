#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level support facilities for matplotlib, a mandatory runtime dependency.

Backends
----------
matplotlib supports numerous **interactive backends** (i.e., bindings to
external GUI-specific widget toolkits), only one of which will be imported by
matplotlib at runtime. If the caller specifies no such backend, a default
backend specific to the current system will be imported. However, all such
backends including such defaults are fairly fragile and hence prone to raising
exceptions under common contexts.

The following table summarizes our current findings:

=========  ========  ======  ========  ======  ========  ======  ========
Backend    Footnote  Is Supported Under?
                     ----------------------------------------------------
                     Linux             OS X              Windows
                     ----------------  ----------------  ----------------
                     Frozen  Unfrozen  Frozen  Unfrozen  Frozen  Unfrozen
=========  ========  ======  ========  ======  ========  ======  ========
CocoaAgg   1         No      No        No      No        No      No
Gtk3Agg    2         No      No        No      No        No      No
Gtk3Cairo  3         No      No        ???     ???       ???     ???
MacOSX               No      No        Yes     Yes       No      No
Qt4Agg     3         No      No        ???     ???       ???     ???
Qt5Agg               ???     ???       ???     ???       ???     ???
TkAgg                Yes     Yes       No      Yes       ???     ???
=========  ========  ======  ========  ======  ========  ======  ========

Footnote descriptions are as follows:

1. Backend "CocoaAgg" is officially deprecated and known to be broken.
2. Backend "Gtk3Agg" is known to be broken under Python 3.
3. These backends do *not* support our current animation method, despite
   otherwise working (e.g., to display static plots).
'''

#FIXME: Refactor backend_names() to discover backend names via the standard
#module "pkg_utils" rather than by manually delving through the filesystem,
#which fails under frozen executables.

#FIXME: It'd be great to raise human-readable exceptions on the specified
#backends *NOT* being available. This is certainly feasible, as the
#following stackoverflow answer demonstrates -- if somewhat involved:
#    https://stackoverflow.com/questions/5091993/list-of-all-available-matplotlib-backends
#That said, we really want to do this *ANYWAY* to print this list when running
#"betse info". So, let's just get this done, please.

#FIXME: Consider contributing most or all of this submodule back to matplotlib.

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To permit matplotlib's default verbosity and backend to be replaced
# by BETSE-specific values, no matplotlib package or module may be imported
# until *AFTER* the MplConfig.init() method has been called -- including
# here at the top-level.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import sys
from betse.exceptions import BetseMatplotlibException
from betse.util.io import ioexceptions
from betse.util.io.log import logconfig, logs
from betse.util.io.log.logenum import LogLevel
from betse.util.os import displays, kernels, oses
from betse.util.path import dirs, pathnames
from betse.util.py import pyfreeze
from betse.util.type import iterables, modules
from betse.util.type.decorator.decmemo import property_cached
from betse.util.type.mapping.mapcls import OrderedArgsDict
from betse.util.type.numeric import versions
from betse.util.type.text import regexes, strs
from betse.util.type.types import (
    type_check, MappingType, SequenceTypes, SetType, StrOrNoneTypes,)
from contextlib import contextmanager

# ....................{ GLOBALS ~ str                      }....................
_BACKEND_NAME_HEADLESS = 'Agg'
'''
Name of the non-GUI-based matplotlib backend to fallback to in the event that
*NO* GUI-based matplotlib backend is usable on this system. For portability,
this backend is guaranteed to be usable on all platforms and systems regardless
of matplotlib version.
'''

# ....................{ GLOBALS ~ dict                     }....................
_LOG_LEVEL_TO_VERBOSITY_LEVEL_NAME = {
    LogLevel.NONE:     'silent',
    LogLevel.CRITICAL: 'silent',
    LogLevel.ERROR:    'silent',
    LogLevel.WARNING:  'silent',
    LogLevel.INFO:     'helpful',
    LogLevel.DEBUG:    'debug',
    LogLevel.ALL:      'debug-annoying',
}
'''
Dictionary mapping from each standard logging level expected by newer matplotlib
versions (i.e., matplotlib >= 2.2.0) to the corresponding non-standard verbosity
level name expected by older matplotlib versions (i.e., matplotlib < 2.2.0).
'''


_RC_PARAMS = {
    #FIXME: This doesn't appear to do anything anymore. Excise, please. *sigh*

    # Unconditionally print terse messages. By default, matplotlib prints *NO*
    # messages. Valid values include: "silent", "helpful", "debug", and
    # "debug-annoying".
    # 'verbose.level': 'helpful',
    # 'verbose.level': 'debug',
}
'''
Dictionary mapping ``matplotlibrc`` option names to corresponding values.

The :meth:`MplConfig.init` method subsequently updates
the :attr:`matplotlib.rcParams` dictionary of default ``matplotlibrc`` options
deserialized from the current ``matplotlibrc`` file with this dictionary. Ergo,
the custom options specified by this dictionary override the default options
defined by that file.
'''

# ....................{ CLASSES                            }....................
class MplConfig(object):
    '''
    High-level wrapper simplifying low-level configuration and introspection of
    matplotlib.
    '''

    # ..................{ INITIALIZERS                       }..................
    def init(self, backend_name: StrOrNoneTypes) -> None:
        '''
        Reconfigure matplotlib with sane defaults specific to the current
        system.

        On first importation, matplotlib configures itself by loading the
        contents of the first ``matplotlibrc`` file found in any of several
        candidate directories. Technically, BETSE _could_ supply an
        application-specific version of this file to force matplotlib to adopt
        application-specific configuration settings. Since synchronizing this
        local copy with remote changes is an onerous (if not ultimately
        infeasible) chore, we elect instead to reconfigure matplotlib _after_
        this file has already been loaded at application startup. While this
        does increase startup costs, the alternatives are all absurd at best.

        Parameters
        ----------
        backend_name: optional[str]
            Name of the matplotlib backend to explicitly enable. Defaults to
            `None`, in which case this method implicitly enables the first
            importable backend known to be both usable and supported by
            application requirements (_in descending order of preference_).

        See Also
        ----------
        http://matplotlib.org/users/customizing.html
            Matplotlib configuration details.
        '''

        # Avoid circular import dependencies.
        from betse.lib.matplotlib import mplcolormap

        # Initialize matplotlib in a safe manner.
        self._init_matplotlib()

        # Establish the matplotlib backend *AFTER* initializing matplotlib.
        self._init_backend(backend_name=backend_name)

        # Register all custom colormaps *AFTER* initializing matplotlib.
        mplcolormap.init()

        # Import all animation writer classes *AFTER* establishing the
        # matplotlib backend, thus implicitly registering these classes with
        # matplotlib. This permits other code in the codebase to conveniently
        # refer to these classes by name without having to manually instantiate
        # these classes. (While it may or may not be necessary to import these
        # classes after establishing the backend, it is only prudent to do so.)
        from betse.lib.matplotlib.writer import mplcls
        if False: mplcls    # silence contemptible IDE warning messages


    def _init_matplotlib(self) -> None:
        '''
        Initialize matplotlib in a safe manner, guaranteeing all subsequent
        importations of matplotlib elsewhere in the codebase to be safe.

        This method imports matplotlib in a safe manner and should be the first
        call to do so for the active Python process. Unfortunately, the
        :mod:`matplotlib.__init__` submodule implicitly imported on the first
        matplotlib importation performs the following unsafe logic:

        * The current command-line argument list :data:`sys.argv` is
          iteratively searched for arguments matching patterns, including:
          * Arguments prefixed by `-d` of length greater than or equal to 3
            (e.g., `-dtkagg` but _not_ simply `-d`). For each such argument,
            :mod:`matplotlib.__init__` attempts to enable the backend whose
            name is given by such argument excluding the prefixing `-d`,
            silently ignoring exceptions.
          * Arguments prefixed by `--verbose-` matching a matplotlib-specific
            verbosity level (e.g., `--verbose-debug`). For each such argument,
            :mod:`matplotlib.__init__` coerces the global verbosity to that
            level.

        This is utterly horrible. Since enabling arbitrary backends can have
        non-negligible side effects, the :mod:`matplotlib.__init__` submodule
        *must* be prevented from performing this logic. Since this submodule is
        imported only on the first importation of a matplotlib module, it
        performing this preventation *only* on the first importation of the
        top-level :mod:`matplotlib` package by this method suffices to sanitize
        matplotlib behaviour.
        '''

        # Copy the current argument list into a temporary list.
        _sys_argv_old = sys.argv[:]

        # Remove all arguments following the mandatory basename of the current
        # process from the current argument list, preventing matplotlib from
        # inspecting arguments and hence enabling arbitrary backends.
        del(sys.argv[1:])

        # Import the "matplotlib.__init__" module. Since this module prints
        # verbose messages if and only if the matplotlib-specific CLI option
        # "--verbose-debug" is in this argument list, this option must be
        # explicitly added to this list *BEFORE* importing this module. Awful!
        try:
            # Logging configuration singleton.
            log_config = logconfig.get()

            # Logging level of the application-wide standard output handler,
            # representing the default logging level for this Python process.
            betse_log_level = log_config.handler_stdout.level

            # Logging level of matplotlib itself, applicable *ONLY* to
            # matplotlib >= 2.2.0. See comments below for contemptible details.
            matplotlib_log_level = None

            # Name of the matplotlib-specific verbosity level corresponding to
            # this numeric level to the corresponding, applicable *ONLY* to
            # matplotlib < 2.2.0.
            #
            # Since there exist fewer matplotlib-specific verbosity levels than
            # standard logging lovels, this conversion is necessarily lossy. A
            # hearty thanks for failing to conform to well-established Python
            # standards yet again, matplotlib!
            #
            # Note that matplotlib-specific verbosity levels have been
            # obsoleted as of matplotlib 2.2.0 by the standard "logging" API.
            # Ergo, all of the logic below pertaining to this variable is
            # required *ONLY* to ensure backward compatibility with older
            # matplotlib versions (i.e., matplotlib < 2.2.0).
            #
            # Ideally, this logic would be conditionally performed *ONLY* if the
            # current version of matplotlib is sufficiently old. Unfortunately,
            # attempting to do so raises chicken-and-egg issues; notably, the
            # current version of matplotlib can only be determined by importing
            # the top-level "matplotlib" module, but doing so requires *ALL*
            # matplotlib-specific CLI options (including those setting
            # verbosity) to have already been injected into the "sys.argv" list.
            # Since these two constraints are mutually exclusive *AND* since
            # squelching verbosity in older versions of matplotlib takes
            # precedence over doing so more intelligently, we do so
            # unintelligently (i.e., unconditionally) -- regardless of whether
            # the current version of matplotlib even supports these options.
            matplotlib_log_level_name = None

            # If BETSE-specific debugging is enabled, enable the mildest form
            # of matplotlib-specific logging. While matplotlib also supports
            # the "debug" and "debug-annoying" verbosity levels, both produce
            # far more effluvia than is generally useful.
            if betse_log_level <= LogLevel.DEBUG:
                matplotlib_log_level = LogLevel.INFO
                matplotlib_log_level_name = 'helpful'
            # Else, squelch all matplotlib-specific logging. Since matplotlib
            # supports no fine-grained verbosity levels between "helpful" and
            # "silent", reduce all non-debugging levels to merely "silent".
            else:
                matplotlib_log_level = LogLevel.WARNING
                matplotlib_log_level_name = 'silent'

            # Convert this name into a matplotlib-specific CLI option.
            # print('matplotlib verbosity: ' + matplotlib_log_level_name)
            # matplotlib_log_level_name='debug'
            sys.argv.append('--verbose-' + matplotlib_log_level_name)

            #FIXME: We should additionally set the "ffmpeg"-specific CLI option
            #"-loglevel" based on the above "betse_log_level" as well -- perhaps by
            #setting "rcParams['animation.ffmpeg_args'] = '-report'" or some
            #such *AFTER* importing matplotlib below.

            # Log this initialization.
            logs.log_debug('Initializing matplotlib with options: %s', sys.argv)

            # Set matplotlib's logging level *BEFORE* importing from matplotlib.
            # Since matplotlib now defaults to the sane WARNING level, this
            # could technically also be performed *AFTER* importing from
            # matplotlib; doing so would, however, ignore all messages logged
            # during this importation with a lower level. (That would be bad.)
            #
            # Note that this logic effectively reduces to a noop under older
            # matplotlib versions (i.e., matplotlib < 2.2.0). Although
            # non-ideal, this appears to have *NO* harmful side effects.
            #
            # Ideally, the following two lines would reduce to this one-liner:
            #
            #     self.log_level = matplotlib_log_level
            #
            # Unfortunately, the property setter implicitly invoked by such an
            # assignment assumes matplotlib to already have been initialized --
            # which, of course, it hasn't. (Chicken-and-egg issues, people.)
            matplotlib_logger = logs.get('matplotlib')
            matplotlib_logger.setLevel(matplotlib_log_level)

            # Import matplotlib submodules *AFTER* setting all matplotlib CLI
            # options, which this importation internally parses in older
            # versions of matplotlib as a non-idempotent side effect.
            import matplotlib      # this is horrible
            from matplotlib import rcParams  # , font_manager

            # If the "matplotlib.verbose" object exists in this version,
            # monkeypatch this object to prevent the
            # matplotlib.verbose.set_level() method from reducing to a noop, as
            # occurs when this private attribute is *NOT* nullified. waat?
            #
            # For reasons, matplotlib 2.2.0 critically broke backward
            # compatibility in an apocalyptic manner by removing this object.
            # For discussion, see the following matplotlib issue:
            #
            #     https://github.com/matplotlib/matplotlib/issues/10716
            if hasattr(matplotlib, 'verbose'):  # this is horrible, too!
                matplotlib.verbose._commandLineVerbose = None
        # Guarantee the prior argument list to be restored from this temporary
        # list even in the event of exceptions.
        finally:
            sys.argv = _sys_argv_old
            del(_sys_argv_old)

        # Unconditionally enable all settings defined by this global.
        rcParams.update(_RC_PARAMS)

        #FIXME: Sadly, the "font_manager.USE_FONTCONFIG" global is currently
        #only modifiable by physically modifying the contents of the
        #"matplotlib/font_manager.py" file at installation time. Yes, this is
        #patently insane -- but such are the vagaries of a pale life.

        # If the external "fc-match" command is in the current ${PATH}, the
        # POSIX-compliant fontconfig subsystem for portably locating system
        # fonts is available. Since fontconfig is the canonical means of
        # performing font lookups under most POSIX-compliant platforms (e.g.,
        # Linux, macOS), instruct matplotlib to do so. For unknown reasons,
        # matplotlib's fontconfig support is labelled "experimental" and hence
        # disabled by default. Since this support appears to be sufficiently
        # robust for our limited use case, manually enable this support.
        # if pathables.is_pathable('fc-match'):
        #     font_manager.USE_FONTCONFIG = True


    #FIXME: Revise docstring in accordance with the iterative platform-specific
    #fallback process now employed for deciding the default Matplotlib backend.
    def _init_backend(self, backend_name: StrOrNoneTypes) -> None:
        '''
        Set the default matplotlib backend to be implicitly used for subsequent
        plotting, which matplotlib requires to be configured _before_ the first
        importation of any the following modules: :mod:`matplotlib.pyplot` or
        :mod:`matplotlib.backends`.

        Specifically:

        * If this method is called by `py.test`-based testing and hence
          possibly by headless continuous integration (CI) with no access to a
          window manager, the non-interactive `Agg` backend is used.
        * Else, if the current platform is:
          * Either Linux or Windows, the `TkAgg` backend is used. This is the
            only backend currently known to survive freezing into executables.
            Alternatives include:
            * `Qt4Agg`, an (_arguably_) aesthetically inferior backend with
              (_inarguably_) significant performance concerns.
          * OS X, the `MacOSX` backend is used. This is the only backend
            currently known to survive freezing into executables. Alternatives
            include:
            * `CocoaAgg`, a non-native backend leveraging the cross-platform
              C++ library AGG (Anti-grain Geometry). Yes! This backend is
              officially deprecated and fundametally broken, however. No!

        Parameters
        ----------
        backend_name: optional[str]
            Name of the matplotlib backend to explicitly enable. Defaults to
            `None`, in which case this method implicitly enables the first
            importable backend known to be both usable and supported by
            application requirements (_in descending order of preference_).
        '''

        # If no specific backend is requested *AND* the active Python process is
        # headless and hence supports only headless backends...
        if backend_name is None and displays.is_headless():
            # Log this observation.
            logs.log_info(
                'Headless display environment detected. '
                'Defaulting to headless backend "%s"...',
                _BACKEND_NAME_HEADLESS)

            # Default to a headless backend.
            backend_name = _BACKEND_NAME_HEADLESS

        # If a specific backend is requested...
        if backend_name is not None:
            # If this backend is usable, enable this backend and return. Note
            # that, as the is_backend_usable() method internally enables this
            # backend as a necessary side effect of testing this backend's
            # usability, this backend need not (and, indeed, should not) be
            # explicitly re-enabled here.
            if self.is_backend_usable(backend_name):
                return
            # Else, log a non-fatal error and continue.
            else:
                logs.log_error('Preferred backend "%s" unusable.', backend_name)
        # Else, no specific backend is requested. Since the default headless
        # backend (e.g., "Agg") should *ALWAYS* be usable, this typically
        # implies this process to be headfull and hence support GUI backends.
        # Default to the first backend usable under the current system.

        # Name of the current platform (e.g., "Linux", "Darwin", "Windows").
        kernel_name = kernels.get_name()

        # Tuple of the names of all matplotlib backends to iteratively
        # fallback to on this platform if supported or None otherwise.
        backend_names = self._kernel_name_to_backend_names_prefer.get(
            kernel_name, None)

        # If this platform is unsupported, raise an exception.
        if backend_names is None:
            raise BetseMatplotlibException(
                'Platform "{}" unsupported.'.format(kernel_name))

        # Log this iteration.
        logs.log_debug(
            'Finding usable backend in: %r', backend_names)

        #FIXME: If this is Matplotlib >= 2.0.0 and the only available backend is
        #"TkAgg", a non-fatal warning should be logged instructing the user to
        #install a more stable backend (e.g., "Qt5Agg", "Qt4Agg", "WxAgg"). To
        #leverage the local "backend_names" list, this warning should probably
        #be emitted here rather than in the is_backend_usable() method.

        # For each such backend (in descending order of preference)...
        for backend_name in backend_names:
            # If this backend is usable, this tester has already implicitly
            # enabled this backend. Our work is done here.
            if self.is_backend_usable(backend_name):
                return
        # Else, no preferred GUI-based backend is usable on the current system
        # (e.g., due to no external widget library being installed).

        # If the fallback non-GUI-based backend is usable...
        if self.is_backend_usable(_BACKEND_NAME_HEADLESS):
            # Log a non-fatal warning.
            logs.log_warning(
                'No usable GUI-based matplotlib backend found. '
                'Defaulting to usable CLI-based backend "%s". '
                'Consider installing support for GUI-based backends %s.',
                _BACKEND_NAME_HEADLESS,
                strs.join_as_disconjunction_double_quoted(*backend_names))

            # Default to this backend and return.
            self.backend_name = _BACKEND_NAME_HEADLESS
            return

        # Else, no backends appear to be usable on the current system. Due to
        # the ubiquity of the headless fallback backend (e.g., "Agg"), this
        # should only occur in extremely uncommon edge cases (e.g., manual
        # recompilation of the entire system). If this does occur, this is
        # typically the result of internal issues in this codebase. To assist in
        # debugging these issues, a simple exception message is raised.
        raise BetseMatplotlibException(
            'No usable matplotlib backend found. '
            '{}-supported backends tested include (in order): {}.'.format(
                kernel_name,
                strs.join_as_conjunction_double_quoted(*backend_names)))

    # ..................{ PROPERTIES                         }..................
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # NOTE: To avoid desynchronization issues between low-level matplotlib
    # internals and the following high-level properties wrapping these
    # internals, most of these properties leverage the less efficient but safer
    # @property decorator rather than the more efficient but less safe
    # @property_cached decorator.
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    @property_cached
    def version(self) -> str:
        '''
        Currently installed version of matplotlib as a ``.``-delimited version
        specifier (e.g., ``2.1.0``).
        '''

        # Delay importation of the "matplotlib.__init__" module.
        import matplotlib

        # Currently installed version of matplotlib to be returned.
        version_return = matplotlib.__version__

        # Either:
        #
        # * If this version exhibits formatting discrepancies introduced in
        #   recent matplotlib versions violating Python conventions for version
        #   specifiers (e.g., the non-compliant version "2.1.0-python3_4" rather
        #   than the compliant version ""2.1.0"), a sequence whose single item
        #   is this version stripped of these discrepencies.
        # * Else, "None".
        #
        # Note that this regular expression need *NOT* be compiled for
        # efficiency, as this property is cached and hence called only once.
        version_munged_or_none = regexes.get_match_groups_numbered_or_none(
            text=version_return, regex=r'^([^-]+)-')

        # Else, this version is non-compliant. In this case, reduce this to a
        # compliant version.
        if version_munged_or_none is not None:
            logs.log_debug(
                'Sanitizing matplotlib version from "%s" to "%s"...',
                version_return, version_munged_or_none[0])
            version_return = version_munged_or_none[0]

        # Return this possibly munged version.
        return version_return

    # ..................{ PROPERTIES ~ log level             }..................
    @property
    def log_level(self) -> LogLevel:
        '''
        Matplotlib-specific logging level.

        If this is a sufficiently old version of matplotlib *not* supporting the
        :mod:`logging` API (i.e., matplotlib < 2.2.0), this method attempts to
        internally convert from a matplotlib-specific verbosity level name to a
        :mod:`logging` level. If this name is:

        * ``silent``, hiding all Matplotlib output *except* non-fatal warnings,
          :attr:`LogLevel.WARNING` is returned. Not exactly silent, is it?
        * ``helpful``, emitting terse Matplotlib debugging output,
          :attr:`LogLevel.INFO` is returned.
        * `debug`, emitting verbose Matplotlib debugging output,
          :attr:`LogLevel.DEBUG` is returned.
        * `debug-annoying`, emitting *very* verbose Matplotlib debugging output.
          :attr:`LogLevel.ALL` is returned. You *never* want this. Trust us.
        '''

        # Matplotlib-specific logger.
        #
        # Since this object unconditionally sets this logger's level regardless
        # of whether the current matplotlib version actually logs to this
        # logger, this logging level is returned as is. Ergo, the docstring
        # above is an innocent lie. (Ignore it, please.)
        matplotlib_logger = logs.get('matplotlib')

        # Return this logger's level.
        return matplotlib_logger.level


    @log_level.setter
    @type_check
    def log_level(self, log_level: LogLevel) -> None:
        '''
        Set the current matplotlib-specific logging level to the passed level.

        Parameters
        -----------
        log_level: LogLevel
            Matplotlib-specific logging level to be set.
        '''

        # Matplotlib-specific logger.
        #
        # Note that this logic effectively reduces to a noop under older
        # matplotlib versions (i.e., matplotlib < 2.2.0). Although
        # non-ideal, this appears to have *NO* harmful side effects.
        matplotlib_logger = logs.get('matplotlib')

        # Set this logger's level to the passed level.
        matplotlib_logger.setLevel(log_level)

        # Delay importation of the "matplotlib.__init__" module. To ensure
        # messages logged by this importation are logged with the passed level,
        # we do so *AFTER* setting this level.
        import matplotlib

        # If the "matplotlib.verbose" object exists in this version, internally
        # convert the passed level to an obsolete verbosity level name required
        # by this obsolete version of matplotlib. *SIGH*
        matplotlib_verbose = getattr(matplotlib, 'verbose', None)
        if matplotlib_verbose is not None:  # this is horrible
            # Non-standard verbosity level name corresponding to the passed
            # level expected by older matplotlib versions (i.e., matplotlib <
            # 2.2.0) *OR* "None" if no such name corresponds.
            verbosity_level_name = _LOG_LEVEL_TO_VERBOSITY_LEVEL_NAME.get(
                log_level, None)

            # If this name is unrecognized, raise an exception. For unknown
            # reasons, the Verbose.set_level() method called below unsafely
            # emits non-fatal warnings rather than raising fatal exceptions on
            # receiving an unrecognized level name.
            if (
                verbosity_level_name is None or
                verbosity_level_name not in matplotlib_verbose.levels
            ):
                raise BetseMatplotlibException(
                    'Matplotlib verbosity level "{}" unrecognized.'.format(
                        verbosity_level_name))

            # Set this verbosity level.
            matplotlib_verbose.set_level(verbosity_level_name)

    # ..................{ PROPERTIES ~ path                  }..................
    @property
    def cache_dirname(self) -> str:
        '''
        Absolute pathname of the platform- and typically user-specific directory
        to which matplotlib caches metadata (e.g., about fonts).
        '''

        # Delay importation of the "matplotlib.__init__" module.
        import matplotlib

        # Return this path.
        return matplotlib.get_cachedir()


    @property
    def rc_filename(self) -> str:
        '''
        Absolute pathname of the currently selected ``matplotlibrc`` file
        establishing default matplotlib options.
        '''

        # Delay importation of the "matplotlib.__init__" module.
        import matplotlib

        # Return this path.
        return matplotlib.matplotlib_fname()

    # ..................{ PROPERTIES ~ backend               }..................
    @property
    def backend(self) -> type(sys):
        '''
        In-memory submodule object corresponding to the current backend.

        Raises
        ----------
        BetseMatplotlibException
            If no backend has been set yet.
        '''

        # If no backend has been set yet, raise an exception.
        if not self.is_backend():
            raise BetseMatplotlibException('Matplotlib backend not yet set.')

        # Fully-qualified name of this backend's submodule.
        backend_module_name = (
            'matplotlib.backends.backend_' + self.backend_name)

        # This backend's submodule. Since this backend has been set, this
        # submodule *SHOULD* still be cached in-memory. Let's be sure.
        backend_module = sys.modules.get(backend_module_name, None)
        if backend_module is None:
            raise BetseMatplotlibException(
                'Current matplotlib backend module "{}" not imported.'.format(
                    backend_module_name))

        # Return this submodule.
        return backend_module


    @property
    def backend_canvas_class(self) -> type(object):
        '''
        `FigureCanvas` subclass corresponding to the current backend (e.g.,
        `FigureCanvasQt5Agg` for the `qt5agg` backend).

        Note that, whereas each backend implements a backend-specific canvas,
        each backend shares the same `Figure` class (i.e., there exist no
        corresponding backend-specific `Figure` subclasses).

        If no backend has been set yet, an exception is raised.
        '''

        # This backend's module.
        backend = self.backend

        # This subclass. While all standard backends *SHOULD* set the module
        # attribute "FigureCanvas" to the "FigureCanvas" subclass specific to
        # that backend, let's be sure.
        backend_canvas_class = getattr(backend, 'FigureCanvas', None)
        if backend_canvas_class is None:
            raise BetseMatplotlibException(
                'Matplotlib backend canvas class '
                '"{}.FigureCanvas" not found.'.format(backend.__name__))
        return backend_canvas_class


    @property
    def backend_figure_filetypes(self) -> list:
        '''
        List of all figure filetypes supported by the current backend (e.g.,
        `['bmp', 'gif', 'jpg', 'png', 'tiff']`).

        If no backend has been set yet, an exception is raised.
        '''

        # There exist two means of acquiring this metadata:
        #
        # 1. Creating an empty canvas specific to the current backend.
        # 2. Creating an empty figure specific to the current backend.
        #
        # Since the latter also implicitly requires creating an empty canvas
        # and is hence both more expensive and more complex, we prefer the
        # former.
        backend_canvas_class = self.backend_canvas_class

        # The "FigureCanvasBase" superclass of all canvas subclasses defines
        # the public class attribute "filetypes" to be the dictionary of all
        # filetypes supported by that canvas. Moreover, all subclasses reliably
        # redefine this attribute in the expected way. While the same list is
        # also obtainable by creating an instance of this subclass and calling
        # that subclass' get_supported_filetypes() method, the current approach
        # is substantially more efficient.
        #
        # Magic is magic. Do not question magic, for it is magical.
        return list(backend_canvas_class.filetypes.keys())

    # ..................{ PROPERTIES ~ backend : name        }..................
    @property
    def backend_name(self) -> str:
        '''
        Human-readable name (e.g., ``Qt5Agg``) of the current backend.

        This name is *not* guaranteed to be lowercase. This name is typically a
        mix of lower- and uppercase alphanumeric characters.

        This backend is *not* guaranteed to be usable. If no backend has been
        enabled yet (either explicitly or implicitly), this is the name of the
        default backend implicitly enabled by importing the
        :mod:`motplotlib.pyplot` submodule.
        '''

        # Delay importation of the "matplotlib.__init__" module.
        import matplotlib

        # Return this backend's name.
        return matplotlib.get_backend()


    @backend_name.setter
    @type_check
    def backend_name(self, backend_name: str) -> None:
        '''
        Set the current backend to the usable backend with the passed name.

        This name is interpreted case-insensitively and hence may be in any
        case including lower- and uppercase (e.g., `tkagg`, `TKAGG`, `TkAgg`).

        This backend is assumed to be usable (i.e., the
        :meth:`is_backend_usable` method is assumed to succeed when passed this
        name). If this is *not* the case, matplotlib will internally raise an
        exception of indeterminate origin and readability. Welcome to hell.
        '''

        # Log this attempt.
        logs.log_debug('Enabling matplotlib backend "%s"...', backend_name)

        # If enabling the "TkAgg" backend *AND* this is matplotlib >= 2.0.0,
        # log a prominent warning. This backend is known to behave unstably
        # under all newer versions of matplotlib. See the
        # "_kernel_name_to_backend_names_prefer" poperty for further details.
        if (
            backend_name == 'TkAgg' and
            versions.is_at_least(self.version, '2.0.0')
        ):
            logs.log_warning(
                'Matplotlib backend "TkAgg" known to be unstable '
                'under matplotlib >= 2.0.0.')

        # Attempt to enable this backend.
        try:
            self._enable_backend(backend_name)
        # If a presumably non-human-readable exception was raised...
        except Exception:
            # Log a human-readable error first.
            logs.log_error(
                'Matplotlib backend "%s" not found or unusable.', backend_name)

            # Re-raise this exception.
            raise

    # ..................{ PROPERTIES ~ backend : names       }..................
    @property_cached
    def backend_names(self) -> SequenceTypes:
        '''
        Sequence of the strictly lowercase names of all currently available
        matplotlib-specific backends (e.g., `['gtk3agg', 'tkagg', 'qt4agg']`).

        While matplotlib provides the canonical lists
        :attr:`matplotlib.rcsetup.interactive_bk`,
        :attr:`matplotlib.rcsetup.non_interactive_bk`, and
        :attr:`matplotlib.rcsetup.all_backends`, even the latter typically
        fails to list all possible backends (e.g., `mixed` tends to be
        missing). For completeness, this function instead iteratively inspects
        the current filesystem.

        For efficiency, this property is created and cached on the first access.
        '''

        # Importing such module has side effects and hence is deferred.
        from matplotlib import backends

        # Absolute path of the directory containing all backends for the
        # currently imported "matplotlib".
        backends_dir = modules.get_dirname(backends)

        # If this directory exists, find all backends in this directory.
        if dirs.is_dir(backends_dir):
            # String prefixing the basenames of backend-specific modules.
            BACKEND_BASENAME_PREFIX = 'backend_'

            # Return and cache a list of these names, discovered by:
            #
            # * Filtering all basenames in this directory for modules.
            # * Converting the remaining basenames to backend names.
            # * Sorting these names in ascending lexicographic order for
            #   readability (e.g., in the "info" subcommand).
            return iterables.sort_ascending([
                pathnames.get_pathname_sans_filetype(
                    strs.remove_prefix_if_found(
                        backend_basename, BACKEND_BASENAME_PREFIX))
                for backend_basename in dirs.iter_basenames(backends_dir)
                if strs.is_prefix(
                    backend_basename, BACKEND_BASENAME_PREFIX) and
                   pathnames.is_filetype_equals(backend_basename, 'py')
            ])
        # Else, this directory does *NOT* exist.
        else:
            # If the active Python interpreter is frozen, this is expected
            # and hence ignorable; else, this is unexpected, in which case a
            # non-fatal warning is logged and such list is cleared.
            if not pyfreeze.is_frozen():
                logs.log_warning(
                    'Directory "{}" not found. '
                    'Matplotlib backends not queryable.'.format(
                        backends_dir))

            # In either case, return and cache the empty list.
            return []

    # ..................{ PROPERTIES ~ backend : names : pri }..................
    @property_cached
    def _backend_names_blacklist(self) -> SetType:
        '''
        Unordered set of the names of all **blacklisted backends** (i.e.,
        backends known to be dangerously unusable on the current system).

        The :meth:`is_backend_usable` method directly reports these backends to
        be unusable *without* unsafely testing this to be the case. Attempting
        to programmatically test the unusability of these backends is known to
        induce extreme fatal errors in the worst case, including segnmentation
        faults with no human-readable exception messages.

        These backends include:

        * `Gtk3Agg`, which emits the following non-fatal warning when enabled:

            UserWarning: The Gtk3Agg backend is not known to work on Python
            3.x.

        * `Gtk3cairo`, which, despite claiming to be a GTK+ 3.x-specific
          backend, appears to attempt to dynamically load GTK+ 2.x-specific
          shared libraries -- inducing the fatal segmentation fault above.
        * All GTK+ 2.x-specific backends (e.g., `Gtk`, `Gtkagg`), conflicting
          with GTK+ 3.x-specific backends (e.g., `Gtk3`, `Gtk3agg`). Attempting
          to switch to the latter after having already switched to the former
          typically induces the following fatal segmentation fault immediately
          halting the current process:

            (betse:5000): Gtk-ERROR **: GTK+ 2.x symbols detected. Using GTK+
            2.x and GTK+ 3 in the same process is not supported
            zsh: trace trap (core dumped)  betse info

          This is *not* a high-level Python exception and hence cannot be
          caught from within Python. This is a POSIX-level process signal.

        Sometimes, the only winning move is not to play at all.
        '''

        # Blacklist:
        return {
            # All GTK+ 2.x-specific backends.
            'gtk', 'gtkagg', 'gtkcairo',

            # All GTK+ 2.x-specific backends.
            'gtk3', 'gtk3agg', 'gtk3cairo',
        }


    @property_cached
    def _kernel_name_to_backend_names_prefer(self) -> MappingType:
        '''
        Dictionary mapping from the name of each supported platform to a
        sequence of the names of all matplotlib backends that the
        :meth:`_init_backend` method iteratively defers to on that platform (in
        descending order of preference) in the event the end user fails to
        explicitly set this name (e.g., via the `--matplotlib-backend` option).

        In that case, that method defers to the first matplotlib backend usable
        on the current system whose name is in this tuple.

        For efficiency, this property is created and cached on the first access.
        '''

        # List of the names of all backends known to be unconditionally
        # supported by this application across *ALL* supported platforms.
        #
        # These are (in descending order of preference):
        #
        # * "Qt5Agg", a GUI backend with (arguably) superior aesthetics but
        #   (inarguably) significant performance concerns by compare to more
        #   preferable backends. Also, note that enabling the experimental
        #   non-blocking behaviour via "pyplot.show(block=False)" reliably
        #   induces segmentation faults with *NO* exception traceback under
        #   at least Windows. Something is better than nothing, though.
        # * "Qt4Agg", a GUI backend with (arguably) inferior aesthetics and
        #   (inarguably) significant performance concerns by compare to more
        #   preferable backends.
        # * "WxAgg", a GUI backend with (arguably) inferior aesthetics and
        #   (inarguably) significant performance *AND* stability concerns by
        #   compare to more preferable backends. As of this writing, the
        #   much-vaunted Python 3.x version of WxPython (i.e., Phoenix) has
        #   only recently been released and supported by matplotlib. Issues
        #   are likely to plague this backend for decades.
        backend_names_prefer = ['Qt5Agg', 'Qt4Agg', 'WxAgg',]

        # If this is matplotlib >= 2.0.0, deprioritize the "TkAgg" backend.
        #
        # For unknown reasons, all recent versions of matplotlib have
        # fundamentally broken this backend with respect to non-blocking
        # animations -- either with or without experimental non-blocking
        # behaviour. This has been extensively tested and isolated to matplotlib
        # itself, which... is frustrating. Let me tell you: we are displeased.
        if versions.is_at_least(self.version, '2.0.0'):
            # Append this backend to the end of this list.
            backend_names_prefer.append('TkAgg')
        # Else, prioritize the "TkAgg" backend. Why? Because:
        #
        # * "TkAgg", a GUI backend with adequate (albeit not particularly
        #   impressive) aesthetics and superior performance by compare to
        #   less preferable backends. (Tcl/Tk: who would've thought?)
        else:
            # Prepend this backend to the beginning of this list.
            backend_names_prefer.insert(0, 'TkAgg')

        # Darwin-specific list of such names, prioritizing the only truly usable
        # Darwin-specific backend before the platform-agnostic backends. Since
        # Darwin and Linux are both POSIX-compatible, cross-platform backends
        # (e.g., "Qt5Agg") tend to behave similarly under both platforms.
        backend_names_prefer_darwin = ['MacOSX',] + backend_names_prefer

        # Under Windows, prefer the exact same backends as preferred under
        # Linux. While POSIX-incompatible and hence irregular, Windows still
        # supports the same backends preferred under Linux in the same order.
        backend_names_prefer_windows = backend_names_prefer

        # Dictionary to be returned and cached.
        kernel_name_to_backend_names_prefer = {
            'Linux':   backend_names_prefer,
            'Darwin':  backend_names_prefer_darwin,
            'Windows': backend_names_prefer_windows,
        }

        # Return and cache this dictionary
        return kernel_name_to_backend_names_prefer

    # ..................{ TESTERS ~ backend                  }..................
    def is_backend(self) -> bool:
        '''
        ``True`` only if a backend has been successfully enabled.

        This function returns ``True`` only if a backend has been either:

        * Explicitly enabled by setting the :meth:`backend_name` property
          internally calling the :func:`matplotlib.use` function.
        * Implicitly enabled by importing either of the following submodules:
          * :mod:`matplotlib.pyplot`.
          * :mod:`matplotlib.backends`.
        '''

        # While ridiculous, this test corresponds exactly to the test performed
        # by the matplotlib.use() method itself to detect repetitious calls.
        # Since the "matplotlib.pyplot" subpackage internally imports this
        # subpackage, testing only this subpackage suffices.
        return modules.is_imported('matplotlib.backends')


    @type_check
    def is_backend_usable(self, backend_name: str) -> bool:
        '''
        ``True`` only if the backend with the passed name is usable on the
        current system.

        Specifically, this method returns ``True`` only if this backend:

        * May be switched to *without* raising exceptions.
        * May both create and destroy a hidden (i.e., invisible) empty figure
          *without* raising exceptions.

        If this backend is usable, this method switches the current backend to
        this backend *without* restoring the previously set backend. This is the
        unavoidable price of robust, reproducible test results. Callers
        requiring the previously set backend to be restored must do so manually
        (e.g., by setting the :func:`property` attribute to the name of that
        backend) *after* calling this method.
        '''

        # Log this attempt.
        logs.log_debug(
            'Testing matplotlib backend "%s" usability...', backend_name)

        # Lowercase this name, ensuring case-insensitive backend names.
        backend_name = backend_name.lower()

        #FIXME: Overkill and uninformative. To blacklist backends properly, what
        #we *REALLY* want to do instead is to:
        #
        #* Remove the "_backend_names_blacklist" property.
        #* Remove this conditional block entirely.
        #* Improve the _enable_backend() method to initially:
        #  * Define a local dictionary "BLACKLIST_BACKEND_NAME_TO_REASON"
        #    mapping from the names of all blacklisted backend names to
        #    human-readable exception messages describing exactly why those
        #    backends have been blacklisted.
        #  * Test whether the passed backend name is a key of this dictionary.
        #  * If so, raise an exception embedding the corresponding value of this
        #    dictionary.
        #
        #Yeah. That's blatantly heaps better. Hindsight is always 20-20, neh?

        # If this backend is unconditionally blacklisted, report this backend to
        # be unusable *WITHOUT* unsafely attempting to switch to this backend.
        if backend_name in self._backend_names_blacklist:
            logs.log_debug(
                'Matplotlib backend "%s" unusable, '
                'due to being internally blacklisted.', backend_name)
            return False

        # Attempt to...
        try:
            # Enable this backend.
            self._enable_backend(backend_name)

            # If doing so succeeds, this backend is usable.
            logs.log_debug('Matplotlib backend "%s" usable.', backend_name)
            return True
        # If doing so raises an exception, this backend is unusable.
        except Exception as exception:
            # Traceback for this exception.
            traceback = ioexceptions.get_traceback(exception)

            # Log this traceback for debuggability.
            logs.log_debug(
                'Matplotlib backend "%s" unusable, '
                'due to raising the following exception:\n%s',
                backend_name, traceback)

            # Report this backend to be unusable.
            return False

    # ..................{ TESTERS ~ block                    }..................
    def is_backend_current_nonblockable(self) -> bool:
        '''
        ``True`` only if the current backend supports true non-blocking display
        on the current system.

        See Also
        ----------
        :meth:`is_backend_nonblockable`
            Further details.
        '''

        return self.is_backend_nonblockable(self.backend_name)


    @type_check
    def is_backend_nonblockable(self, backend_name: str) -> bool:
        '''
        ``True`` only if the backend with the passed name supports **true
        non-blocking display** (i.e., in a backend-specific GUI eventloop in
        a separate thread of the current process) on the current system.

        If this method returns:

        * ``True``, this backend may reliably display plots and animations in a
          non-blocking manner by passing the experimental ``block=False``
          parameter to the :func:`matplotlib.pyplot.show` method.
        * ``False``, the experimental ``block=False`` parameter must *not* be
          passed to the :func:`matplotlib.pyplot.show` function. Doing so will
          either silently hide plots and animations with no warnings or errors
          *or* raise a low-level segmentation fault with no high-level exception
          traceback. Instead, consider falling back to the following approach:

          .. code:: python

              # Setup the current plot or animation for non-blocking display.
              matplotlib.interactive(True)
              matplotlib.pyplot.show()

              # Display the current frame of this plot or animation.
              matplotlib.pyplot.draw()
              matplotlib.pyplot.pause(0.0001)

              # Teardown the current plot or animation.
              matplotlib.interactive(False)

        Caveats
        ----------
        The fallback approach detailed above typically emits the following
        non-fatal warning on the first ``matplotlib.interactive(True)`` call:

        .. code::

           [betse] /usr/lib64/python3.4/site-packages/matplotlib/backend_bases.py:2437: MatplotlibDeprecationWarning: Using default event loop until function specific to this GUI is implemented
             warnings.warn(str, mplDeprecation)

        While verbose, this warning appears to be safely ignorable.
        '''

        # For safety, a whitelist approach is preferred. That is, only backends
        # explicitly known to support experimental "pyplot.show(block=False)"
        # behavior are assumed to do so; all other backends are assumed to *NOT*
        # support this behavior.
        return backend_name == 'TkAgg'

    # ..................{ GETTERS                            }..................
    def get_rc_param(self, param_name) -> object:
        '''
        Value of the parameter with the passed ``.``-delimited name (e.g.,
        ``savefile.dpi``) in the external ``matplotlibrc`` file.

        If no such parameter exists, an exception is raised.

        Raises
        ----------
        KeyError
            If no such parameter exists.
        '''

        # Delay importation of the "matplotlib.__init__" module.
        from matplotlib import rcParams

        # Return this parameter's value.
        return rcParams[param_name]


    def get_metadata(self) -> OrderedArgsDict:
        '''
        Ordered dictionary synopsizing the current matplotlib installation.
        '''

        # This dictionary.
        metadata = OrderedArgsDict(
            'rc file', self.rc_filename,
            'cache dir', self.cache_dirname,
            'current backend', self.backend_name,
        )

        # For each available backend, add metadata synopsizing that backend.
        for backend_name in self.backend_names:
            metadata[
                'backend {} usable'.format(backend_name.capitalize())] = str(
                    self.is_backend_usable(backend_name)).lower()

        # Get this dictionary.
        return metadata

    # ..................{ CONTEXTS                           }..................
    @contextmanager
    def reducing_log_level_to_debug_if_info(self):
        '''
        Context manager setting the matplotlib-specific verbosity level to
        :attr:`LogLevel.DEBUG` if currently :attr:`LogLevel.INFO` for the
        duration of this context.

        This context manager temporarily increases this level by one level.
        Although the :attr:`LogLevel.INFO` and :attr:`LogLevel.DEBUG` levels
        *both* produce debug output, only the latter produces debug output for
        external commands invoked by matplotlib (e.g., for encoding video via
        ``ffmpeg``); the former produces *no* such output. Since the latter is
        overly verbose for general use and hence useful only for specific cases,
        consider instead:

        * Defaulting to the :attr:`LogLevel.INFO` level.
        * Escalating to the :attr:`LogLevel.DEBUG` level by explicitly entering
          this context manager for the duration of special-case work requiring
          verbosity.

        This context manager guaranteeably reverts this level to the prior
        level even when fatal exceptions are raised. If this level is *not*
        currently :attr:`LogLevel.INFO`, this context manager is a noop.

        Returns
        -----------
        contextlib._GeneratorContextManager
            Context manager setting the level as described above.

        Yields
        -----------
        None
            Since this context manager yields no value, the caller's ``with``
            statement must be suffixed by *no* ``as`` clause.
        '''

        # If the current level is "helpful"...
        if self.log_level is LogLevel.INFO:
            # Escalate to the "debug" level temporarily.
            self.log_level = LogLevel.DEBUG

            # Yield control to the body of the caller's "with" block.
            try:
                yield
            # Revert to the prior level even if that block raised an exception.
            finally:
                self.log_level = LogLevel.INFO
        # Else, the current level is *NOT* "helpful". Reduce to a noop.
        else:
            yield

    # ..................{ MAKERS                             }..................
    #FIXME: Define a new make_backend_figure_manager() method as well. This
    #method is fine, for now. We might very well want to call this sometime.
    def make_backend_figure(self, *args, **kwargs):
        '''
        Create and return a new :class:`Figure` instance associated with a new
        :class:`FigureCanvasBase` instance corresponding to the current backend
        (e.g., :class:`FigureCanvasQt5Agg` for the `qt5agg` backend), passed
        the passed arguments.

        This figure will _not_ be associated with a new or existing
        :class:`FigureManagerBase` instance corresponding to the current
        backend (e.g., :class:`FigureManagerQt5Agg` for the `qt5agg` backend)
        _or_ with the :mod:`matplotlib.pyplot` interface. In theory, this
        implies that this figure should be implicitly garbage collected on:

        * If either _not_ displayed or displayed in a blocking manner, leaving
          scope.
        * If displayed in a non-blocking manner, this figure's window being
          closed.

        If no backend has been set yet, an exception is raised.

        Returns
        ----------
        matplotlib.figure.Figure
            Such figure.
        '''

        # Delay importation of the "matplotlib.__init__" module.
        from matplotlib.figure import Figure

        # Create this figure.
        figure = Figure(*args, **kwargs)
        self.backend_figure_canvas(figure)

        # Return this figure.
        return figure

    # ..................{ PRIVATE                            }..................
    @type_check
    def _enable_backend(self, backend_name: str) -> None:
        '''
        Set the current backend to the backend with the passed name.

        This name is interpreted case-insensitively and hence may be in any
        case including lower- and uppercase (e.g., `tkagg`, `TKAGG`, `TkAgg`).

        This low-level method is principally intended to be called by
        higher-level methods (e.g., :meth:`backend_name`,
        :meth:`is_backend_usable`) wrapping this method with additional
        caller-friendly logic.
        '''

        # Delay importation of the "matplotlib.__init__" module.
        import matplotlib

        # Lowercase this name, as backend names are case-insensitive.
        backend_name = backend_name.lower()

        # True only if a backend has already been enabled.
        is_backend = self.is_backend()

        # Wrap all attempts to import from unsafe matplotlib subpackages (e.g.,
        # "matplotlib.backends", "matplotlib.pyplot") with logic safely undoing
        # these imports when *NO* backend has already been set. (See below.)
        try:
            # If no backend has been enabled...
            if not is_backend:
                # Log this attempt.
                logs.log_debug(
                    'Enabling matplotlib backend "%s" via use()...',
                    backend_name)

                # Enable this backend by calling the use() function, which is
                # expected to behave in a stable manner.
                matplotlib.use(backend_name)
            # Else, a backend has already been enabled.
            #
            # If this backend is tbe currently enabled backend, noop. This is
            # substantially safer than calling the switch_backend() function in
            # the subsequent conditional branch, which is known to be unstable.
            elif backend_name == self.backend_name:
                # Log this attempt.
                logs.log_debug(
                    'Ignoring already enabled matplotlib backend "%s".',
                    backend_name)

                # Reduce to a noop.
                return
            # In this unfortunate case, we have no recourse but to call the
            # switch_backend() function known to be unstable.
            else:
                # Log this attempt.
                logs.log_debug(
                    'Enabling matplotlib backend "%s" via switch_backend()...',
                    backend_name)

                # Unfortunately, if the current platform is macOS *AND* the new
                # backend to be enabled is "TkAgg", enabling this backend is *NOT*
                # safe and must absolutely be prohibited by raising an exception.
                # Attempting to enable this backend under this edge case commonly
                # results in a segmentation fault, terminating the active Python
                # process in a non-human-readable manner resembling:
                #
                #     [betse] Testing matplotlib backend "tkagg"...
                #     backend TkAgg version 8.5
                #     2016-12-31 01:53:08.886 Python[19521:163945] -[NSApplication _setup:]: unrecognized selector sent to instance 0x7fc278fbec60
                #     2016-12-31 01:53:08.890 Python[19521:163945] An uncaught exception was raised
                #     2016-12-31 01:53:08.890 Python[19521:163945] -[NSApplication _setup:]: unrecognized selector sent to instance 0x7fc278fbec60
                #     2016-12-31 01:53:08.891 Python[19521:163945] (
                #              0   CoreFoundation                      0x00007fff936ea452 __exceptionPreprocess + 178
                #              ...
                #              97  ???                                 0x0000000000000004 0x0 + 4
                #     )
                #     2016-12-31 01:53:08.892 Python[19521:163945] *** Terminating app due to uncaught exception 'NSInvalidArgumentException', reason: '-[NSApplication _setup:]: unrecognized selector sent to instance 0x7fc278fbec60'
                #     *** First throw call stack:
                #     (
                #             0   CoreFoundation 0x00007fff936ea452 __exceptionPreprocess + 178
                #             ...
                #             97  ???  0x0000000000000004 0x0 + 4
                #     )
                #     libc++abi.dylib: terminating with uncaught exception of type NSException
                #     Abort trap: 6
                #
                # This is a common issue afflicting numerous matplotlib users under
                # all versions of macOS. The core issue appears to be that the
                # Xcode-bundled installation of Tcl/Tk is *NOT* usable by Python.
                # There exist two solutions (in no particular order):
                #
                # 1. Reinstall Python to use a Homebrew- or MacPorts-compiled
                #    installation of Tcl/Tk instead. While this does constitute a
                #    valid solution for end users, BETSE itself has no means of
                #    enforcing this dictate and *MUST* thus assume the current
                #    installation of Tcl/Tk to be the Xcode-bundled version.
                # 2. Call "matplotlib.pyplot.use('TkAgg')" *BEFORE* the first
                #    import of the "matplotlib.pyplot" submodule. While BETSE
                #    itself could technically attempt to enforce this by
                #    preferentially detecting the "TkAgg" backend *BEFORE* all
                #    other backends, the "MacOS" backend is always preferable under
                #    macOS and should thus always be detected first.
                #
                # In short, no sane solution exists. The only sane solution is to
                # refuse to play the game at all.
                if backend_name == 'tkagg' and oses.is_macos():
                    raise BetseMatplotlibException(
                        'Matplotlib backend "TkAgg" not '
                        'safely switchable to under macOS.')

                # Delay importation of this submodule until *ABSOLUTELY* necessary.
                # Importing this submodule implicitly enables the default matplotlib
                # backend defined by the "backend" RC parameter in the current
                # matplotlib configuration if no backend has been enabled, which
                # can exhibit harmful side effects in common edge cases.
                from matplotlib import pyplot

                # Switch from the current to the passed backend.
                pyplot.switch_backend(backend_name)

            # In either case, this backend has now been set. Since this does
            # *NOT* imply this backend to be usable, further testing is needed.

            # Delay importation of this submodule until *ABSOLUTELY* necessary.
            # See pertinent commentary above.
            from matplotlib import pyplot

            # Validate the usability of this backend by attempting to create and
            # destroy a hidden empty figure. Technically, doing so could exhibit
            # harmful side effects on uncooperative platforms (e.g., Windows) in
            # common edge cases.
            #
            # Unfortunately, doing so is also essential. The success of
            # switching to this backend above is a necessary but insufficient
            # condition of this backend's usability. While the success of
            # switching to some backends (e.g., "TkAgg") does reliably imply
            # those backends to be usable, the success of switching to other
            # backends (e.g., "Qt5Agg") only implies that the corresponding
            # packages (e.g., "PyQt5") are successfully importable; it does
            # *NOT* imply these packages and thus these backends to actually be
            # usable in any sense. Notably, the PyQt family of backends are
            # infamous for raising non-human-readable exceptions on attempting
            # to display the first figure. "Why, PyQt? WHY!?!?"
            pyplot.figure()
            pyplot.close()
        # If an exception is raised, this backend is unusable,
        except:
            # If no backend was set before the current call to this method,
            # unimport all top-level matplotlib subpackages internally
            # imported by the above call. Doing so ensures that the next
            # call to the is_backend() method reports False and hence that
            # this same conditional branch is re-entered on the next call to
            # this method.
            #
            # This is *NOT* merely a caller convenience. This is critical
            # functionality required for sane testing of backend usability.
            # In particular, if this is *NOT* done, the next call to this
            # method is guaranteed to raise a non-human-readable exception.
            # Why? Because, in that case:
            #
            # * The modules unimported here will still have been imported.
            # * Importing this module sets the current matplotlib backend,
            #   either to the backend requested by the first call to the
            #   matplotlib.use() function if called or to the default
            #   matplotlib backend set by "matplotlibrc" otherwise. In
            #   either case, a backend has now been set.
            # * Since that call to the matplotlib.use() function raised an
            #   exception, that backend is unusable. If that backend is the
            #   default matplotlib backend (e.g., "Qt5Agg"), then the
            #   default matplotlib backend is also unusable.
            # * The "matplotlib.pyplot" subpackage uses the current
            #   matplotlib backend if the "matplotlib.backends" subpackage
            #   has been imported.
            # * Since the latter has indeed been imported, the subsequent
            #   attempt to import the former will re-raise the same
            #   exception as initially raised here.
            # * Ergo, this method may only be safely called once.
            #
            # Since the above behaviour is insane, this subpackage is
            # unimported instead to sanitize life.
            if not is_backend:
                modules.unimport_module_if_imported(
                    'matplotlib.backends', 'matplotlib.pyplot')

            # Re-raise this exception.
            raise

# ....................{ SINGLETONS                         }....................
mpl_config = MplConfig()
'''
Singleton matplotlib configuration wrapper.
'''
