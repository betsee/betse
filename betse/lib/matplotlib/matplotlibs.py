#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

#FIXME: Redirect Matplotlib output through our logging interface. Currently, we
#configure Matplotlib via the "RC_PARAMS" global constant to simply output all
#non-debug messages to stdout, thereby preventing these messages from being
#logged. While Matplotlib can be reconfigured to redirect output to a file
#instead, doing so then prevents that output from also being redirected to
#stdout. (See the matplotlib.Verbose.set_fileo() method for details.)
#
#Ultimately, this is what needs to happen. Either:
#
#* Monkey-patch the matplotlib.verbose.report() method to leverage the
#  standard "logging" framework rather than simply printing to the currently
#  configured "fileo" object. Note that "matplotlib.verbose" is an instance of
#  the "matplotlib.Verbose" class. Sensible, that.
#* The non-monkey-patch approach:
#  1. Define a new BETSE-specific "LoggingVerbose" subclass of the
#     "matplotlib.Verbose" class. This subclass should leverage the standard
#     "logging" framework (e.g., by redefining the Verbose.report() method).
#  2. Assign "matplotlib.verbose = LoggingVerbose()".
#
#Actually, I'd be quite satisfied by the latter approach. Let's do this right.
#That should allow us to ignore all "verbose.*" rc parameters (e.g.,
#"verbose.level"), simplifying logic below. We probably still want the same
#"RC_PARAMS" global dictionary for use with future rc parameters -- but at least
#we'd be able to eliminate this bit of "verbose.level" hackiness.

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

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To permit matplotlib's default verbosity and backend to be replaced
# by BETSE-specific values, no matplotlib package or module may be imported
# until *AFTER* the MplConfig.init() method has been called -- including
# here at the top-level.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import sys
from betse.exceptions import BetseMatplotlibException
from betse.util.io.log import logconfig, logs
from betse.util.io.log.logenum import LogLevel
from betse.util.os import displays, kernels, oses
from betse.util.path import dirs, paths
from betse.util.py import freezers
from betse.util.type import iterables, regexes, strs, modules
from betse.util.type.call.memoizers import property_cached
from betse.util.type.mappings import OrderedArgsDict
from betse.util.type.types import type_check, StrOrNoneTypes
from contextlib import contextmanager

# ....................{ CONSTANTS                          }....................
RC_PARAMS = {
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

# ....................{ CONSTANTS ~ private                }....................
_BACKEND_NAME_HEADLESS = 'Agg'
'''
Name of the non-GUI-based matplotlib backend to fallback to in the event that
*NO* GUI-based matplotlib backend is usable on this system. For portability,
this backend is guaranteed to be usable on all platforms and systems regardless
of matplotlib version.
'''

# ....................{ CONSTANTS ~ private : prefer       }....................
_KERNEL_NAME_TO_BACKEND_NAMES_PREFERRED = {
    #FIXME: Inject the "WxAgg" backend somewhere into this list after
    #shown to be working under Python 3.x. Sadly, since a stable
    #version of WxPython Phoenix has yet to be released, this may take
    #considerably longer than first assumed. The vapourware: it burns!

    # Under Linux, the following backends are preferred:
    #
    # * "TkAgg", a GUI backend with adequate (albeit not particularly
    #   impressive) aesthetics and superior performance by compare to
    #   less preferable backends. Tcl/Tk: who would have ever thought?
    # * "Qt4Agg", a GUI backend with (arguably) inferior aesthetics and
    #   (inarguably) significant performance concerns by compare to
    #   more preferable backends. Unlike all less preferable backends (e.g.,
    #   "Qt5Agg"), however, this backend has no remaining stability concerns.
    #   Since stability trumps aesthetics and performance, this backend is
    #   preferable to all remaining backends. It could be worse.
    # * "Qt5Agg", a GUI backend with (arguably) superior aesthetics but
    #   (inarguably) significant performance *AND* stability concerns by compare
    #   to more preferable backends. In particular, enabling the experimental
    #   non-blocking behaviour with "pyplot.show(block=False)" reliably induces
    #   low-level segmentation faults with no high-level exception traceback
    #   under at least Windows. Something is better than nothing, though.
    'Linux': ('TkAgg', 'Qt4Agg', 'Qt5Agg',),

    # Preferred backends for the following platforms reuse the the
    # preferred backends for Linux defined above.
    'Darwin': None,
    'Windows': None,
}
'''
Dictionary mapping from the name of each supported platform to a tuple of the
names of all matplotlib backends to iteratively fallback to on that platform
(in descending order of preference) in the event the the end user fails to
explicitly specify such a name (e.g., via the `--matplotlib-backend` option).
In this case, this method subsequently falls back to the first matplotlib
backend usable on the current system whose name is in this tuple.
'''

# Under OS X and iOS, prefer the only genuinely usable Darwin-specific
# backend before the same backends as preferred under Linux. Since
# Darwin and Linux are both POSIX-compatible, cross-platform backends
# (e.g., "Qt5Agg") tend to behave similarly under both platforms.
_KERNEL_NAME_TO_BACKEND_NAMES_PREFERRED['Darwin'] = (
    ('MacOSX',) + _KERNEL_NAME_TO_BACKEND_NAMES_PREFERRED['Linux'])

# Under Windows, prefer the exact same backends as preferred under
# Linux. While POSIX-incompatible and hence irregular, Windows still
# supports the same backends preferred under Linux in the same order.
_KERNEL_NAME_TO_BACKEND_NAMES_PREFERRED['Windows'] = (
    _KERNEL_NAME_TO_BACKEND_NAMES_PREFERRED['Linux'])

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

        # Initialize matplotlib in a safe manner.
        self._init_matplotlib()

        # Set the default matplotlib backend *AFTER* initializing matplotlib.
        self._init_backend(backend_name=backend_name)

        # Import all animation writer classes *AFTER* establishing the backend,
        # thus implicitly registering these classes with matplotlib. This, in
        # turn, permits other portions of the codebase to conveniently refer to
        # these classes by name without having to manually instantiate them.
        # (While it may or may not be necessary to import these classes after
        # establishing the backend, it only seems prudent to do so.)
        from betse.lib.matplotlib.writer import mplclass
        if False: mplclass    # silence contemptible IDE warning messages


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

            # Logging level of the stdout handler, representing the default
            # logging level for the active Python process.
            log_level = log_config.handler_stdout.level

            # Convert this numeric level to the corresponding name of the
            # matplotlib-specific verbosity level. Since there exist fewer of
            # the latter than the former, this conversion is necessarily lossy.
            # A hearty thanks for failing to conform to well-established Python
            # standards yet again, matplotlib.
            verbosity_level_name = None

            # If BETSE-specific debugging is enabled, enable the mildest form
            # of matplotlib-specific logging. While matplotlib also supports
            # the "debug" and "debug-annoying" verbosity levels, both produce
            # far more effluvia than is helpful.
            if log_level <= LogLevel.DEBUG:
                verbosity_level_name = 'helpful'
            # Else, squelch all matplotlib-specific logging. Since matplotlib
            # supports no fine-grained verbosity levels between "helpful" and
            # "silent", reduce all non-debugging levels to merely "silent".
            else:
                verbosity_level_name = 'silent'

            # Convert this name into a matplotlib-specific CLI option.
            # print('matplotlib verbosity: ' + verbosity_level_name)
            # verbosity_level_name='debug'
            sys.argv.append('--verbose-' + verbosity_level_name)

            #FIXME: We should additionally set the "ffmpeg"-specific CLI option
            #"-loglevel" based on the above "log_level" as well -- perhaps by
            #setting "rcParams['animation.ffmpeg_args'] = '-report'" or some
            #such *AFTER* importing matplotlib below.

            # Log this initialization.
            logs.log_debug('Initializing matplotlib with options: %s', sys.argv)

            # Import matplotlib submodules *AFTER* setting all matplotlib CLI
            # options above, which this importation parses.
            from matplotlib import rcParams, verbose  # , font_manager

            # Prevent the verbose.set_level() method from reducing to a noop,
            # as occurs when this private attribute is *NOT* nullified. waat?
            verbose._commandLineVerbose = None     # yes, this is horrible
        # Guarantee the prior argument list to be restored from this temporary
        # list even in the event of exceptions.
        finally:
            sys.argv = _sys_argv_old
            del(_sys_argv_old)

        # Unconditionally enable settings defined by the "RC_PARAMS" global.
        rcParams.update(RC_PARAMS)

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

        # If a specific backend is requested, enable this backend and return.
        if backend_name is not None:
            self.backend_name = backend_name
            return
        # Else, no specific backend is requested.

        # If this process is headless and hence supports only CLIs...
        if displays.is_headless():
            # Log this observation.
            logs.log_info(
                'Headless display environment detected. '
                'Defaulting to non-GUI backend "%s". ',
                _BACKEND_NAME_HEADLESS)

            # Default to a headless backend and return.
            self.backend_name = _BACKEND_NAME_HEADLESS
            return
        # Else, this process is headfull and hence supports GUIs. Default to
        # the first backend usable under the current system.

        # Name of the current platform (e.g., "Linux", "Darwin", "Windows").
        kernel_name = kernels.get_name()

        # Tuple of the names of all matplotlib backends to iteratively
        # fallback to on this platform if supported or None otherwise.
        backend_names = _KERNEL_NAME_TO_BACKEND_NAMES_PREFERRED.get(
            kernel_name, None)

        # If this platform is unsupported, raise an exception.
        if backend_names is None:
            raise BetseMatplotlibException(
                'Platform "{}" unsupported.'.format(kernel_name))

        # For each such backend (in descending order of preference)...
        for backend_name in backend_names:
            # If this backend is usable, this tester has already implicitly
            # enabled this backend. Our work is done here.
            if self.is_backend_usable(backend_name):
                break
        # Else, no preferred GUI-based backend is usable on the current system
        # (e.g., due to no external widget library being installed).
        else:
            # If the fallback non-GUI-based backend is usable, log a non-fatal
            # warning and default to this backend.
            if self.is_backend_usable(_BACKEND_NAME_HEADLESS):
                logs.log_warning(
                    'Usable matplotlib GUI backend not found. '
                    'Defaulting to non-GUI backend "%s". '
                    'Consider installing support for GUI backends %s.',
                    _BACKEND_NAME_HEADLESS,
                    strs.join_as_disconjunction_double_quoted(*backend_names),
                )
            # Else, no backends appear to be usable on the current system. For
            # safety, raise an exception. Due to the ubiquity of the fallback
            # backend, this should *NEVER* happen. Damn you, Murphy!
            else:
                raise BetseMatplotlibException(
                    'Usable matplotlib backend not found.')

    # ..................{ PROPERTIES ~ path                  }..................
    @property
    def cache_dirname(self) -> str:
        '''
        Absolute path of the platform- (and typically user-) specific directory
        to which matplotlib caches metadata (e.g., on fonts).
        '''

        # Delay importation of the "matplotlib.__init__" module.
        import matplotlib

        # Return this path.
        return matplotlib.get_cachedir()


    @property
    def rc_filename(self) -> str:
        '''
        Absolute path of the currently selected ``matplotlibrc`` file
        establishing default matplotlib options.
        '''

        # Delay importation of the "matplotlib.__init__" module.
        import matplotlib

        # Return this path.
        return matplotlib.matplotlib_fname()

    # ..................{ PROPERTIES ~ verbosity             }..................
    @property
    def verbosity_level_name(self) -> str:
        '''
        Name of the current matplotlib-specific verbosity level.

        This name is guaranteed to be one of the following human-readable
        lowercase alphabetic words (in order of increasing verbosity):

        * `silent`, hiding all Matplotlib output _except_ non-fatal warnings.
        * `helpful`, emitting terse Matplotlib debugging output.
        * `debug`, emitting verbose Matplotlib debugging output.
        * `debug-annoying`, emitting _very_ verbose Matplotlib debugging output.
          (You never want this.)
        '''

        # Delay importation of the "matplotlib.__init__" module.
        from matplotlib import verbose

        # Return this name.
        return verbose.level


    @verbosity_level_name.setter
    @type_check
    def verbosity_level_name(self, verbosity_level_name: str) -> None:
        '''
        Set the current matplotlib-specific verbosity level to the level with
        the passed name.

        If this name is _not_ a recognized level name, an exception is raised.

        Parameters
        -----------
        verbosity_level_name : str
            Name of this level.

        See Also
        -----------
        verbosity_level_name
            This property's corresponding getter documents the set of all
            recognized level names.
        '''

        # Delay importation of the "matplotlib.__init__" module.
        from matplotlib import verbose

        # If this name is unrecognized, raise an exception. For unknown
        # reasons, the Verbose.set_level() method called below unsafely emits
        # non-fatal warnings rather than raising fatal exceptions on receiving
        # an unrecognized level name.
        if verbosity_level_name not in verbose.levels:
            raise BetseMatplotlibException(
                'Matplotlib verbosity level "{}" unrecognized.'.format(
                    verbosity_level_name))

        # Set this verbosity level.
        verbose.set_level(verbosity_level_name)

    # ..................{ PROPERTIES ~ backend               }..................
    @property
    def backend(self) -> type(sys):
        '''
        In-memory module object corresponding to the current backend.

        If no backend has been set yet, an exception is raised.
        '''

        # If no backend has been set yet, raise an exception.
        if not self.is_backend():
            raise BetseMatplotlibException('Matplotlib backend undefined.')

        # Name of this backend's module.
        backend_module_name = (
            'matplotlib.backends.backend_' + self.backend_name)

        # This backend's module. Since this backend has been set, this module
        # *SHOULD* still be cached in-memory. Let's be sure.
        backend_module = sys.modules.get(backend_module_name, None)
        if backend_module is None:
            raise BetseMatplotlibException(
                'Matplotlib backend module "{}" not found.'.format(
                    backend_module_name))

        # Return this module.
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

    # ..................{ PROPERTIES ~ backend : names       }..................
    @property_cached
    def backend_names(self) -> list:
        '''
        List of the strictly lowercase names of all currently available
        matplotlib-specific backends (e.g., `['gtk3agg', 'tkagg', 'qt4agg']`).

        While matplotlib provides the canonical lists
        :attr:`matplotlib.rcsetup.interactive_bk`,
        :attr:`matplotlib.rcsetup.non_interactive_bk`, and
        :attr:`matplotlib.rcsetup.all_backends`, even the latter typically
        fails to list all possible backends (e.g., `mixed` tends to be
        missing). For completeness, this function instead iteratively inspects
        the current filesystem.
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
                paths.get_pathname_sans_filetype(
                    strs.remove_prefix_if_found(
                        backend_basename, BACKEND_BASENAME_PREFIX))
                for backend_basename in dirs.list_basenames(backends_dir)
                if strs.is_prefix(
                    backend_basename, BACKEND_BASENAME_PREFIX) and
                    paths.is_filetype_equals(backend_basename, 'py')
            ])
        # Else, this directory does *NOT* exist.
        else:
            # If the active Python interpreter is frozen, this is expected
            # and hence ignorable; else, this is unexpected, in which case a
            # non-fatal warning is logged and such list is cleared.
            if not freezers.is_frozen():
                logs.log_warning(
                    'Directory "{}" not found. '
                    'Matplotlib backends not queryable.'.format(
                        backends_dir))

            # In either case, return and cache the empty list.
            return []

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
        Set the current backend to the backend with the passed name.

        This name is interpreted case-insensitively and hence may be in any
        case including lower- and uppercase (e.g., `tkagg`, `TKAGG`, `TkAgg`).
        '''

        # Log this attempt.
        logs.log_debug('Enabling matplotlib backend "%s"...', backend_name)

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

    # ..................{ TESTERS                            }..................
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
        return 'matplotlib.backends' in sys.modules


    @type_check
    def is_backend_usable(self, backend_name: str) -> bool:
        '''
        ``True`` only if the backend with the passed name is **usable** (i.e.,
        safely switchable to without raising exceptions) on the current system.

        If this backend is usable, this method switches the current backend
        to this backend *without* restoring the previously set backend. This is
        the unavoidable price of robust, reproducible test results. Callers
        requiring the previously set backend to be restored must do so manually
        (e.g., by setting the :func:`property` attribute to the name of that
        backend) *after* calling this method.
        '''

        # Log this attempt.
        logs.log_debug('Testing matplotlib backend "%s"...', backend_name)

        # Lowercase this name, ensuring case-insensitive backend names.
        backend_name = backend_name.lower()

        # Unconditionally report the following backends to be unusable without
        # actually testing these backends:
        #
        # * "Gtk3Agg", emitting the following non-fatal warning when enabled:
        #     UserWarning: The Gtk3Agg backend is not known to work on Python
        #     3.x.
        # * All GTK+ 2.x-specific backends (e.g., "Gtk", "Gtkagg"), conflicting
        #   with GTK+ 3.x-specific backends (e.g., "Gtk3", "Gtk3agg").
        #   Specifically, attempting to switch to the latter after having
        #   already switched to the former typically induces the following
        #   fatal segmentation fault immediately halting the current process:
        #
        #     (betse:5000): Gtk-ERROR **: GTK+ 2.x symbols detected. Using GTK+
        #     2.x and GTK+ 3 in the same process is not supported
        #     zsh: trace trap (core dumped)  betse info
        #
        #   This is *NOT* a high-level Python exception exception and hence
        #   cannot be caught from within Python. This is a low-level process
        #   signal. The only winning move is not to play at all.
        # * "Gtk3cairo", which, despite claiming to be a GTK+ 3.x-specific
        #   backend, appears to attempt to dynamically load GTK+ 2.x-specific
        #   shared libraries -- inducing the fatal segmentation fault above.

        #FIXME: For efficiency, compile this regex for reuse during iteration.

        if regexes.is_match(
            text=backend_name, regex=r'^gtk(?:3(?:agg|cairo)|[^3]|$)'):
            # print("bad gtk backend: " + backend_name)
            return False

        # Attempt to...
        try:
            # Enable this backend.
            self._enable_backend(backend_name)

            # If doing so succeeds, this backend is usable.
            return True
        # Else, this backend is unusable.
        except:
            return False

    # ..................{ TESTERS                            }..................
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
    def verbosity_debug_if_helpful(self):
        '''
        Context manager setting the matplotlib-specific verbosity level to
        `debug` if currently `helpful` for the duration of this context.

        This context manager temporarily increases this level by one level.
        Although the `helpful` and `debug` levels _both_ produce debug output,
        only the latter produces debug output for external commands invoked by
        matplotlib (e.g., for encoding video via `ffmpeg`); the former produces
        _no_ such output. However, the latter is overly verbose for general use
        and hence useful only for specific cases. Consider instead:

        * Defaulting to the `helpful` level.
        * Escalating to the `debug` level by explicitly entering this context
          manager for the duration of special-case work requiring verbosity.

        This context manager guaranteeably reverts this level to the prior
        level even when fatal exceptions are raised. If this level is _not_
        currently `helpful`, this context manager is a noop.

        Returns
        -----------
        contextlib._GeneratorContextManager
            Context manager setting the level as described above.

        Yields
        -----------
        None
            Since this context manager yields no value, the caller's `with`
            statement must be suffixed by _no_ `as` clause.
        '''

        # If the current level is "helpful"...
        if self.verbosity_level_name == 'helpful':
            # Escalate to the "debug" level temporarily.
            self.verbosity_level_name = 'debug'

            # Yield control to the body of the caller's "with" block.
            try:
                yield
            # Revert to the prior level even if that block raised an exception.
            finally:
                self.verbosity_level_name = 'helpful'
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
        high-level methods (e.g., :meth:`backend_name`,
        :meth:`is_backend_usable`) wrapping this method with additional
        caller-friendly logic.
        '''

        # Delay importation of the "matplotlib.__init__" module.
        import matplotlib

        # Lowercase this name, as backend names are case-insensitive.
        backend_name = backend_name.lower()

        # If no backend has been enabled, enable this backend by calling
        # the known-to-be-stable use() function.
        if not self.is_backend():
            matplotlib.use(backend_name)
        # Else, a backend has already been enabled. In this unfortunate case,
        # we have no recourse but to call the known-to-be-unstable
        # switch_backend() function.
        else:
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
                    'Matplotlib backend "TkAgg" unsafe on macOS.')

            # Delay importation of this submodule until *ABSOLUTELY*
            # necessary. Importing this submodule implicitly enables the
            # default matplotlib backend (defined by the "backend" RC
            # parameter in the current matplotlib configuration) if no
            # backend has been enabled, which can have harmful side effects
            # in edge cases.
            from matplotlib import pyplot

            # Switch from the current to the passed backend.
            pyplot.switch_backend(backend_name)

        # Delay importation of this submodule until *ABSOLUTELY* necessary.
        # See pertinent commentary above.
        from matplotlib import pyplot

        # validate the usability of this backend by attempting to create
        # and destroy a hidden empty figure. Technically, doing so could
        # incur unintended side effects on uncooperative platforms (e.g.,
        # Windows) in edge cases.
        #
        # Unfortunately, doing so is also essential. The success of
        # switching to this backend above is a necessary but insufficient
        # condition of this backend's usability. While the success of
        # switching to some backends (e.g., "TkAgg") does reliably imply
        # those backends to be usable, the success of switching to other
        # backends (e.g., "Qt5Agg") only implies that the corresponding
        # packages (e.g., "PyQt5") are successfully importable; this does
        # *NOT* imply these packages and hence these backends to actually
        # be usable in real-world use cases. Notably, the PyQt family of
        # backends are infamous for raising non-human-readable exceptions
        # on attempting to create the first figure. (Why, PyQt? WHY!?!?)
        pyplot.figure()
        pyplot.close()

# ....................{ SINGLETONS                         }....................
mpl_config = MplConfig()
'''
Singleton matplotlib configuration wrapper.
'''
