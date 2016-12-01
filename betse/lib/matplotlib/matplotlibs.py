#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
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
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To permit Matplotlib's default verbosity and backend to be
# initialized to BETSE-specific values, no matplotlib package or module may be
# imported until *AFTER* the MatplotlibConfig.init() method has been called --
# including here at the top-level.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

import sys
from betse.exceptions import BetseMatplotlibException
from betse.util.io.log import logconfig, logs
from betse.util.os import kernels
from betse.util.path import dirs, paths
from betse.util.py import freezers, pys
from betse.util.type import iterables, strs, modules
from betse.util.type.objects import property_cached
from betse.util.type.types import type_check
from collections import OrderedDict
from contextlib import contextmanager
from matplotlib.colors import Colormap

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
Dictionary mapping `matplotlibrc`-specific option names to corresponding values.

The `init()` function subsequently updates `matplotlib.rcParams` (i.e., the
in-memory dictionary of default `matplotlibrc` options deserialized from the
current `matplotlibrc` file) with this dictionary. Hence, the custom options
specified by this dictionary override the default options specified by the
on-disk `matplotlibrc` file.
'''

# ....................{ CONSTANTS ~ z-order                }....................
# These values derive from the following canonical z-order example:
#     http://matplotlib.org/examples/pylab_examples/zorder_demo.html

#FIXME: Shift these constants into a new "betse.lib.matplotlib.zorder" module.

ZORDER_PATCH = 1
'''
Default **z-order** (i.e., positive integer ordering artist drawing, such that
artists with larger z-orders are drawn over artists with smaller z-orders) for
patch artists (e.g., `Patch`, `PatchCollection`).

This is the lowest default z-order, thus drawing patch artists under all other
artists by default.
'''

ZORDER_LINE = 2
'''
Default **z-order** (i.e., positive integer ordering artist drawing, such that
artists with larger z-orders are drawn over artists with smaller z-orders) for
line artists (e.g., `Line2D`, `LineCollection`, `StreamplotSet`).

This is the middle default z-order, thus drawing line artists over all patch
artists but under all text artists by default.
'''

ZORDER_TEXT = 3
'''
Default **z-order** (i.e., positive integer ordering artist drawing, such that
artists with larger z-orders are drawn over artists with smaller z-orders) for
text artists (e.g., `Text`).

This is the highest default z-order, thus drawing text artists over all other
artists by default.
'''

ZORDER_STREAM = (ZORDER_LINE + ZORDER_TEXT) / 2
'''
BETSE-specific **z-order** (i.e., positive integer ordering artist drawing, such
that artists with larger z-orders are drawn over artists with smaller z-orders)
for streamplots (e.g., `StreamplotSet`).

This magic number has been chosen such that streamplots with this z-order will
be drawn over all line and patch artists but under all text artists by default.
Streamplots are technically a variant of line artists but sufficiently non-
linear (and visually busy) to warrant separate handling.

It may also be pertinent to note that recent Matplotlib releases as of this
writing (1.40) accidentally broke backward compatibility with respect to default
streamplot z-order. Specifically, streamplots were assigned a default z-order of
`ZORDER_PATCH` rather than `ZORDER_LINE`:

    https://github.com/matplotlib/matplotlib/pull/5567
'''

# ....................{ GETTERS                            }....................
@type_check
def get_colormap(colormap_name: str) -> Colormap:
    '''
    Matplotlib colormap uniquely identified by the passed name.

    Parameters
    ----------
    colormap_name : str
        Name of the attribute in the `matplotlib.cm` module corresponding to the
        desired colormap (e.g., `Blues`, `Greens`, `jet`, `rainbow).

    Returns
    ----------
    matplotlib.colors.Colormap
        Such colormap.

    See Also
    ----------
    http://matplotlib.org/examples/color/colormaps_reference.html
        List of supported colormaps.
    '''

    # Delay importation of the "matplotlib.__init__" module.
    from matplotlib import cm as colormaps
    from matplotlib.colors import Colormap

    # Colormap with the passed name if any or "None" otherwise.
    colormap = getattr(colormaps, colormap_name, None)
    if not isinstance(colormap, Colormap):
        raise BetseMatplotlibException(
            'Matplotlib colormap "{}" not found.'.format(colormap_name))

    # Return this colormap.
    return colormap

# ....................{ CLASSES                            }....................
class MatplotlibConfig(object):
    '''
    High-level wrapper simplifying low-level configuration and introspection of
    Matplotlib.
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self):

        # Initialize our superclass.
        super().__init__()


    def init(self) -> None:
        '''
        Reconfigure matplotlib with sane defaults specific to the current
        system.

        On first importation, matplotlib configures itself by loading the
        contents of the first `matplotlibrc` file found in any of several
        candidate directories. Technically, BETSE _could_ supply an
        application-specific version of this file to force matplotlib to adopt
        application-specific configuration settings. Since synchronizing this
        local copy with remote changes is an onerous (if not ultimately
        infeasible) chore, we elect instead to reconfigure matplotlib _after_
        this file has already been loaded at application startup. While this
        does increase startup costs, the alternatives are all absurd at best.

        See Also
        ----------
        http://matplotlib.org/users/customizing.html
            Matplotlib configuration details.
        '''

        # Initialize matplotlib in a safe manner.
        self._init_matplotlib()

        # Set the default matplotlib backend *AFTER* initializing matplotlib.
        self._init_backend()

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

        * The current command-line argument list :data:`sys.argv` is iteratively
          searched for arguments matching patterns, including:
          * Arguments prefixed by `-d` of length greater than or equal to 3
            (e.g., `-dtkagg` but _not_ simply `-d`). For each such argument,
            :mod:`matplotlib.__init__` attempts to enable the backend whose name
            is given by such argument excluding the prefixing `-d`, silently
            ignoring exceptions.
          * Arguments prefixed by `--verbose-` matching a matplotlib-specific
            verbosity level (e.g., `--verbose-debug`). For each such argument,
            :mod:`matplotlib.__init__` coerces the global verbosity to that
            level.

        This is utterly horrible. Since enabling arbitrary backends can have
        non-negligible side effects, the :mod:`matplotlib.__init__` submodule
        _must_ be prevented from performing this logic. Since this submodule is
        imported only on the first importation of a matplotlib module, it
        performing this preventation _only_ on the first importation of the
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

            # If BETSE-specific debugging is enabled, enable the mildest form of
            # matplotlib-specific logging. While matplotlib also supports the
            # "debug" and "debug-annoying" verbosity levels, both produce far
            # more effluvia than is helpful.
            if log_level <= logconfig.DEBUG:
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
            logs.log_debug(
                'Initializing matplotlib with options: %s', sys.argv)

            # Import matplotlib *AFTER* setting all matplotlib-specific CLI
            # options above, which this importation parses.
            from matplotlib import rcParams, verbose

            # Prevent the verbose.set_level() method from reducing to a noop, as
            # occurs when this private attribute is *NOT* nullified. (waaaat?)
            verbose._commandLineVerbose = None     # yes, this is horrible
        # Guarantee the prior argument list to be restored from this temporary
        # list even in the event of exceptions.
        finally:
            sys.argv = _sys_argv_old
            del(_sys_argv_old)

        # Unconditionally enable the settings defined by the "RC_PARAMS" global.
        rcParams.update(RC_PARAMS)


    #FIXME: Revise docstring in accordance with the iterative platform-specific
    #fallback process now employed for deciding the default Matplotlib backend.
    def _init_backend(self) -> None:
        '''
        Set the default matplotlib backend to be implicitly used for subsequent
        plotting, which matplotlib requires to be configured _before_ the first
        importation of any the following modules: :mod:`matplotlib.pyplot` or
        :mod:`matplotlib.backends`.

        Specifically:

        * If this method is called by `py.test`-driven testing and hence
          possibly by headless remote continuous integration (CI) with no access
          to a window manager, the non-interactive `Agg` backend is used.
        * Else, if the current platform is:
          * Either Linux or Windows, the `TkAgg` backend is used. This is the
            only backend currently known to survive freezing into executables.
            Alternatives include:
            * `Qt4Agg`, an (_arguably_) aesthetically inferior backend with
              (_inarguably_) significant performance concerns.
          * OS X, the `MacOSX` backend is used. This is the only backend
            currently known to survive freezing into executables. Alternatives
            include:
            * `CocoaAgg`, a non-native backend leveraging the cross-platform C++
              library AGG (Anti-grain Geometry). That's good. Unfortunately,
              this backend is officially deprecated and fundametally broken.
              That's bad.
        '''

        #FIXME: Add transparent support for headless environments here, for
        #which we'll want to default to a headless backend as we currently do
        #when testing (e.g., "Agg"). To do so, add OS-specific logic testing for
        #the presence of a windowing manager as follows:
        #
        #* Under Linux, test for (in order):
        #  * X11 connectivity.
        #  * Wayland connectivity.
        #  * Mir connectivity.
        #* Under OS X and Windows, don't bother testing anything. Windows
        #  explicitly fails to support headless operation. OS X technically
        #  supports headless operation via an obscure (albeit well-documented)
        #  hack where one enters ">console" as the login username, but basically
        #  fails to support headless operation for most intents and purposes.
        #
        #This is critical under Linux, as it would render BETSE amenable to
        #scripted, remote, and (hopefully) parallelized usage.
        #FIXME: Actually, is the above comment relevant anymore in light of the
        #inclusion of the headless "Agg" backend listed below? Specifically, are
        #non-headless backends (e.g., "TkAgg") actually importable under
        #headless environments? Hopefully, this is *NOT* the case. If this is
        #indeed not the case, then no further work needs be done; the logic
        #below already implicitly handles headless environments.

        # Name of the non-GUI-based matplotlib backend to fallback to in the
        # event that *NO* GUI-based matplotlib backend is usable on this system.
        # For portability, this backend is guaranteed to be usable on all
        # platforms and systems regardless of matplotlib version.
        _BACKEND_NAME_FALLBACK = 'Agg'

        # Dictionary mapping from the name of each supported platform to a tuple
        # of the names of all matplotlib backends to iteratively fallback to on
        # that platform (in descending order of preference) in the event the the
        # end user fails to explicitly specify such a name (e.g., via the
        # `--matplotlib-backend` option). In this case, this method subsequently
        # falls back to the first matplotlib backend usable on the current
        # system whose name is in this tuple.
        _KERNEL_NAME_TO_BACKEND_NAMES_PREFERRED = {
            #FIXME: Inject the "Qt5Agg" backend somewhere into this list after
            #shown to be working. Sadly, this backend appears to be overly
            #fragile and hence unusable. Attempting to use this backend on
            #Gentoo Linux with matplotlib 1.5.1 yields the following fatal
            #exception on attempting to plot or animate:
            #
            #    QObject::connect: Cannot connect NavigationToolbar2QT::message(QString) to (null)::_show_message()
            #    TypeError: connect() failed between NavigationToolbar2QT.message[str] and _show_message()
            #FIXME: Inject the "WxAgg" backend somewhere into this list after
            #shown to be working under Python 3.x. Sadly, since a stable version
            #of WxPython Phoenix has yet to be released, this may take
            #considerably longer than first assumed. The vapourware: it burns!

            # Under Linux, the following backends are preferred:
            #
            # * "TkAgg", a GUI backend with adequate (albeit not particularly
            #   impressive) aesthetics and superior performance by compare to
            #   less preferable backends. Tcl/Tk: who would have ever thought?
            # * "Qt4Agg", a GUI backend with (arguably) inferior aesthetics and
            #   (inarguably) significant performance concerns by compare to more
            #   preferable backends. Something is better than nothing.
            'Linux': ('TkAgg', 'Qt4Agg',),
            # 'Linux': (),

            # Under OS X and iOS, the preferred backends are defined below in
            # terms of the preferred Linux-specific backends defined above.
            'Darwin': None,

            # Under Windows, the following backends are preferred:
            #
            # * "Qt4Agg". While typically less preferable than "TkAgg", this
            #   backend's implementation is somewhat stabler under many Windows
            #   systems than that of "TkAgg". Hence, stability wins.
            # * "TkAgg". (See the prior item.)
            'Windows': ('Qt4Agg', 'TkAgg',),
        }

        # Under OS X and iOS, prefer the only genuinely usable Darwin-specific
        # backend before the same backends as preferred under Linux. Since
        # Darwin and Linux are both POSIX-compatible, cross-platform backends
        # (e.g., "Qt5Agg") tend to behave similarly under both platforms.
        _KERNEL_NAME_TO_BACKEND_NAMES_PREFERRED['Darwin'] = (
            ('MacOSX',) + _KERNEL_NAME_TO_BACKEND_NAMES_PREFERRED['Linux'])

        #FIXME: Refactor the test suite to:
        #
        #* In functional tests, explicitly pass the "--matplotlib-backend=agg"
        #  CLI option to all invocations of the "betse" command.
        #* In unit tests, explicitly set the default backend to be used...
        #  somehow. But how, exactly? Perhaps this method could accept the name
        #  of the desired matplotlib backend (if any) as a parameter; this
        #  method could then be called by a unit test fixture to default to the
        #  "Agg" backend.
        #
        #After doing so, remove the "if pys.is_testing():" branch below.

        # If tests are being run, default to the non-interactive "Agg" backend
        if pys.is_testing():
            self.backend_name = _BACKEND_NAME_FALLBACK
            return

        # Else, tests are *NOT* being run. Default to the first backend usable
        # under the current system.
        #
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
            if self.is_backend_usable(_BACKEND_NAME_FALLBACK):
                logs.log_warning(
                    'Usable matplotlib GUI backend not found. '
                    'Defaulting to non-GUI backend "%s". '
                    'Consider installing support for GUI backends %s.',
                    _BACKEND_NAME_FALLBACK,
                    strs.join_as_disconjunction_double_quoted(*backend_names),
                )
            # Else, no backends appear to be usable on the current system. For
            # safety, raise an exception. Due to the ubiquity of the fallback
            # backend, this should *NEVER* happen. Damn you, Murphy!
            else:
                raise BetseMatplotlibException(
                    'Usable matplotlib backend not found.')

    # ..................{ PROPERTIES ~ rc                    }..................
    @property
    def rc_filename(self) -> str:
        '''
        Absolute path of the current `matplotlibrc` file establishing default
        matplotlib options.
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

        # If this name is unrecognized, raise an exception. For unknown reasons,
        # the Verbose.set_level() method called below unsafely emits non-fatal
        # warnings rather than raising fatal exceptions on receiving an
        # unrecognized level name.
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
        # Since the latter also implicitly requires creating an empty canvas and
        # is hence both more expensive and more complex, we adopt the former.
        backend_canvas_class = self.backend_canvas_class

        # The "FigureCanvasBase" superclass of all canvas subclasses defines the
        # public class attribute "filetypes" to be the dictionary of all
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
        `matplotlib.rcsetup.interactive_bk`,
        `matplotlib.rcsetup.non_interactive_bk`, and
        `matplotlib.rcsetup.all_backends`, even the latter such list typically
        fails to list all possible backends (e.g., `mixed` tends to be missing).
        Instead, this function iteratively inspects the current filesystem.
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
        Human-readable name (e.g., `Qt5Agg`) of the current backend.

        This name is _not_ guaranteed to be lowercase and, in fact, is typically
        a mix of upper- and lowercase alphanumeric characters.
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

        This name is interpreted case-insensitively and hence may be in any case
        including mixed lower and uppercase (e.g., `tkagg`, `TKAGG`, `TkAgg`).
        '''

        # Delay importation of the "matplotlib.__init__" module.
        import matplotlib

        # Lowercase this name, as backend names are case-insensitive.
        backend_name = backend_name.lower()

        try:
            # If no backend has been set yet, set this backend by calling the
            # non-experimental and hence safer use() function.
            if not self.is_backend():
                matplotlib.use(backend_name)
            # Else, we have no recourse but to call the experimental
            # switch_backend() function.
            else:
                # Avoid importing this module at the top-level, for safety.
                from matplotlib import pyplot
                pyplot.switch_backend(backend_name)
        # If presumably non-human-readable exception was raised...
        except Exception:
            # Log a human-readable error.
            logs.log_error(
                'Matplotlib backend "{}" not found or not loadable.'.format(
                    backend_name))

            # Reraise this exception up the callstack.
            raise

        # Log this setting *AFTER* successfully setting this backend.
        logs.log_debug(
            'Enabled matplotlib backend "{}".'.format(backend_name))

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

        This context manager guaranteeably reverts this level to the prior level
        even when fatal exceptions are raised. If this level is _not_ currently
        `helpful`, this context manager is a noop.

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

    # ..................{ TESTERS                            }..................
    def is_backend(self) -> bool:
        '''
        `True` if a backend has been set (i.e., if the `matplotlib.use()` method
        has been called for the current Python session) or `False` otherwise.
        '''

        # While ridiculous, this test corresponds exactly to the test performed
        # by the matplotlib.use() method itself to detect repetitious calls.
        return 'matplotlib.backends' in sys.modules


    @type_check
    def is_backend_usable(self, backend_name: str) -> bool:
        '''
        `True` if the backend with the passed name is **usable** (i.e., safely
        switchable to without raising exceptions) on the current system.

        If this backend is usable, this method also switches the current backend
        to this backend _without_ restoring the previously set backend. This is
        as the unavoidable result of performing this test in a reliable manner.
        Caller requiring the previously set backend to be restored must do so
        manually _after_ calling this method.
        '''

        # Delay importation of the "matplotlib.__init__" module.
        from matplotlib import pyplot

        # Lowercase this name, ensuring case-insensitive backend names.
        backend_name = backend_name.lower()

        # If this backend is "Gtk3Agg", immediately report this backend to be
        # unusable without testing such report. Attempting to use this backend
        # currently emits the following warning, which only serves to confound
        # the issue for end users:
        #
        #     UserWarning: The Gtk3Agg backend is not known to work on Python 3.x.
        if backend_name == 'gtk3agg':
            return False

        # If this backend is successfully enabled and hence usable, return True.
        try:
            pyplot.switch_backend(backend_name)
            return True
        # Else, return False.
        except:
            return False

    # ..................{ GETTERS                            }..................
    def get_rc_param(self, param_name) -> object:
        '''
        Value of the parameter with the passed `.`-delimited name (e.g.,
        `savefile.dpi`) in the external `matplotlibrc` file.

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


    def get_metadata(self) -> OrderedDict:
        '''
        Ordered dictionary synopsizing the current matplotlib installation.
        '''

        # This dictionary.
        metadata = OrderedDict((
            ('rc file', self.rc_filename),
            ('current backend', self.backend_name),
            # ('[backend] current', self.backend_name),
        ))

        # For each available backend, add metadata synopsizing that backend.
        for backend_name in self.backend_names:
            metadata[
                'backend {} usable'.format(backend_name.capitalize())] = str(
                    self.is_backend_usable(backend_name)).lower()

        # Get this dictionary.
        return metadata

    # ..................{ MAKERS                             }..................
    #FIXME: Define a new make_backend_figure_manager() method as well. This
    #method is fine, for now. We might very well want to call this sometime.
    def make_backend_figure(self, *args, **kwargs):
        '''
        Create and return a new `Figure` instance passed the passed arguments
        _and_ associated with a new `FigureCanvasBase` instance corresponding to
        the current backend (e.g., `FigureCanvasQt5Agg` for the `qt5agg`
        backend).

        This figure will _not_ be associated with a new or existing
        `FigureManagerBase` instance corresponding to the current backend (e.g.,
        `FigureManagerQt5Agg` for the `qt5agg` backend) _or_ with the `pyplot`
        interface. In theory, this implies that this figure should be implicitly
        garbage collected on:

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

# ....................{ SINGLETONS                         }....................
mpl_config = MatplotlibConfig()
'''
Singleton matplotlib configuration wrapper.
'''
