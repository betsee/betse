#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level support facilities for `matplotlib`, a mandatory runtime dependency.

Backends
----------
`matplotlib` supports numerous **interactive backends** (i.e., bindings to
external GUI-specific widget toolkits), only one of which will be imported by
`matplotlib` at runtime. If the caller specifies no such backend, a default
backend specific to the current system will be imported. However, all such
backends including such defaults are fairly fragile and hence prone to raising
exceptions under common contexts.

The following table summarizes our current findings:

=========  ========   ======  ========  ======  ========  ======  ========
Backend    Footnote   Is Supported Under?
                      ----------------------------------------------------
                      Linux             OS X              Windows
                      ----------------  ----------------  ----------------
                      Frozen  Unfrozen  Frozen  Unfrozen  Frozen  Unfrozen
=========  ========   ======  ========  ======  ========  ======  ========
CocoaAgg   1          No      No        No      No        No      No
Gtk3Agg    2          No      No        No      No        No      No
Gtk3Cairo  3          No      No        ???     ???       ???     ???
MacOSX                No      No        Yes     Yes       No      No
Qt4Agg     3          No      No        ???     ???       ???     ???
TkAgg      4          Yes     Yes       No      Yes       ???     ???
=========  ========   ======  ========  ======  ========  ======  ========

Footnote descriptions are as follows:

1. Backend "CocoaAgg" is officially deprecated and known to be broken.
2. Backend "Gtk3Agg" is known to be broken under Python 3.
3. These backends do *not* support our current animation method, despite
   otherwise working (e.g., to display static plots).
4. Backend `TkAgg` is only freezable under `matplotlib` < 1.4.0. `matplotlib`
   switched Python compatibility layers from `2to3` to `six`, the latter of
   which currently obstructs freezing. This is provisionally correctable under
   `matplotlib` >= 1.4.0 by manually patching the
   `site-packages/matplotlib/backends/backend_tkagg.py` script under the current
   Python 3 interpreter as follows:

    # Replace these header imports...
    from six.moves import tkinter as Tk
    from six.moves import tkinter_filedialog as FileDialog

    # ...with this header import.
    import tkinter as Tk
'''

#FIXME: O.K.; so, basically, Tcl/Tk is absolutely terrible and fundamentally
#does *NOT* work anywhere. At the very least, we need to revert back the change
#that conditionally set the "macosx" backend for that operating system. I'm not
#convinced that the frozen Linux version will actually work in the absence of
#Tcl/Tk data files in the expected system locations. We'll want to test this on
#a pristine laptop -- and soon.

#FIXME: "tkinter" support is, frankly, bizarre. It works *ONLY* under
#"matplotlib-1.3.0". It fails both under "matplotlib-1.4.0" and newer *AND* when
#imported directly (e.g., via "import tkinter") with an inscrutable "ValueError"
#Unicode exception of the sort described at:
#
#    https://github.com/pyinstaller/pyinstaller/issues/1164
#
#This makes absolutely no sense. However, it doesn't have to. Here's what we
#factually know, from which the truth should be deducible:
#
#* The versions of Tcl, Tk, tkinter, and Python are the exact same in all three
#  scenarious. We have personally verified this by recompiling Python and hence
#  tkinter.
#* tkinter works under matplotlib 1.3 and *NOT* under 1.4, both of which are
#  readily installable under Gentoo.
#
#Given that, all we have to do is take the difference of the
#"matplotlib/backends/backend_tkagg.py" files between the two versions and
#figure out what changed. The 1.3 version is given locally, of course; the 1.4
#version is available remotely at:
#
#    https://github.com/matplotlib/matplotlib/tree/master/lib/matplotlib/backends
#
#Our gut intuition is that this has something to do with the "tkinter"-specific
#shared library "_tkagg". The 1.3 version does the following, in order:
#
#* Directly imports "tkinter" which directly imports the shared library
#  "_tkinter".
#* Directly imports "tkagg" which directly imports the shared library
#  "_tkagg".
#
#The 1.4 version similarly does the following, in order:
#
#* Indirectly imports "tkinter" which directly imports the shared library
#  "_tkinter".
#* Directly imports "tkagg" which:
#  * Indirectly imports "tkinter" which directly imports the shared library
#    "_tkinter".
#  * Directly imports the shared library "_tkagg".
#
#It's a small difference, but it's probably enough. The key here is that, since
#PyInstaller fails to detect the indirect importation of "tkinter" and hence
#"_tkinter", "tkagg" and hence "_tkagg" is imported *BEFORE* "_tkinter". This
#absolutely has to be an order-of-shared-library-loading thing. Well, it doesn't
#*HAVE* to be, but that's all we've got. This could be determined by, at some
#point:
#
#* Reinstalling matplotlib 1.4.
#* Iteratively reverting the "backends/backend_tkagg.py" and "backends/tkagg.py"
#  files installed with 1.4 to their 1.3 counterparts *UNTIL* BETSE is
#  freezable. This is pretty much guaranteed to work, albeit annoyingly. But
#  there's no alternative, given the blatant inscrutability of the "tkinter"
#  exception message.

#FIXME: Refactor backend_names() to discover backend names via the standard
#module "pkg_utils" rather than by manually delving through the filesystem,
#which fails under frozen executables.

#FIXME: On attempting to enable the "TkAgg" backend under OS X, we receive the
#following runtime exception from PyInstaller-frozen executables despite such
#freezing appearing to have succeeded:
#
#    Traceback (most recent call last):
#      File "/Users/osxguest/py/betse/betse/cli/cli.py", line 82, in run
#        self._run()
#      File "/Users/osxguest/py/betse/betse/cli/clicli.py", line 114, in _run
#        subcommand_method()
#      File "/Users/osxguest/py/betse/betse/cli/clicli.py", line 329, in _run_sim
#        subcommand_method()
#      File "/Users/osxguest/py/betse/betse/cli/clicli.py", line 356, in _run_sim_plot_init
#        self._get_sim_runner().plotInit()
#      File "/Users/osxguest/py/betse/betse/cli/clicli.py", line 371, in _get_sim_runner
#        from betse.science.simrunner import SimRunner
#      File "/opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/PyInstaller-3.0dev_bfa3b79-py3.4.egg/PyInstaller/loader/pyi_importers.py", line 302, in load_module
#        exec(bytecode, module.__dict__)
#      File "/Users/osxguest/py/betse/betse/science/simrunner.py", line 11, in <module>
#        from betse.science import visualize as viz
#      File "/opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/PyInstaller-3.0dev_bfa3b79-py3.4.egg/PyInstaller/loader/pyi_importers.py", line 302, in load_module
#        exec(bytecode, module.__dict__)
#      File "/Users/osxguest/py/betse/betse/science/visualize.py", line 46, in <module>
#        import matplotlib.pyplot as plt
#      File "/opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/PyInstaller-3.0dev_bfa3b79-py3.4.egg/PyInstaller/loader/pyi_importers.py", line 302, in load_module
#        exec(bytecode, module.__dict__)
#      File "/opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/matplotlib/pyplot.py", line 109, in <module>
#        _backend_mod, new_figure_manager, draw_if_interactive, _show = pylab_setup()
#      File "/opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/matplotlib/backends/__init__.py", line 32, in pylab_setup
#        globals(),locals(),[backend_name],0)
#      File "/opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/PyInstaller-3.0dev_bfa3b79-py3.4.egg/PyInstaller/loader/pyi_importers.py", line 302, in load_module
#        exec(bytecode, module.__dict__)
#      File "/opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/matplotlib/backends/backend_tkagg.py", line 7, in <module>
#        from six.moves import tkinter_filedialog as FileDialog
#      File "/opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/six.py", line 89, in __get__
#        result = self._resolve()
#      File "/opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/six.py", line 108, in _resolve
#        return _import_module(self.mod)
#      File "/opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/six.py", line 79, in _import_module
#        __import__(name)
#    ImportError: No module named 'tkinter.filedialog'
#
#The crux of the issue is the hidden import
#"from six.moves import tkinter_filedialog as FileDialog". Given the popularity
#of the "six" module for purposes of cross-Python[23]-portability, the long-
#term solution is to patch PyInstaller's imports detection to also detect "six"-
#style imports. Since there are only a relatively small number of possible
#imports that "six" supports in this manner, this should be *VERY* feasible.
#Indeed, "py2exe" recently added support for this very functionality via the
#following commit:
#
#    http://sourceforge.net/p/py2exe/svn/764
#
#This is incredibly important. matplotlib officially dropped 2to3 in favor of
#six starting with version 1.4, implying that PyInstaller now basically requires
#six support to also support newer matplotlib versions. To quote: "As of
#matplotlib 1.4, the six library is used to support Python 2 and 3 from a
#single code base."
#
#A short-term solution, alternately, is to patch the "matplotlib_backends()"
#function in "hookutils.py" to append such hidden imports to the list of
#returned hidden imports. This is significantly easier and hence the way to go,
#for now.

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To permit the current backend to be subsequently specified by the
# init() function, the *ONLY* matplotlib module importable at the top-level of
# this module is the eponymous top-level module "matplotlib". See init().
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
from betse.util.io import loggers
from betse.util.path import dirs, paths
from betse.util.python import modules
from betse.util.system import oses
from betse.util.type import containers, strs
from collections import OrderedDict
import sys

# ....................{ IMPORTS ~ matplotlib               }....................
# Import matplotlib in a safe manner. Unfortunately, the "matplotlib.__init__"
# module implicitly imported on the first matplotlib importation performs the
# following unsafe logic:
#
# * The current command-line argument list "sys.argv" is iteratively searched
#   for arguments prefixed by "-d" of length >= 3 (e.g., "-dtkagg" but *NOT*
#   simply "-d").
# * For each such argument, "matplotlib.__init__" attempts to enable the backend
#   whose name is given by such argument excluding the prefixing "-d", silently
#   ignoring exceptions.
#
# This is utterly horrible. Since enabling arbitrary backends can have non-
# negligible side effects, "matplotlib.__init__" *MUST* be prevented from
# performing such logic. Since such module is imported *ONLY* on the first
# importation of a matplotlib module, it suffices to perform such preventation
# *ONLY* for the following importation of the top-level matplotlib package.

# Copy the current argument list into a temporary list.
_sys_argv_old = sys.argv[:]

# Remove all arguments following the mandatory basename of the current process
# from the current argument list, preventing matplotlib from inspecting
# arguments and hence enabling arbitrary backends.
del(sys.argv[1:])

# Import matplotlib safely.
import matplotlib

# Restore the prior argument list from such temporary list.
sys.argv = _sys_argv_old
del(_sys_argv_old)

# ....................{ CONSTANTS                          }....................
RCPARAMS = {
    # Print terse messages. By default, *NO* messages are printed. Valid
    # values include: "silent", "helpful", "debug", and "debug-annoying".
    'verbose.level': 'helpful',
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

# ....................{ CLASSES                            }....................
class MatplotlibConfig(object):
    '''
    `matplotlib` wrapper simplifying configuration and introspection.

    Attributes
    ----------
    _backend_names : list
        List of the strictly lowercase names of all currently available
        `matplotlib`-specific backends (e.g., `['gtk3agg', 'tkagg', 'qt4agg']`).
    '''
    def __init__(self):
        super().__init__()
        self._backend_names = None

    # ..................{ INITIALIZERS                       }..................
    def init(self) -> None:
        '''
        Reconfigure `matplotlib` with sane defaults specific to the current
        system.

        On first importation, `matplotlib` configures itself by loading the
        contents of the first `matplotlibrc` file found in any of several
        candidate directories. Technically, `betse` *could* supply an
        application-specific version of such file to force `matplotlib` to adopt
        application-specific configuration settings. Since synchronizing such
        local copy with remote changes is an onerous (if not ultimately
        infeasible) chore, we elect instead to reconfigure `matplotlib` *after*
        such file has already been loaded at application startup. While this
        slightly increases the cost of such startup, the alternatives are
        impractical at best.

        See Also
        ----------
        http://matplotlib.org/users/customizing.html
            `matplotlib` configuration details.
        '''
        # Reconfigure the following settings, whose keys are the names of settings
        # provided by the official "matplotlibrc" file.
        matplotlib.rcParams.update(RCPARAMS)

        #FIXME: Excise commentary.
        # Configure the backend to be implicitly used for subsequent plotting.
        # Such backend *MUST* be configured prior to the first importation of
        # matplotlib's "pyplot", "matplotlib", or "backends" modules.
        #
        # If the current operating system is Linux or Windows, the "TkAgg"
        # backend is currently set. Alternatives include:
        #
        # * "Qt4Agg", an aesthetically inferior backend *NOT* appearing to
        #   support animation out of the box.
        #
        # If the current operating system is Apple OS X, any non-native backend
        # leveraging the cross-platform C++ library AGG (Anti-grain Geometry) is
        # preferred to native backends. There exist two of the latter:
        #
        # * "CocoaAgg", which leverages AGG but is officially deprecated and
        #   fundametally broken. (Really.)
        # * "MacOSX", which does *NOT* leverage AGG and is known to have
        #   outstanding issues (e.g., the show() method refusing to block).
        #
        # There remain numerous non-native, AGG-based backends. We currently
        # prefer "TkAgg" for the following reasons:
        #
        # * Parity (both visual and functional) with the Linux and Windows
        #   versions.
        # * Stability. While once well-known to sporadically crash, such backend
        #   has proved stable throughout the development of our CLI interface.
        # * Installability. When installing "matplotlib", MacPorts enables by
        #   default the "tkinter" variant and hence such backend but *NO* other
        #   AGG-based backends.

        # If the current operating system is OS X, enable the only backend known
        # to survive freezing: the typical default for such system, "MacOSX".
        if oses.is_os_x():
            self.backend_name = 'MacOSX'
        # Else, the current operating system is Linux or Windows. In such case,
        # enable the only backend known to survive freezing: again, the typical
        # default for such systems, "TkAgg". Unlike backend "MacOSX", however,
        # backend "TkAgg" is fragile and currently requires matplotlib < 1.4.0.
        else:
            # self.backend_name = 'Gtk3Cairo'
            # self.backend_name = 'Qt4Agg'
            self.backend_name = 'TkAgg'

    # ..................{ TESTERS                            }..................
    def is_backend_usable(self, backend_name: str) -> bool:
        '''
        True if the `matplotlib` backend with the passed name is switchable to
        and hence usable on the current system.

        This method modifies but does *not* restore the previously set backend.
        If required, the caller *must* manually save and restore such backend.
        '''
        assert isinstance(backend_name, str),\
            '"{}" not a string.'.format(backend_name)

        # Importing such module has side effects and hence is deferred.
        from matplotlib import pyplot

        # Since backend names are case-insensitive, lowercase such name.
        backend_name = backend_name.lower()

        # If such backend is "Gtk3Agg", immediately report such backend to be
        # unusable without testing such report. Attempting to use such backend
        # currently emits the following warning, which only serves to confound
        # the issue for end users:
        #     UserWarning: The Gtk3Agg backend is not known to work on Python 3.x.
        if backend_name == 'gtk3agg':
            return False

        # Test such backend.
        try:
            pyplot.switch_backend(backend_name)
            return True
        except:
            return False

    # ..................{ GETTERS                            }..................
    def get_metadata(self) -> OrderedDict:
        '''
        Get an ordered dictionary synopsizing the current `matplotlib`
        installation.
        '''
        # Such dictionary.
        metadata = OrderedDict((
            ('rc file', self.rc_filename),
            ('current backend', self.backend_name),
            # ('[backend] current', self.backend_name),
        ))

        # For each available backend, add metadata synopsizing such backend.
        for backend_name in self.backend_names:
            metadata[
                'backend {} usable'.format(backend_name.capitalize())] =\
                str(self.is_backend_usable(backend_name)).lower()

        # Get such dictionary.
        return metadata

    # ..................{ PROPERTIES ~ rc file               }..................
    @property
    def rc_filename(self) -> str:
        '''
        Get the absolute path of the current `matplotlibrc` file establishing
        default `matplotlib` options.
        '''
        return matplotlib.matplotlib_fname()

    # ..................{ PROPERTIES ~ backend name          }..................
    @property
    def backend_name(self) -> str:
        '''
        Get the strictly lowercase name of the current `matplotlib` backend.
        '''
        return matplotlib.get_backend()

    @backend_name.setter
    def backend_name(self, backend_name: str) -> None:
        '''
        Set the current `matplotlib` backend to the backend with the passed
        name.

        Such name is interpreted case-insensitively and hence may be in any case
        including mixed lower and uppercase (e.g., `tkagg`, `TKAGG`, `TkAgg`).
        '''
        assert isinstance(backend_name, str),\
            '"{}" not a string.'.format(backend_name)

        # Since backend names are case-insensitive, lowercase such name.
        backend_name = backend_name.lower()

        #FIXME: Excise this after adding "six" support to PyInstaller.

        # If such backend is "TkAgg", manually import "tkinter". matplotlib 1.4
        # dropped the "2to3" compatibility layer in favor of "six", which hides
        # imports in the form of:
        #
        #     from six.moves import tkinter_filedialog as FileDialog
        #
        # Since PyInstaller fails to detect such imports, manually import the
        # offending modules to notify PyInstaller of such requirements.
        if backend_name == 'tkagg':
            pass
            # import os
            # print('tcl library: ' + str(os.environ.get('TCL_LIBRARY')))
            # print('tk library: ' + str(os.environ.get('TK_LIBRARY')))
            # import tkinter
            # import _tkinter
            # import tkinter.filedialog
            # import tkinter.messagebox

        try:
            # If neither the "matplotlib.pyplot" nor "matplotlib.pylab" modules
            # have been imported yet, prefer setting such backend by calling the
            # non-experimental and hence safer use() function.
            if not modules.is_imported(
                'matplotlib.pyplot', 'matplotlib.pylab'):
                matplotlib.use(backend_name)
            # Else, we have no recourse but to call the experimental
            # switch_backend() function.
            else:
                # Importing such module has side effects and hence is deferred.
                from matplotlib import pyplot
                pyplot.switch_backend(backend_name)
        # Since such functions tend to raise non-human-readable exceptions, log
        # a human-readable error on catching such exception before reraising
        # such exception.
        except Exception:
            loggers.log_error(
                'Matplotlib backend "{}" not found or not loadable.'.format(
                    backend_name))
            raise

        # Log such setting *AFTER* succeeding.
        loggers.log_debug(
            'Enabled matplotlib backend "{}".'.format(backend_name))

    # ..................{ PROPERTIES ~ backend names         }..................
    @property
    def backend_names(self) -> str:
        '''
        Get a list of the strictly lowercase names of all currently available
        `matplotlib`-specific backends (e.g., `['gtk3agg', 'tkagg', 'qt4agg']`).

        While `matplotlib` provides the canonical lists
        `matplotlib.rcsetup.interactive_bk`,
        `matplotlib.rcsetup.non_interactive_bk`, and
        `matplotlib.rcsetup.all_backends`, even the latter such list typically
        fails to list all possible backends (e.g., `mixed` tends to be missing).
        Instead, this function iteratively inspects the current filesystem.
        '''
        # If this funuction has *NOT* been called at least once, create and
        # cache such list.
        if self._backend_names is None:
            # Importing such module has side effects and hence is deferred.
            from matplotlib import backends

            # Absolute path of the directory containing all backends for the
            # currently imported "matplotlib".
            backends_dir = modules.get_dirname(backends)

            # If such directory exists, find all backends in such directory.
            if dirs.is_dir(backends_dir):
                # String prefixing the basenames of backend-specific modules.
                BACKEND_BASENAME_PREFIX = 'backend_'

                # Such names, discovered by:
                #
                # * Filtering all basenames in such directory for modules.
                # * Converting the remaining basenames to backend names.
                # * Sorting such names in ascending lexicographic order. While
                #   *NOT* strictly necessary, of course, such sorting improves
                #   output (e.g., from the "info" subcommand).
                self._backend_names = containers.sort_as_lexicographic_ascending([
                    paths.get_pathname_sans_filetype(
                        strs.remove_prefix_if_found(
                            backend_basename, BACKEND_BASENAME_PREFIX))
                    for backend_basename in dirs.list_basenames(backends_dir)
                    if strs.is_prefix(
                        backend_basename, BACKEND_BASENAME_PREFIX) and
                        paths.is_filetype(backend_basename, 'py')
                ])
            # Else, such directory does *NOT* exist. Log a non-fatal warning and
            # clear such list.
            else:
                loggers.log_warning(
                    'Directory "{}" not found. Matplotlib backends not inspectable.'.format(
                        backends_dir))
                self._backend_names = []

        # Get the cached list.
        return self._backend_names

# ....................{ SINGLETONS                         }....................
config = MatplotlibConfig()
'''
Singleton `matplotlib` configuration wrapper.
'''

# --------------------( WASTELANDS                         )--------------------
            # self.backend_name = 'Cairo'
# In the above table, `(`- and `)`-delimited numbers signify the following
# footnotes of the same number:
# (e.g., under Windows or when frozen
# (e.g., on operating systems, freeze),
# , exactly one of which will be used to
# . If `matplotlib` is *not* explicitly
#  of a the caller fails to notify
#  of its preferred  backend is explicitly specified `matplotlib`
#  *and* , the following

#FUXME: This appears to be required due to a PyInstaller bug. Research.
#FUXME: Actually, even this appears to fail under OS X. We've temporarily
#shifted this back to "betse.ignition", which is hardly ideal.
# import tkinter.filedialog

            #FUXME: Uncomment after correcting this issue.
            # pass
            # for backend_basename in dirs.list_basenames(backends_dir):
            #     if strs.is_prefix(backend_basename, BACKEND_BASENAME_PREFIX) and\
            #        paths.is_filetype(backend_basename, 'py'):
            #         print('backend: ' + strs.remove_prefix_if_found(
            #             backend_basename, BACKEND_BASENAME_PREFIX))

        # if oses.is_os_x():
        #     self.backend_name = 'TkAgg'
        # else:
        #     self.backend_name = 'TkAgg'
        # '''
        # Initialize the `backend_names` attribute.
        # While "self.backend_name" *COULD* be returned instead, such attribute
        # could be incorrect in the event of an external call to matplotlib's
        # use() or switch_backend() functions. For safety, defer to matplotlib.
        # # Nullify attributes for safety.
        # self.backend_name = None

    # backend_name : str
    #     Strictly lowercase name of the current `matplotlib` backend.
        # # If this is the first call to this setter, call the safer use()
        # # function. Such function may *ALWAYS* be safely called
        # if not is_imported(self.backend_name:
        #     matplotlib.use(backend_name)
        #     self.backend_name = backend_name

        # return self.backend if self.backend else matplotlib.get_backend()
# ....................{ SETTERS                            }....................
# def _set_matplotlib_backend_tkagg() -> None:
#     '''
#     Set the current `matplotlib` backend to `TkAgg`, the customary default
#     backend leveraging the `tkinter` GUI toolkit and AGG rendering framework.
#     '''
#     # matplotlib 1.4 dropped the "2to3" compatibility layer in favor of "six",
#     # which hides imports in the form of
#     # "from six.moves import tkinter_filedialog as FileDialog". Since
#     # PyInstaller currently fails to detect such imports, we manually import the
#     # offending modules to notify PyInstaller of these requirements.
#     import tkinter.filedialog
#
#     # Set such backend.
#     matplotlib.use('tkagg')

    # Inspector analyzing currently available `matplotlib` backends.
        # Basenames of all backend-specific Python modules in such directory.
        # backend_basenames = [
    #FUXME: If we continue to be plagued by OS X plotting issues, consider
    #configuring the "verbose.level" parameter at runtime to "debug" below.

# from betse.util.io import loggers
