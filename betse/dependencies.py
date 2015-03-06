#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Mandatory runtime dependency checking.

This module defines functions intended to be called by high-level interface
modules (e.g., `betse.cli.cli`) *before* attempting to import such dependencies.
'''

#FIXME: Refactor as follows:
#
#* Create a new subpackage "betse.util.dependency".
#* Split this module into two modules:
#  * "betse.util.dependency.dependencies", providing the exception handling
#    function.
#  * "betse.util.dependency.matplotlib", providing the matplotlib-specific
#    functions.
#
#Hence, this module survives, albeit in *EXTREMELY* limited form. That's fine,
#for now. (Better minimal than overkill, we should think.)

#FIXME: It'd be great to raise human-readable exceptions on the specified
#backends *NOT* being available. This is certainly feasible, as the
#following stackoverflow answer demonstrates -- if somewhat involved:
#    https://stackoverflow.com/questions/5091993/list-of-all-available-matplotlib-backends
#That said, we really want to do this *ANYWAY* to print such list when running
#"betse info". So, let's just get this done, please.

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
# WARNING: To raise human-readable exceptions on missing mandatory dependencies,
# the top-level of this module may import *ONLY* from packages guaranteed to
# exist at installation time (e.g., stock Python packages).
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
from betse import metadata
from betse.util.python import modules
from betse.util.system import oses

# ....................{ INITIALIZERS                       }....................
def init() -> None:
    '''
    Initialize all mandatory runtime dependencies of `betse`.

    Specifically (in order):

    * Raise an exception unless all such dependencies are currently satisfiable.
    * Reconfigure `matplotlib` with sane defaults specific to the current
      system.
    '''
    # Ensure that all mandatory dependencies exist *BEFORE* subsequent logic
    # (possibly) importing such dependencies.
    die_unless_satisfiable()

    # Configure such dependencies.
    configure_matplotlib()

# ....................{ EXCEPTIONS                         }....................
def die_unless_satisfiable() -> None:
    '''
    Raise an exception unless mandatory runtime dependencies of `betse` are
    **satisfiable** (i.e., importable and of a satisfactory version).

    Equivalently, an exception is raised if at least one such dependency is
    unsatisfied.

    Specifically, this function unconditionally validates the existence of all
    such dependencies. For such dependencies with `setuptools`-specific metadata
    (e.g., `.egg-info/`-suffixed subdirectories of the `site-packages/`
    directory for the active Python 3 interpreter, typically created by
    `setuptools` at install time), this function additionally validates the
    versions of such dependencies to satisfy `betse` requirements.
    '''
    # Template for exception messages raised on missing dependencies.
    exception_template = 'Mandatory dependency "{}" not found.'

    # If the setuptools-specific "pkg_resources" dependency is missing, fail
    # *BEFORE* attempting to import such dependency below.
    modules.die_unless(
        module_name = 'pkg_resources',
        exception_message = exception_template.format('pkg_resources')
    )

    # Import such dependency and all required classes in such dependency.
    from pkg_resources import (Environment, VersionConflict, WorkingSet)
    import pkg_resources

    # Set of all betse-specific dependencies in a manner usable below.
    requirements = pkg_resources.parse_requirements(
        metadata.DEPENDENCIES_RUNTIME)

    # Set of all setuptools-installed Python packages available under the active
    # Python 3 interpreter.
    working_set = WorkingSet()

    # Helper for finding the former in the latter.
    environment = Environment(working_set.entries)

    # Validate each such dependency.
    for requirement_required in requirements:
        # Name of the top-level importable module provided by such project.
        module_name = requirement_required.project_name

        # If such dependency is missing, fail.
        if not modules.is_module(module_name):
            # If such dependency is "setuptools", reduce such dependency to
            # merely "pkg_resources" and try again. "pkg_resources" is a single
            # module installed with (but inexplicably outside of the package
            # tree of) setuptools. BETSE requires setuptools and hence
            # "pkg_resources" at install time but *ONLY* the latter at runtime.
            if module_name == 'setuptools':
                module_name = 'pkg_resources'

            # Try again.
            modules.die_unless(
                module_name, exception_template.format(
                    str(requirement_required)))

        # Else, such dependency exists.
        #
        # Best match for such dependency under the active Python 3 interpreter
        # if any or None otherwise.
        requirement_provided = environment.best_match(
            requirement_required, working_set)

        # If a match was found with version conflicting with that required,
        # fail.
        if requirement_provided is not None and\
           requirement_provided not in requirement_required:
               raise VersionConflict(
                   'Mandatory dependency {} required but only {} found.'.format(
                       requirement_required, requirement_provided))

# ....................{ CONFIGURERS                        }....................
def configure_matplotlib() -> None:
    '''
    Reconfigure `matplotlib` with sane defaults specific to the current system.

    On first importation, `matplotlib` configures itself by loading the contents
    of the first `matplotlibrc` file found in any of several candidate
    directories. Technically, `betse` *could* supply an application-specific
    version of such file to force `matplotlib` to adopt application-specific
    configuration settings. Since synchronizing such local copy with remote
    changes is an onerous (if not ultimately infeasible) chore, we elect instead
    to reconfigure `matplotlib` *after* such file has already been loaded at
    application startup. While this slightly increases the cost of such startup,
    the alternatives are impractical at best.

    See Also
    ----------
    http://matplotlib.org/users/customizing.html
        Further details on `matplotlib` configuration.
    '''
    # Avoid circular import dependencies, as well as optional dependency
    # imports prohibited at the top-level of this module.
    from betse.util.io import loggers

    #FIXME: Nonsense. After shifting these functions to a new module, shift this
    #importation to the top of such module.
    import matplotlib

    # Reconfigure the following settings, whose keys are the names of settings
    # provided by the official "matplotlibrc" file.
    matplotlib.rcParams.update({
        #FIXME: If we continue to be plagued by OS X plotting issues, consider
        #configuring the "verbose.level" parameter at runtime to "debug" below.

        # Print terse messages. By default, *NO* messages are printed. Valid
        # values include: "silent", "helpful", "debug", and "debug-annoying".
        'verbose.level': 'helpful',
        # 'verbose.level': 'debug',
    })

    # Configure the backend to be implicitly used for *ALL* subsequent plotting.
    # Such backend *MUST* be configured prior to the first importation of either
    # the "matplotlib.pyplot" or "matplotlib.pylab" modules.  Since backend
    # names are case-insensitive, lowercase names are preferred below.
    #
    # If the current operating system is Apple OS X, prefer any non-native
    # backend leveraging the cross-platform C++ library AGG (Anti-grain
    # Geometry) to native backends. There exist two native backends for OS X:
    #
    # * "CocoaAgg", which leverages AGG but is officially deprecated and
    #   fundametally broken. (Really.)
    # * "MacOSX", which does *NOT* leverage AGG and is known to have outstanding
    #   issues (e.g., the show() method refusing to block).
    #
    # There remain numerous non-native, AGG-based backends. We currently prefer
    # "TkAgg" for the following reasons:
    #
    # * Parity (both visual and functional) with the Linux and Windows versions.
    # * Stability. While once well-known to sporadically crash, such backend
    #   has proved stable throughout the development of our CLI interface.
    # * Installability. When installing "matplotlib", MacPorts enables by
    #   default the "tkinter" variant and hence such backend but *NO* other AGG-
    #   based backends.
    if oses.is_os_x():
        _set_matplotlib_backend_tkagg()
    # Else, prefer the "TkAgg" backend. Alternatives include:
    #
    # * "Qt4Agg", an aesthetically inferior backend *NOT* appearing to support
    #   animation out of the box. (That's interesting, in the bad way.)
    else:
        _set_matplotlib_backend_tkagg()

def _set_matplotlib_backend_tkagg() -> None:
    '''
    Set `TkAgg` to be the current `matplotlib` backend.
    '''
    #FIXME: Nonsense. See above.
    import matplotlib

    # matplotlib 1.4 dropped the "2to3" compatibility layer in favor of "six",
    # which hides imports in the form of
    # "from six.moves import tkinter_filedialog as FileDialog". Since
    # PyInstaller currently fails to detect such imports, we manually import the
    # offending modules to notify PyInstaller of these requirements.
    import tkinter.filedialog

    # Set such backend.
    matplotlib.use('tkagg')

# --------------------( WASTELANDS                         )--------------------
    # If the current operating system is Apple OS X, prefer the "CocoaAgg"
    # backend to the "MacOSX" backend. The former leverages the cross-platform
    # C++ library AGG (Anti-grain Geometry) and hence tends to be better
    # supported; the latter does not.
        #FIXME: Extract into a new packages.is_package_PyObjC() function.
        # If PyObjC is installed, enable the "CocoaAgg" backend, which
        # internally requires such dependency.
        #if modules.is_module('PyObjCTools'):
#           loggers.log_warning(
#               'Optional dependency "PyObjC" not found. '
#               'Falling back from matplotlib backend "CocoaAgg" to "MacOSX".'
#           )
#           matplotlib.use('macosx')

        #FUXME: Alternately, perhaps we want to redirect
            # exception_template.format(metadata.DEPENDENCY_SETUPTOOLS)
    # List of setuptools requirements strings signifying all safely testable
    # mandatory dependencies. Sadly, the following mandatory dependencies are
    # *NOT* safely testable:
    #
    # * PySide. Under numerous Linux distributions, PySide is installed with
    #   cmake rather than setuptools
    # 'pyside >= 1.2.0',

#FUXME: Implement setuptools-based version checking as well, for dependencies
#providing setuptools-specific egg metadata. stackoverflow offered the following
#useful example code:
#
#    from pkg_resources import (WorkingSet, DistributionNotFound)
#    working_set = WorkingSet()
#
#    # Printing all installed modules
#    print tuple(working_set)
#
#    # Detecting if module is installed
#    try:
#        dep = working_set.require('paramiko>=1.0')
#    except DistributionNotFound:
#        pass
#
#Given that, consider the following approach:
#
#    import metadata
#    import pkg_resources
#    from pkg_resources import (
#        DistributionNotFound, Environment, VersionConflict, WorkingSet)
#    working_set = WorkingSet()
#    environment = Environment(working_set.entries)
#    requirements = pkg_resources.parse_requirements(
#        metadata.REQUIREMENTS)
#    for requirement in requirements:
#        if not importlib.find_loader(requirement.project_name):
#            raise DistributionNotFound(
#                "Mandatory dependency {} not found.".format(requirement))
#            )
#
#        requirement_distribution = environment.best_match(
#            requirement, working_set)
#
#        # Oops, the "best" so far conflicts with a dependency.
#        if requirement_distribution and\
#           requirement_distribution not in requirements:
#               raise VersionConflict(
#                   "Mandatory dependency {} found but {} required.".format(
#                       requirement_distribution, requirement))
#
#The above *SHOULD* work, but assumes use of a new
#"metadata.install_requirements" module dict. Not terribly arduous, of course;
#merely shift such dict from its current use in "setup.py".

        # if not importlib.find_loader(requirement.project_name):
        #     raise DistributionNotFound(
        #         'Mandatory dependency {} not found.'.format(requirement))
    # for module_name in {
    #     'scipy'
    # }:
    # For each such dependency, this function also attempts to validate such
    # dependency's version. Specifically, if such dependency both exists *and*,
    # an exception is raised if

    # Caveats
    # ----------
    # Dependency versions are *not* validated. This is subject to change
