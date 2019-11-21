#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Root test configuration** (i.e., early-time configuration guaranteed to be
run by :mod:`pytest` *before* passed command-line arguments are parsed) for
this test suite.

Caveats
----------
For safety, this configuration should contain *only* early-time hooks
absolutely required by :mod:`pytest` design to be defined in this
configuration. Hooks for which this is the case (e.g.,
:func:`pytest_addoption`) are explicitly annotated as such in official
:mod:`pytest` documentation with a note resembling:

    Note

    This function should be implemented only in plugins or ``conftest.py``
    files situated at the tests root directory due to how pytest discovers
    plugins during startup.

This file is the aforementioned ``conftest.py`` file "...situated at the tests
root directory."

See Also
----------
:mod:`betse_test.conftest`
    Global test configuration applied after this configuration.
'''

# ....................{ IMPORTS                           }....................
import os, sys

# ....................{ HOOKS ~ option                    }....................
def pytest_addoption(parser: '_pytest.config.Parser') -> None:
    '''
    Hook run immediately on :mod:`pytest` startup *before* parsing command-line
    arguments (and hence performing test collection), typically registering
    application-specific :mod:`argparse`-style options and ini-style config
    values.

    Options
    ----------
    After :mod:`pytest` parses these options, the
    ``pytestconfig.getoption({option_var_name})`` method of the
    ``pytestconfig`` fixture provides the value of the argument accepted by
    each option (if any), where ``{option_var_name}`` is the value of the
    ``dest`` keyword argument passed to the :meth:`parser.add_option` method in
    the body of this hook.

    Caveats
    ----------
    This hook should be implemented *only* in plugins or ``conftest.py`` files
    situated at the top-level tests directory for this application (e.g., like
    the current file), due to plugin discovery by :mod:`pytest` at startup.

    Parameters
    ----------
    parser : _pytest.config.Parser
        :mod:`pytest`-specific command-line argument parser, inspired by the
        :mod:`argparse` API.
    '''

    #FIXME: Sample option specification preserved entirely for posterity.
    # # String argument options (i.e., options requiring a string argument),
    # # disabled unless explicitly passed.
    # parser.addoption(
    #     '--export-sim-conf-dir',
    #     dest='export_sim_conf_dirname',
    #     default=None,
    #     help=(
    #         'target directory into which all '
    #         'source simulation configuration directories produced by '
    #         '"@skip_unless_export_sim_conf"-marked tests are to be copied'
    #     ),
    #     metavar='DIRNAME',
    # )

    pass

# ....................{ HOOKS ~ session : start           }....................
def pytest_sessionstart(session: '_pytest.main.Session') -> None:
    '''
    Hook run immediately *before* starting the current test session (i.e.,
    calling the :func:`pytest.session.main` function).

    Parameters
    ----------
    session: _pytest.main.Session
        :mod:`pytest`-specific test session object.
    '''

    # Sanitize import directories *BEFORE* the first module importation.
    _clean_imports()

    # Print test-specific metadata *AFTER* sanitizing these directories.
    _print_metadata()


def _clean_imports() -> None:
    '''
    Sanitize and validate import directories (i.e., the global :attr:`sys.list`
    of the absolute and relative dirnames of all directories to search for
    modules and packages to be imported from).

    Specifically, this function:

    * If this low-level :mod:`pytest` test harness is *not* isolated to a venv
      (e.g., due to being exercised by the low-level ``pytest`` command),
      reduce to a noop.
    * Else, this harness is isolated to a venv (e.g., due to being exercised by
      the high-level ``tox`` command):

      #. If the top-level directory for this project is listed in the global
         list of all import directories (i.e., :attr:`sys.path`), remove this
         directory from this list. Doing so prevents this test session from
         accidentally importing from modules *not* isolated to this venv,
         including this project being tested.
      #. If the first directory on this list is *not* isolated to this venv,
         raise an exception. This condition implies that modules will be
         imported from outside this venv, which entirely defeats the purpose of
         isolating tests with :mod:`tox` to a venv in the first place.
      #. If the top-level :mod:`betse` package is *not* isolated to this venv,
         raise an exception. This condition implies that this project has been
         imported from outside this venv -- again defeating the purpose.

    Raises
    ----------
    ValueError
        If either the first directory on :attr:`sys.path` or the top-level
        :mod:`betse` package are *not* isolated to this venv.
    '''

    # True only if tests are isolated to a venv produced by either...
    #
    # See the betse.util.py.pvenv.is_venv() function, whose implementation is
    # inlined below. While calling that function directly would (of course) be
    # preferable, doing so invites chicken-and-egg issues by importing *BEFORE*
    # sanitizing import directories.
    is_venv = (
        # "virtualenv", which uniquely defines the "sys.real_prefix"
        # attribute to the absolute dirname of the top-level directory
        # containing the system-wide Python interpreter *OR*...
        hasattr(sys, 'real_prefix') or

        # "venv", which (possibly non-uniquely) sets:
        #
        # * The "sys.base_prefix" attribute to the absolute dirname of the
        #   top-level directory containing the system-wide Python interpreter.
        # * The "sys.prefix" attribute to the absolute dirname of the
        #   top-level directory containing the venv-specific Python interpreter
        #   if any *OR* the system-wide Python interpreter otherwise.
        #
        # Note that, as Python >= 3.3 *ALWAYS* defines the "sys.base_prefix"
        # attribute, testing this attribute's existence is unnecessary.
        sys.prefix != sys.base_prefix
    )

    # Print a header for disambiguity.
    print('------[ venv ]------')

    # Print whether tests are isolated to a venv.
    print('venv test isolation: {}'.format(is_venv))

    # If tests are *NOT* isolated to a venv, silently reduce to a noop.
    if not is_venv:
        return
    # ELse, tests are isolated to a venv.

    # Print the absolute dirname of this venv's top-level directory.
    print('venv dir: {}'.format(sys.prefix))

    # List of the absolute dirnames of all directories to search for
    # modules and packages to be imported from, guaranteed to be isolated
    # to this venv.
    sys_path_new = []

    # Absolute dirname of this project's top-level directory.
    PROJECT_DIRNAME = os.path.dirname(__file__)

    # Absolute dirname of this venv's top-level directory, suffixed by a
    # directory separator for disambiguity when calling str.startswith() below.
    VENV_DIRNAME = sys.prefix + os.path.sep

    # For the dirname of each directory to search for Python imports...
    for import_dirname in sys.path:
        # If this dirname is the empty string implying this project's top-level
        # directory, omit this dirname from this list and warn the user.
        if not import_dirname:
            print(
                'WARNING: Ignoring non-isolated empty import directory...',
                file=sys.stderr)
        # Else if this dirname is that of this project's top-level directory,
        # omit this dirname from this list and warn the user.
        elif import_dirname == PROJECT_DIRNAME:
            print(
                'WARNING: '
                'Ignoring non-isolated import directory "{}"...'.format(
                    import_dirname),
                file=sys.stderr)
        # Else, preserve this dirname in this list as is.
        else:
            sys_path_new.append(import_dirname)
    # This list has now been purged of offending dirnames.

    #FIXME: Iteration stripping all non-isolated dirnames from "sys.path".
    #While currently broken due to venv packages failing to adequately isolate
    #venvs from system-wide directories, we hope to reenable this... sometime.
    #
    #Specifically, both "venv" and "virtualenv" appear to create insufficient
    #and arguably broken virtual environments whose
    #"lib/python${PYTHON_VERSION}/" subdirectories contain only a proper subset
    #of all requisite stdlib files -- thus requiring that the equivalent
    #system-wide Python dirnames remain on "sys.path". Removing these dirnames
    #induces the following fatal exception on attempting to import the stdlib
    #"pkgutil" submodule from within a purportedly isolated "tox" test venv:
    #
    #    INTERNALERROR> Traceback (most recent call last):
    #    INTERNALERROR>   File "/home/leycec/py/betse/.tox/py36/lib/python3.6/site-packages/_pytest/main.py", line 194, in wrap_session
    #    INTERNALERROR>     config.hook.pytest_sessionstart(session=session)
    #    INTERNALERROR>   File "/home/leycec/py/betse/.tox/py36/lib/python3.6/site-packages/pluggy/hooks.py", line 286, in __call__
    #    INTERNALERROR>     return self._hookexec(self, self.get_hookimpls(), kwargs)
    #    INTERNALERROR>   File "/home/leycec/py/betse/.tox/py36/lib/python3.6/site-packages/pluggy/manager.py", line 92, in _hookexec
    #    INTERNALERROR>     return self._inner_hookexec(hook, methods, kwargs)
    #    INTERNALERROR>   File "/home/leycec/py/betse/.tox/py36/lib/python3.6/site-packages/pluggy/manager.py", line 86, in <lambda>
    #    INTERNALERROR>     firstresult=hook.spec.opts.get("firstresult") if hook.spec else False,
    #    INTERNALERROR>   File "/home/leycec/py/betse/.tox/py36/lib/python3.6/site-packages/pluggy/callers.py", line 208, in _multicall
    #    INTERNALERROR>     return outcome.get_result()
    #    INTERNALERROR>   File "/home/leycec/py/betse/.tox/py36/lib/python3.6/site-packages/pluggy/callers.py", line 80, in get_result
    #    INTERNALERROR>     raise ex[1].with_traceback(ex[2])
    #    INTERNALERROR>   File "/home/leycec/py/betse/.tox/py36/lib/python3.6/site-packages/pluggy/callers.py", line 187, in _multicall
    #    INTERNALERROR>     res = hook_impl.function(*args)
    #    INTERNALERROR>   File "/home/leycec/py/betse/conftest.py", line 106, in pytest_sessionstart
    #    INTERNALERROR>     _print_metadata()
    #    INTERNALERROR>   File "/home/leycec/py/betse/conftest.py", line 260, in _print_metadata
    #    INTERNALERROR>     from betse.util.py.module import pymodule
    #    INTERNALERROR>   File "/home/leycec/py/betse/.tox/py36/lib/python3.6/site-packages/betse/util/py/module/pymodule.py", line 31, in <module>
    #    INTERNALERROR>     from betse.util.io.log import logs
    #    INTERNALERROR>   File "/home/leycec/py/betse/.tox/py36/lib/python3.6/site-packages/betse/util/io/log/logs.py", line 54, in <module>
    #    INTERNALERROR>     from betse.util.type import types
    #    INTERNALERROR>   File "/home/leycec/py/betse/.tox/py36/lib/python3.6/site-packages/betse/util/type/types.py", line 22, in <module>
    #    INTERNALERROR>     import functools, inspect, logging, pkg_resources, re
    #    INTERNALERROR>   File "/home/leycec/py/betse/.tox/py36/lib/python3.6/site-packages/pkg_resources/__init__.py", line 31, in <module>
    #    INTERNALERROR>     import pkgutil
    #    INTERNALERROR> ModuleNotFoundError: No module named 'pkgutil'
    #    ERROR: InvocationError for command /home/leycec/py/betse/.tox/py36/bin/pytest /home/leycec/py/betse (exited with code 3)
    #
    #Until resolved, the following *MUST* be temporarily disabled:
    # # For the dirname of each directory to search for imports...
    # for import_dirname in sys.path:
    #     # If this directory resides inside this venv, preserve this
    #     # directory in this list as is.
    #     if import_dirname.startswith(VENV_DIRNAME):
    #         sys_path_new.append(import_dirname)
    #     # Else, this directory resides outside this venv. In this case,
    #     # narn that this directory will *NOT* be importable from.
    #     else:
    #         print(
    #             'WARNING: '
    #             'Ignoring non-isolated import directory "{}"...'.format(
    #                 import_dirname),
    #             file=sys.stderr)

    # Replace the original such list with this redacted list.
    sys.path = sys_path_new
    # print('import paths: ' + str(sys.path))

    # First dirname on this list *AFTER* replacing this list.
    import_first_dirname = sys.path[0]

    # If this dirname is *NOT* suffixed by a directory separator, do so.
    if import_first_dirname[-1] != os.path.sep:
        import_first_dirname += os.path.sep

    # If this dirname is *NOT* isolated to this venv, raise an exception.
    if not import_first_dirname.startswith(VENV_DIRNAME):
        raise ValueError(
            'Leading import directory "{}" not isolated to '
            'venv directory "{}".'.format(import_first_dirname, VENV_DIRNAME))

    # BETSE, imported *AFTER* performing sanity checks above.
    import betse

    # Absolute dirname of the directory containing the top-level
    # "betse.__init__" submodule.
    BETSE_DIRNAME = os.path.dirname(betse.__file__)

    # If this dirname is *NOT* suffixed by a directory separator, do so.
    if BETSE_DIRNAME[-1] != os.path.sep:
        BETSE_DIRNAME += os.path.sep

    # If this directory is *NOT* isolated to this venv, raise an exception.
    if not BETSE_DIRNAME.startswith(VENV_DIRNAME):
        raise ValueError(
            'Code directory "{}" not isolated to '
            'venv directory "{}".'.format(BETSE_DIRNAME, VENV_DIRNAME))


def _print_metadata() -> None:
    '''
    Print test-specific metadata for debuggability and quality assurance (QA).
    '''

    # Print a header for disambiguity.
    print('------[ paths ]------')

    # Print the absolute dirname of the system-wide Python prefix and
    # current Python prefix, which differs from the former under venvs.
    print('python prefix (system [base]): ' + sys.base_prefix)
    print('python prefix (system [real]): ' + getattr(sys, 'real_prefix', ''))
    print('python prefix (current): ' + sys.prefix)

    # Print the current list of the (absolute or relative) dirnames of all
    # directories to be iteratively searched for importable modules and
    # packages, initialized from the "${PYTHONPATH}" environment variable and
    # subsequently extended by pytest. Since Python searches this list in
    # descending order, directories listed earlier assume precedence over
    # directories listed later.
    print('import paths: ' + str(sys.path))

    # Defer heavyweight imports until *AFTER* printing the above metadata.
    import betse
    from betse.util.py.module import pymodule

    # Print the absolute dirname of the top-level "betse" package.
    print('project path: ' + pymodule.get_dirname_canonical(betse))

    # Print all imported module names for debugging purposes.
    # from betse.util.py.module import pyimport
    # print(
    #     '------[ imported modules ]------\n' +
    #     pyimport.to_str_modules_imported_name())

    # Print all environment variables for debugging purposes.
    # from betse.util.os.shell import shellenv
    # print('------[ environment variables ]------\n' + shellenv.to_str())

    # Print a POSIX-specific process tree for debugging purposes.
    # from betse.util.path.command import cmdrun
    # cmdrun.run_or_die(command_words=('pstree',))

    # # Print all attributes of the "site" module for debugging purposes.
    # import site
    # from betse.util.type.obj import objiter
    # print('------[ "site" attributes ]------')
    # for attr_name, attr_value in objiter.iter_attrs(site):
    #     print('{}: {}'.format(attr_name, attr_value))
    #
    # # Print all attributes of the "sys" module for debugging purposes.
    # print('------[ "sys" attributes ]------')
    # for attr_name, attr_value in objiter.iter_attrs(sys):
    #     print('{}: {}'.format(attr_name, attr_value))

# ....................{ HOOKS ~ session : stop            }....................
def pytest_sessionfinish(session, exitstatus) -> None:
    '''
    Hook run immediately *after* completing the current test session (i.e.,
    calling the :func:`pytest.session.main` function).
    '''

    pass
