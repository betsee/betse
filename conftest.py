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
#FIXME: Strip all non-"tox"-isolated directories from "sys.path" when running
#under "tox". To do so, see the informative question at:
#    https://stackoverflow.com/questions/55737714/how-does-a-tox-environment-set-its-sys-path
#FIXME: Raise an exception if running under "tox" but the project path
#(i.e., "pymodule.get_dirname_canonical(betse)") is *NOT* isolated to a
#"tox"-isolated venv.
def pytest_sessionstart(session: '_pytest.main.Session') -> None:
    '''
    Hook run immediately *before* starting the current test session (i.e.,
    calling the :func:`pytest.session.main` function).

    Parameters
    ----------
    session: _pytest.main.Session
        :mod:`pytest`-specific test session object.
    '''

    # Sanitize import paths *BEFORE* the first module importation.
    _clean_imports()

    # Print test-specific metadata *AFTER* sanitizing import paths.
    _print_metadata()


def _clean_imports() -> None:
    '''
    Sanitize import directories (i.e., the global :attr:`sys.list` of the
    absolute and relative dirnames of all directories to search for modules and
    packages to be imported from).

    Specifically, this function:

    * If this low-level :mod:`pytest` test harness is isolated to a venv (e.g.,
      due to being exercised by a higher-level :mod:`tox` wrapper), remove all
      import directories *not* isolated to this venv from the global list of
      all import directories (i.e., :attr:`sys.path`). Doing so prevents this
      test session from accidentally importing from modules and packages *not*
      isolated to this venv, including this application being tested.
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

    # If tests are isolated to a venv...
    if is_venv:
        # Print the absolute dirname of this venv's top-level directory.
        print('venv dir: {}'.format(sys.prefix))

        # List of the absolute dirnames of all directories to search for
        # modules and packages to be imported from, guaranteed to be isolated
        # to this venv.
        sys_path_new = []

        # Absolute dirname of this venv's top-level directory, suffixed by a
        # directory separator for disambiguity.
        sys_prefix = sys.prefix + os.path.sep

        # For the dirname of each directory to search for imports...
        for import_dirname in sys.path:
            # If this directory resides inside this venv, preserve this
            # directory in this list as is.
            if import_dirname.startswith(sys_prefix):
                sys_path_new.append(import_dirname)
            # Else, this directory resides outside this venv. In this case,
            # narn that this directory will *NOT* be importable from.
            else:
                print(
                    'WARNING: '
                    'Ignoring non-isolated import directory "{}"...'.format(
                        import_dirname),
                    file=sys.stderr)

        # Replace the original such list with this redacted list.
        sys.path = sys_path_new


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
