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

# ....................{ HOOKS ~ session                   }....................
#FIXME: Strip all non-"tox"-isolated directories from "sys.path" when running
#under "tox". To do so, see the informative question at:
#    https://stackoverflow.com/questions/55737714/how-does-a-tox-environment-set-its-sys-path
#FIXME: Raise an exception if running under "tox" but the project path
#(i.e., "pymodule.get_dirname_canonical(betse)") is *NOT* isolated to a
#"tox"-isolated venv.
def pytest_sessionstart(session):
    '''
    Hook run immediately *before* starting the current test session (i.e.,
    calling the :func:`pytest.session.main` function).
    '''

    # Defer heavyweight imports.
    import betse, sys
    from betse.util.py.module import pymodule

    # Print a header for disambiguity.
    print('------[ python paths ]------')

    # Print the absolute dirname of the system-wide Python prefix and
    # current Python prefix, which differs from the former under venvs.
    print('python prefix (system):  ' + sys.base_prefix)
    print('python prefix (current): ' + sys.prefix)

    # Print the absolute dirname of the top-level "betse" package.
    print('project path: ' + pymodule.get_dirname_canonical(betse))

    # Print the current list of the (absolute or relative) dirnames of all
    # directories to be iteratively searched for importable modules and
    # packages, initialized from the "${PYTHONPATH}" environment variable and
    # subsequently extended by pytest. Since Python searches this list in
    # descending order, directories listed earlier assume precedence over
    # directories listed later.
    print('import paths: ' + str(sys.path))

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


def pytest_sessionfinish(session, exitstatus):
    '''
    Hook run immediately *after* completing the current test session (i.e.,
    calling the :func:`pytest.session.main` function).
    '''

    pass
