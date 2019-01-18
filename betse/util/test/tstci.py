#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level **testing-specific continuous integration** (i.e., remote host
exercising all functional and unit tests within this application's test suite
across one or more Python versions and one or more platforms) facilities.

See Also
----------
https://media.readthedocs.org/pdf/python-semantic-release/latest/python-semantic-release.pdf
    Section 3.2.1, *Environment checks,* from the official documentation for
    the completely unrelated albeit significantly more popular
    ``python-semantic-version`` project. This section reliably documents the
    set of all environment variables uniquely identifying continuous
    integration (CI) hosts, thus serving as the principal inspiration for this
    whole submodule. Thanks for all the CI, ``python-semantic-version``.
'''

# ....................{ IMPORTS                           }....................
# from betse.util.io.log import logs
from betse.util.type.types import type_check

# ....................{ TESTERS                           }....................
def is_ci_circle() -> bool:
    '''
    ``True`` only if the active Python interpreter is currently running tests
    under the remote CircleCI continuous integration (CI) host.
    '''

    return _is_ci_host(env_var_name='CIRCLECI')


def is_ci_frigg() -> bool:
    '''
    ``True`` only if the active Python interpreter is currently running tests
    under the remote Frigg continuous integration (CI) host.
    '''

    return _is_ci_host(env_var_name='FRIGG')


def is_ci_gitlab() -> bool:
    '''
    ``True`` only if the active Python interpreter is currently running tests
    under a remote GitLab-CI continuous integration (CI) host.
    '''

    return _is_ci_host(env_var_name='GITLAB_CI')


def is_ci_semaphore() -> bool:
    '''
    ``True`` only if the active Python interpreter is currently running tests
    under the remote Semphare continuous integration (CI) host.
    '''

    return _is_ci_host(env_var_name='SEMAPHORE')


def is_ci_travis() -> bool:
    '''
    ``True`` only if the active Python interpreter is currently running tests
    under the remote Travis-CI continuous integration (CI) host.
    '''

    return _is_ci_host(env_var_name='TRAVIS')

# ....................{ TESTERS ~ private                 }....................
@type_check
def _is_ci_host(env_var_name: str) -> bool:
    '''
    ``True`` only if the active Python interpreter is currently running tests
    (e.g., with the :mod:`pytest` test harness) *and* an environment variable
    with the passed name uniquely associated with one remote continuous
    integration (CI) host exists in the shell environment encapsulating this
    interpreter, implying this host to be running this interpreter and hence
    these tests.

    Parameters
    ----------
    env_var_name : str
        Name of the environment variable uniquely identifying this interpreter
        to be running under one remote continuous integration (CI) host.

    Returns
    ----------
    bool
        ``True`` only if this interpreter is running tests under the CI host
        uniquely identified by this variable.
    '''

    # Avoid circular import dependencies.
    from betse.util.os.shell import shellenv

    # Return true only if this interpreter is running tests *AND* an
    # environment variable with the passed name exists, implying these tests to
    # be run under the CI host uniquely identified with that variable.
    #
    # Ideally, this implementation would also call the
    # betse.util.test.tests.is_calling() function to ensure that tests are
    # running as expected. However, that function is *NOT* safely callable from
    # module scope. Calling that function here would thus prevent this function
    # and functions calling this function from also being safely callable from
    # module scope -- which, in turn, would prevent test decorators (e.g.,
    # @betse.util.test.pytest.mark.pytskip.skip_if_ci_gitlab) from being safely
    # used... *AT ALL.* Since that would rather defeat the purpose of defining
    # any of this functionality in the first place, we have little choice but
    # to defer to a simple environment variable existence test.
    return shellenv.is_var(env_var_name)
