#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Autouse fixtures** (i.e., fixtures unconditionally applicable to all tests
within a given scope of the current test session).
'''

# ....................{ IMPORTS                           }....................
from betse.util.type.types import GeneratorType
from contextlib import contextmanager
from pytest import fixture

# ....................{ FIXTURES ~ package                }....................
@fixture(scope='package', autouse=True)
def betse_init_package() -> None:
    '''
    **Autouse package initialization fixture** (i.e., fixture unconditionally
    applicable to all tests within the test package importing this fixture,
    which :mod:`pytest` automatically invokes once at session startup *after*
    collecting these tests but *before* running these tests).

    This fixture initializes the application metadata singleton -- a mandatory
    prerequisite before exercising *any* application functionality. Since
    functional (but *not* unit) tests already implicitly initialize this
    singleton, this fixture is typically only imported by the top-level unit
    testing subpackage (e.g., :mod:`betse_test.unit`).

    See Also
    ----------
    :func:`_run_app_tests`
        Further details.
    '''

    # Make it so... or else.
    yield from _run_app_tests()

# ....................{ FIXTURES ~ session                }....................
#FIXME: This fixture is currently unused until we inevitably require this yet
#again. Ergo, this should be preserved for the sad future that comes.
# @fixture(scope='session', autouse=True)
# def betse_init_session() -> None:
#     '''
#     **Autouse session initialization fixture** (i.e., fixture unconditionally
#     applicable to all tests for the current test session, which :mod:`pytest`
#     automatically invokes once at session startup *after* collecting these
#     tests but *before* running these tests).
#
#     This fixture initializes the application metadata singleton -- a mandatory
#     prerequisite before exercising *any* application functionality. Since a
#     proper subset of testing fixtures and decorators (e.g., the
#     :func:`betse.util.test.pytest.mark.pytskip.skip_if_requirement` decorator)
#     assume this singleton to have already been initialized, this fixture is
#     typically imported by the top-level plugin (e.g.,
#     :mod:`betse_test.conftest`).
#     '''
#
#     pass

# ....................{ PRIVATE ~ runners                 }....................
def _run_app_tests() -> GeneratorType:
    '''
    Generator exercising all application-specific tests scoped by the calling
    fixture by ensuring that the application metadata singleton is initialized
    *before* the first such test and deinitialized *after* the last such test.

    Specifically, this generator (in order):

    #. Initializes the application metadata singleton.
    #. Yields control to the caller.
    #. Deinitializes this singleton, regardless of whether the caller raised an
       exception or not.

    See Also
    ----------
    :func:`init_app`
    :func:`deinit_app`
        Further details.
    '''

    # Initialize the application metadata singleton.
    init_app()

    # Yield control to the body of the caller's "with" block.
    try:
        yield
    # Deinitialize this singleton even if that block raised an exception.
    finally:
        deinit_app()

# ....................{ [DE]INITIALIZERS                  }....................
def init_app() -> None:
    '''
    Initialize the application metadata singleton (e.g., for unit tests, which
    fail to implicitly initialize this singleton ala functional tests).

    Specifically, this fixture (in order):

    #. Coerces the active Python interpreter into running **headless** (i.e.,
       with *no* access to a GUI display). Allowing headfull operation would
       would allow tests erroneously attempting to connect to an X11 server to
       locally succeed but remotely fail, as headless continuous integration
       (CI) pipelines typically have no access to an X11 server. Coercing
       headlessness ensures orthogonality between these cases by coercing the
       former to fail as well.
    #. Initializes the application core.
    #. Initializes all third-party dependencies thereof.
    #. Enables the default non-interactive matplotlib backend ``Agg``,
       *guaranteed* to be usable on all platforms. By default, matplotlib
       enables an interactive backend (e.g., ``Qt5Agg``) unsuitable for use
       under typically headless test automation.

    Motivation
    ----------
    This function performs early test-specific initialization of this
    application and dependencies thereof. Doing so prevents the magic
    :func:`betse.science.__init__` function from attempting to perform a
    subsequent test-agnostic initialization of either this application or
    dependencies on the first importation of the :mod:`betse.science`
    subpackage -- as in fixtures importing from that subpackage (e.g., the
    :mod:`betse_test.fixture.simconf.simconfer` fixture importing the
    :mod:`betse_test.fixture.simconf.simconfwrapper` submodule importing the
    :med:`betse.science.config.confwrap` submodule).
    '''

    # Defer heavyweight imports.
    from betse.util.app.meta import appmetaone
    from betse.util.os import displays
    from betse.util.test.pytest import pytests

    # If an application metadata singleton already exists, reduce to a noop.
    if appmetaone.is_app_meta():
        return
    # Else, no application metadata singleton currently exists.

    # Prepend a leading newline, which py.test curiously neglects to do itself.
    print('\n')

    # Print this initialization.
    pytests.output('Initializing application metadata singleton...')

    # Initialize a BETSE-specific application metadata singleton if the
    # appmetaone.set_app_meta() function has yet to be called.
    app_meta = appmetaone.set_app_meta_betse_if_unset()

    # Coerce the active Python interpreter into running headless *AFTER*
    # initializing this singleton, which enables the default logging
    # configuration to which this setter logs this operation.
    #
    # Note that this operation technically needs to be performed:
    #
    # * Only once for the entire test suite when py.test is *NOT* parallelized
    #   with "xdist", in which case all tests run in the same process and hence
    #   share the same global variables.
    # * Once for each test when py.test is parallelized with "xdist", in which
    #   case each test is run in a distinct subprocess and hence does *NOT*
    #   share the same global variables.
    #
    # Since setting global variables is fast, doing so here transparently
    # supports both use cases detailed above with no discernable downside. See
    # the docstring for additional commentary.
    displays.set_headless(True)

    # Initialize all mandatory third-party dependencies with a standard
    # non-interactive matplotlib backend guaranteed to exist *AFTER* coercing
    # the active Python interpreter into running headless. Why? Because
    # dependencies typically detect headless environments.
    app_meta.init_libs(matplotlib_backend_name='Agg')

    # Print the completion of this initialization.
    pytests.output('Initialized application metadata singleton.')


def deinit_app() -> None:
    '''
    Deinitialize the application metadata singleton.
    '''

    # Defer heavyweight imports.
    from betse.util.app.meta import appmetaone
    from betse.util.test.pytest import pytests

    # If no application metadata singleton currently exists, reduce to a noop.
    if not appmetaone.is_app_meta():
        return
    # Else, an application metadata singleton currently exists.

    # Prepend a leading newline, which py.test curiously neglects to do itself.
    print('\n')

    # Print this deinitialization.
    pytests.output('Deinitializing application metadata singleton...')

    # Deinitialize the application metadata singleton.
    appmetaone.deinit()

    # Print this deinitialization.
    pytests.output('Deinitialized application metadata singleton.')
