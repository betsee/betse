#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Autouse session fixtures** (i.e., fixtures generically applicable to the
current :mod:`pytest` session as a whole, which :mod:`pytest` automatically
invokes exactly once at session startup *after* collecting all functional and
unit tests to be subsequently run under this session).
'''

# ....................{ IMPORTS                           }....................
from pytest import fixture

# ....................{ FIXTURES                          }....................
@fixture(scope='session', autouse=True)
def betse_autouse() -> None:
    '''
    **Autouse session fixture** (i.e., fixture generically applicable to the
    current :mod:`pytest` session as a whole, which :mod:`pytest` automatically
    invokes exactly once at session startup *after* collecting all functional
    and unit tests to be subsequently run under this session).

    Specifically, this fixture:

    * Temporarily unsets the external ``${DISPLAY}`` environment variable if
      currently set (e.g., to the X11-specific socket to be connected to
      display GUI components) for the duration of this session. Allowing this
      variable to remain set would allow tests erroneously attempting to
      connect to an X11 server to locally succeed but remotely fail. Why?
      Because headless continuous integration (CI) typically has no access to
      an X11 server. Unsetting this variable ensures orthogonality between
      these cases by coercing the former to fail as well.
    * Instantiates and initializes the application metadata singleton in a
      manner suitable for unit testing. Specifically (in order):

      #. Coerces the active Python interpreter into running **headless** (i.e.,
         with *no* access to a GUI display). Allowing headfull operation would
         would allow tests erroneously attempting to connect to an X11 server
         to locally succeed but remotely fail, as headless continuous
         integration (CI) pipelines typically have no access to an X11 server.
         Coercing headlessness ensures orthogonality between these cases by
         coercing the former to fail as well.
      #. Initializes the BETSE core.
      #. Initializes all third-party dependencies thereof.
      #. Enables the default non-interactive matplotlib backend ``Agg``,
         *guaranteed* to be usable on all platforms. By default, matplotlib
         enables an interactive backend (e.g., ``Qt5Agg``) unsuitable for use
         under typically headless test automation.

    Motivation
    ----------
    This fixture is automatically requested by *all* tests (functional, unit,
    or otherwise) without needing to be explicitly requested. Moreover, this
    fixture is imported into the top-level :mod:`betse_test.conftest` plugin
    and hence *guaranteed* to be run prior to all other application-specific
    fixtures. Doing so avoids spurious issues in other fixtures and tests.

    Notably, the early test-specific initialization of both BETSE and
    dependencies prevents the magic :func:`betse.science._ignite` function from
    attempting to perform a subsequent test-agnostic initialization of either
    BETSE or dependencies on the first importation of the :mod:`betse.science`
    subpackage, as occurs in fixtures importing from that subpackage (e.g., the
    :mod:`betse_test.fixture.simconf.simconfer` fixture importing the
    :mod:`betse_test.fixture.simconf.simconfwrapper` submodule importing the
    :med:`betse.science.config.confwrap` submodule).
    '''

    # Defer heavyweight imports.
    from betse.util.os import displays
    from betse.util.app.meta import metaappton

    # Prepend a leading newline, which py.test curiously neglects to do itself.
    print('\n')

    # Inform callers of application initialization.
    print('[py.test] Initializing BETSE for testing...')

    # Instantiate and set a BETSE-specific application metadata singleton if
    # the metaappton.set_app_meta() function has yet to be called.
    app_meta = metaappton.make_app_meta_betse()

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

    # Inform callers of the completion of this initialization.
    print('[py.test] Initialized BETSE for testing.')
