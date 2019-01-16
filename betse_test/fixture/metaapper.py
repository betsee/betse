#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Fixtures initializing the application metadata singleton in a general-purpose
manner suitable for *all* functional and unit tests.
'''

# ....................{ IMPORTS                           }....................
from pytest import fixture

# ....................{ FIXTURES                          }....................
@fixture(scope='session', autouse=True)
def betse_app_meta() -> 'betse.util.app.meta.metaappabc.MetaAppABC':
    '''
    Automatically run per-session fixture instantiating and initializing the
    application metadata singleton in a manner suitable for unit testing,
    returning this singleton for caller convenience.

    Specifically, this fixture:

    * Initializes the BETSE core.
    * Initializes all third-party dependencies thereof.
    * Enables the default non-interactive matplotlib backend ``Agg``,
      *guaranteed* to be usable on all platforms. By default, matplotlib
      enables an interactive backend (e.g., ``Qt5Agg``) unsuitable for use
      under typically headless test automation.

    Motivation
    ----------
    This fixture is automatically requested by *all* tests (functional, unit,
    or otherwise) without needing to be explicitly requested. Moreover, this
    fixture is imported into the top-level :mod:`betse_test.conftest` plugin
    and hence *guaranteed* to be run prior to all other BETSE-specific
    fixtures. Doing so avoids spurious issues in other fixtures and tests.

    Notably, the early test-specific initialization of both BETSE and
    dependencies prevents the magic :func:`betse.science._ignite` function from
    attempting to perform a subsequent test-agnostic initialization of either
    BETSE or dependencies on the first importation of the :mod:`betse.science`
    subpackage, as occurs in fixtures importing from that subpackage (e.g., the
    :mod:`betse_test.fixture.simconf.simconfer` fixture importing the
    :mod:`betse_test.fixture.simconf.simconfwrapper` submodule importing the
    :med:`betse.science.config.confwrap` submodule).

    Returns
    ----------
    betse.util.app.meta.metaappabc.MetaAppABC
        Application metadata singleton.
    '''

    # Defer heavyweight imports.
    from betse.util.app.meta import metaappton

    # Inform callers of this initialization.
    print('\n[py.test] Initializing BETSE for testing...')

    # Instantiate and set a BETSE-specific application metadata singleton if
    # the metaappton.set_app_meta() function has yet to be called *AND*, in
    # either case, initialize all mandatory third-party dependencies with a
    # standard non-interactive matplotlib backend guaranteed to exist.
    app_meta = metaappton.make_app_meta_betse()
    app_meta.init_libs(matplotlib_backend_name='Agg')

    # Inform callers of the completion of this initialization.
    print('[py.test] Initialized BETSE for testing.')

    # Return this singleton for caller convenience.
    return app_meta
