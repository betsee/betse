#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Fixtures initializing BETSE in a manner suitable for _all_ testing.
'''

# ....................{ IMPORTS                            }....................
from pytest import fixture

# ....................{ FIXTURES                           }....................
@fixture(scope='session', autouse=True)
def betse_init() -> None:
    '''
    Automatically run per-session fixture initializing BETSE in a manner
    suitable for unit testing.

    Specifically, this fixture:

    * Initializes the BETSE core.
    * Initializes all third-party dependencies thereof.
    * Enables the default non-interactive matplotlib backend `Agg`, _guaranteed_
      to be usable on all platforms. By default, matplotlib enables an
      interactive backend (e.g., `Qt5Agg`) unsuitable for use under typically
      headless test automation.

    Motivation
    ----------
    This fixture is automatically requested by _all_ tests (functional, unit,
    or otherwise) without needing to be explicitly requested. Moreover, this
    fixture is imported into the top-level `betse_test.conftest` plugin and
    hence _guaranteed_ to be run prior to all other BETSE-specific fixtures.
    Doing so avoids spurious issues in other fixtures and tests.

    Notably, the early test-specific initialization of both BETSE and
    dependencies prevents the magic :func:`betse.science._ignite` function from
    attempting to perform a subsequent test-agnostic initialization of either
    BETSE or dependencies on the first importation of the :mod:`betse.science`
    subpackage, as occurs in fixtures importing from that subpackage (e.g., the
    :mod:`betse_test.fixture.simconfig.simconfer` fixture importing the
    :mod:`betse_test.fixture.simconfig.simconfwrapper` submodule importing the
    :med:`betse.science.config.confwrap` submodule).
    '''

    # Defer heavyweight imports.
    from betse import ignition
    from betse.lib import libs

    # Inform users of this initialization.
    print('\n[py.test] Initializing BETSE...')

    # Initialize the core application. Note that the higher-level
    # ignition.ignote() function is intentionally *NOT* called here, as doing
    # so could erroneously enable a headfull matplotlib backend.
    ignition.init()

    # Initialize all dependencies *AFTER* the core application.
    libs.init(matplotlib_backend_name='Agg')
