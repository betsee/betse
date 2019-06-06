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

#FIXME: This fixture is currently unused until we inevitably require this yet
#again. Ergo, this should be preserved for the sad future that comes. See also
#the betse_test.conftest._init_app() function for similar functionality.

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

    Specifically, this fixture...

    Motivation
    ----------
    This fixture is automatically requested by *all* tests (functional, unit,
    or otherwise) without needing to be explicitly requested. Moreover, this
    fixture is imported into the top-level :mod:`betse_test.conftest` plugin
    and hence *guaranteed* to be run prior to all other application-specific
    fixtures. Doing so avoids spurious issues in other fixtures and tests.
    '''

    pass
