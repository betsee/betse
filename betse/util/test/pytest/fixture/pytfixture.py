#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level :mod:`pytest` fixtures.
'''

# ....................{ IMPORTS                           }....................
from pytest import fixture

# ....................{ FIXTURES                          }....................
@fixture(scope='session')
def monkeypatch_session() -> '_pytest.monkeypatch.MonkeyPatch':
    '''
    Session-scoped ``monkeypatch`` fixture.

    Caveats
    ----------
    **The ``monkeypatch`` fixture should never be referenced from
    session-scoped fixtures.** Doing so raises the following fatal error on
    session startup:

        ScopeMismatch: You tried to access the 'function' scoped fixture
        'monkeypatch' with a 'session' scoped request object, involved
        factories

    This equivalent fixture suffers no such deficits and hence should be
    referenced from session-scoped fixtures instead.

    Returns
    ----------
    _pytest.monkeypatch.MonkeyPatch
        Object returned by the `monkeypatch` fixture keeping a record of
        :func:`setattr`, item, environment, and :attr:`sys.path` changes.

    See Also
    ----------
    https://github.com/pytest-dev/pytest/issues/363#issuecomment-406536200
        GitHub comment strongly inspiring this implementation.
    '''

    # Defer violations of privacy encapsulation for safety. Tragically, this
    # class is inaccessible via any other means.
    from _pytest.monkeypatch import MonkeyPatch

    # Create and yield to the calling session-scoped fixture a new
    # monkey-patcher instance.
    monkey_patch = MonkeyPatch()
    yield monkey_patch

    # Revert all temporary changes applied to this instance by the calling
    # session-scoped fixture.
    monkey_patch.undo()
