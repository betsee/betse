#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :func:`betse.util.py.pymodule.is_c_extension` function.
'''

# ....................{ IMPORTS                           }....................
# from pytest import fixture
# import pytest, sys

# ....................{ TESTS                             }....................
def test_is_c_extension_unmonkeypatched() -> None:
    '''
    Unit test the :func:`betse.util.py.pymodule.is_c_extension` function
    *without* monkeypatching the PEP 302-specific `__loader__` attribute of
    tested modules.
    '''

    # Imports deferred for safety.
    from betse.util.py import pymodule
    from numpy.core import multiarray

    # Ensure that a C extension is correctly detected as such. Testing this
    # particular C extension is critical, as the core codebase explicitly
    # performs the same call at application startup.
    assert pymodule.is_c_extension(multiarray) is True

    # Ensure that a pure-Python submodule is *NOT* detected to be a C extension.
    assert pymodule.is_c_extension(pymodule) is False


def test_is_c_extension_monkeypatched(monkeypatch) -> None:
    '''
    Test all calls of the :func:`betse.util.py.pymodule.is_c_extension`
    function by monkeypatching the PEP 302-specific `__loader__` attribute of
    tested modules.

    This monkeypatch improves code coverage by exercising an infrequently used
    fallback codepath in this function.

    Parameters
    ----------
    monkeypatch : MonkeyPatch
        Builtin fixture object permitting object attributes to be safely
        modified for the duration of this unit test.
    '''

    # Imports deferred for safety.
    from betse.util.py import pymodule
    from numpy.core import multiarray

    # Remove the the PEP 302-specific "__loader__" attribute of each such
    # module *BEFORE* performing testing. Since this attribute need *NOT* exist,
    # prevent exceptions from being raised when this attribute does *NOT* exist.
    for module_object in (pymodule, multiarray):
        monkeypatch.delattr(module_object, name='__loader__', raising=False)

    # Defer to the existing base unit test.
    test_is_c_extension_unmonkeypatched()
