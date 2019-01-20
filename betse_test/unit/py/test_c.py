#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests for the :func:`betse.util.py.pymodule.is_c` function.
'''

# ....................{ IMPORTS                           }....................

# ....................{ TESTS                             }....................
def test_is_c_unmonkeypatched() -> None:
    '''
    Unit test the :func:`betse.util.py.module.pymodule.is_c` function *without*
    monkeypatching the special PEP 302-specific ``__loader__`` attribute of
    tested modules.
    '''

    # Defer heavyweight imports.
    from betse.lib.numpy import numpys
    from betse.util.py.module import pymodule

    # Arbitrary Numpy submodule guaranteed to be implemented as a C extension.
    # Note that testing this particular C module is critical, as the core
    # codebase explicitly performs the same call at application startup.
    numpy_c_extension = numpys.get_c_extension()

    # Ensure that a C module is correctly detected as such.
    assert pymodule.is_c(numpy_c_extension) is True

    # Ensure that a pure-Python submodule is *NOT* detected to be a C
    # extension.
    assert pymodule.is_c(pymodule) is False


def test_is_c_monkeypatched(monkeypatch) -> None:
    '''
    Test all calls of the :func:`betse.util.py.module.pymodule.is_c` function
    by monkeypatching the special PEP 302-specific ``__loader__`` attribute of
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
    from betse.lib.numpy import numpys
    from betse.util.py.module import pymodule

    # Arbitrary Numpy submodule guaranteed to be implemented as a C extension.
    numpy_c_extension = numpys.get_c_extension()

    # Remove the the PEP 302-specific "__loader__" attribute of each such
    # module before performing testing. Since this attribute need *NOT* exist,
    # prevent exceptions from being raised if this attribute does *NOT* exist.
    for module_object in (pymodule, numpy_c_extension):
        monkeypatch.delattr(module_object, name='__loader__', raising=False)

    # Defer to the existing base unit test.
    test_is_c_unmonkeypatched()
