#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising \*ML (e.g., HTML, SGML, XML) functionality.
'''

# ....................{ IMPORTS                            }....................
# import pytest

# ....................{ TESTS                              }....................
def test_is_ml() -> None:
    '''
    Unit test the :func:`betse.util.type.text.mls.is_ml` tester.
    '''

    # Defer heavyweight imports.
    from betse.util.type.text import mls

    # Assert this tester to behave as expected.
    assert mls.is_ml(
        '<nausicaa kaze no="tane">And that one shall come to you</nausicaa>')
    assert not mls.is_ml(
        '<garbed in "raiment" of blue and descending upon a field of gold>...')
