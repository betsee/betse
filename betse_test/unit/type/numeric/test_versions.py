#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Unit tests exercising the :mod`betse.util.type.numeric.versions` submodule.

Unit testing this submodule as critical, both because of its centrality to
sane dependency resolution *and* because of the inability of the
:mod:`setuptools`-bundled :mod:`pkg_resources` package underlying this
submodule to maintain adequate (or, *any*) backward compatibility guarantees.
'''

# ....................{ IMPORTS                           }....................
# import pytest

# ....................{ TESTS ~ base 10                   }....................
def test_version_comparators() -> None:
    '''
    Unit test the incestous family of
    ``betse.util.type.numeric.version.is_*_than*`` comparators (e.g.,
    :func:`betse.util.type.numeric.version.is_less_than_or_equal_to`).
    '''

    # Defer heavyweight imports.
    from betse.util.type.numeric import versions

    # Assert each comparator to behave as expected on common edge cases.
    assert versions.is_greater_than('0.6.43.410', (0, 6, 43, 409))
    assert versions.is_greater_than_or_equal_to(('3', '300', '33'), '3.300.33')
    assert versions.is_less_than(('0', '412', '454'), (0, 412, 455))
    assert versions.is_less_than_or_equal_to('23.140.692.63', '23.140.692.63')


def test_version_converters() -> None:
    '''
    Unit test the :func:`betse.util.type.numeric.version.to_comparable`
    function.
    '''

    # Defer heavyweight imports.
    from betse.util.type.numeric import versions

    # Assert this converter to behave as expected on common edge cases --
    # notably, pre-release (i.e., "alpha candidate") version numbers for the
    # PySide2 dependency required by the BETSEE downstream consumer. For
    # unknown reasons:
    #
    # * These numbers are suffixed by "~"-prefixed labels violating PEP 440
    #   compliance. The pkg_resources package underlying this converter creates
    #   and returns a unique class of objects to represent such versions.
    # * These objects are sorted incorrectly. Versions containing tilde
    #   characters are incorrectly sorted as strictly less than versions *NOT*
    #   containing tilde characters (e.g., "5.9.0~a1" < "5.7.0").
    version_bad_new = versions.to_comparable('5.9.0~a1')
    version_bad_old = versions.to_comparable('5.6.0~a1')

    # Assert these uncompliant versions to be sorted in the correct manner.
    assert not versions.is_less_than(version_bad_new, '5.7.0')
    assert     versions.is_less_than(version_bad_old, '5.7.0')
