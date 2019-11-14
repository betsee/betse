#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Low-level :mod:`pytest` functionality.

This submodule principally defines utility functions simplifying access to
:mod:`pytest` subpackages, submodules, and attributes unclassifiable under a
more specific submodule.
'''

# ....................{ IMPORTS                           }....................
from betse.util.type.types import ClassType, ModuleType

# ....................{ GETTERS                           }....................
def get_pytest_fixtures_submodule() -> ModuleType:
    '''
    Version-specific private :mod:`pytest` submodule providing fixture classes.

    Although private, this submodule publishes classes of public interest.
    Unfortunately, as this submodule *is* private, its fully-qualified name has
    changed across :mod:`pytest` versions. This function guaranteeably
    retrieves this submodule regardless of version.

    Specifically, under:

    * :mod:`pytest` >= 3.0.0, the private :mod:`_pytest.fixtures` submodule is
      imported and returned.
    * :mod:`pytest` < 3.0.0, the private :mod:`_pytest.python` submodule is
      imported and returned.
    '''

    # Attempt to import and return the newer "pytest" fixture module.
    try:
        from _pytest import fixtures as pytest_fixtures
    # Failing that, fallback to importing and returning the older such module.
    except ImportError:
        from _pytest import python as pytest_fixtures

    # In either case, return this module.
    return pytest_fixtures


def get_pytest_fixture_lookup_error_type() -> ClassType:
    '''
    Type of all :class:`FixtureLookupError` exceptions raised by :mod:`pytest`.

    Specifically, under:

    * :mod:`pytest` >= 3.0.0, the private
      :class:`_pytest.fixtures.FixtureLookupError` class is imported and
      returned.
    * :mod:`pytest` < 3.0.0, the private
      :class:`_pytest.python.FixtureLookupError` class is imported and
      returned.
    '''

    # Version-specific private `pytest` submodule providing fixture classes.
    pytest_fixtures_submodule = get_pytest_fixtures_submodule()

    # Return this exception class defined by this submodule.
    return pytest_fixtures_submodule.FixtureLookupError

# ....................{ OUTPUTTERS                        }....................
def output(*objs) -> None:
    '''
    Print all passed objects as is to standard output in a format mimicking
    that of standard :mod:`pytest` messages *without* logging these objects.

    This function is intended to be called *only* by :mod:`pytest`-specific
    fixtures, decorators, and helpers.

    This function is intentionally *not* named ``print()`` to avoid conflict
    with the builtin function of the same name.

    Examples
    ----------
        >>> from betse.util.test.pytest import pytests
        >>> pytests.output('Ego, ergo simulare.')
        [py.test] Ego, ergo simulare.
    '''

    # Avoid circular import dependencies.
    from betse.util.os.command import cmds

    # Print these messages in a "py.test"-friendly format.
    print('[{}] {}'.format(cmds.get_current_basename(), ''.join(objs)))
