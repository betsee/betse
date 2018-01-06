#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Utility functions simplifying access to various :mod:`pytest` subpackages,
submodules, and attributes unclassifiable with a more specific designation.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import ClassType, ModuleType

# ....................{ GETTERS                            }....................
def get_pytest_fixtures_submodule() -> ModuleType:
    '''
    Version-specific private :mod:`pytest` submodule providing fixture classes.

    Although private, this submodule publishes classes of public interest.
    Unfortunately, as this submodule *is* private, its fully-qualified name has
    changed across :mod:`pytest` versions. This function guaranteeably retrieves
    this submodule regardless of version.

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

# ....................{ GETTERS ~ type                     }....................
def get_fixture_lookup_error_type() -> ClassType:
    '''
    Class of all :class:`FixtureLookupError` exceptions raised by :mod:`pytest`.

    Specifically, under:

    * :mod:`pytest` >= 3.0.0, the private
      :class:`_pytest.fixtures.FixtureLookupError` class is imported and
      returned.
    * :mod:`pytest` < 3.0.0, the private
      :class:`_pytest.python.FixtureLookupError` class is imported and returned.
    '''

    # Version-specific private `pytest` submodule providing fixture classes.
    pytest_fixtures_submodule = get_pytest_fixtures_submodule()

    # Return this exception class defined by this submodule.
    return pytest_fixtures_submodule.FixtureLookupError
