#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Utility functions simplifying access to various py.test subpackages, submodules,
and attributes unclassifiable with a more specific designation.
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.types import ClassType, ModuleType

# ....................{ GETTERS                            }....................
def get_pytest_fixtures_submodule() -> ModuleType:
    '''
    Version-specific private `pytest` submodule providing fixture classes.

    Although private, this submodule publishes classes of public interest.
    Unfortunately, as this submodule _is_ private, its fully-qualified name has
    changed across `pytest` versions. This function guaranteeably retrieves this
    submodule regardless of version.

    Specifically, under:

    * `pytest` >= 3.0.0, the private `_pytest.fixtures` submodule is imported
      and returned.
    * `pytest` < 3.0.0, the private `_pytest.python` submodule is imported and
      returned.
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
    Class of all `FixtureLookupError` exceptions raised by `pytest`.

    Specifically, under:

    * `pytest` >= 3.0.0, the private `_pytest.fixtures.FixtureLookupError` class
      is imported and returned.
    * `pytest` < 3.0.0, the private `_pytest.python.FixtureLookupError` class is
      imported and returned.
    '''

    # Version-specific private `pytest` submodule providing fixture classes.
    pytest_fixtures_submodule = get_pytest_fixtures_submodule()

    # Return this exception class defined by this submodule.
    return pytest_fixtures_submodule.FixtureLookupError
