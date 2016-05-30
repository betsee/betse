#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Decorators skipping tests.

These decorators conditionally mark their decorated tests as skipped depending
on whether the conditions signified by the passed parameters are satisfied
(e.g., the importability of the passed module name).
'''

# ....................{ IMPORTS                            }....................
import pytest
from betse.util.type import types

# ....................{ IMPORTS ~ private                  }....................
# Sadly, the following imports require private modules and packages.
from _pytest.runner import Skipped

# ....................{ SKIP ~ alias                       }....................
skip_if = pytest.mark.skipif
'''
Conditionally skip the decorated test with the passed human-readable
justification if the passed boolean is `False`.

Parameters
----------
boolean : bool
    Boolean to be tested.
reason : str
    Human-readable message justifying the skipping of this test.
'''

# ....................{ SKIP ~ module                      }....................
def skip(reason: str) -> callable:
    '''
    Unconditionally skip the decorated test with the passed human-readable
    justification.

    This decorator is intended to be called both directly as a function _and_
    indirectly as a decorator, which differs from both:

    * `pytest.skip()`, intended to be called only directly as a function.
      Attempting to call this function indirectly as a decorator produces
      extraneous ignorable messages on standard output resembling
      `SKIP [1] betse_test/unit/test_import.py:66: could not import 'xdist'`,
      for unknown (and probably uninteresting) reasons.
    * `pytest.mark.skip()`, intended to be called only indirectly as a
      decorator. Attempting to call this decorator directly as a function
      reduces to a noop, for unknown (and probably uninteresting) reasons.

    Parameters
    ----------
    reason : str
        Human-readable message justifying the skipping of this test.
    '''
    assert types.is_str_nonempty(reason), (
        types.assert_not_str_nonempty(reason, 'Reason'))

    return skip_if(True, reason=reason)

# ....................{ SKIP ~ module                      }....................
def skip_unless_module(module_name: str, minimum_version: str = None):
    '''
    Skip the decorated test if the module with the passed name is unimportable
    _or_ importable but of a version less than the passed minimum version if
    non-`None`.

    This decorator is typically applied to tests requiring optional third-party
    Python dependencies.

    Parameters
    ----------
    module_name : str
        Fully-qualified name of the module to be tested for.
    minversion : optional[str]
        Optional minimum version of this module as a dot-delimited string (e.g.,
        `0.4.0`) to be tested for if any _or_ `None` otherwise, in which case
        any version is acceptable. Defaults to `None`.

    Returns
    ----------
    pytest.skipif
        Decorator describing these requirements if unmet _or_ the identity
        decorator reducing to a noop otherwise.
    '''
    assert types.is_str_nonempty(module_name), (
        types.assert_not_str_nonempty(module_name, 'Module name'))
    assert types.is_str_nonempty_or_none(minimum_version), (
        types.assert_not_str_nonempty_or_none(
            minimum_version, 'Minimum version'))

    # Attempt to import this module and module version.
    try:
        pytest.importorskip(module_name, minimum_version)
    # If unimportable, skip this test.
    except Skipped as exc:
        return skip(str(exc))
    # If an unexpected exception was raised, skip this test with a stack trace.
    except Exception as exc:
        # Defer heavyweight imports to their point of use.
        from betse.util.io import stderrs

        # Print this exception's stack trace to stderr.
        stderrs.output_exception(heading=(
            'skip_unless_module({}, {}) '
            'raised unexpected exception:\n'.format(
                module_name, minimum_version)))

        # Skip this test with this exception's message.
        return skip(str(exc))
    # Else, this module is importable. Reduce this decoration to a noop.
    else:
        return identity


def identity(obj: object) -> object:
    '''
    Return the passed object unmodified.

    This function is typically used as the **identity decorator** (i.e.,
    decorator returning the decorated object unmodified).
    '''

    return obj

# ....................{ SKIP ~ plugin                      }....................
skip_unless_plugin_xdist = skip_unless_module('xdist')
'''
Skip the decorated test if the `pytest-xdist` plugin is _not_ installed.

This decorator is typically applied to tests requiring **process isolation**
(i.e., isolating tests to dedicated subprocesses of the current test session).
While this plugin provides such isolation out-of-the-box, vanilla `py.test` does
_not_. Hence, these tests _must_ be skipped in the absence of this plugin.

Examples of tests requiring process isolation include:

* Unit tests testing importability. Since Python caches imports performed by the
  active Python interpreter (e.g., via `sys.modules`) _and_ since the order in
  which `py.test` runs tests should be assumed to be non- deterministic,
  importability _cannot_ be reliably tested within a single Python process.
'''
