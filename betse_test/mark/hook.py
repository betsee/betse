#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Decorators marking tests with custom keywords inspected by custom py.test hooks.

These decorators unconditionally mark their decorated tests with custom (and
hence unofficial) keywords inspected _only_ by custom (and hence unofficial)
py.test hooks defined by `conftest` plugins.
'''

# ....................{ IMPORTS                            }....................
import pytest

# ....................{ FAIL ~ alias                       }....................
serial_parametrized = pytest.mark.serial_parametrized
'''
Unconditionally mark the decorated parametrized test as **serial**.

If the decorated test is _not_ parametrized, an exception is raised.  Else, each
set of parameter values passed to this test is tested serially. The success of
each subsequently passed set of parameter values depends on the success of all
previously passed sets of parameter values. On the first failure of this test
induced by a passed set of parameter values, each subsequent call to this test
passed a subsequent set of parameter values is automatically marked as an
`XFAIL` and hence fails _without_ being run.

Implementation
----------
The majority of the black magic required by this decoration is implemented as
low-level py.test hooks in the top-level `betse_func.conftest` plugin. To
preserve state between parametrized calls to the same test, these hooks
dynamically add the following attributes to this test's underlying function or
method object:

* `_betse_first_failing_param_id`, the unique identifier of the first set of
  parameter values passed to this test raising an exception for the current test
  session if any _or_ `None` otherwise (i.e., if this test has yet to raise an
  exception for any parameters).
'''
