#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Decorators failing tests.

These decorators conditionally mark their decorated tests as failing depending
on whether the conditions signified by the passed parameters are satisfied
(e.g., the importability of the passed module name).
'''

# ....................{ IMPORTS                            }....................
import pytest
from betse.util.type.types import type_check

# ....................{ FAIL ~ alias                       }....................
xfail_if = pytest.mark.xfail
'''
Conditionally mark the decorated test as ignorably known to fail with the
passed human-readable justification if the passed boolean is ``False``.

Parameters
----------
boolean : bool
    Boolean to be tested.
reason : str
    Human-readable message justifying the failure of this test.

See Also
----------
:func:`pytest.mark.xfail`
    Further details on the ``XFAIL`` test state.
'''


@type_check
def xfail(reason: str):
    '''
    Unconditionally mark the decorated test as ignorably known to fail with the
    passed human-readable justification.

    py.test will neither run this test nor accumulate this failure. While
    superficially similar to tests unconditionally skipped via the `@skip()`
    decorator, this failure will be collected as an `XFAIL` by py.test reporting.

    Parameters
    ----------
    reason : str
        Human-readable message justifying the failure of this test.
    '''

    return xfail_if(True, reason=reason)
