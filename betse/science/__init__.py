#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ CONVENIENCES                       }....................
#FIXME: Revise docstring. This function no longer contains test-specific logic.
def _ignite() -> None:
    '''
    Initialize both the current application _and_ all mandatory third-party
    dependencies of this application with sane defaults if this application
    both has not already been initialized _and_ is not currently being tested.

    This function is a convenience intended to streamline interactive use (e.g.,
    from web-based IPython notebooks or CLI-based Python consoles) by implicitly
    initializing this application on the first import of this subpackage.

    If this application either has already been manually initialized _or_ is
    currently being tested by a test suite that that may subsequently do so,this
    call silently reduces to a noop -- which is a good thing.
    '''

    # Localize *ALL* imports to this function, preventing the package namespace
    # from being polluted with imports requiring deletion via the del() builtin.
    from betse import ignition
    # from betse.util.py import pys

    # If tests are being run, avoid pre-initializing this application.
    # if pys.is_testing():
    #     return

    # Initialize both this application and all dependencies thereof.
    ignition.ignite()

#FIXME: Overly cumbersome. Just inline the two statements above here; then,
#excise this function.
# Initialize both this application and all dependencies thereof.
_ignite()

# ....................{ CLEANUP                            }....................
# Delete *ALL* attributes (including callables) defined above, preventing the
# package namespace from being polluted with these attributes.
del _ignite
