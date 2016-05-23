#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS                            }....................
from betse import ignition

# ....................{ CONVENIENCES                       }....................
# As a convenience for interactive use (e.g., from web-based IPython notebooks
# or CLI-based Python consoles), implicitly initialize BETSE on the first import
# of this subpackage. If BETSE has already been initialized (as typically occurs
# if this subpackage has been imported by a non-interactive CLI or GUI process),
# this call silently reduces to a noop -- which is a good thing.
ignition.init()
