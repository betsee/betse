#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Top-level application namespace.

For PEP 8 compliance, this namespace exposes a subset of the metadata constants
provided by the :mod:`betse.metadata` module commonly inspected by external
automation.
'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from modules guaranteed to exist at the start of
# installation. This includes all standard Python and application modules but
# *NOT* third-party dependencies, which if currently uninstalled will only be
# installed at some later time in the installation.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ IMPORTS                            }....................
# For PEP 8 compliance, versions constants expected by external automation are
# imported under their PEP 8-mandated names.
from betse.metadata import VERSION as __version__
from betse.metadata import VERSION_PARTS as __version_info__

# ....................{ GLOBALS                            }....................
# Document all global variables imported into this namespace above.

__version__
'''
Human-readable application version as a ``.``-delimited string.

For PEP 8 compliance, this specifier has the canonical name ``__version__``
rather than that of a typical global (e.g., ``VERSION_STR``).
'''


__version_info__
'''
Machine-readable application version as a tuple of integers.

For PEP 8 compliance, this specifier has the canonical name ``__version_info__``
rather than that of a typical global (e.g., ``VERSION_PARTS``).
'''
