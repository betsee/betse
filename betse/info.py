#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''Metadata constants describing high-level `betse` behaviour.'''

# ....................{ IMPORTS                            }....................
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# WARNING: To avoid race conditions during setuptools-based installation, this
# module may import *ONLY* from packages guaranteed to exist at installation
# time (e.g., stock Python packages). In particular, this implies that
# packages possibly installed by such installation (e.g., mandatory
# dependencies) must *NOT* be imported from.
#
# Yes, this is terrible. No, we didn't make the ad-hoc rules.
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ....................{ CONSTANTS                          }....................
__version__ = '0.0.1'
'''
Version specifier.

For PEP 8 compliance, such specifier is exposed as the canonical variable
`__variable__` rather than a constant `VARIABLE`.
'''

NAME = 'BETSE'
'''Program name.'''

DESCRIPTION = ''.join((
    NAME, ' ',
    '(Bioelectric Tissue Simulation Environment) simulates propagation of ',
    'electrical phenomena in biological tissue (e.g., ion channel-gated ',
    'current flow).',
))
'''Program description.'''

AUTHORS = 'Alexis Pietak, Cecil Curry, et al.'
'''Program authors as a comma-delimited list.'''

# ....................{ PYTHON ~ version                   }....................
# Validate the version of the active Python interpreter *BEFORE* subsequent code
# possibly depending on such version. Since such version should be validated
# both at setuptools-based install time and post-install runtime *AND* since
# this module is imported sufficiently early by both, stash such validation here
# to avoid duplication of such logic and hence the hardcoded Python version.
#
# The "sys" module exposes three version-related constants for such purpose:
#
# * "hexversion", an integer intended to be used in hexadecimal format: e.g.,
#    >>> sys.hexversion
#    33883376
#    >>> '%x' % sys.hexversion
#    '20504f0'
# * "version", a human-readable string: e.g.,
#    >>> sys.version
#    2.5.2 (r252:60911, Jul 31 2008, 17:28:52)
#    [GCC 4.2.3 (Ubuntu 4.2.3-2ubuntu7)]
# * "version_info", a tuple of three or more integers or strings: e.g.,
#    >>> sys.version_info
#    (2, 5, 2, 'final', 0)
#
# Only the first lends itself to reliable comparison.
#
# Note that the nearly decade-old and officially accepted PEP 345 proposed a new
# field "requires_python" configured via a key-value pair passed to the call to
# setup() in "setup.py" (e.g., "requires_python = ['>=2.2.1'],"), such field has
# yet to be integrated into either disutils or setuputils. Hence, such field is
# validated manually in the typical way. Behead the infidel setuptools!
import sys
if sys.hexversion < 0x03030000:
    raise RuntimeError(''.join((
        NAME, ' ',
        'requires at least Python 3.3. However, the active Python ',
        'interpreter is of version:\n\n',
        sys.version, '\n\n',
        'We feel sadness for you.',
    )))

# --------------------( WASTELANDS                         )--------------------
# Such type is validated by comparison of
# And yes: this *IS* the typical way
# for validating Python versions. Note that
