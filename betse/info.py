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
# Version specifier. For PEP 8 compliance, such specifier is stored in the
# canonical variable "__variable__" rather than as a constant "VARIABLE".
__version__ = '0.0.1'

# Program description.
DESCRIPTION = (
    'betse (Bioelectric Tissue Simulation Environment) simulates '
    'propagation of electrical phenomena in biological tissue (e.g., ion '
    'channel-gated current flow).')

# --------------------( WASTELANDS                         )--------------------
