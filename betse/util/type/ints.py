#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Low-level integer facilities.
'''

# ....................{ IMPORTS                            }....................
# import sys

# ....................{ CONSTANTS                          }....................
# Number of bits per byte.
BITS_PER_BYTE = 8

# Largest unsigned integer value representable with one byte.
BYTE = 255

# Size denominations in base 2 rather than base 10, for mild efficiency.
KB = 1 << 10
MB = 1 << 20
GB = 1 << 30
TB = 1 << 40

# --------------------( WASTELANDS                         )--------------------
