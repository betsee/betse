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
BITS_PER_BYTE = 8
'''
Number of bits per byte.
'''

# Size denominations in base 2 rather than base 10, for mild efficiency.
KB = 1 << 10
MB = 1 << 20
GB = 1 << 30
TB = 1 << 40

# ....................{ CONSTANTS ~ max                    }....................
BYTE_VALUE_MAX = 255
'''
Maximum value of unsigned bytes.
'''

INT_VALUE_MAX_32_BIT = 1 << 32
'''
Maximum value of variables of internal type `Py_ssize_t` on 32-bit systems.

Such value is suitable for comparison against `sys.maxsize`, the maximum value
for such variables on the current system.
'''

# --------------------( WASTELANDS                         )--------------------
