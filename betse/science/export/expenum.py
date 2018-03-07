#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Export-specific enumeration types** (i.e., :class:`enum.Enum` subclasses
expressing mutually exclusive choices specific to post-simulation exports).
'''

# ....................{ IMPORTS                            }....................
from betse.util.type.enums import make_enum

# ....................{ ENUMS                              }....................
SimExportType = make_enum(
    class_name='SimExportType', member_names=('ANIM', 'CSV', 'PLOT',))
'''
Enumeration of all supported types of **post-simulation exports** (e.g., plots,
animations, comma-separated value (CSV) files).

Attributes
----------
ANIM : enum
    Post-simulation animations, visualizing the cell cluster across one or more
    sampled time steps of an initialization or simulation phase as the same
    number of frames of an animated video.
CSV : enum
    Post-simulation comma-separated value (CSV) files, aggregating raw data for
    the cell cluster across one or more sampled time steps of an initialization
    or simulation phase.
PLOT : enum
    Post-simulation plots, visualizing the cell cluster across one or more
    sampled time steps of an initialization or simulation phase as a single
    non-animated image.
'''
