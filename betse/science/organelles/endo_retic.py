#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

"""

Creates an ER (endoplasmic reticulum) class, which includes ER-specific pumps, channels, and specific methods
relating to calcium dynamics including calcium induced calcium release controlled by inositol-triphosphate.
This class also contains the facilities to initialize, define the core computations for a simulation loop,
remove ER during a cutting event, save and report on data, and plot.

"""

import os
import os.path
import numpy as np
from betse.science import toolbox as tb
from betse.science import sim_toolbox as stb
from betse.util.io.log import logs
import matplotlib.pyplot as plt
from betse.exceptions import BetseExceptionParameters
from betse.science.plot import plot as viz
from betse.science.plot.anim.anim import AnimCellsTimeSeries, AnimEnvTimeSeries