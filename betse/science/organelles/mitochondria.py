#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

"""

Creates a mitochondria class, which includes a suite of mitochondria-specific molecules, along with
mitochondria-specific pumps (i.e. electron transport chain), channels, and specific methods.
This class also contains the facilities to initialize, define the core computations for a simulation loop,
remove mitochondria during a cutting event, save and report on data, and plot.

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