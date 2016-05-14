#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import copy
import os
import os.path
import time
from random import shuffle

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate as interp
from scipy.ndimage.filters import gaussian_filter

from betse.exceptions import BetseExceptionSimulation
from betse.util.io.log import logs

from betse.science import filehandling as fh
from betse.science import finitediff as fd
from betse.science import toolbox as tb
from betse.science.plot.anim.anim import AnimCellsWhileSolving
from betse.science import sim_toolbox as stb
from betse.science.tissue.channels_o import Gap_Junction
from betse.science.tissue.handler import TissueHandler