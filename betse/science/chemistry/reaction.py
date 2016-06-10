#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import copy
import os
import os.path
import time
import numpy as np
from random import shuffle
from scipy import interpolate as interp
from scipy.ndimage.filters import gaussian_filter
from betse.util.io.log import logs
from betse.science import filehandling as fh
from betse.science import finitediff as fd
from betse.science import toolbox as tb

def enzyme(reactants, products, Km_reactants, Km_products, rate):

    pass