#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

"""

Create and electrodiffuses a suite of customizable general molecule in the BETSE ecosystem,
including functionality to pump the molecule, use it as a gating ligand, produce and consume it,
and use it an enzyme to facilitate another reaction. The molecule is assumed to be at low
concentrations and to not have a significant effect on system voltages or currents. This
module creates a structure containing all user-defined molecules, along with the facilities
to initialize, define the core computations for a simulation loop, save and report on data,
and plot.

"""

import copy
import os
import os.path
import time
import numpy as np
from betse.science import filehandling as fh
from betse.science import toolbox as tb
from betse.science import sim_toolbox as stb
from betse.util.io.log import logs
import matplotlib.pyplot as plt
from betse.exceptions import BetseExceptionParameters, BetseExceptionSimulation
from betse.science.plot import plot as viz
from betse.science.plot.anim.anim import AnimCellsTimeSeries, AnimEnvTimeSeries
from betse.science.chemistry.molecule import MasterOfMolecules, Molecule
from betse.science.config import sim_config

class MasterOfMetabolism(object):

    def __init__(self, sim, cells, p):

        pass

    def read_config(self, sim, cells, p):

        if p.metabolism_enabled:

            # create the path to read the metabolism config file:
            self.configPath = os.path.expanduser(p.metabo_config_filename)

            # read the config file into a dictionary:
            self.config_dic = sim_config.read_metabo(self.configPath)

            # obtain specific sub-dictionaries from the config file:
            substances_config = self.config_dic['biomolecules']
            reactions_config = self.config_dic['reactions']

            # initialize the substances of metabolism in a core field encapsulating
            # Master of Molecules:
            self.core = MasterOfMolecules(sim, substances_config, p)

            # initialize the reactions of metabolism:
            self.core.read_reactions(reactions_config, sim, cells, p)









