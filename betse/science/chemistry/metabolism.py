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

    def __init__(self, p):

        # Make the BETSE-specific cache directory if not found.
        betse_cache_dir = os.path.expanduser(p.init_path)
        os.makedirs(betse_cache_dir, exist_ok=True)

        # Define data paths for saving an initialization and simulation run:
        self.savedMoM = os.path.join(betse_cache_dir, 'MetabolicNetwork.betse')

    def read_metabo_config(self, sim, cells, p):  # can this be implemented in MasterofMolecules instead? To cut down on classes?

        # create the path to read the metabolism config file:

        self.configPath = os.path.join(p.config_dirname, p.metabo_config_filename)

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

    def run_core_sim(self, sim, cells, p):

        # specify a time vector
        loop_time_step_max = p.init_tsteps
        # Maximum number of seconds simulated by the current run.
        loop_seconds_max = loop_time_step_max * p.dt
        # Time-steps vector appropriate for the current run.
        tt = np.linspace(0, loop_seconds_max, loop_time_step_max)

        # create a time-samples vector
        tsamples = set()
        i = 0
        while i < len(tt) - p.t_resample:
            i += p.t_resample
            tsamples.add(tt[i])

        self.core.clear_cache()
        self.time = []

        for t in tt:

            self.core.run_loop_reactions(t, sim, cells, p)

            if t in tsamples:

                logs.log_info('------------------' + str(np.round(t,3)) +' s --------------------')
                self.time.append(t)
                self.core.write_data(sim, p)
                self.core.report(sim, p)

        logs.log_info('Saving simulation...')
        datadump = [self, cells, p]
        fh.saveSim(self.savedMoM, datadump)
        message = 'Metabolic network simulation saved to' + ' ' + self.savedMoM
        logs.log_info(message)

        logs.log_info('-------------------Simulation Complete!-----------------------')


















