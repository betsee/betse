#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

"""
Controls a gene regulatory network.

Creates and electrodiffuses a suite of customizable general gene products in the
BETSE ecosystem, where the gene products are assumed to activate and/or inhibit
the expression of other genes (and therefore the production of other gene
products) in the gene regulatory network (GRN).
"""

import numpy as np
from betse.science import filehandling as fh
from betse.util.io.log import logs
from betse.util.path import pathnames
from betse.science.chemistry.networks import MasterOfNetworks
from betse.science.chemistry.netplot import set_net_opts
# from betse.science.math.systems import SimMaster


class MasterOfMolecules(object):

    def __init__(self, p):

        #FIXME: Extract the "GeneralNetwork.betse" basename into a new
        #configuration option, perhaps in the "init file saving" section.
        #FIXME: Replace "GeneralNetwork.betse" with "GeneralNetwork.betse.gz" to
        #compress this pickled file.

        # Define data paths for saving an initialization and simulation run:
        self.savedMoM = pathnames.join(p.init_pickle_dirname, 'GeneralNetwork.betse')

    def read_mol_config(self, sim, cells, p):

        # read the config file into a dictionary:
        config_dic = p.network_config

        # obtain specific sub-dictionaries from the config file:
        substances_config = config_dic['biomolecules']
        reactions_config = config_dic.get('reactions', None)
        transporters_config = config_dic.get('transporters', None)
        channels_config = config_dic.get('channels', None)
        modulators_config = config_dic.get('modulators', None)

        # initialize the substances of metabolism in a core field encapsulating
        # Master of Molecules:
        self.core = MasterOfNetworks(sim, cells, substances_config, p)

        # Time dilation:
        self.core.time_dila = float(config_dic.get('time dilation factor', 1.0))

        # read in substance properties from the config file, and initialize basic properties:
        self.core.read_substances(sim, cells, substances_config, p)
        self.core.tissue_init(sim, cells, substances_config, p)


        if reactions_config is not None:
            # initialize the reactions of metabolism:
            self.core.read_reactions(reactions_config, sim, cells, p)
            self.core.write_reactions()
            self.core.create_reaction_matrix()
            self.core.write_reactions_env()
            self.core.create_reaction_matrix_env()

            self.reactions = True

        else:
            self.core.create_reaction_matrix()
            self.core.create_reaction_matrix_env()
            self.reactions = False

        # initialize transporters, if defined:
        if transporters_config is not None:
            self.core.read_transporters(transporters_config, sim, cells, p)
            self.core.write_transporters(sim, cells, p)

            self.transporters = True

        else:
            self.transporters = False

        # initialize channels, if desired:
        if channels_config is not None:
            self.core.read_channels(channels_config, sim, cells, p)
            self.channels = True

        else:
            self.channels = False

        # initialize modulators, if desired:

        if modulators_config is not None:
            self.core.read_modulators(modulators_config, sim, cells, p)
            self.modulators = True

        else:
            self.modulators = False

        # read in network plotting options:
        self.core.net_plot_opts = config_dic.get('network plotting', None)

        # set plotting options for the network:
        set_net_opts(self.core, self.core.net_plot_opts, p)

        # plot up direction surfaces:  # FIXME move this to its own command in Simrunner
        # dyna_surf = SimMaster(config_dic, p)

        # after primary initialization, check and see if optimization required:
        opti = config_dic['optimization']['optimize network']
        self.core.opti_N = config_dic['optimization']['optimization steps']
        self.core.opti_method = config_dic['optimization']['optimization method']
        self.core.target_vmem = float(config_dic['optimization']['target Vmem'])
        self.core.opti_T = float(config_dic['optimization']['optimization T'])
        self.core.opti_step = float(config_dic['optimization']['optimization step'])

        if opti is True:
            logs.log_info("The General Network is being analyzed for optimal rates...")
            self.core.optimizer(sim, cells, p)
            self.reinitialize(sim, cells, p)

    #FIXME: Most of this method appears to have been copy-and-pasted
    #from the read_mol_config() method above. That's... not the best. Let's
    #extract the code shared in common between these two methods into a new
    #_init_mols() method internally called by these two methods.
    def reinitialize(self, sim, cells, p):

        # read the config file into a dictionary:
        config_dic = p.network_config

        # Time dilation:
        self.core.time_dila = float(config_dic.get('time dilation factor', 1.0))

        # obtain specific sub-dictionaries from the config file:
        substances_config = config_dic['biomolecules']
        reactions_config = config_dic.get('reactions', None)
        transporters_config = config_dic.get('transporters', None)
        channels_config = config_dic.get('channels', None)
        modulators_config = config_dic.get('modulators', None)

        self.core.tissue_init(sim, cells, substances_config, p)

        if reactions_config is not None:
            # initialize the reactions of metabolism:
            self.core.read_reactions(reactions_config, sim, cells, p)
            self.core.write_reactions()
            self.core.create_reaction_matrix()
            self.core.write_reactions_env()
            self.core.create_reaction_matrix_env()

            self.reactions = True

        else:
            self.core.create_reaction_matrix()
            self.core.create_reaction_matrix_env()
            self.reactions = False

        # initialize transporters, if defined:
        if transporters_config is not None:
            self.core.read_transporters(transporters_config, sim, cells, p)
            self.core.write_transporters(sim, cells, p)

            self.transporters = True

        else:
            self.transporters = False

        # initialize channels, if desired:
        if channels_config is not None:
            self.core.read_channels(channels_config, sim, cells, p)
            self.channels = True

        else:
            self.channels = False

        # initialize modulators, if desired:
        if modulators_config is not None:
            self.core.read_modulators(modulators_config, sim, cells, p)
            self.modulators = True

        else:
            self.modulators = False

    def run_core_sim(self, sim, cells, p):

        sim.vm = -50e-3*np.ones(sim.mdl)

        sim.time = []


        # initialize key fields of simulator required to interface (dummy init)
        sim.rho_pump = 1.0
        sim.rho_channel = 1.0

        # set molecules to not affect charge for sim-grn test-drives:
        p.substances_affect_charge = False

        # specify a time vector
        loop_time_step_max = p.init_tsteps
        # Maximum number of seconds simulated by the current run.
        loop_seconds_max = loop_time_step_max * p.dt
        # Time-steps vector appropriate for the current run.
        tt = np.linspace(0, loop_seconds_max, loop_time_step_max)

        tsamples = set()
        i = 0
        while i < len(tt) - p.t_resample:
            i = int(i + p.t_resample)
            tsamples.add(tt[i])

        self.core.clear_cache()
        self.time = []

        for t in tt:

            if self.transporters:
                self.core.run_loop_transporters(t, sim, cells, p)

            # if self.modulators:
            #     self.core.run_loop_modulators(sim, self.core, cells, p)

            self.core.run_loop(t, sim, cells, p)




            if t in tsamples:
                sim.time.append(t)

                logs.log_info('------------------' + str(np.round(t,3)) +' s --------------------')
                self.time.append(t)
                self.core.write_data(sim, cells, p)
                self.core.report(sim, p)

        logs.log_info('Saving simulation...')
        datadump = [self, cells, p]
        fh.saveSim(p.grn_pickle_filename, datadump)
        self.core.init_saving(cells, p, plot_type='init', nested_folder_name='Gene')
        self.core.export_eval_strings(p)
        self.core.export_equations(p)
        message = 'Gene regulatory network simulation saved to' + ' ' + p.grn_pickle_filename
        logs.log_info(message)

        logs.log_info('-------------------Simulation Complete!-----------------------')
