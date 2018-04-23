#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Create and electrodiffuses a suite of customizable general molecule in the BETSE
ecosystem, including functionality to pump the molecule, use it as a gating
ligand, produce and consume it, and use it an enzyme to facilitate another
reaction.

The molecule is assumed to be at low concentrations and to not have a
significant effect on system voltages or currents. This module creates a
structure containing all user-defined molecules, along with the facilities to
initialize, define the core computations for a simulation loop, save and report
on data, and plot.
'''

#FIXME: This submodule no longer appears to be used anywhere in the codebase.
#Can we safely remove this now? If not, no worries. Lumbering lumberjacks!

import numpy as np
from betse.science import filehandling as fh
from betse.util.io.log import logs
from betse.util.path import pathnames
from betse.science.chemistry.networks import MasterOfNetworks
from betse.science.config import confio
from betse.exceptions import BetseSimConfException
from betse.science.chemistry.netplot import set_net_opts


class MasterOfMetabolism(object):

    def __init__(self, p):

        #FIXME: Extract the "MetabolicNetwork.betse" basename into a new
        #configuration option, perhaps in the "init file saving" section.
        #FIXME: Replace "MetabolicNetwork.betse" with
        #"MetabolicNetwork.betse.gz" to compress this pickled file.

        # Define data paths for saving an initialization and simulation run:
        self.savedMoM = pathnames.join(p.init_pickle_dirname, 'MetabolicNetwork.betse')

    def read_metabo_config(self, sim, cells, p):

        # create the path to read the metabolism config file:
        self.configPath = pathnames.join(
            p.conf_dirname, p.metabo_config_filename)

        # read the config file into a dictionary:
        self.config_dic = confio.read_metabo(self.configPath)

        # determine if mitochondria are enabled:
        self.mit_enabled = self.config_dic['enable mitochondria']

        # obtain specific sub-dictionaries from the config file:
        substances_config = self.config_dic['biomolecules']
        reactions_config = self.config_dic.get('reactions', None)
        transporters_config = self.config_dic.get('transporters', None)
        channels_config = self.config_dic.get('channels', None)
        modulators_config = self.config_dic.get('modulators', None)

        # initialize the substances of metabolism in a core field encapsulating
        # Master of Molecules:
        self.core = MasterOfNetworks(sim, cells, substances_config, p, mit_enabled=self.mit_enabled)

        if reactions_config is not None:

            # initialize the reactions of metabolism:
            self.core.read_reactions(reactions_config, sim, cells, p)
            self.core.write_reactions()
            self.core.create_reaction_matrix()
            self.core.write_reactions_env()
            self.core.create_reaction_matrix_env()

            if self.mit_enabled is True:
                self.core.write_reactions_mit()
                self.core.create_reaction_matrix_mit()

            self.reactions = True

        else:
            self.core.create_reaction_matrix()
            self.core.create_reaction_matrix_env()
            self.reactions = False

        # initialize transporters, if defined:
        if transporters_config is not None and len(transporters_config) >0:
            self.core.read_transporters(transporters_config, sim, cells, p)
            self.core.write_transporters(sim, cells, p)
            self.transporters = True

        else:
            self.transporters = False

        # initialize any custom channels:-------------

        if channels_config is not None and len(channels_config) > 0:
            self.core.read_channels(channels_config, sim, cells, p)
            self.channels = True

        else:
            self.channels = False

        # initialize any modulators------------------

        if modulators_config is not None and len(modulators_config) > 0:
            self.core.read_modulators(modulators_config, sim, cells, p)
            self.modulators = True

        else:
            self.modulators = False

        # test to make sure the metabolic simulation includes core components:
        if 'ATP' not in self.core.molecules or 'ADP' not in self.core.molecules or 'Pi' not in self.core.molecules:

            raise BetseSimConfException("This metabolic simulation does not contain key substances."
                                           "Please define 'ATP', 'ADP' and 'Pi' biomolecules in your "
                                           "metabolism configuration file and try again.")

        # read in network plotting options:
        self.core.net_plot_opts = self.config_dic.get('network plotting', None)

        # set plotting options for the network:
        set_net_opts(self.core, self.core.net_plot_opts, p)

            # after primary initialization, check and see if optimization required:
        opti = self.config_dic['optimization']['optimize network']
        self.core.opti_N = self.config_dic['optimization']['optimization steps']
        self.core.opti_method = self.config_dic['optimization']['optimization method']
        self.core.target_vmem = float(self.config_dic['optimization']['target Vmem'])
        self.core.opti_T = float(self.config_dic['optimization']['optimization T'])
        self.core.opti_step = float(self.config_dic['optimization']['optimization step'])
    #     self.core.opti_run = self.config_dic['optimization']['run from optimization']
    #
        if opti is True:
            logs.log_info("The Metabolic Network is being analyzed for optimal rates...")
            self.core.optimizer(sim, cells, p)
            self.reinitialize(sim, cells, p)
    #
    #     if self.core.opti_run is True:
    #         self.run_from_init(self, sim, cells, p)
    #
    # def run_from_init(self, sim, cells, p):
    #
    #     pass

    #FIXME: Oh, boy. Most of this method appears to have been copy-and-pasted
    #from the read_metabo_config() method above. That's... not the best. Let's
    #extract the code shared in common between these two methods into a new
    #_init_metabo() method internally called by these two methods.
    def reinitialize(self, sim, cells, p):

        # create the path to read the metabolism config file:
        self.configPath = pathnames.join(
            p.conf_dirname, p.metabo_config_filename)

        # read the config file into a dictionary:
        self.config_dic = confio.read_metabo(self.configPath)

        # determine if mitochondria are enabled:
        self.mit_enabled = self.config_dic['enable mitochondria']

        # obtain specific sub-dictionaries from the config file:
        substances_config = self.config_dic['biomolecules']
        reactions_config = self.config_dic.get('reactions', None)
        transporters_config = self.config_dic.get('transporters', None)
        channels_config = self.config_dic.get('channels', None)
        modulators_config = self.config_dic.get('modulators', None)

        # initialize the substances of metabolism in a core field encapsulating
        # Master of Molecules:
        self.core.tissue_init(sim, cells, substances_config, p)

        if reactions_config is not None and len(reactions_config):

            # initialize the reactions of metabolism:
            self.core.read_reactions(reactions_config, sim, cells, p)
            self.core.write_reactions()
            self.core.create_reaction_matrix()
            self.core.write_reactions_env()
            self.core.create_reaction_matrix_env()

            if self.mit_enabled is True:
                self.core.write_reactions_mit()
                self.core.create_reaction_matrix_mit()

            self.reactions = True

        else:
            self.core.create_reaction_matrix()
            self.core.create_reaction_matrix_env()
            self.reactions = False

        # initialize transporters, if defined:
        if transporters_config is not None and len(transporters_config) > 0:
            self.core.read_transporters(transporters_config, sim, cells, p)
            self.core.write_transporters(sim, cells, p)
            self.transporters = True

        else:
            self.transporters = False

        # initialize any custom channels:-------------

        if channels_config is not None and len(channels_config) > 0:
            self.core.read_channels(channels_config, sim, cells, p)
            self.channels = True

        else:
            self.channels = False

        # initialize any modulators------------------

        if modulators_config is not None and len(modulators_config) > 0:
            self.core.read_modulators(modulators_config, sim, cells, p)
            self.modulators = True

        else:
            self.modulators = False

        # test to make sure the metabolic simulation includes core components:
        if 'ATP' not in self.core.molecules or 'ADP' not in self.core.molecules or 'Pi' not in self.core.molecules:
            raise BetseSimConfException("This metabolic simulation does not contain key substances."
                                           "Please define 'ATP', 'ADP' and 'Pi' biomolecules in your "
                                           "metabolism configuration file and try again.")

    def run_core_sim(self, sim, cells, p):
        """
        Runs a simulation of the biochemical reaction network only, with a dummy sim and dyna module.
        This allows the user to test the reaction network without the influence of bioelectrical dynamics.

        This method is called in the BETSE CLI command betse sim-brn my_yaml.yaml and data plotted via
        betse plot sim-brn my_yaml.yaml

        """

        # # point sim.metabo to this object
        # sim.metabo = self

        # create a dictionary pointing to key metabolic molecules used in sim: ATP, ADP and Pi:
        sim.met_concs = {'cATP': self.core.cell_concs['ATP'][cells.mem_to_cells],
            'cADP': self.core.cell_concs['ADP'][cells.mem_to_cells],
            'cPi': self.core.cell_concs['Pi'][cells.mem_to_cells]}

        # initialize Vmem to an initial value common to many cell types:
        sim.vm = -50e-3*np.ones(sim.mdl)

        sim.conc_J_x = np.zeros(sim.edl)
        sim.conc_J_y = np.zeros(sim.edl)

        # initialize key fields of simulator required to interface (dummy init)
        sim.rho_pump = 1.0
        sim.rho_channel = 1.0

        sim.J_env_x = np.zeros(sim.edl)
        sim.J_env_y = np.zeros(sim.edl)

        # p.substances_affect_charge = False

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
        sim.time = []

        for t in tt:

            if self.transporters:
                self.core.run_loop_transporters(t, sim, cells, p)

            # if self.channels:
            #     self.core.run_loop_channels(sim, self.core, cells, p)
            #
            # if self.modulators:
            #     self.core.run_loop_modulators(sim, self.core, cells, p)

            self.core.run_loop(t, sim, cells, p)

            # # update core ions in sim:
            # for i in sim.movingIons:
            #
            #     # update the ion concentration intra-cellularly:
            #     sim.cc_mems[i][:], sim.cc_cells[i][:], _ = \
            #         stb.update_intra(sim, cells, sim.cc_mems[i][:],
            #             sim.cc_cells[i][:],
            #             sim.D_free[i],
            #             sim.zs[i], p)


            if t in tsamples:

                sim.time.append(t)

                logs.log_info('------------------' + str(np.round(t,3)) +' s --------------------')
                self.time.append(t)
                self.core.write_data(sim, cells, p)
                self.core.report(sim, p)

        logs.log_info('Saving simulation...')
        datadump = [self, cells, p]
        fh.saveSim(self.savedMoM, datadump)
        self.core.init_saving(cells, p, plot_type='init', nested_folder_name='Metabolism')
        self.core.export_eval_strings(p)
        self.core.export_equations(p)
        message = 'Metabolic network simulation saved to' + ' ' + self.savedMoM
        logs.log_info(message)

        logs.log_info('-------------------Simulation Complete!-----------------------')

    def update_ATP(self, flux, sim, cells, p):

        """
        Update ATP, ADP and Pi concentrations using a
        concentration change defined on membranes or cell
        centres.

        This method is typically called in sim after ATPase pumps are run.

        flux        concentration flux [mol/m2 s], negative consumes ATP; array must be of length sim.mdl

        """

        cATP = self.core.cell_concs['ATP']
        cADP = self.core.cell_concs['ADP']
        cPi = self.core.cell_concs['Pi']

        deltac = (np.dot(cells.M_sum_mems, flux*cells.mem_sa)/cells.cell_vol)*p.dt

        self.core.cell_concs['ATP'] = cATP + deltac
        self.core.cell_concs['ADP'] = cADP - deltac
        self.core.cell_concs['Pi'] = cPi - deltac

        sim.met_concs = {'cATP': self.core.cell_concs['ATP'][cells.mem_to_cells],
            'cADP': self.core.cell_concs['ADP'][cells.mem_to_cells],
            'cPi': self.core.cell_concs['Pi'][cells.mem_to_cells]}
