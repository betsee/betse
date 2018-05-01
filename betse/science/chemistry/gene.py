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
from betse.science.config import confio
from betse.science.chemistry.netplot import set_net_opts
# from betse.science import sim_toolbox as stb
from betse.science.visual.plot import plotutil as viz
import matplotlib.pyplot as plt
from betse.science.organelles.microtubules import Mtubes


class MasterOfGenes(object):

    def __init__(self, p):

        pass


    def read_gene_config(self, sim, cells, p):

        # Previously loaded GRN-specific configuration file as a dictionary.
        config_dic = p.grn.conf

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

        # reset microtubules?
        self.reset_MT = config_dic.get('reset microtubules', False)

        # recalculate fluid flow?
        self.recalc_fluid = config_dic.get('recalculate fluid', False)

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

        # after primary initialization, check and see if optimization required:
            # after primary initialization, check and see if optimization required:
        opti = config_dic['optimization']['optimize network']
        self.core.opti_N = config_dic['optimization']['optimization steps']
        self.core.opti_method = config_dic['optimization']['optimization method']
        self.core.target_vmem = float(config_dic['optimization']['target Vmem'])
        self.core.opti_T = float(config_dic['optimization']['optimization T'])
        self.core.opti_step = float(config_dic['optimization']['optimization step'])
        # self.core.opti_run = config_dic['optimization']['run from optimization']

        if opti is True:
            logs.log_info("The Gene Network is being analyzed for optimal rates...")
            self.core.optimizer(sim, cells, p)
            self.reinitialize(sim, cells, p)

    def reinitialize(self, sim, cells, p):

        # Previously loaded GRN-specific configuration file as a dictionary.
        config_dic = p.grn.conf

        # Time dilation:
        self.core.time_dila = float(config_dic.get('time dilation factor', 1.0))

        # reset microtubules?
        self.reset_MT = config_dic.get('reset microtubules', False)

        # recalculate fluid flow?
        self.recalc_fluid = config_dic.get('recalculate fluid flow', False)

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

        # set molecules to not affect charge for sim-grn test-drives:
        p.substances_affect_charge = False

        # Set sim-grn specific time step and sampling:
        p.dt = p.grn_dt
        p.init_time_total = p.grn_total_time
        p.init_time_sampling = p.grn_tsample

        p.init_tsteps = int(p.init_time_total / p.dt)
        p.resample = p.init_time_sampling
        p.total_time = p.init_time_total

        # Number of time steps (including sampled and unsampled) between each
        # unsampled time step, including that unsampled time step itself.
        p.t_resample = p.resample / p.dt

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
            # logs.log_debug('Time sample i: {!r}'.format(i))
            tsamples.add(tt[i])

        # if p.grn_runmodesim:
        self.reinitialize(sim, cells, p)

        self.core.clear_cache()
        self.time = []

        self.mod_after_cut = False # set this to false

        if self.recalc_fluid:  # If user requests the GRN recalculate/calculate fluid:

            logs.log_info("Calculating fluid in terms of endogenous currents...")

            p.fluid_flow = True  # turn fluid flow on (in case it was off)

            # calculate or re-calculate fluid flow in terms of curl-component of current:
            cc = sim.cc_env.mean(axis=0).reshape(cells.X.shape)
            zz = sim.zs.mean()

            sim.u_env_x = -sim.J_env_x / (p.F * cc * zz)
            sim.u_env_y = -sim.J_env_y / (p.F * cc * zz)

        if p.use_microtubules:

            sim.mtubes.reinit(cells, p)

            self.mtubes_x_time = []
            self.mtubes_y_time = []

            if self.reset_MT:
                logs.log_info("Resetting microtubules for sim-grn simulation...")
                sim.mtubes = Mtubes(sim, cells, p)

        for t in tt:

            if self.transporters:
                self.core.run_loop_transporters(t, sim, cells, p)

            self.core.run_loop(t, sim, cells, p)


            if p.use_microtubules: # update the microtubules:
                sim.mtubes.update_mtubes(cells, sim, p)

            if p.grn_runmodesim is True and t > p.cut_time and self.mod_after_cut is False:

                sim.dyna.runAllDynamics(sim, cells, p, t)

                if sim.dyna.event_cut.is_fired and self.mod_after_cut is False: # if a cutting event has just been run:

                    self.core.mod_after_cut_event(sim.target_inds_cell_o, sim.target_inds_mem_o, sim, cells, p)
                    logs.log_info(
                        "Redefining dynamic dictionaries to point to the new sim...")
                    self.core.redefine_dynamic_dics(sim, cells, p)

                    logs.log_info(
                        "Reinitializing the gene regulatory network for simulation...")
                    self.reinitialize(sim, cells, p)
                    # self.core.clear_cache()
                    sim.uxmt, self.uymt = sim.mtubes.mtubes_to_cell(cells, p)

                    self.mod_after_cut = True  # set the boolean to avoid repeat action


            if t in tsamples:
                sim.time.append(t)

                logs.log_info('------------------' + str(np.round(t,3)) +' s --------------------')
                self.time.append(t)
                self.core.write_data(sim, cells, p)
                self.core.report(sim, p)

                if p.use_microtubules:
                    # microtubules:
                    self.mtubes_x_time.append(sim.mtubes.mtubes_x * 1)
                    self.mtubes_y_time.append(sim.mtubes.mtubes_y * 1)



        logs.log_info('Saving simulation...')
        datadump = [self, cells, p]
        fh.saveSim(p.grn_pickle_filename, datadump)
        self.core.init_saving(cells, p, plot_type='grn', nested_folder_name='RESULTS')

        # microtubules plot------------------------------------------------------------------------
        if p.use_microtubules:

            logs.log_info("Plotting microtubules used in GRN simulation...")

            plt.figure()
            ax = plt.subplot(111)

            umtx, umty = sim.mtubes.mtubes_to_cell(cells, p)

            plt.figure()
            ax = plt.subplot(111)

            viz.plotVectField(
                umtx,
                umty,
                cells,
                p,
                plot_ecm=False,
                title='Final Microtubule Alignment Field',
                cb_title='Aligned MT Fraction',
                colorAutoscale=False,
                minColor=0.0,
                maxColor=1.0,
            )

            # viz.mem_quiver(
            #     sim.mtubes.mtubes_x,
            #     sim.mtubes.mtubes_y,
            #     ax,
            #     cells,
            #     p,
            # )

            ax.set_xlabel('X-Distance [um]')
            ax.set_ylabel('Y-Distance [um]')
            ax.set_title('Microtubule arrangement in cells')

            if p.autosave is True:
                savename = self.core.imagePath + 'Microtubules' + '.png'
                plt.savefig(savename, format='png', transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

        if self.recalc_fluid:

            logs.log_info("Plotting fluid flow used in GRN simulation")

            plt.figure()
            ax = plt.subplot(111)

            viz.plotStreamField(
                1e9 * sim.u_env_x,
                1e9 * sim.u_env_y,
                cells, p,
                plot_ecm=True,
                title='Final Fluid Velocity in Environment',
                cb_title='Velocity [nm/s]'
            )

            if p.autosave is True:
                savename = self.core.imagePath + 'Fluid_ECM' + '.png'
                plt.savefig(savename, format='png', transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)


        self.core.export_eval_strings(p)
        self.core.export_equations(p)
        message = 'Gene regulatory network simulation saved to' + ' ' + p.grn_pickle_filename
        logs.log_info(message)

        logs.log_info('-------------------Simulation Complete!-----------------------')
