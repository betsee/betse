#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Gene regulatory network (GRN).

This submodule creates and electrodiffuses a suite of customizable general gene
products in the BETSE ecosystem, where the gene products are assumed to
activate and/or inhibit the expression of other genes (and therefore the
production of other gene products) in the gene regulatory network (GRN).
'''

#FIXME: Unify the large amount of code shared in common between this and the
#"gene" submodule as follows:
#
#* Define a new "grnabc" submodule defining a new "GrnABC" abstract base class,
#  defining at least:
#  * run_core_sim().
#  * All other methods shared in common between the "MasterOfGenes" and
#    "MasterOfMolecules" classes.
#* Refactor the "MasterOfGenes" and "MasterOfMolecules" classes to subclass the
#  "GrnABC" abstract base class.

# ....................{ IMPORTS                           }....................
import matplotlib.pyplot as plt
import numpy as np
from betse.science import filehandling as fh
from betse.science.chemistry.netplot import set_net_opts
from betse.science.chemistry.networks import MasterOfNetworks
from betse.science.organelles.microtubules import Mtubes
from betse.science.phase.phasecls import SimPhase
from betse.science.visual.plot import plotutil as viz
from betse.util.io.log import logs
from betse.util.type.types import type_check

# ....................{ CLASSES                           }....................
class MasterOfGenes(object):
    '''
    Gene regulatory network (GRN).
    '''

    # ..................{ INITIALIZERS                      }..................
    #FIXME: Let's try to remove this method altogether, if we can. Tree frogs!
    def __init__(self, p):
        pass


    @type_check
    def reinitialize(self, phase: SimPhase) -> None:

        # Localize high-level phase objects for convenience.
        cells = phase.cells
        p     = phase.p
        sim   = phase.sim

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

        self.core.tissue_init(phase, substances_config)

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
            self.core.read_transporters(transporters_config, phase)
            self.core.write_transporters(sim, cells, p)

            self.transporters = True

        else:
            self.transporters = False

        # initialize channels, if desired:
        if channels_config is not None:
            self.core.read_channels(channels_config, phase)
            self.channels = True

        else:
            self.channels = False

        # initialize modulators, if desired:
        if modulators_config is not None:
            self.core.read_modulators(modulators_config, sim, cells, p)
            self.modulators = True

        else:
            self.modulators = False

    # ..................{ READERS                           }..................
    @type_check
    def read_gene_config(self, phase: SimPhase) -> None:

        # Localize high-level phase objects for convenience.
        cells = phase.cells
        p     = phase.p
        sim   = phase.sim

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
        self.core.tissue_init(phase, substances_config)

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
            self.core.read_transporters(transporters_config, phase)
            self.core.write_transporters(sim, cells, p)

            self.transporters = True
        else:
            self.transporters = False

        # initialize channels, if desired:
        if channels_config is not None:
            self.core.read_channels(channels_config, phase)
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

        optim_exists = config_dic.get('optimization', None)
        if optim_exists is not None:

            # after primary initialization, check and see if optimization required:
                # after primary initialization, check and see if optimization required:
            opti = config_dic['optimization']['optimize network']
            self.core.opti_N = config_dic['optimization']['optimization steps']
            self.core.opti_method = config_dic['optimization']['optimization method']
            self.core.target_vmem = float(config_dic['optimization']['target Vmem'])
            self.core.opti_T = float(config_dic['optimization']['optimization T'])
            self.core.opti_step = float(config_dic['optimization']['optimization step'])
            # self.core.opti_run = config_dic['optimization']['run from optimization']

            if opti:
                logs.log_info('Analyzing gene network for optimal rates...')
                self.core.optimizer(sim, cells, p)
                self.reinitialize(phase)

    # ..................{ RUNNERS                           }..................
    @type_check
    def run_core_sim(self, phase: SimPhase) -> None:
        '''
        Simulate this gene regulatory network (GRN) in an isolated manner
        ignoring all other biophysicality (e.g., bioelectricity, fluid flow).

        Parameters
        ----------
        phase : SimPhase
            Current simulation phase.
        '''

        # Localize high-level phase objects for convenience.
        cells = phase.cells
        p     = phase.p
        sim   = phase.sim

        # set molecules to not affect charge for sim-grn test-drives:
        p.substances_affect_charge = False

        #FIXME: This... is quite unfortunate. Ideally, core parameters in the
        #"Parameters" object should *NEVER* be modified, as doing so actually
        #modifies the underlying YAML dictionary structure. If the current
        #simulation configuration file is saved to disk, the user-defined
        #values for these parameters will be overwritten by the new hacky
        #values defined below. Which would be bad. Is there no alternative?

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
        self.reinitialize(phase)

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

        # if p.use_microtubules:
        #     sim.mtubes.reinit(cells, p)
        #
        #     self.mtubes_x_time = []
        #     self.mtubes_y_time = []
        #
        #     if self.reset_MT:
        #         logs.log_info(
        #             'Resetting microtubules for sim-grn simulation...')
        #         sim.mtubes = Mtubes(sim, cells, p)

        for t in tt:
            if self.transporters:
                self.core.run_loop_transporters(t, sim, cells, p)

            self.core.run_loop(phase=phase, t=t)

            # if p.use_microtubules: # update the microtubules:
            #     sim.mtubes.update_mtubes(cells, sim, p)

            # If...
            if (
                # The simulation phase is being run...
                p.grn_runmodesim and
                # The cutting event has yet to be performed here...
                not self.mod_after_cut and
                # The current time step is at least the time step at which the
                # cutting event is scheduled to occur...
                t >= p.event_cut_time
            ):
                phase.dyna.fire_events(phase=phase, t=t)

                #FIXME: This logic should ideally reside somewhere in the
                #"betse.science.tissue" subpackage. For example, similar logic
                #already exists in the
                #betse.science.tissue.tishandler.TissueHandler._sim_events_tissue()
                #method called by the above phase.dyna.fire_events() call. It
                #should thus be possible to shift the body of the following if
                #conditional into the body of the similar if conditional in the
                #TissueHandler._sim_events_tissue() method. After doing so,
                #this outer if conditional and the "mod_after_cut" variable
                #would then be safely removable. Overlord, unite your power!

                # If a cutting event has just been run...
                if (phase.dyna.event_cut.is_fired and
                    not self.mod_after_cut):
                    self.core.mod_after_cut_event(
                        phase,
                        sim.target_inds_cell_o,
                        sim.target_inds_mem_o,
                    )
                    self.core.redefine_dynamic_dics(sim, cells, p)

                    logs.log_info(
                        "Reinitializing the gene regulatory network for simulation...")
                    self.reinitialize(phase)
                    # self.core.clear_cache()
                    # sim.uxmt, self.uymt = sim.mtubes.mtubes_to_cell(cells, p)

                    self.mod_after_cut = True  # set the boolean to avoid repeat action

            if t in tsamples:
                sim.time.append(t)

                logs.log_info('------------------' + str(np.round(t,3)) +' s --------------------')
                self.time.append(t)
                self.core.write_data(sim, cells, p)
                self.core.report(sim, p)

                # if p.use_microtubules:
                #     # microtubules:
                #     self.mtubes_x_time.append(sim.mtubes.mtubes_x * 1)
                #     self.mtubes_y_time.append(sim.mtubes.mtubes_y * 1)

        logs.log_info('Saving simulation...')
        datadump = [self, cells, p]
        fh.saveSim(p.grn_pickle_filename, datadump)
        self.core.init_saving(cells, p, plot_type='grn', nested_folder_name='RESULTS')

        # # microtubules plot------------------------------------------------------------------------
        # if p.use_microtubules:
        #
        #     logs.log_info("Plotting microtubules used in GRN simulation...")
        #
        #     plt.figure()
        #     ax = plt.subplot(111)
        #
        #     umtx, umty = sim.mtubes.mtubes_to_cell(cells, p)
        #
        #     plt.figure()
        #     ax = plt.subplot(111)
        #
        #     viz.plotVectField(
        #         umtx,
        #         umty,
        #         cells,
        #         p,
        #         plot_ecm=False,
        #         title='Final Microtubule Alignment Field',
        #         cb_title='Aligned MT Fraction',
        #         colorAutoscale=False,
        #         minColor=0.0,
        #         maxColor=1.0,
        #     )
        #
        #     # viz.mem_quiver(
        #     #     sim.mtubes.mtubes_x,
        #     #     sim.mtubes.mtubes_y,
        #     #     ax,
        #     #     cells,
        #     #     p,
        #     # )
        #
        #     ax.set_xlabel('X-Distance [um]')
        #     ax.set_ylabel('Y-Distance [um]')
        #     ax.set_title('Microtubule arrangement in cells')
        #
        #     if p.autosave is True:
        #         savename = self.core.imagePath + 'Microtubules' + '.png'
        #         plt.savefig(savename, format='png', transparent=True)
        #
        #     if p.turn_all_plots_off is False:
        #         plt.show(block=False)

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
