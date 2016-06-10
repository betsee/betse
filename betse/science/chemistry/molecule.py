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

class MasterOfMolecules(object):

    def __init__(self, sim, config_settings, p):

        """
         Initializes the MasterOfMolecules object.

        sim:                An instance of simulator
        config_settings     List of dictionaries storing key settings (p.molecules_config)
        p                   An instance of params
        """

        self.init(sim, config_settings, p)

    def init(self, sim, config_settings, p):
        """
            Initializes all of the main variables for all molecules included in the simulation.
            'config settings' is typically = p.molecules_config

        """
        # Initialize a list that will keep track of molecule names in the simulation
        self.molecule_names = []

        for mol_dic in config_settings:
            # get each user-defined name-filed in the dictionary:
            name = str(mol_dic['name'])

            # add the name to the name catalogue:
            self.molecule_names.append(name)

            # add a field to the MasterOfMolecules corresponding to Molecule object with that name
            setattr(self, name, Molecule())

            # now set the attributes of that Molecule object with the cornucopia of variables:

            # get MasterOfMolecules.name
            obj = getattr(self, name)

            # assign general properties
            obj.name = name   # let object know who it is
            obj.Dm = mol_dic['Dm']  # membrane diffusion coefficient [m2/s]
            obj.Do = mol_dic['Do']  # free diffusion constant in extra and intracellular spaces [m2/s]
            obj.z = mol_dic['z']  # charge (oxidation state)
            obj.c_envo = mol_dic['env conc']  # initial concentration in the environment [mmol/L]
            obj.c_cello = mol_dic['cell conc']  # initial concentration in the cytoplasm [mmol/L]

            # factors involving auto-catalytic growth and decay in the cytoplasm
            gad = mol_dic['growth and decay']

            obj.r_production = gad['production rate']
            obj.r_decay = gad['decay rate']
            obj.Kgd = gad['Km']

            # assign ion channel gating properties
            icg = mol_dic['ion channel gating']

            gating_ion_o = icg['ion channel target']  # get a target ion label to gate membrane to (or 'None')

            if gating_ion_o != 'None':
                obj.use_gating_ligand = True
                self.gating_ion = sim.get_ion(gating_ion_o)

            else:
                obj.use_gating_ligand = False
                self.gating_ion = []

            obj.gating_Hill_K = icg['target Hill coefficient']
            obj.gating_Hill_n = icg['target Hill exponent']
            obj.gating_max_val = icg['peak channel opening']
            obj.gating_extracell = icg['acts extracellularly']

            # assign active pumping properties
            ap = mol_dic['active pumping']
            obj.use_pumping = ap['turn on']
            obj.pump_to_cell = ap['pump to cell']
            obj.pump_max_val = ap['maximum rate']
            obj.pump_Km = ap['pump Km']

            # assign boundary change event properties
            cab = mol_dic['change at bounds']
            obj.change_at_bounds = cab['event happens']
            obj.change_bounds_start = cab['change start']
            obj.change_bounds_end = cab['change finish']
            obj.change_bounds_rate = cab['change rate']
            obj.change_bounds_target = cab['concentration']

            # assign plotting properties
            pd = mol_dic['plotting']
            obj.make_plots = pd['plot 2D']
            obj.make_ani = pd['animate']
            obj.plot_autoscale = pd['autoscale colorbar']
            obj.plot_max = pd['max val']
            obj.plot_min = pd['min val']

            # create data structures to use with sim --------
            # initialize concentrations in cells:
            obj.c_cells = np.ones(sim.cdl)*obj.c_cello
            obj.c_mems = np.ones(sim.mdl) * obj.c_cello

            # initialize concentration in the environment:
            if p.sim_ECM is False:
                obj.c_env = np.ones(sim.mdl) * obj.c_envo
            else:
                obj.c_env = np.ones(sim.edl) * obj.c_envo

            # initialize concentration at the boundary
            obj.c_bound = obj.c_envo

    def init_saving(self, cells, p, plot_type = 'init'):

            # init files
            if p.autosave is True:

                if plot_type == 'sim':
                    results_path = os.path.join(p.sim_results, 'Molecules')
                    p.plot_type = 'sim'

                elif plot_type == 'init':
                    results_path = os.path.join(p.init_results, 'Molecules')
                    p.plot_type = 'init'

                self.resultsPath = os.path.expanduser(results_path)
                os.makedirs(self.resultsPath, exist_ok=True)

                self.imagePath = os.path.join(self.resultsPath, 'fig_')

            # check that the plot cell is in range of the available cell indices:
            if p.plot_cell not in cells.cell_i:
                raise BetseExceptionParameters(
                    'The "plot cell" defined in the "results" section of your '
                    'configuration file does not exist in your cluster. '
                    'Choose a plot cell number smaller than the maximum cell number.')

    def run_loop(self, t, sim, cells, p):
        """
        Runs the main simulation loop steps for each of the molecules included in the simulation.

        """

        # get the name of the specific substance:
        for name in self.molecule_names:

            obj = getattr(self,name)

            obj.pump(sim, cells, p)
            obj.transport(sim, cells, p)
            obj.updateIntra(sim, cells, p)

            if p.run_sim is True:
                obj.gating(sim, cells, p)
                obj.update_boundary(t, p)

    def mod_after_cut_event(self,target_inds_cell, target_inds_mem, sim, cells, p):

        # get the name of the specific substance:
        for name in self.molecule_names:
            obj = getattr(self, name)

            obj.remove_cells(target_inds_cell, target_inds_mem, sim, cells, p)

    def clear_cache(self):
        """
        Initializes or clears the time-storage vectors at the begining of init and sim runs.

        """

        # get the name of the specific substance:
        for name in self.molecule_names:
            obj = getattr(self, name)

            obj.c_mems_time = []
            obj.c_cells_time = []
            obj.c_env_time = []

    def write_data(self, sim, p):
        """
        Writes concentration data from a time-step to time-storage vectors.

        """

        # get the name of the specific substance:
        for name in self.molecule_names:
            obj = getattr(self, name)

            obj.c_mems_time.append(obj.c_mems)
            obj.c_cells_time.append(obj.c_cells)
            obj.c_env_time.append(obj.c_env)

    def report(self, sim, p):
        """
        At the end of the simulation, tell user about mean, final concentrations of each molecule.

        """

        # get the name of the specific substance:
        for name in self.molecule_names:

            obj = getattr(self, name)

            logs.log_info('Final average concentration of ' + str(name) + ' in the cell: ' +
                                           str(np.round(1.0e3*obj.c_cells.mean(), 4)) + ' umol/L')

            logs.log_info('Final average concentration of ' + str(name) + ' in the environment: ' +
                                          str(np.round(1.0e3*obj.c_env.mean(), 4)) + ' umol/L')

    def export_all_data(self, sim, cells, p):

        """

        Exports concentration data from each molecule to a file for a single cell
        (plot cell defined in params) as a function of time.

        """
        logs.log_info('Exporting raw data for auxiliary molecules...')
        # get the name of the specific substance:
        for name in self.molecule_names:
            obj = getattr(self, name)

            obj.export_data(sim, cells, p, self.resultsPath)

    def plot(self, sim, cells, p):
        """
        Creates plots for each molecule included in the simulation.

        """

        logs.log_info('Plotting data for auxiliary molecules...')

        # get the name of the specific substance:
        for name in self.molecule_names:
            obj = getattr(self, name)

            if p.plot_single_cell_graphs:
            # create line graphs for the substance
                obj.plot_1D(sim, p, self.imagePath)

            # create 2D maps for the substance
            obj.plot_cells(sim, cells, p, self.imagePath)

            # if there's a real environment, plot 2D concentration in the environment
            if p.sim_ECM:

                obj.plot_env(sim, cells, p, self.imagePath)

    def anim(self, sim, cells, p):
        """
        Animates 2D data for each molecule in the simulation.

        """

        logs.log_info('Animating data for auxiliary molecules...')
        # get the name of the specific substance:
        for name in self.molecule_names:
            obj = getattr(self, name)

            if p.createAnimations:

                # create 2D animations for the substance in cells
                obj.anim_cells(sim, cells, p)

                # create 2D animations for the substance in the environment
                if p.sim_ECM:

                    obj.anim_env(sim, cells, p)

class Molecule(object):

    def transport(self, sim, cells, p):
        """
        Transports the molecule across the membrane,
        through gap junctions, and if p.sim_ECM is true,
        through extracellular spaces and the environment.

        """

        self.c_mems, self.c_env, _, _, _, _ = stb.molecule_mover(sim, self.c_mems,
                                                                self.c_env, cells, p, z=self.z, Dm = self.Dm,
                                                                Do = self.Do, c_bound = self.c_bound,
                                                                ignoreECM = True)

    def updateC(self, flux, sim, cells, p):
        """

        General updater for a flux defined on membranes and updating concentrations in
        cells and environment.

        """
        self.c_mems, self.c_env = stb.update_Co(sim, self.c_mems, self.c_cells, flux, cells, p, ignoreECM=False)

    def updateIntra(self, sim, cells, p):

        self.c_mems, self.c_cells, _ = stb.update_intra(sim, cells, self.c_mems, self.c_cells, self.Do, self.z, p)

    def pump(self, sim, cells, p):

        """
        Defines a generic active transport pump that can be used to move
        a general molecule (such as serotonin or glutamate)
        into or out of the cell by active transport.

        Works on the basic premise of enzymatic pumps defined elsewhere:

        pump_out is True:

        cX_cell + cATP  -------> cX_env + cADP + cPi

        pump_out is False:

        cX_env + cATP  <------- cX_cell + cADP + cPi

        """

        if self.use_pumping:

            self.c_mems, self.c_env, _ = stb.molecule_pump(sim, self.c_mems, self.c_env,
                                                                 cells, p, Df=self.Do, z=self.z,
                                                                 pump_into_cell=self.pump_to_cell,
                                                                 alpha_max=self.pump_max_val, Km_X=self.pump_Km,
                                                                 Km_ATP=1.0)

    def gating(self, sim, cells, p):
        """
        Uses the molecule concentration to open an ion channel in the cell membranes.

        """

        # update membrane permeability if dye targets an ion channel:
        if self.use_gating_ligand:

            if self.gating_extracell is False:

                sim.Dm_mod_dye = sim.rho_channel*self.gating_max_val*tb.hill(self.c_mems,
                                                                        self.gating_Hill_K,self.gating_Hill_n)

                sim.Dm_morpho[self.gating_ion] = sim.Dm_mod_dye


            elif self.gating_extracell is True and p.sim_ECM is True:

                sim.Dm_mod_dye = self.gating_max_val*tb.hill(self.c_env,self.gating_Hill_K,self.gating_Hill_n)

                sim.Dm_morpho[self.gating_ion] = sim.rho_channel*sim.Dm_mod_dye[cells.map_mem2ecm]

    def reaction(self, sim, cells, p):
        """
        Grows and/or decays the molecule concentration in the cell cytosol using a simple rate equation
        representing saturating autocatalytic production/decay via Michaelis-Menten kinetics.

        """
        cc = self.c_cells/self.Kgd

        delta_cells = self.r_production*(cc/(1 + cc)) - self.r_decay*(cc/(1 + cc))

        self.c_cells = self.c_cells + delta_cells*p.dt

        # make sure the concs inside the cell are evenly mixed after production/decay:
        self.updateIntra(sim, cells, p)

    def remove_cells(self, target_inds_cell, target_inds_mem, sim, cells, p):
        """
        During a cutting event, removes the right cells from the simulation network,
        while preserving additional information.

        """

        # remove cells from the cell concentration list:
        ccells2 = np.delete(self.c_cells, target_inds_cell)
        # reassign the new data vector to the object:
        self.c_cells = ccells2

        # remove cells from the mems concentration list:
        cmems2 = np.delete(self.c_mems, target_inds_mem)
        # reassign the new data vector to the object:
        self.c_mems = cmems2

    def update_boundary(self, t, p):
        """
        Run a dynamic event in the sim, which alters concentration at the global boundary.

        t:          simulation world time
        p:          parameters instance

        """

        if self.change_at_bounds:

            effector_MorphEnv = tb.pulse(t,self.change_bounds_start,self.change_bounds_end,self.change_bounds_rate)

            if p.sim_ECM is False:
                self.c_env[:] = self.conc_MorphEnv*effector_MorphEnv + self.c_envo*(1-effector_MorphEnv)

            elif p.sim_ECM is True: # simulate addition of counter salt to maintain charge neutrality:
                self.c_bound = self.conc_MorphEnv*effector_MorphEnv + self.c_envo*(1-effector_MorphEnv)

    def export_data(self, sim, cells, p, savePath):

        saveName = 'ExportData_' + self.name + '.csv'

        saveData = os.path.join(savePath, saveName)

        ci = p.plot_cell  # index of cell to get time-dependent data for

        # create the header, first entry will be time:
        headr = 'time_s' + ','

        ccell = [arr[ci] for arr in self.c_cells_time]

        headr = headr + 'Cell_Conc_' + self.name + '_mmol/L' + ','

        if p.sim_ECM is True:

            cenv = [obj_cenv[cells.map_cell2ecm][ci] for obj_cenv in self.c_env_time]

        else:

            cenv = [np.dot(cells.M_sum_mems, obj_cenv) / cells.num_mems for obj_cenv in self.c_env_time]

        headr = headr + 'Env_Conc_' + self.name + '_mmol/L' + ','

        time = np.asarray(sim.time)
        ccell = np.asarray(ccell)
        cenv = np.asarray(cenv)

        dataM = np.column_stack((time, ccell, cenv))

        np.savetxt(saveData, dataM, delimiter=',', header=headr)

    def plot_1D(self, sim, p, saveImagePath):
        """
        Create 1D plot of concentration in cell and environment for a single cell (params plot cell)
        as a function of simulation time.

        """

        c_cells = [1.0e3*arr[p.plot_cell] for arr in self.c_cells_time]
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(sim.time, c_cells)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Concentration [umol/L]')
        ax.set_title('Concentration of ' + self.name + ' in cell ' + str(p.plot_cell))

        if p.autosave is True:
            savename = saveImagePath + 'CellConcentration_' + self.name + '.png'
            plt.savefig(savename, format='png', transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    def plot_cells(self, sim, cells, p, saveImagePath):
        """
        Create 2D plot of cell concentrations.

        """

        fig, ax, cb = viz.plotPrettyPolyData(self.c_mems*1e3,
            sim, cells, p,
            number_cells=p.enumerate_cells,
            clrAutoscale=self.plot_autoscale,
            clrMin=self.plot_min,
            clrMax=self.plot_max,
            clrmap=p.default_cm)

        ax.set_title('Final ' + self.name + ' Concentration in Cells')
        ax.set_xlabel('Spatial distance [um]')
        ax.set_ylabel('Spatial distance [um]')
        cb.set_label('Concentration umol/L')

        if p.autosave is True:
            savename = saveImagePath + 'cell_conc_' + self.name + '.png'
            plt.savefig(savename,format='png', transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    def plot_env(self, sim, cells, p, saveImagePath):
        """
        Create 2D plot of environmental concentration.

        """

        fig = plt.figure()
        ax = plt.subplot(111)

        dyeEnv = (self.c_env*1e3).reshape(cells.X.shape)

        xmin = cells.xmin*p.um
        xmax = cells.xmax*p.um
        ymin = cells.ymin*p.um
        ymax = cells.ymax*p.um

        bkgPlot = ax.imshow(dyeEnv,origin='lower',extent=[xmin,xmax,ymin,ymax], cmap=p.default_cm)


        if self.plot_autoscale is False:
            bkgPlot.set_clim(self.plot_min, self.plot_max)

        cb = fig.colorbar(bkgPlot)

        ax.axis('equal')

        ax.axis([xmin,xmax,ymin,ymax])

        ax.set_title('Final ' + self.name + ' Concentration in Environment')
        ax.set_xlabel('Spatial distance [um]')
        ax.set_ylabel('Spatial distance [um]')
        cb.set_label('Concentration umol/L')

        if p.autosave is True:
            savename = saveImagePath + 'env_conc_' + self.name + '.png'
            plt.savefig(savename,format='png',dpi = 300.0, transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    def anim_cells(self, sim, cells, p):

        """
        Create 2D animation of cell concentration.

        """

        type_name = self.name + '_cells'
        fig_tit = 'Cytosolic ' + self.name

        AnimCellsTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=[1e3*arr for arr in self.c_mems_time],
            type=type_name,
            figure_title=fig_tit,
            colorbar_title='Concentration [umol/L]',
            is_color_autoscaled=self.plot_autoscale,
            color_min=self.plot_min,
            color_max=self.plot_max)

    def anim_env(self, sim, cells, p):

        """
        Create 2D animation of env concentration.

        """

        type_name = self.name + '_env'
        fig_tit = 'Environmental ' + self.name

        env_time_series = [env.reshape(cells.X.shape)*1e3 for env in self.c_env_time]
        AnimEnvTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=env_time_series,
            type=type_name,
            figure_title=fig_tit,
            colorbar_title='Concentration [umol/L]',
            is_color_autoscaled=self.plot_autoscale,
            color_min=self.plot_min,
            color_max=self.plot_max)

