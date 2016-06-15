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


# FIXME we will need methods to update mitochondria concs and Vmit, as well as ER concs and Vmit defined
# in both MasterOfMolecules and Molecule...perhaps 'update_organelle()'

# FIXME see if we need rates of all reactions...

class MasterOfMolecules(object):

    def __init__(self, sim, config_substances, p):

        """
         Initializes the MasterOfMolecules object.

        sim:                An instance of simulator
        config_settings     List of dictionaries storing key settings (p.molecules_config)
        p                   An instance of params
        """
        # all metabolic simulations require ATP, ADP and Pi. Initialize these fields to None so that we can test
        # for their pressence in metabolic sims:
        self.ATP = None
        self.ADP = None
        self.Pi = None

        self.read_substances(sim, config_substances, p)

    def read_substances(self, sim, config_substances, p):
        """
            Initializes all of the main variables for all molecules included in the simulation.

            config_substances:  dictionary containing BETSE biomolecule template fields

        """
        # Initialize a list that will keep track of molecule names in the simulation
        self.molecule_names = []
        # Initialize a list that keeps track of molecule index
        self.molecule_index = []

        for q, mol_dic in enumerate(config_substances):
            # get each user-defined name-filed in the dictionary:
            name = str(mol_dic['name'])

            # add the name to the name catalogue:
            self.molecule_names.append(name)
            self.molecule_index.append(q)

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

            obj.c_mito = mol_dic.get('mit conc',None)  # initialized to None if optional fields not present
            obj.c_ero = mol_dic.get('er conc', None)

            # factors involving auto-catalytic growth and decay in the cytoplasm
            gad = mol_dic.get('growth and decay', None)

            if gad != 'None' and gad is not None:

                obj.simple_growth = True

                obj.r_production = gad['production rate']
                obj.r_decay = gad['decay rate']
                obj.Kgd = gad['Km']
                obj.n_production = gad['n']

                obj.growth_activators_list = gad['activators']
                obj.growth_activators_k = gad['k activators']
                obj.growth_activators_Km = gad['Km activators']

                obj.growth_activators_n = gad['k activators']
                obj.growth_inhibitors_list = gad['inhibitors']
                obj.growth_inhibitors_k = gad['k inhibitors']
                obj.growth_inhibitors_Km = gad['Km inhibitors']
                obj.growth_inhibitors_n = gad['k inhibitors']

            else:
                obj.simple_growth = False

            # assign ion channel gating properties
            icg = mol_dic.get('ion channel gating', None)

            if icg is not None:

                obj.ion_channel_gating = True

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

            else:

                obj.ion_channel_gating = False

            # assign active pumping properties
            ap = mol_dic.get('active pumping', None)

            if ap is not None:

                obj.active_pumping = True
                obj.use_pumping = ap['turn on']
                obj.pump_to_cell = ap['pump to cell']
                obj.pump_max_val = ap['maximum rate']
                obj.pump_Km = ap['pump Km']
                obj.pumps_use_ATP = ap['uses ATP']

            else:
                obj.active_pumping = False


            # assign boundary change event properties
            cab = mol_dic.get('change at bounds', None)

            if cab is not None:
                obj.change_bounds = True
                obj.change_at_bounds = cab['event happens']
                obj.change_bounds_start = cab['change start']
                obj.change_bounds_end = cab['change finish']
                obj.change_bounds_rate = cab['change rate']
                obj.change_bounds_target = cab['concentration']

            else:
                obj.change_bounds = False

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

            # if there is an initial concentration for mitochondria, define a conc vector for it:
            if obj.c_mito is not None:

                obj.c_mit = np.ones(sim.cdl)*obj.c_mito

            # if there is an initial concentration for endo retic, define a conc vector for it:
            if obj.c_ero is not None:

                obj.c_er = np.ones(sim.cdl)*obj.c_ero

            # initialize concentration in the environment:
            if p.sim_ECM is False:
                obj.c_env = np.ones(sim.mdl) * obj.c_envo
            else:
                obj.c_env = np.ones(sim.edl) * obj.c_envo

            # initialize concentration at the boundary
            obj.c_bound = obj.c_envo

    def read_reactions(self, config_reactions, sim, cells, p):

        """
          Read in and initialize parameters for all user-defined reactions.

          config_options:  dictionary containing BETSE reaction template fields

          """

        # Initialize a list that will keep track of reaction names in the simulation
        self.reaction_names = []

        # Initialize a list that keeps the index of the reaction:
        self.reaction_index = []

        for q, react_dic in enumerate(config_reactions):
            # get each user-defined name-filed in the dictionary:
            name = str(react_dic['name'])

            # add the name to the name catalogue:
            self.reaction_names.append(name)
            self.reaction_index.append(q)

            # add a field to the MasterOfReactions corresponding to Reaction object with that name:
            #self._reactions[name] = Reaction()
            setattr(self, name, Reaction())

            # now set the attributes of that Reaction object with the cornucopia of variables:

            # get MasterOfMolecules.name
            obj = getattr(self, name)

            # assign general properties
            obj.name = name  # let object know who it is

            # list where the reaction takes place; if field not specified default to 'cell':
            obj.reaction_zone = react_dic.get('reaction zone', 'cell')

            # set the main fields of the reaction object
            obj.reactants_list = react_dic['reactants']
            obj.reactants_coeff = react_dic['reactant multipliers']

            obj.products_list = react_dic['products']
            obj.products_coeff = react_dic['product multipliers']

            obj.Km_reactants_list = react_dic['Km reactants']
            obj.Km_products_list = react_dic['Km products']
            obj.vmax = float(react_dic['max rate'])

            obj.delta_Go = react_dic['standard free energy']

            if obj.delta_Go == 'None':
                obj.delta_Go = None     # make the field a proper None variable

            else:
                obj.delta_Go = float(obj.delta_Go)

            obj.reaction_activators_list = react_dic.get('reaction activators', None)
            obj.reaction_activators_Km = react_dic.get('activator Km', None)
            obj.reaction_activators_n = react_dic.get('activator n', None)
            obj.reaction_inhibitors_list = react_dic.get('reaction inhibitors', None)
            obj.reaction_inhibitors_Km = react_dic.get('inhibitor Km', None)
            obj.reaction_inhibitors_n = react_dic.get('inhibitor n', None)

            # now we want to load the right concentration data arrays for reactants into the Reaction object:

            # check the reaction states to make sure all reagents/products are named biomolecules or ions:
            # self.check_reactions(obj, sim, cells, p)

    def read_transporters(self, config_transporters, sim, cells, p):

        """
            Read in and initialize parameters for all user-defined transporters.

            config_options:  dictionary containing BETSE transporter template fields

            """

        # Initialize a list that will keep track of reaction names in the simulation
        self.transporter_names = []

        # Initialize a list that keeps the index of the reaction:
        self.transporter_index = []

        for q, trans_dic in enumerate(config_transporters):
            # get each user-defined name-filed in the dictionary:
            name = str(trans_dic['name'])

            # add the name to the name catalogue:
            self.transporter_names.append(name)
            self.transporter_index.append(q)

            # add a field to the MasterOfReactions corresponding to Transporter object with that name:

            setattr(self, name, Transporter())

            # now set the attributes of that Transporter object with the cornucopia of variables:

            # get MasterOfMolecules.name
            obj = getattr(self, name)

            # assign general properties
            obj.name = name  # let object know who it is

            # list where the reaction takes place; if field not specified default to 'cell':
            obj.reaction_zone = trans_dic.get('reaction zone', 'cell')

            # set the main fields of the reaction object
            obj.reactants_list = trans_dic['reactants']
            obj.reactants_coeff = trans_dic['reactant multipliers']

            obj.products_list = trans_dic['products']
            obj.products_coeff = trans_dic['product multipliers']

            obj.Km_reactants_list = trans_dic['Km reactants']
            obj.Km_products_list = trans_dic['Km products']

            obj.transport_out_list = trans_dic['transfered out of cell']
            obj.transport_in_list = trans_dic['transfered into cell']

            obj.vmax = float(trans_dic['max rate'])

            obj.delta_Go = trans_dic['standard free energy']

            if obj.delta_Go == 'None':
                obj.delta_Go = None  # make the field a proper None variable

            else:
                obj.delta_Go = float(obj.delta_Go)

            obj.transporter_activators_list = trans_dic.get('transporter activators', None)
            obj.transporter_activators_Km = trans_dic.get('activator Km', None)
            obj.transporter_activators_n = trans_dic.get('activator n', None)
            obj.transporter_inhibitors_list = trans_dic.get('transporter inhibitors', None)
            obj.transporter_inhibitors_Km = trans_dic.get('inhibitor Km', None)
            obj.transporter_inhibitors_n = trans_dic.get('inhibitor n', None)

            # now we want to load the right concentration data arrays for reactants into the Reaction object:

            # check the reaction states to make sure all reagents/products are named biomolecules or ions:
            # self.check_reactions(obj, sim, cells, p)

    def check_reactions(self, obj, sim, cells, p):

        # FIXME complete this

        if obj.reaction_zone == 'cell':

            # Case 1: no transfer across the membrane -- these are in-cytoplasm reactions -------------------------------
            for i, reactant_name in enumerate(obj.reactants_list):

                self.set_react_concs(obj, sim, cells, p, reactant_name, i, 'c_cells', 'cc_cells')

            # Now load the right concentration data arrays for products into the Reaction object:
            for j, product_name in enumerate(obj.products_list):

                self.set_prod_concs(obj, sim, cells, p, product_name, j, 'c_cells', 'cc_cells')


        elif obj.reaction_zone == 'mitochondria':

            # Case 1: no transfer across the membrane -- these are in-cytoplasm reactions -------------------------------
            for i, reactant_name in enumerate(obj.reactants_list):

                self.set_react_concs(obj, sim, cells, p, reactant_name, i, 'c_mit', 'cc_mit')

            # Now load the right concentration data arrays for products into the Reaction object:
            for j, product_name in enumerate(obj.products_list):

                self.set_prod_concs(obj, sim, cells, p, product_name, j, 'c_mit', 'cc_mit')

        else:
            raise BetseExceptionParameters("You have requested a reaction zone that does not exist."
                                           "Valid options include: 'cell' and 'mitochondria'."
                                           "Please check your config(s) settings defining reactions"
                                           " and try again. ")

    def set_react_concs(self, obj, sim, cells, p, reactant_name, i, reactant_type_self, reactant_type_sim):
        # FIXME put in a clause for c_mit or cc_mit undefined!

        """

        obj:                    Molecule field of MasterOfMolecules
        sim:                    Instance of BETSE sim
        reactant_type_self:     Concentration field of Molecule: 'c_cells', 'c_env', 'c_mit' or 'c_er'
        reactant_type_sim:      Concentration field of sim: 'cc_cells', 'cc_env',  'cc_mit' or 'cc_er'

        """

        label = 'i' + reactant_name
        ion_check = getattr(sim, label, None)

        if ion_check is None:

            obj.reactant_source_object.append(id(self))
            obj.reactant_source_type.append(reactant_type_self)

        else:

            obj.reactant_source_object.append(id(sim))
            obj.reactant_source_type.append(reactant_type_sim)

        # add the index of the reaction to the list so we can access modifiers like reaction coefficients
        # and Km values at a later point:
        obj.inds_react.append(i)

    def set_prod_concs(self, obj, sim, cells, p, product_name, j, product_type_self, product_type_sim):

        # Now load the right concentration data arrays for products into the Reaction object:
        # see if the name is an ion defined in sim:
        label = 'i' + product_name
        ion_check = getattr(sim, label, None)

        if ion_check is None:

            obj.product_source_object.append(id(self))
            obj.product_source_type.append(product_type_self)

        else:

            obj.product_source_object.append(id(sim))
            obj.product_source_type.append(product_type_sim)

        # add the index of the reaction to the list so we can access modifiers like reaction coefficients
        # and Km values at a later point:
        obj.inds_prod.append(j)

    def assign_new_concentrations(self, obj, delta_c, new_reactants, new_products, sim, cells, p):

        for i, react_name in enumerate(obj.reactants_list):

            if obj.reactant_source_object[i] == id(self):

                source_obj = getattr(self, react_name)

                setattr(source_obj, obj.reactant_source_type[i], new_reactants[i])

            elif obj.reactant_source_object[i] == id(sim):

                source_obj = getattr(sim, obj.reactant_source_type[i])
                ion_label = 'i' + react_name
                ion_type = getattr(sim, ion_label)

                source_obj[ion_type] = new_reactants[i]

        for j, prod_name in enumerate(obj.products_list):

                if obj.product_source_object[j] == id(self):

                    source_obj = getattr(self, prod_name)

                    setattr(source_obj, obj.product_source_type[j], new_products[j])

                elif obj.product_source_object[j] == id(sim):

                    source_obj = getattr(sim, obj.product_source_type[j])
                    ion_label = 'i' + prod_name
                    ion_type = getattr(sim, ion_label)

                    source_obj[ion_type] = new_products[j]

    def init_saving(self, cells, p, plot_type = 'init', nested_folder_name = 'Molecules'):

            # init files
            if p.autosave is True:

                if plot_type == 'sim':
                    results_path = os.path.join(p.sim_results, nested_folder_name)
                    p.plot_type = 'sim'

                elif plot_type == 'init':
                    results_path = os.path.join(p.init_results, nested_folder_name)
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

            # if pumping is enabled:
            if obj.active_pumping:
                obj.pump(sim, cells, p)

            # update the production-decay regulatory network (if defined):
            if obj.simple_growth is True:
                obj.growth_and_decay(self, sim, cells, p)

            if p.run_sim is True:
                # use the substance as a gating ligand (if desired)
                if obj.ion_channel_gating:
                    obj.gating(sim, cells, p)

                # update the global boundary (if desired)
                if obj.change_bounds:
                    obj.update_boundary(t, p)

            # transport the molecule through gap junctions and environment:
            obj.transport(sim, cells, p)

            # update the substance on the inside of the cell:
            obj.updateIntra(sim, cells, p)

    def updateInside(self, sim, cells, p):
        """
        Runs the main simulation loop steps for each of the molecules included in the simulation.

        """

        # get the name of the specific substance:
        for name in self.molecule_names:

            obj = getattr(self, name)
            obj.updateIntra(sim, cells, p)

    def run_loop_reactions(self, t, sim, sim_metabo, cells, p):

        # get the object corresponding to the specific reaction:
        for i, name in enumerate(self.reaction_names):

            # get the Reaction object:
            obj = getattr(self, name)

            # compute the new reactants and products
            rate = obj.compute_reaction(sim, sim_metabo, cells, p)

    def run_loop_transporters(self, t, sim, sim_metabo, cells, p):

        # get the object corresponding to the specific transporter:
        for i, name in enumerate(self.transporter_names):

            # get the Reaction object:
            obj = getattr(self, name)

            # compute the new reactants and products
            rate = obj.compute_reaction(sim, sim_metabo, cells, p)

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

    def export_all_data(self, sim, cells, p, message = 'for auxiliary molecules...'):

        """

        Exports concentration data from each molecule to a file for a single cell
        (plot cell defined in params) as a function of time.

        """
        logs.log_info('Exporting raw data for ' + message)
        # get the name of the specific substance:
        for name in self.molecule_names:
            obj = getattr(self, name)

            obj.export_data(sim, cells, p, self.resultsPath)

    def plot(self, sim, cells, p, message = 'for auxiliary molecules...'):
        """
        Creates plots for each molecule included in the simulation.

        """

        logs.log_info('Plotting 1D and 2D data for ' + message)

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


        data_all1D = []
        fig_all1D = plt.figure()
        ax_all1D = plt.subplot(111)

        for name in self.molecule_names:
            obj = getattr(self, name)

            c_cells = [arr[p.plot_cell] for arr in obj.c_cells_time]

            ax_all1D.plot(sim.time, c_cells, linewidth = 2.0, label=name)

        legend = ax_all1D.legend(loc = 'upper right', shadow = False, frameon = False)

        ax_all1D.set_xlabel('Time [s]')
        ax_all1D.set_ylabel('Concentration [mmol/L]')
        ax_all1D.set_title('Concentration of all substances in cell ' + str(p.plot_cell))

        if p.autosave is True:
            savename = self.imagePath + 'AllCellConcentrations' + '.png'
            plt.savefig(savename, format='png', transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    def anim(self, sim, cells, p, message = 'for auxilary molecules...'):
        """
        Animates 2D data for each molecule in the simulation.

        """

        logs.log_info('Animating data for ' + message)
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

    def __init__(self):

        # Set all fields to None -- these will be dynamically set by MasterOfMolecules

        self.c_cello = None
        self.c_memo = None
        self.c_env = None
        self.z = None
        self.Dm = None
        self.Do = None
        self.c_bound = None
        self.Kgd = None
        self.use_pumping = None
        self.pumps_use_ATP = None
        self.pump_to_cell = None
        self.pump_max_val = None
        self.pump_Km = None
        self.use_gating_ligand = None
        self.gating_extracell = None
        self.gating_max_val = None
        self.gating_Hill_K = None
        self.gating_Hill_n = None
        self.gating_ion = None

        self.change_at_bounds = None
        self.change_bounds_start = None
        self.change_bounds_end = None
        self.change_bounds_rate = None
        self.change_bounds_target = None
        self.make_plots = None
        self.make_ani = None
        self.plot_autoscale = None
        self.plot_max = None
        self.plot_min = None

        self.r_production = None
        self.r_decay = None
        self.n_production = None

        self.growth_activators_list = None
        self.growth_activators_k = None
        self.growth_activators_Km = None
        self.growth_activators_n = None
        self.growth_inhibitors_list = None
        self.growth_inhibitors_k = None
        self.growth_inhibitors_Km = None
        self.growth_inhibitors_n = None

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
        self.c_mems, self.c_env = stb.update_Co(sim, self.c_mems, self.c_cells, flux, cells, p, ignoreECM=True)

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

            if self.pumps_use_ATP:

                if p.metabolism_enabled:
                    met_vect = sim.met_concs
                else:
                    met_vect = None

                self.c_mems, self.c_env, flux = stb.molecule_pump(sim, self.c_mems, self.c_env,
                                                                     cells, p, Df=self.Do, z=self.z,
                                                                     pump_into_cell=self.pump_to_cell,
                                                                     alpha_max=self.pump_max_val, Km_X=self.pump_Km,
                                                                     Km_ATP=1.0, met = met_vect)
                if p.metabolism_enabled:
                    # update ATP concentrations after pump action:
                    sim.metabo.update_ATP(flux, sim, cells, p)

            else:

                self.c_mems, self.c_env, flux = stb.molecule_transporter(sim, self.c_mems, self.c_env,
                    cells, p, Df=self.Do, z=self.z, pump_into_cell=self.pump_to_cell, alpha_max=self.pump_max_val,
                    Km_X=self.pump_Km, Keq= 1.0)

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

    def growth_and_decay(self, super_self, sim, cells, p):
        """
        Grows and/or decays the molecule concentration in the cell cytosol using a simple rate equation
        representing saturating autocatalytic production/decay via Michaelis-Menten kinetics.

        """

        cc = self.c_cells/self.Kgd

        activator_alpha, inhibitor_alpha = get_influencers_grn(sim, super_self, self.growth_activators_list,
            self.growth_activators_k, self.growth_activators_Km, self.growth_activators_n,
            self.growth_inhibitors_list, self.growth_inhibitors_k, self.growth_inhibitors_Km,
            self.growth_inhibitors_n, reaction_zone='cell')


        # RK4 method (appears identical to faster Euler method)
        # delta_cells = tb.RK4(lambda cc: self.r_production*inhibitor_alpha*activator_alpha - self.r_decay*cc)
        # c_change = delta_cells(self.c_cells, p.dt)
        # self.c_cells = self.c_cells + c_change
        # self.c_mems = self.c_mems + c_change[cells.mem_to_cells]

        delta_cells = self.r_production*inhibitor_alpha*activator_alpha - self.r_decay*cc
        self.c_cells = self.c_cells + delta_cells*p.dt

        self.c_mems = self.c_mems + delta_cells[cells.mem_to_cells]*p.dt

        # make sure the concs inside the cell are evenly mixed after production/decay:
        # self.updateIntra(sim, cells, p)

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

class Reaction(object):

    def __init__(self):

        # pre-populate the object with fields that will be assigned by MasterOfMolecules

        self.name = None
        self.reactants_list = None
        self.reactants_coeff = None
        self.products_list = None
        self.products_coeff = None
        self.Km_reactants_list = None
        self.Km_products_list = None
        self.vmax = None
        self.delta_Go = None

        self.reaction_zone = None

        self.c_reactants = []  # concentrations of reactions, defined on cell-centres
        self.z_reactants = []  # charge of reactants
        self.inds_react = []  # indices to each reactant, to keep their order
        self.c_products = []  # concentrations of reactions, defined on cell-centres
        self.z_products = []  # charge of products
        self.inds_prod = []  # indices to each reactant, to keep their order

        self.product_source_object = []  # object from which product concentrations come from
        self.product_source_type = []    # type of product concentrations sourced from object
        self.reactant_source_object = []  # object from which reactants concentrations come from
        self.reactant_source_type = []  # type of reactant concentrations sourced from object

        # activator and inhibitors of the reaction
        self.reaction_activators_list = None
        self.reaction_activators_Km = None
        self.reaction_activators_n = None
        self.reaction_inhibitors_list = None
        self.reaction_inhibitors_Km = None
        self.reaction_inhibitors_n = None

    def get_reactants(self, sim, sim_metabo, reactant_type_self, reactant_type_sim):

        """
        Get updated concentrations direct from sources for chemical reaction's reactants.

        sim:                    An instance of simulator
        sim_metabo:             Sim instance of reaction block (sim.metabo or sim.reacto for general reactions)
        reactant_type_self:     Data type to retrieve: c_cells, c_mems, c_mit
        reactant_type_sim:      Data typpe to retrieve: cc_cells, cc_mems, cc_mit

        """

        self.c_reactants = []

        for reactant_name in self.reactants_list:

            label = 'i' + reactant_name
            ion_check = getattr(sim, label, None)

            if ion_check is None:

                try:
                    obj_reactant = getattr(sim_metabo, reactant_name)
                    c_react = getattr(obj_reactant, reactant_type_self)
                    self.c_reactants.append(c_react)

                except KeyError:

                    raise BetseExceptionParameters('Name of product is not a defined chemical, or is not'
                                                   'an ion currently included in the ion profile being used.'
                                                   'Please check biomolecule definitions and ion profile'
                                                   'settings of your config(s) and try again.')

            else:

                # define the reactant as the ion concentration from the cell concentrations object in sim:
                sim_conco = getattr(sim, reactant_type_sim)

                sim_conc = sim_conco[ion_check]

                self.c_reactants.append(sim_conc)

    def get_products(self, sim, sim_metabo, product_type_self, product_type_sim):

        """
        Get updated concentrations direct from sources for chemical reaction's products.

        sim:                    An instance of simulator
        sim_metabo:             Sim instance of reaction block (sim.metabo or sim.reacto for general reactions)
        product_type_self:     Data type to retrieve: c_cells, c_mems, c_mit
        product_type_sim:      Data typpe to retrieve: cc_cells, cc_mems, cc_mit

        """

        self.c_products = []

        for product_name in self.products_list:

            label = 'i' + product_name
            ion_check = getattr(sim, label, None)

            if ion_check is None:

                try:
                    obj_prod = getattr(sim_metabo, product_name)
                    c_prod = getattr(obj_prod, product_type_self)
                    self.c_products.append(c_prod)

                except KeyError:

                    raise BetseExceptionParameters('Name of product is not a defined chemical, or is not'
                                                   'an ion currently included in the ion profile being used.'
                                                   'Please check biomolecule definitions and ion profile'
                                                   'settings of your config(s) and try again.')

            else:

                # define the reactant as the ion concentration from the cell concentrations object in sim:
                sim_conco = getattr(sim, product_type_sim)

                sim_conc = sim_conco[ion_check]

                self.c_products.append(sim_conc)

    def set_reactant_c(self, deltac, sim, sim_metabo, reactant_type_self, reactant_type_sim):

        for i, reactant_name in enumerate(self.reactants_list):

            conc = self.c_reactants[i] - deltac*self.reactants_coeff[i]

            label = 'i' + reactant_name
            ion_check = getattr(sim, label, None)

            if ion_check is None:

                obj_reactant = getattr(sim_metabo, reactant_name)
                setattr(obj_reactant, reactant_type_self, conc)

            else:
                # define the reactant as the ion concentration from the cell concentrations object in sim:
                sim_conc = getattr(sim, reactant_type_sim)

                sim_conc[ion_check] = conc

    def set_product_c(self, deltac, sim, sim_metabo, product_type_self, product_type_sim):

        for i, product_name in enumerate(self.products_list):

            conc = self.c_products[i] + deltac*self.products_coeff[i]

            label = 'i' + product_name
            ion_check = getattr(sim, label, None)

            if ion_check is None:

                obj_reactant = getattr(sim_metabo, product_name)
                setattr(obj_reactant, product_type_self, conc)

            else:
                # define the reactant as the ion concentration from the cell concentrations object in sim:
                sim_conc = getattr(sim, product_type_sim)

                sim_conc[ion_check] = conc

    def compute_reaction(self, sim, sim_metabo, cells, p):


        # get up-to-date concentration data for the reaction:
        if self.reaction_zone == 'cell':

            self.get_reactants(sim, sim_metabo, 'c_cells', 'cc_cells')
            self.get_products(sim, sim_metabo, 'c_cells', 'cc_cells')

        elif self.reaction_zone == 'mitochondria':

            self.get_reactants(sim, sim_metabo, 'c_mit', 'cc_mit')
            self.get_products(sim, sim_metabo, 'c_mit', 'cc_mit')

        # if reaction is reversible, calculate reaction quotient Q, the equilibrium constant, and backwards rate
        if self.delta_Go is not None:

            # define the reaction equilibrium coefficient:
            Keqm = np.exp(-self.delta_Go/(p.R*sim.T))

            # calculate the MM rate coefficient for the backwards reaction direction and Q:
            backwards_term = []
            Q_deno_list = []
            Q_numo_list = []

            for i, prod in enumerate(self.c_products):
                # calculate a factor in the backwards rate term:
                coeff = self.products_coeff[i]
                Km = self.Km_products_list[i]
                cs = (prod / Km) ** coeff
                term = cs / (1 + cs)
                backwards_term.append(term)

                # calculate a factor in the reaction quotient numerator term:
                ci = prod**coeff
                Q_numo_list.append(ci)

            backwards_term = np.asarray(backwards_term)

            backwards_rate = self.vmax * np.prod(backwards_term, axis=0)

            for j, react in enumerate(self.c_reactants):
                # calculate a factor in the reaction quotient denominator term:
                coeff = self.reactants_coeff[j]
                cj = react**coeff
                Q_deno_list.append(cj)

            # get the numerator and denomenator of the reaction quotient:
            Q_deno = np.prod(Q_deno_list, axis=0)
            Q_numo = np.prod(Q_numo_list, axis=0)

            # finally, *the* reaction quotient:
            Q = Q_numo/Q_deno

        else:

            Q  = 0
            backwards_rate = 0
            Keqm = 1

        # Next, calculate the MM rate coefficient for the forward reaction direction:
        forward_term = []

        for i, react in enumerate(self.c_reactants):

            coeff = self.reactants_coeff[i]
            Km = self.Km_reactants_list[i]
            cs = (react/Km)**coeff

            term = cs/(1 + cs)

            forward_term.append(term)

        forward_term = np.asarray(forward_term)

        forward_rate = self.vmax*np.prod(forward_term, axis=0)

        reaction_rate = forward_rate - (Q/Keqm)*backwards_rate

        if self.reaction_zone == 'cells':

            tag = 'mems'

        elif self.reaction_zone == 'mitochondria':

            tag = 'mitochondria'

        # get net effect of any activators or inhibitors of the reaction:
        activator_alpha, inhibitor_alpha = get_influencers(sim, sim_metabo, self.reaction_activators_list,
            self.reaction_activators_Km, self.reaction_activators_n, self.reaction_inhibitors_list,
            self.reaction_inhibitors_Km, self.reaction_inhibitors_n, reaction_zone=tag)

        # final degree of change, returned for use elsewhere (?):
        flux = activator_alpha*inhibitor_alpha*reaction_rate
        deltaC = activator_alpha*inhibitor_alpha*reaction_rate*p.dt

        if self.reaction_zone == 'cell':

            self.set_reactant_c(deltaC, sim, sim_metabo,'c_cells', 'cc_cells')
            self.set_product_c(deltaC, sim, sim_metabo, 'c_cells', 'cc_cells')

        if self.reaction_zone == 'mitochondria':

            self.set_reactant_c(deltaC, sim, sim_metabo,'c_mit', 'cc_mit')
            self.set_product_c(deltaC, sim, sim_metabo, 'c_mit', 'cc_mit')

        return flux

class Transporter(object):

    def __init__(self):

        # pre-populate the object with fields that will be assigned by MasterOfMolecules

        self.name = None
        self.reactants_list = None
        self.reactants_coeff = None
        self.products_list = None
        self.products_coeff = None
        self.Km_reactants_list = None
        self.Km_products_list = None
        self.vmax = None
        self.delta_Go = None

        self.reaction_zone = None

        self.c_reactants = []  # concentrations of reactions, defined on cell-centres
        self.z_reactants = []  # charge of reactants
        self.inds_react = []  # indices to each reactant, to keep their order
        self.c_products = []  # concentrations of reactions, defined on cell-centres
        self.z_products = []  # charge of products
        self.inds_prod = []  # indices to each reactant, to keep their order

        self.product_source_object = []  # object from which product concentrations come from
        self.product_source_type = []    # type of product concentrations sourced from object
        self.reactant_source_object = []  # object from which reactants concentrations come from
        self.reactant_source_type = []  # type of reactant concentrations sourced from object

        # activator and inhibitors of the reaction
        self.transporter_activators_list = None
        self.transporter_activators_Km = None
        self.transporter_activators_n = None
        self.transporter_inhibitors_list = None
        self.transporter_inhibitors_Km = None
        self.transporter_inhibitors_n = None

        self.transport_out_list = None
        self.transport_in_list = None

    def get_reactants(self, sim, sim_metabo, reactant_type_self, reactant_type_sim):

        """
        Get updated concentrations direct from sources for chemical reaction's reactants.

        sim:                    An instance of simulator
        sim_metabo:             Sim instance of reaction block (sim.metabo or sim.reacto for general reactions)
        reactant_type_self:     Data type to retrieve: c_cells, c_mems, c_mit
        reactant_type_sim:      Data typpe to retrieve: cc_cells, cc_mems, cc_mit

        """

        self.c_reactants = []

        for reactant_name in self.reactants_list:

            label = 'i' + reactant_name
            ion_check = getattr(sim, label, None)

            if ion_check is None:

                try:
                    obj_reactant = getattr(sim_metabo, reactant_name)
                    c_react = getattr(obj_reactant, reactant_type_self)
                    self.c_reactants.append(c_react)

                except KeyError:

                    raise BetseExceptionParameters('Name of product is not a defined chemical, or is not'
                                                   'an ion currently included in the ion profile being used.'
                                                   'Please check biomolecule definitions and ion profile'
                                                   'settings of your config(s) and try again.')

            else:

                # define the reactant as the ion concentration from the cell concentrations object in sim:
                sim_conco = getattr(sim, reactant_type_sim)

                sim_conc = sim_conco[ion_check]

                self.c_reactants.append(sim_conc)

    def get_products(self, sim, sim_metabo, product_type_self, product_type_sim):

        """
        Get updated concentrations direct from sources for chemical reaction's products.

        sim:                    An instance of simulator
        sim_metabo:             Sim instance of reaction block (sim.metabo or sim.reacto for general reactions)
        product_type_self:     Data type to retrieve: c_cells, c_mems, c_mit
        product_type_sim:      Data typpe to retrieve: cc_cells, cc_mems, cc_mit

        """

        self.c_products = []

        for product_name in self.products_list:

            label = 'i' + product_name
            ion_check = getattr(sim, label, None)

            if ion_check is None:

                try:
                    obj_prod = getattr(sim_metabo, product_name)
                    c_prod = getattr(obj_prod, product_type_self)
                    self.c_products.append(c_prod)

                except KeyError:

                    raise BetseExceptionParameters('Name of product is not a defined chemical, or is not'
                                                   'an ion currently included in the ion profile being used.'
                                                   'Please check biomolecule definitions and ion profile'
                                                   'settings of your config(s) and try again.')

            else:

                # define the reactant as the ion concentration from the cell concentrations object in sim:
                sim_conco = getattr(sim, product_type_sim)

                sim_conc = sim_conco[ion_check]

                self.c_products.append(sim_conc)

    def set_reactant_c(self, deltaMoles, sim, sim_metabo, reactant_tags, cells, p):

        for i, reactant_name in enumerate(self.reactants_list):

            type_tag = reactant_tags[i]

            if type_tag != 'c_env':

                if type_tag == 'c_mems':

                    deltaC = deltaMoles/cells.mem_vol

                    conc = self.c_reactants[i] - deltaC*self.reactants_coeff[i]

                    label = 'i' + reactant_name
                    ion_check = getattr(sim, label, None)

                    if ion_check is None:

                        obj_reactant = getattr(sim_metabo, reactant_name)
                        setattr(obj_reactant, 'c_mems', conc)

                    else:
                        # define the reactant as the ion concentration from the cell concentrations object in sim:
                        sim_conc = getattr(sim, 'cc_mems')

                        sim_conc[ion_check] = conc

                elif type_tag == 'c_mit':

                    deltaC = deltaMoles / cells.mit_vol

                    conc = self.c_reactants[i] - deltaC * self.reactants_coeff[i]

                    label = 'i' + reactant_name
                    ion_check = getattr(sim, label, None)

                    if ion_check is None:

                        obj_reactant = getattr(sim_metabo, reactant_name)
                        setattr(obj_reactant, 'c_mit', conc)

                    else:
                        # define the reactant as the ion concentration from the cell concentrations object in sim:
                        sim_conc = getattr(sim, 'cc_mit')

                        sim_conc[ion_check] = conc

            elif type_tag == 'c_env':

                if p.sim_ECM is False:

                    deltaC = deltaMoles/cells.mem_sa

                    conc = self.c_reactants[i]

                    conc = conc - deltaC * self.reactants_coeff[i]

                    mean_conc = conc.mean()
                    conc[:] = mean_conc

                    label = 'i' + reactant_name
                    ion_check = getattr(sim, label, None)

                    if ion_check is None:

                        obj_reactant = getattr(sim_metabo, reactant_name)
                        setattr(obj_reactant, 'c_env', conc)

                    else:
                        # define the reactant as the ion concentration from the cell concentrations object in sim:
                        sim_conc = getattr(sim, 'cc_env')
                        sim_conc[ion_check] = conc


                elif p.sim_ECM is True:

                    flux = deltaMoles/cells.mem_sa

                    flux_env = np.zeros(sim.edl)
                    flux_env[cells.map_mem2ecm] = flux

                    # save values at the cluster boundary:
                    bound_vals = flux_env[cells.ecm_bound_k]

                    # set the values of the global environment to zero:
                    flux_env[cells.inds_env] = 0

                    # finally, ensure that the boundary values are restored:
                    flux_env[cells.ecm_bound_k] = bound_vals

                    delta_env = (flux_env * cells.memSa_per_envSquare) / cells.ecm_vol

                    label = 'i' + reactant_name
                    ion_check = getattr(sim, label, None)

                    if ion_check is None:

                        obj_reactant = getattr(sim_metabo, reactant_name)

                        conco = obj_reactant.c_env

                        conco = conco - delta_env*self.reactants_coeff[i]

                        setattr(obj_reactant, 'c_env', conco)

                    else:
                        # define the reactant as the ion concentration from the cell concentrations object in sim:
                        sim_conc = getattr(sim, 'cc_env')

                        conco = sim_conc[ion_check]

                        conco = conco - delta_env*self.reactants_coeff[i]

                        sim_conc[ion_check] = conco[:]

    def set_product_c(self, deltaMoles, sim, sim_metabo, product_tags, cells, p):
        # FIXME set_reactant_c and set_product_c must take into account transfer between environment and cell!

        for i, product_name in enumerate(self.products_list):

            type_tag = product_tags[i]

            if type_tag != 'c_env':

                if type_tag == 'c_mems':

                    deltaC = deltaMoles / cells.mem_vol

                    conc = self.c_products[i] + deltaC * self.products_coeff[i]

                    label = 'i' + product_name
                    ion_check = getattr(sim, label, None)

                    if ion_check is None:

                        obj_product = getattr(sim_metabo, product_name)
                        setattr(obj_product, 'c_mems', conc)

                    else:
                        # define the reactant as the ion concentration from the cell concentrations object in sim:
                        sim_conc = getattr(sim, 'cc_mems')

                        sim_conc[ion_check] = conc

                elif type_tag == 'c_mit':

                    deltaC = deltaMoles / cells.mit_vol

                    conc = self.c_products[i] + deltaC * self.products_coeff[i]

                    label = 'i' + product_name
                    ion_check = getattr(sim, label, None)

                    if ion_check is None:

                        obj_product = getattr(sim_metabo, product_name)
                        setattr(obj_product, 'c_mit', conc)

                    else:
                        # define the reactant as the ion concentration from the cell concentrations object in sim:
                        sim_conc = getattr(sim, 'cc_mit')

                        sim_conc[ion_check] = conc

            elif type_tag == 'c_env':

                if p.sim_ECM is False:

                    deltaC = deltaMoles / cells.mem_sa

                    conc = self.c_products[i]

                    conc = conc + deltaC * self.products_coeff[i]

                    mean_conc = conc.mean()
                    conc[:] = mean_conc

                    label = 'i' + product_name
                    ion_check = getattr(sim, label, None)

                    if ion_check is None:

                        obj_product = getattr(sim_metabo, product_name)
                        setattr(obj_product, 'c_env', conc)

                    else:
                        # define the reactant as the ion concentration from the cell concentrations object in sim:
                        sim_conc = getattr(sim, 'cc_env')
                        sim_conc[ion_check] = conc


                elif p.sim_ECM is True:

                    flux = deltaMoles / cells.mem_sa

                    flux_env = np.zeros(sim.edl)
                    flux_env[cells.map_mem2ecm] = flux

                    # save values at the cluster boundary:
                    bound_vals = flux_env[cells.ecm_bound_k]

                    # set the values of the global environment to zero:
                    flux_env[cells.inds_env] = 0

                    # finally, ensure that the boundary values are restored:
                    flux_env[cells.ecm_bound_k] = bound_vals

                    delta_env = (flux_env * cells.memSa_per_envSquare) / cells.ecm_vol

                    label = 'i' + product_name
                    ion_check = getattr(sim, label, None)

                    if ion_check is None:

                        obj_product = getattr(sim_metabo, product_name)

                        conco = obj_product.c_env

                        conco = conco + delta_env*self.products_coeff[i]

                        setattr(obj_product, 'c_env', conco)

                    else:
                        # define the reactant as the ion concentration from the cell concentrations object in sim:
                        sim_conc = getattr(sim, 'cc_env')

                        conco = sim_conc[ion_check]

                        conco = conco + delta_env * self.products_coeff[i]

                        sim_conc[ion_check] = conco[:]

    def compute_reaction(self, sim, sim_metabo, cells, p):

        if self.reaction_zone == 'cell':
            type_self_out = 'c_env'
            type_sim_out = 'cc_env'

            type_self_in = 'c_mems'
            type_sim_in = 'cc_mems'

            type_self = 'c_mems'
            type_sim = 'cc_mems'

        elif self.reaction_zone == 'mitochondria':
            type_self_out = 'c_cells'
            type_sim_out = 'cc_cells'

            type_self_in = 'c_mit'
            type_sim_in = 'cc_mit'

            type_self = 'c_mit'
            type_sim = 'cc_mit'

        echem_terms = []
        c_reactants_trans = []   # substances transferred across membrane -- start state concs
        c_products_trans = []    # substances transferred across membrane -- end state concs
        trans_react_index = []
        trans_prod_index = []
        react_transfer_tag = []
        prod_transfer_tag = []

        if self.transport_out_list != None:

            for out_name in self.transport_out_list:

                _, c_in = get_conc(sim, sim_metabo, out_name, type_self_in, type_sim_in, cells, p)
                z_out, c_out = get_conc(sim, sim_metabo, out_name, type_self_out, type_sim_out, cells, p)

                c_reactants_trans.append(c_in)
                react_transfer_tag.append(type_self_in)

                c_products_trans.append(c_out)
                prod_transfer_tag.append(type_self_out)

                ind_r = self.reactants_list.index(out_name)
                coeff = self.reactants_coeff[ind_r]

                ind_p = self.products_list.index(out_name)

                trans_react_index.append(ind_r)
                trans_prod_index.append(ind_p)

                # get the electrochemical potential term for this reagent
                # it's negative because it starts inside the cell
                out_term = -z_out*coeff*sim.vm*p.F

                echem_terms.append(out_term)


        if self.transport_in_list != None:

            for in_name in self.transport_in_list:

                z_in, c_in = get_conc(sim, sim_metabo, in_name, type_self_in, type_sim_in, cells, p)
                _, c_out = get_conc(sim, sim_metabo, in_name, type_self_out, type_sim_out, cells, p)

                c_reactants_trans.append(c_out)
                react_transfer_tag.append(type_self_out)

                c_products_trans.append(c_in)
                prod_transfer_tag.append(type_self_in)

                ind_r = self.products_list.index(in_name)
                coeff = self.products_coeff[ind_r]

                ind_p = self.products_list.index(in_name)

                trans_react_index.append(ind_r)
                trans_prod_index.append(ind_p)

                # get the electrochemical potential term for this reagent
                # it's positive because it ends inside the cell
                in_term = z_in*coeff*sim.vm*p.F

                echem_terms.append(in_term)

        echem_terms = np.asarray(echem_terms)
        vmem_term = np.sum(echem_terms, axis = 0)

        # modification factor for the equilibrium constant due to transfer of charged item across
        #  transmembrane voltage:

        Kmod = np.exp(-vmem_term/(p.R*sim.T))

        # get up-to-date concentration data for the reaction:

        self.get_reactants(sim, sim_metabo, type_self, type_sim)
        self.get_products(sim, sim_metabo, type_self, type_sim)

        self.reactant_transfer_tag = ['c_cell' for x in range(0, len(self.c_reactants))]
        self.product_transfer_tag = ['c_cell' for x in range(0, len(self.c_products))]

        for ind_r in trans_react_index:

            self.c_reactants[ind_r] = c_reactants_trans[ind_r]
            self.reactant_transfer_tag[ind_r] = react_transfer_tag[ind_r]


        for ind_p in trans_prod_index:
            self.c_products[ind_p] = c_products_trans[ind_p]
            self.product_transfer_tag[ind_p] = prod_transfer_tag[ind_p]

        # define the reaction equilibrium coefficient:
        Keqm = np.exp(-self.delta_Go/(p.R*sim.T))*Kmod

        # calculate the MM rate coefficient for the backwards reaction direction and Q:
        backwards_term = []
        Q_deno_list = []
        Q_numo_list = []

        for i, prod in enumerate(self.c_products):
            # calculate a factor in the backwards rate term:
            coeff = self.products_coeff[i]
            Km = self.Km_products_list[i]
            cs = (prod / Km) ** coeff
            term = cs / (1 + cs)
            backwards_term.append(term)

            # calculate a factor in the reaction quotient numerator term:
            ci = prod**coeff
            Q_numo_list.append(ci)

        backwards_term = np.asarray(backwards_term)

        backwards_rate = self.vmax * np.prod(backwards_term, axis=0)

        for j, react in enumerate(self.c_reactants):
            # calculate a factor in the reaction quotient denominator term:
            coeff = self.reactants_coeff[j]
            cj = react**coeff
            Q_deno_list.append(cj)

        # get the numerator and denomenator of the reaction quotient:
        Q_deno = np.prod(Q_deno_list, axis=0)
        Q_numo = np.prod(Q_numo_list, axis=0)

        # finally, *the* reaction quotient:
        Q = Q_numo/Q_deno

        # Next, calculate the MM rate coefficient for the forward reaction direction:
        forward_term = []

        for i, react in enumerate(self.c_reactants):

            coeff = self.reactants_coeff[i]
            Km = self.Km_reactants_list[i]
            cs = (react/Km)**coeff

            term = cs/(1 + cs)

            forward_term.append(term)

        forward_term = np.asarray(forward_term)

        forward_rate = self.vmax*np.prod(forward_term, axis=0)

        reaction_rate = forward_rate - (Q/Keqm)*backwards_rate

        # get net effect of any activators or inhibitors of the reaction:

        if self.reaction_zone == 'cells':

            tag = 'mems'

        elif self.reaction_zone == 'mitochondria':

            tag = 'mitochondria'




        activator_alpha, inhibitor_alpha = get_influencers(sim, sim_metabo, self.transporter_activators_list,
            self.transporter_activators_Km, self.transporter_activators_n, self.transporter_inhibitors_list,
            self.transporter_inhibitors_Km, self.transporter_inhibitors_n, reaction_zone=tag)



        deltaC = activator_alpha*inhibitor_alpha*reaction_rate*p.dt

        # multiply by compartment volume to get a mole transfer:
        if self.reaction_zone == 'cell':
            deltaMoles = deltaC*cells.mem_vol

        elif self.reaction_zone == 'mitochondria':
            deltaMoles = deltaC*cells.mit_vol

        self.set_reactant_c(deltaMoles, sim, sim_metabo,self.reactant_transfer_tag, cells, p)
        self.set_product_c(deltaMoles, sim, sim_metabo, self.product_transfer_tag, cells, p)

        return deltaC

def get_influencers(sim, sim_metabo, a_list, Km_a_list, n_a_list, i_list, Km_i_list,
                    n_i_list, reaction_zone='cell'):
    """
    Get coefficients representing the net effect of all activators and inhibitors on a particular reaction.

    Parameters
    ------------
    sim                 Instance of BETSE simulator
    sim_metabo:         Instance of MasterOfMetabolism
    a_list:             activator names list
    Km_a_list:          activator half-max constants
    n_a_list:           activator Hill exponents
    i_list:             inhibitor names list
    Km_i_list:          inhibitor half-max constants
    n_i_list:           inhibitor Hill exponents
    reaction_zone:      Reaction occurring in 'cell' or 'mitochondria'

    Returns
    ------------
    activator_alpha         Coefficient of net effect of activators
    inhibitor_alpha         Coefficient of net effect of inhibitors
    """

    if reaction_zone == 'cell':
        type_self = 'c_cells'
        type_sim = 'cc_cells'

    if reaction_zone == 'mems':
        type_self = 'c_mems'
        type_sim = 'cc_mems'

    elif reaction_zone == 'mitochondria':

        type_self = 'c_mit'
        type_sim = 'cc_mit'

    # initialize a blank list
    activator_terms = []

    if a_list is not None and a_list != 'None' and len(a_list) > 0:  # if user specified activators for growth/decay

        # get reaction zone for data type:

        # get the activator concentration for the substance, and
        # create a term based on Hill form:
        for i, activator_name in enumerate(a_list):

            label = 'i' + activator_name
            ion_check = getattr(sim, label, None)

            if ion_check is None:

                try:
                    obj_activator = getattr(sim_metabo, activator_name)
                    c_act = getattr(obj_activator, type_self)

                except KeyError:

                    raise BetseExceptionParameters('Name of reaction activator is not a defined chemical, '
                                                   'or is not an ion currently included in the ion profile '
                                                   'being used.'
                                                   'Please check biomolecule definitions and ion profile'
                                                   'settings of your config(s) and try again.')

            else:
                # define the reactant as the ion concentration from the cell concentrations object in sim:
                sim_conco = getattr(sim, type_sim)
                c_act = sim_conco[ion_check]

            Km_act = Km_a_list[i]
            n_act = n_a_list[i]

            cs = (c_act / Km_act)**n_act

            act_term = cs / (1 + cs)

            activator_terms.append(act_term)

        activator_terms = np.asarray(activator_terms)

        # calculate the net effect of all activator terms:
        activator_alpha = np.prod(activator_terms, axis=0)


    else:

        activator_alpha = 1


    # initialize a blank list
    inhibitor_terms = []

    if i_list is not None and i_list != 'None' and len(i_list) > 0:  # if user specified inhibitors for growth/decay

        # get the inhibitor concentration for the substance, and
        # create a term based on Hill form:
        for j, inhibitor_name in enumerate(i_list):

            label = 'i' + inhibitor_name
            ion_check = getattr(sim, label, None)

            if ion_check is None:

                try:
                    obj_inhibitor = getattr(sim_metabo, inhibitor_name)
                    c_inh = getattr(obj_inhibitor, type_self)

                except KeyError:

                    raise BetseExceptionParameters('Name of substance is not a defined chemical, '
                                                   'or is not an ion currently included in the ion profile '
                                                   'being used.'
                                                   'Please check biomolecule definitions and ion profile'
                                                   'settings of your config(s) and try again.')

            else:
                # define the reactant as the ion concentration from the cell concentrations object in sim:
                sim_conco = getattr(sim, type_sim)
                c_inh = sim_conco[ion_check]

            Km_inh = Km_i_list[j]
            n_inh = n_i_list[j]

            cs = (c_inh / Km_inh) ** n_inh

            inh_term = 1 / (1 + cs)

            inhibitor_terms.append(inh_term)

        inhibitor_terms = np.asarray(inhibitor_terms)

        # calculate the net effect of all activator terms:
        inhibitor_alpha = np.prod(inhibitor_terms, axis = 0)

    else:
        inhibitor_alpha = 1


    return activator_alpha, inhibitor_alpha

def get_influencers_grn(sim, sim_metabo, a_list, k_a_list, Km_a_list, n_a_list, i_list, k_i_list, Km_i_list,
        n_i_list, reaction_zone='cell'):


    """
    Get coefficients representing the net effect of all activators and inhibitors on a particular reaction.

    Parameters
    ------------
    sim                 Instance of BETSE simulator
    sim_metabo:         Instance of MasterOfMetabolism
    a_list:             activator names list
    Km_a_list:          activator half-max constants
    n_a_list:           activator Hill exponents
    i_list:             inhibitor names list
    Km_i_list:          inhibitor half-max constants
    n_i_list:           inhibitor Hill exponents
    reaction_zone:      Reaction occurring in 'cell' or 'mitochondria'

    Returns
    ------------
    activator_alpha         Coefficient of net effect of activators
    inhibitor_alpha         Coefficient of net effect of inhibitors
    """

    if reaction_zone == 'cell':
        type_self = 'c_cells'
        type_sim = 'cc_cells'

    if reaction_zone == 'mems':
        type_self = 'c_mems'
        type_sim = 'cc_mems'

    elif reaction_zone == 'mitochondria':

        type_self = 'c_mit'
        type_sim = 'cc_mit'

    # initialize a blank list
    activator_terms = []

    if a_list is not None and a_list != 'None' and len(a_list) > 0:  # if user specified activators for growth/decay

        # get reaction zone for data type:

        # get the activator concentration for the substance, and
        # create a term based on Hill form:
        for i, activator_name in enumerate(a_list):

            label = 'i' + activator_name
            ion_check = getattr(sim, label, None)

            if ion_check is None:

                try:
                    obj_activator = getattr(sim_metabo, activator_name)
                    c_act = getattr(obj_activator, type_self)

                except KeyError:

                    raise BetseExceptionParameters('Name of reaction activator is not a defined chemical, '
                                                   'or is not an ion currently included in the ion profile '
                                                   'being used.'
                                                   'Please check biomolecule definitions and ion profile'
                                                   'settings of your config(s) and try again.')

            else:
                # define the reactant as the ion concentration from the cell concentrations object in sim:
                sim_conco = getattr(sim, type_sim)
                c_act = sim_conco[ion_check]

            k_act = k_a_list[i]
            Km_act = Km_a_list[i]
            n_act = n_a_list[i]

            cs = ((c_act*k_act)/Km_act) ** n_act

            act_term = cs / (1 + cs)

            activator_terms.append(act_term)

        activator_terms = np.asarray(activator_terms)

        # calculate the net effect of all activator terms:
        activator_alpha = np.prod(activator_terms, axis=0)


    else:

        activator_alpha = 1

    # initialize a blank list
    inhibitor_terms = []

    if i_list is not None and i_list != 'None' and len(i_list) > 0:  # if user specified inhibitors for growth/decay

        # get the inhibitor concentration for the substance, and
        # create a term based on Hill form:
        for j, inhibitor_name in enumerate(i_list):

            label = 'i' + inhibitor_name
            ion_check = getattr(sim, label, None)

            if ion_check is None:

                try:
                    obj_inhibitor = getattr(sim_metabo, inhibitor_name)
                    c_inh = getattr(obj_inhibitor, type_self)

                except KeyError:

                    raise BetseExceptionParameters('Name of substance is not a defined chemical, '
                                                   'or is not an ion currently included in the ion profile '
                                                   'being used.'
                                                   'Please check biomolecule definitions and ion profile'
                                                   'settings of your config(s) and try again.')

            else:
                # define the reactant as the ion concentration from the cell concentrations object in sim:
                sim_conco = getattr(sim, type_sim)
                c_inh = sim_conco[ion_check]

            k_inh = k_i_list[j]
            Km_inh = Km_i_list[j]
            n_inh = n_i_list[j]

            cs = ((c_inh*k_inh)/Km_inh) ** n_inh

            inh_term = 1 / (1 + cs)

            inhibitor_terms.append(inh_term)

        inhibitor_terms = np.asarray(inhibitor_terms)

        # calculate the net effect of all activator terms:
        inhibitor_alpha = np.prod(inhibitor_terms, axis=0)

    else:
        inhibitor_alpha = 1

    return activator_alpha, inhibitor_alpha

def get_conc(sim, sim_metabo, name, type_self, type_sim, cells, p):

        label = 'i' + name
        ion_check = getattr(sim, label, None)

        if ion_check is None:

            try:
                obj = getattr(sim_metabo, name)
                c = getattr(obj, type_self)
                z = obj.z

            except KeyError:

                raise BetseExceptionParameters('Name of substance is not a defined chemical, '
                                               'or is not an ion currently included in the ion profile '
                                               'being used.'
                                               'Please check biomolecule definitions and ion profile'
                                               'settings of your config(s) and try again.')

        else:
            # define the reactant as the ion concentration from the cell concentrations object in sim:
            sim_conco = getattr(sim, type_sim)
            c = sim_conco[ion_check]
            z = sim.zs[ion_check]

        if type_self == 'c_env' or type_sim == 'cc_env':

            if p.sim_ECM:
                c = c[cells.map_mem2ecm]

        return z, c








