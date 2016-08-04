#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import os
import os.path
import numpy as np
from betse.science import toolbox as tb
from betse.science import sim_toolbox as stb
from betse.science.tissue.handler import TissueHandler
from betse.science.event import modulators as mods
from betse.util.io.log import logs
import matplotlib.pyplot as plt
from betse.exceptions import BetseParametersException
from betse.science.plot import plot as viz
from betse.science.plot.anim.anim import AnimCellsTimeSeries, AnimEnvTimeSeries
from betse.science.organelles.mitochondria import Mito
from betse.util.type.mappings import DynamicValue, DynamicValueDict
from collections import OrderedDict
from matplotlib import colors
from matplotlib import cm

from betse.science.tissue.channels import vg_na as vgna
from betse.science.tissue.channels import vg_nap as vgnap
from betse.science.tissue.channels import vg_k as vgk
from betse.science.tissue.channels import vg_kir as vgkir
from betse.science.tissue.channels import vg_funny as vgfun
from betse.science.tissue.channels import vg_ca as vgca


class MasterOfNetworks(object):

    def __init__(self, sim, cells, config_substances, p, mit_enabled = False):

        """
         Initializes the MasterOfNetworks object.

        sim:                An instance of simulator
        config_settings     List of dictionaries storing key settings (p.molecules_config)
        p                   An instance of params
        """

        # Initialize dictionary mapping from molecule names to Molecule objects
        self.molecules = OrderedDict({})
        # Initialize a dict that keeps the Reaction objects:
        self.reactions = OrderedDict({})
        # Initialize a dict of Transporter objects:
        self.transporters = OrderedDict({})
        # Initialize a dict of Channels:
        self.channels = OrderedDict({})
        # Initialize a dict of modulators:
        self.modulators = OrderedDict({})

        # Initialize reaction rates array to None (filled in later, if applicable):
        self.reaction_rates = None

        # all metabolic simulations require ATP, ADP and Pi. Initialize these fields to None so that we can test
        # for their presence in metabolic sims:
        self.ATP = None
        self.ADP = None
        self.Pi = None

        # set the key controlling presence of mitochondria:
        self.mit_enabled = mit_enabled

        # read in substance properties from the config file, and initialize basic properties:
        self.read_substances(sim, cells, config_substances, p)
        self.tissue_init(sim, cells, config_substances, p)


        # write substance growth and decay equations:
        self.write_growth_and_decay()

        self.ave_cell_vol = cells.cell_vol.mean()  # average cell volume

        # colormap for plotting series of 1D lines:
        self.plot_cmap = 'viridis'

        self.globals = globals()
        self.locals = locals()

    #------------Initializers-------------------------------------------------------------------------------------------

    def read_substances(self, sim, cells, config_substances, p):
        """
            Initializes all core data structures and concentration variables
            for all molecules included in the simulation, as well as any ions present in sim.

            config_substances:  dictionary containing BETSE biomolecule template fields

        """

        logs.log_info("Reading additional substance data...")

        # Initialize a dictionaries that will eventually hold dynamic values for cell, env and mit concentrations:
        cell_concs_mapping = {}
        mem_concs_mapping = {}
        env_concs_mapping = {}
        bound_concs_mapping = {}

        # if mitochondria are enabled:
        if self.mit_enabled:
            self.mit = Mito(sim, cells, p)
            mit_concs_mapping = {}
        else:
            self.mit = None
            mit_concs_mapping = None

        # begin by adding all sim ions to the dynamic concentration mapping structures:
        for k, val in p.ions_dict.items():

            if val == 1: # if the ion is used in the simulation

                # get the numerical index from the short string label (i.e. 'Na', 'K', etc) of the ion:
                ion_index = sim.get_ion(k)

                # add the ion to the conc mappings by its short string name:
                cell_concs_mapping[k] = DynamicValue(
                    lambda ion_index=ion_index: sim.cc_cells[ion_index],
                    lambda value, ion_index=ion_index: setattr(sim.cc_cells.__setindex__(ion_index, value)))

                mem_concs_mapping[k] = DynamicValue(
                    lambda ion_index=ion_index: sim.cc_mems[ion_index],
                    lambda value, ion_index=ion_index: setattr(sim.cc_mems.__setindex__(ion_index, value)))

                env_concs_mapping[k] = DynamicValue(
                    lambda ion_index=ion_index: sim.cc_env[ion_index],
                    lambda value, ion_index=ion_index: setattr(sim.cc_env.__setindex__(ion_index, value)))

                bound_concs_mapping[k] = DynamicValue(
                    lambda ion_index=ion_index: sim.c_env_bound[ion_index],
                    lambda value, ion_index=ion_index: setattr(sim.c_env_bound.__setindex__(ion_index, value)))

                if self.mit_enabled:
                    mit_concs_mapping[k] = DynamicValue(
                    lambda ion_index=ion_index: sim.cc_mit[ion_index],
                    lambda value, ion_index=ion_index: setattr(sim.cc_mit.__setindex__(ion_index, value)))

        # Now move on to building the data structures for the user-defined substances:
        for q, mol_dic in enumerate(config_substances):
            # get each user-defined name-filed in the dictionary:
            name = str(mol_dic['name'])

            # add a field to the MasterOfMolecules corresponding to Molecule object with that name:
            self.molecules[name] = Molecule(sim, cells, p)

            # now set the attributes of that Molecule object with the cornucopia of variables:
            # define an alias so we don't have to keep typing the long version:
            mol = self.molecules[name]

            # assign general properties
            mol.name = name   # let object know who it is

            mol.c_envo = mol_dic['env conc']  # initial concentration in the environment [mmol/L]
            mol.c_cello = mol_dic['cell conc']  # initial concentration in the cytoplasm [mmol/L]

            mol.c_mito = mol_dic.get('mit conc', None)  # initialized to None if optional fields not present

            # create concentration data arrays:
            mol.c_cells = np.ones(sim.cdl) * mol.c_cello
            mol.c_mems = np.ones(sim.mdl) * mol.c_cello

            # create dynamic mappings for the cell and mem conc vectors:
            cell_concs_mapping[name] = DynamicValue(
                lambda name = name: self.molecules[name].c_cells,
                lambda value, name = name: setattr(self.molecules[name], 'c_cells', value))

            mem_concs_mapping[name] = DynamicValue(
                lambda name = name: self.molecules[name].c_mems,
                lambda value, name = name: setattr(self.molecules[name], 'c_mems', value))

            # if there is an initial concentration for mitochondria and mit are enabled, define a conc vector for it:
            if self.mit_enabled:

                if mol.c_mito is not None:

                    mol.c_mit = np.ones(sim.cdl) * mol.c_mito

                elif mol.c_mito is None: # else just set it to zero
                    mol.c_mit = np.zeros(sim.cdl)

                mit_concs_mapping[name] = DynamicValue(
                    lambda name = name: self.molecules[name].c_mit,
                    lambda value, name = name: setattr(self.molecules[name], 'c_mit', value))

            # initialize concentration in the environment:
            if p.sim_ECM is False:
                mol.c_env = np.ones(sim.mdl) * mol.c_envo
            else:
                mol.c_env = np.ones(sim.edl) * mol.c_envo

            env_concs_mapping[name] = DynamicValue(
                lambda name = name: self.molecules[name].c_env,
                lambda value, name = name: setattr(self.molecules[name], 'c_env', value))

            # initialize concentration at the boundary
            mol.c_bound = mol.c_envo

            bound_concs_mapping[name] = DynamicValue(
                lambda name = name: self.molecules[name].c_bound,
                lambda value, name = name: setattr(self.molecules[name], 'c_bound', value))

            if self.mit_enabled:
                mol.mit_enabled = True
            else:
                mol.mit_enabled = False

        self.cell_concs = DynamicValueDict(cell_concs_mapping)
        self.mem_concs = DynamicValueDict(mem_concs_mapping)
        self.env_concs = DynamicValueDict(env_concs_mapping)
        self.bound_concs = DynamicValueDict(bound_concs_mapping)

        if self.mit_enabled:
            self.mit_concs = DynamicValueDict(mit_concs_mapping)

        else:
            self.mit_concs = None

        # # transform self.molecules into an ordered dictionary so that we have guaranteed indices and order:
        # self.molecules = OrderedDict(self.molecules)

    def tissue_init(self, sim, cells, config_substances, p):
        """
        Completes the initialization process of each molecule with additional
        fields, but doesn't touch the concentrations. This is for the case where
        the user changes config file settings after running an init, and ensures
        new parameters are updated.

        """

        logs.log_info("Initializing substances/reaction network...")

        for q, mol_dic in enumerate(config_substances):
            # get each user-defined name-filed in the dictionary:
            name = str(mol_dic['name'])

            if name not in self.molecules:

                logs.log_warning("WARNING: You've added a new substance, which was not in your initialization."
                                 "It will be ignored. Please re-run an init to see all your substances in action.")


            else:
                # otherwise, get the name of the Molecule instance we're looking for as an alias:
                mol = self.molecules[name]

                mol.Dm = mol_dic['Dm']  # membrane diffusion coefficient [m2/s]
                mol.Do = mol_dic['Do']  # free diffusion constant in extra and intracellular spaces [m2/s]
                mol.z = mol_dic['z']  # charge (oxidation state)

                mol.ignore_ECM_pump = mol_dic['ignore ECM']

                mol.ignoreTJ = mol_dic['TJ permeable']  # ignore TJ?
                mol.ignoreGJ = mol_dic['GJ impermeable'] # ignore GJ?

                mol.TJ_factor = float(mol_dic['TJ factor'])

                # factors involving growth and decay (gad) in the cytoplasm
                gad = mol_dic.get('growth and decay', None)

                if gad != 'None' and gad is not None:

                    mol.simple_growth = True

                    mol.r_production = gad['production rate']
                    mol.r_decay = gad['decay rate']
                    mol.Kgd = gad['Km']
                    mol.n_production = gad['n']

                    mol.growth_profiles_list = gad['apply to']

                    modulator_function_name = gad.get('modulator function', None)

                    mol.growth_activators_list = gad.get('activators', None)
                    mol.growth_activators_zone = gad.get('zone activators', None)
                    mol.growth_activators_Km = gad.get('Km activators', None)

                    mol.growth_activators_n = gad.get('k activators', None)
                    mol.growth_inhibitors_list = gad.get('inhibitors', None)
                    mol.growth_inhibitors_zone = gad.get('zone inhibitors', None)
                    mol.growth_inhibitors_Km = gad.get('Km inhibitors', None)
                    mol.growth_inhibitors_n = gad.get('k inhibitors', None)

                    mol.init_growth(cells, p)

                    # the modulator function, if requested, makes gad occur modulated by a function.
                    # Make this happen, if it's requested:

                    if modulator_function_name != 'None' and modulator_function_name is not None:
                        mol.growth_mod_function_mems, _ = getattr(mods, modulator_function_name)(mol.growth_targets_mem,
                                                                                                  cells, p)
                        mol.growth_mod_function_cells, _ = getattr(mods, modulator_function_name)(mol.growth_targets_cell,
                                                                                                   cells, p)

                    else:
                        mol.growth_mod_function_mems = 1
                        mol.growth_mod_function_cells = 1

                else:
                    mol.simple_growth = False
                    mol.growth_profiles_list = None

                # assign ion channel gating properties
                icg = mol_dic.get('ion channel gating', None)

                if icg is not None:

                    mol.ion_channel_gating = True

                    gating_ion_o = icg['ion channel target']  # get a target ion label to gate membrane to (or 'None')

                    if gating_ion_o != 'None':
                        mol.use_gating_ligand = True

                        mol.gating_ion = []

                        for ion_o in gating_ion_o:
                            mol.gating_ion.append(sim.get_ion(ion_o))

                    else:
                        mol.use_gating_ligand = False
                        mol.gating_ion = []

                    mol.gating_Hill_K = float(icg['target Hill coefficient'])
                    mol.gating_Hill_n = float(icg['target Hill exponent'])
                    mol.gating_max_val = float(icg['peak channel opening'])
                    mol.gating_extracell = icg['acts extracellularly']

                    # get any optional activators and inhibitors for the channel:
                    mol.activators_list = icg.get('activators', None)
                    mol.activators_Km = icg.get('Km activators', None)
                    mol.activators_n = icg.get('n activators', None)

                    mol.inhibitors_list = icg.get('inhibitors', None)
                    mol.inhibitors_Km = icg.get('Km inhibitors', None)
                    mol.inhibitors_n = icg.get('n inhibitors', None)

                else:

                    mol.ion_channel_gating = False

                # assign active pumping properties
                ap = mol_dic.get('active pumping', None)

                if ap is not None:

                    mol.active_pumping = True
                    mol.use_pumping = ap['turn on']
                    mol.pump_to_cell = ap['pump to cell']
                    mol.pump_max_val = ap['maximum rate']
                    mol.pump_Km = ap['pump Km']
                    mol.pumps_use_ATP = ap['uses ATP']


                else:
                    mol.active_pumping = False

                # assign boundary change event properties
                cab = mol_dic.get('change at bounds', None)

                if cab is not None:
                    mol.change_bounds = True
                    mol.change_at_bounds = cab['event happens']
                    mol.change_bounds_start = cab['change start']
                    mol.change_bounds_end = cab['change finish']
                    mol.change_bounds_rate = cab['change rate']
                    mol.change_bounds_target = cab['concentration']

                else:
                    mol.change_bounds = False

                # assign plotting properties
                pd = mol_dic['plotting']
                mol.make_plots = pd['plot 2D']
                mol.make_ani = pd['animate']

                mol.plot_autoscale = pd['autoscale colorbar']
                mol.plot_max = pd['max val']
                mol.plot_min = pd['min val']

    def plot_init(self, config_substances):

        for q, mol_dic in enumerate(config_substances):
            # get each user-defined name-filed in the dictionary:
            name = str(mol_dic['name'])

            # assign alias for convenience
            mol = self.molecules[name]

            # assign plotting properties
            pd = mol_dic['plotting']
            mol.make_plots = pd['plot 2D']
            mol.make_ani = pd['animate']

            mol.plot_autoscale = pd['autoscale colorbar']
            mol.plot_max = pd['max val']
            mol.plot_min = pd['min val']

    def init_saving(self, cells, p, plot_type='init', nested_folder_name='Molecules'):

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
            raise BetseParametersException(
                'The "plot cell" defined in the "results" section of your '
                'configuration file does not exist in your cluster. '
                'Choose a plot cell number smaller than the maximum cell number.')

    def read_reactions(self, config_reactions, sim, cells, p):

        """
          Read in and initialize parameters for all user-defined reactions.

          config_options:  dictionary containing BETSE reaction template fields

          """

        logs.log_info("Reading reaction input data...")

        for q, react_dic in enumerate(config_reactions):

            # get each user-defined name-filed in the dictionary:
            name = str(react_dic['name'])

            # add a Reaction object to MasterOfReactions reaction dictionary:
            self.reactions[name] = Reaction(sim, cells, p)

            # now set the attributes of that Reaction object with the cornucopia of variables:
            # assign an alias for convenience:
            obj = self.reactions[name]

            # assign general properties
            obj.name = name  # let object know who it is

            # list where the reaction takes place; if field not specified default to 'cell':
            obj.reaction_zone = react_dic['reaction zone']

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
            obj.reaction_activators_zone = react_dic.get('activator zone', None)

            obj.reaction_inhibitors_list = react_dic.get('reaction inhibitors', None)
            obj.reaction_inhibitors_Km = react_dic.get('inhibitor Km', None)
            obj.reaction_inhibitors_n = react_dic.get('inhibitor n', None)
            obj.reaction_inhibitors_zone = react_dic.get('inhibitor zone', None)


            if self.mit_enabled:
                obj.mit_enabled = True
            else:
                obj.mit_enabled = False

            # # Finally, make the self.reactions into an ordered dictionary, so we can guarantee indices:
            # self.reactions = OrderedDict(self.reactions)

    def read_transporters(self, config_transporters, sim, cells, p):

        """
            Read in and initialize parameters for all user-defined transporters.

            config_options:  dictionary containing BETSE transporter template fields

            """

        logs.log_info("Reading transporter input data...")

        for q, trans_dic in enumerate(config_transporters):
            # get each user-defined name-filed in the dictionary:
            name = str(trans_dic['name'])

            # add a field to the MasterOfReactions corresponding to Transporter object with that name:

            self.transporters[name] = Transporter(sim, cells, p)

            # now set the attributes of that Transporter object with the cornucopia of variables:

            # get MasterOfMolecules.name
            obj = self.transporters[name]

            # assign general properties
            obj.name = name  # let object know who it is

            # list where the reaction takes place; if field not specified default to 'cell':
            obj.reaction_zone = trans_dic['reaction zone']

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
            obj.ignore_ECM_transporter = trans_dic['ignore ECM']

            obj.transporter_profiles_list = trans_dic['apply to']
            obj.init_reaction(cells, p)

            if obj.delta_Go == 'None':
                obj.delta_Go = None  # make the field a proper None variable

            else:
                obj.delta_Go = float(obj.delta_Go)

            obj.transporter_activators_list = trans_dic.get('transporter activators', None)
            obj.transporter_activators_Km = trans_dic.get('activator Km', None)
            obj.transporter_activators_n = trans_dic.get('activator n', None)
            obj.transporter_activators_zone = trans_dic.get('activator zone', None)

            obj.transporter_inhibitors_list = trans_dic.get('transporter inhibitors', None)
            obj.transporter_inhibitors_Km = trans_dic.get('inhibitor Km', None)
            obj.transporter_inhibitors_n = trans_dic.get('inhibitor n', None)
            obj.transporter_inhibitors_zone = trans_dic.get('inhibitor zone', None)

            if self.mit_enabled:
                obj.mit_enabled = True
            else:
                obj.mit_enabled = False

    def read_channels(self, config_channels, sim, cells, p):

        logs.log_info("Reading channel input data...")

        for q, chan_dic in enumerate(config_channels):
            # get each user-defined name-filed in the dictionary:
            name = str(chan_dic['name'])

            # add a field to the MasterOfReactions corresponding to Channel object with that name:
            self.channels[name] = Channel(sim, cells, p)

            # now set the attributes of that channel:

            # get MasterOfMolecules.name
            obj = self.channels[name]

            # assign general properties
            obj.name = name  # let object know who it is

            # list where the reaction takes place; for channels, defaults to 'cell':
            obj.reaction_zone = 'cell'

            obj.channel_class = chan_dic['channel class']
            obj.channel_type = chan_dic['channel type']
            obj.channelMax = chan_dic['max conductivity']
            obj.channel_profiles_list = chan_dic['apply to']

            obj.channel_activators_list = chan_dic.get('channel activators', None)
            obj.channel_activators_Km = chan_dic.get('activator Km', None)
            obj.channel_activators_n = chan_dic.get('activator n', None)
            obj.channel_activators_zone = chan_dic.get('activator zone', None)
            obj.channel_inhibitors_list = chan_dic.get('channel inhibitors', None)
            obj.channel_inhibitors_Km = chan_dic.get('inhibitor Km', None)
            obj.channel_inhibitors_n = chan_dic.get('inhibitor n', None)
            obj.channel_inhibitors_zone = chan_dic.get('inhibitor zone', None)

            obj.init_channel(obj.channel_class, obj.channel_type, obj.channelMax, sim, cells, p)

    def read_modulators(self, config_modulators, sim, cells, p):

        logs.log_info("Reading modulator input data...")

        for q, mod_dic in enumerate(config_modulators):

            name = str(mod_dic['name'])

            # add a field to the MasterOfReactions corresponding to Channel object with that name:
            self.modulators[name] = Modulator()

            # now set the attributes of that channel:

            # get MasterOfMolecules.name
            obj = self.modulators[name]

            obj.target_label = str(mod_dic['target'])
            obj.zone = str(mod_dic['zone'])
            obj.max_val = float(mod_dic['max effect'])
            obj.modulator_activators_list = mod_dic.get('activators', None)
            obj.modulator_activators_Km = mod_dic.get('activator Km', None)
            obj.modulator_activators_n = mod_dic.get('activator n', None)
            obj.modulator_activators_zone = mod_dic.get('activator zone', None)
            obj.modulator_inhibitors_list = mod_dic.get('inhibitors', None)
            obj.modulator_inhibitors_Km = mod_dic.get('inhibitor Km', None)
            obj.modulator_inhibitors_n = mod_dic.get('inhibitor n', None)
            obj.modulator_inhibitors_zone = mod_dic.get('inhibitor zone', None)

            obj.init_modulator(sim, cells, p)

    def write_growth_and_decay(self):

        logs.log_info("Writing substance growth/decay equations...")

        for mol_name in self.molecules:

            if self.molecules[mol_name].simple_growth is True:

                cc = "(self.molecules['{}'].c_cells / self.molecules['{}'].Kgd)".format(mol_name, mol_name)

                # get activators and inhibitors and associated data:

                a_list = self.molecules[mol_name].growth_activators_list
                Km_a_list = self.molecules[mol_name].growth_activators_Km
                n_a_list = self.molecules[mol_name].growth_activators_n
                i_list = self.molecules[mol_name].growth_inhibitors_list
                Km_i_list = self.molecules[mol_name].growth_inhibitors_Km
                n_i_list = self.molecules[mol_name].growth_inhibitors_n

                activator_alpha, inhibitor_alpha = self.get_influencers(a_list, Km_a_list, n_a_list, i_list, Km_i_list,
                    n_i_list, reaction_zone='cell')

                # definine remaining portion of substance's growth and decay expression:

                mod_funk = "(self.molecules['{}'].growth_mod_function_cells)".format(mol_name)
                r_prod =  "(self.molecules['{}'].r_production)".format(mol_name)
                r_decay = "(self.molecules['{}'].r_decay)".format(mol_name)

                gad_eval_string = mod_funk + "*" + r_prod + "*" + inhibitor_alpha + "*" + activator_alpha + "-" + \
                                  r_decay + "*" + cc

                self.molecules[mol_name].gad_eval_string = gad_eval_string

            else:

                self.molecules[mol_name].gad_eval_string = "np.zeros(sim.cdl)"

    def write_reactions(self):
        """
        Reactions are now constructed during the init as strings that are evaluated in eval calls in each time-step.
        This function constructs the evaluation strings for each reaction, given the metadata stored
        in each reaction object (e.g. lists of reactants, products, etc).

        """

        logs.log_info("Writing reaction equations...")

        for reaction_name in self.reactions:

        # define aliases for convenience:

            reactant_names = self.reactions[reaction_name].reactants_list
            reactant_coeff = self.reactions[reaction_name].reactants_coeff
            reactant_Km = self.reactions[reaction_name].Km_reactants_list

            product_names = self.reactions[reaction_name].products_list
            product_coeff = self.reactions[reaction_name].products_coeff
            product_Km = self.reactions[reaction_name].Km_products_list

            # activator/inhibitor lists and associated data:
            a_list = self.reactions[reaction_name].reaction_activators_list
            Km_a_list = self.reactions[reaction_name].reaction_activators_Km
            n_a_list = self.reactions[reaction_name].reaction_activators_n
            i_list = self.reactions[reaction_name].reaction_inhibitors_list
            Km_i_list = self.reactions[reaction_name].reaction_inhibitors_Km
            n_i_list = self.reactions[reaction_name].reaction_inhibitors_n

            # first calculate a reaction coefficient Q, as a string expression
            numo_string_Q = "("
            denomo_string_Q = "("

            for i, (name, coeff) in enumerate(zip(reactant_names, reactant_coeff)):

                denomo_string_Q += "(self.cell_concs['{}']".format(name)
                denomo_string_Q += "**{})".format(coeff)

                if i < len(reactant_names) - 1:

                    denomo_string_Q += "*"

                else:

                    denomo_string_Q += ")"

            for i, (name, coeff) in enumerate(zip(product_names, product_coeff)):

                numo_string_Q += "(self.cell_concs['{}']".format(name)
                numo_string_Q += "**{})".format(coeff)

                if i < len(product_names) - 1:

                    numo_string_Q += "*"

                else:

                    numo_string_Q += ")"

            # define the final reaction quotient string, Q:
            Q = "(" + numo_string_Q + '/' + denomo_string_Q + ")"

            # next calculate the forward and backward reaction rate coefficients:---------------------------------------

            forward_coeff = "("
            backward_coeff = "("

            for i, (name, n, Km) in enumerate(zip(reactant_names, reactant_coeff, reactant_Km)):

                numo_string_r = "((self.cell_concs['{}']/{})**{})".format(name, Km, n)
                denomo_string_r = "(1 + (self.cell_concs['{}']/{})**{})".format(name, Km, n)

                term = "(" + numo_string_r + "/" + denomo_string_r + ")"

                forward_coeff += term

                if i < len(reactant_names) - 1:

                    forward_coeff += "*"

                else:

                    forward_coeff += ")"

            for i, (name, n, Km) in enumerate(zip(product_names, product_coeff, product_Km)):

                numo_string_p = "((self.cell_concs['{}']/{})**{})".format(name, Km, n)
                denomo_string_p = "(1 + (self.cell_concs['{}']/{})**{})".format(name, Km, n)

                term = "(" + numo_string_p + "/" + denomo_string_p + ")"

                backward_coeff += term

                if i < len(product_names) - 1:

                    backward_coeff += "*"

                else:

                    backward_coeff += ")"

            # if reaction is reversible deal calculate an equilibrium constant:
            if self.reactions[reaction_name].delta_Go is not None:

                # define the reaction equilibrium coefficient expression:
                Keqm = "(np.exp(-self.reactions['{}'].delta_Go / (p.R * sim.T)))".format(reaction_name)

            else:

                Q = "0"
                backward_coeff = "0"
                Keqm = "1"

            # Put it all together into a final reaction string (+max rate):

            reversed_term = "(" + Q + "/" + Keqm + ")"

            activator_alpha, inhibitor_alpha = self.get_influencers(a_list, Km_a_list, n_a_list, i_list, Km_i_list,
                n_i_list, reaction_zone='cell')

            vmax = "self.reactions['{}'].vmax".format(reaction_name)

            reaction_eval_string = vmax + "*" + activator_alpha + "*" + inhibitor_alpha + "*" + "(" + \
                                   forward_coeff + "-" + "(" + reversed_term + "*" + backward_coeff + ")" + ")"

            # reaction_eval_string = vmax + "*" + activator_alpha + "*" + inhibitor_alpha + "*" + forward_coeff + \
            #                        "*" + "(" + "1 - " + reversed_term  + ")"

            # add the composite string describing the reaction math to a new field:
            self.reactions[reaction_name].reaction_eval_string = reaction_eval_string

            # call statement to evaluate:
            # eval(self.reactions['consume_ATP'].reaction_eval_string, globals(), locals())

    def write_transporters(self, cells, p):
        """
        Reactions are now constructed during the init as strings that are evaluated in eval calls in each time-step.
        This function constructs the evaluation strings for each reaction, given the metadata stored
        in each reaction object (e.g. lists of reactants, products, etc).

        """

        logs.log_info("Writing transporter equations...")

        for transp_name in self.transporters:

        # define aliases for convenience:

            reactant_names = self.transporters[transp_name].reactants_list
            reactant_coeff = self.transporters[transp_name].reactants_coeff
            reactant_Km = self.transporters[transp_name].Km_reactants_list

            product_names = self.transporters[transp_name].products_list
            product_coeff = self.transporters[transp_name].products_coeff
            product_Km = self.transporters[transp_name].Km_products_list

            # activator/inhibitor lists and associated data:
            a_list = self.transporters[transp_name].transporter_activators_list
            Km_a_list = self.transporters[transp_name].transporter_activators_Km
            n_a_list = self.transporters[transp_name].transporter_activators_n
            i_list = self.transporters[transp_name].transporter_inhibitors_list
            Km_i_list = self.transporters[transp_name].transporter_inhibitors_Km
            n_i_list = self.transporters[transp_name].transporter_inhibitors_n

            reaction_zone = self.transporters[transp_name].reaction_zone

            transport_out_list = self.transporters[transp_name].transport_out_list
            transport_in_list = self.transporters[transp_name].transport_in_list

        # initialize lists to hold the reactants and product transfer tags (initialized to default zone):
            react_transfer_tag = ['mem_concs' for x in reactant_names]
            prod_transfer_tag = ['mem_concs' for x in product_names]

        # list holding terms affecting free energy via transfer of charged item across membrane:
            echem_terms_list = []

            if reaction_zone == 'cell':

                type_out = 'env_concs'

                vmem = "sim.vm"   # get the transmembrane voltage for this category


                in_delta_term_react = "-self.transporters['{}'].flux*(cells.mem_sa/cells.mem_vol)".format(transp_name)
                in_delta_term_prod = "self.transporters['{}'].flux*(cells.mem_sa/cells.mem_vol)".format(transp_name)

                if p.sim_ECM is True:

                    if self.transporters[transp_name].ignore_ECM_transporter is True:

                        out_delta_term_react = "-self.transporters['{}'].flux*(cells.memSa_per_envSquare" \
                                               "[cells.map_mem2ecm]/cells.ecm_vol)".format(transp_name)

                        out_delta_term_prod = "self.transporters['{}'].flux*(cells.memSa_per_envSquare" \
                                              "[cells.map_mem2ecm]/cells.ecm_vol)".format(transp_name)

                    else:
                        out_delta_term_react = "-self.transporters['{}'].flux*(cells.memSa_per_envSquare" \
                                               "[cells.map_mem2ecm]/cells.true_ecm_vol[cells.map_mem2ecm])"\
                                                .format(transp_name)

                        out_delta_term_prod = "self.transporters['{}'].flux*(cells.memSa_per_envSquare" \
                                              "[cells.map_mem2ecm]/cells.true_ecm_vol[cells.map_mem2ecm])"\
                                                .format(transp_name)

                else:
                    out_delta_term_react = "-self.transporters['{}'].flux*(cells.mem_sa/cells.mem_vol)".format(transp_name)

                    out_delta_term_prod = "self.transporters['{}'].flux*(cells.mem_sa/cells.mem_vol)".format(transp_name)

            elif reaction_zone == 'mit' and self.mit_enabled is True:

                # initialize lists to hold the reactants and product transfer tags (initialized to zone):
                react_transfer_tag = ['mit_concs' for x in reactant_names]
                prod_transfer_tag = ['mit_concs' for x in product_names]

                type_out = 'cell_concs'

                vmem = "sim.metabo.mit.Vmit"  # get the transmembrane voltage for this category

                in_delta_term_react = "-self.transporters['{}'].flux*(sim.metabo.mit_sa/sim.metabo.mit_vol)".\
                                                                        format(transp_name)

                in_delta_term_prod = "self.transporters['{}'].flux*(sim.metabo.mit_sa/sim.metabo.mit_vol)".\
                                                                format(transp_name)

                out_delta_term_react = "-self.transporters['{}'].flux*(cells.cell_sa/cells.cell_vol)". \
                                                                        format(transp_name)

                out_delta_term_prod = "self.transporters['{}'].flux*(cells.cell_sa/cells.cell_vol)". \
                                                                     format(transp_name)

        # initialize list that will hold expression for calculating net concentration change
            delta_strings_reactants = [in_delta_term_react for x in reactant_names]
            delta_strings_products = [in_delta_term_prod for x in product_names]

            # calculate a reactants zone tag list and terms affecting transporter free energy via Vmem:
            if transport_out_list != 'None':

                for out_name in transport_out_list:

                    # get the index for the substance in the molecules database
                    prod_i = product_names.index(out_name)

                    prod_transfer_tag[prod_i] = type_out

                    coeff = product_coeff[prod_i]

                    # calculate the effect of transfer on transporter's free energy:
                    eterm = "-{}*self.molecules['{}'].z*p.F*{}".format(coeff, out_name, vmem)

                    echem_terms_list.append(eterm)

                    # update delta string for correct transfer:
                    delta_strings_products[prod_i] = out_delta_term_prod


            if transport_in_list != 'None':

                for in_name in transport_in_list:

                    # get the index for the substance in the molecules database
                    react_i = reactant_names.index(in_name)
                    react_transfer_tag[react_i] = type_out

                    coeff = reactant_coeff[react_i]

                    # calculate the effect of transfer on transporter's free energy:
                    eterm = "{}*self.molecules['{}'].z*p.F*{}".format(coeff, in_name, vmem)

                    echem_terms_list.append(eterm)

                    # update delta string for correct transfer:
                    delta_strings_reactants[react_i] = out_delta_term_react


            # create the eterms string expression describing net effect of trans-membrane fluxes on free energy:
            echem_string = "("
            for i, et in enumerate(echem_terms_list):

                echem_string += et

                if i < len(echem_terms_list) -1:

                    echem_string += "+"

                else:
                    echem_string += ")"

            # first calculate a reaction coefficient Q, as a string expression
            numo_string_Q = "("
            denomo_string_Q = "("

            for i, (name, coeff, tag) in enumerate(zip(reactant_names, reactant_coeff, react_transfer_tag)):

                if tag != 'env_concs' or p.sim_ECM is False:

                    denomo_string_Q += "(self.{}['{}']".format(tag, name)

                # get the concentration from the environment, mapped to respective membranes:
                else:

                    denomo_string_Q += "(self.{}['{}'][cells.map_mem2ecm]".format(tag, name)

                denomo_string_Q += "**{})".format(coeff)

                if i < len(reactant_names) - 1:

                    denomo_string_Q += "*"

                else:

                    denomo_string_Q += ")"

            for i, (name, coeff, tag) in enumerate(zip(product_names, product_coeff, prod_transfer_tag)):

                if tag != 'env_concs' or p.sim_ECM is False:

                    numo_string_Q += "(self.{}['{}']".format(tag, name)

                # get the concentration from the environment mapped to the respective membranes:
                else:
                    numo_string_Q += "(self.{}['{}'][cells.map_mem2ecm]".format(tag, name)

                numo_string_Q += "**{})".format(coeff)

                if i < len(product_names) - 1:

                    numo_string_Q += "*"

                else:

                    numo_string_Q += ")"

            # define the final reaction quotient string, Q:
            Q = "(" + numo_string_Q + '/' + denomo_string_Q + ")"

            # next calculate the forward and backward reaction rate coefficients:---------------------------------------

            forward_coeff = "("
            backward_coeff = "("

            for i, (name, n, Km, tag) in enumerate(zip(reactant_names, reactant_coeff, reactant_Km, react_transfer_tag)):

                if tag != 'env_concs' or p.sim_ECM is False:

                    numo_string_r = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    denomo_string_r = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)

                else:

                    numo_string_r = "((self.{}['{}'][cells.map_mem2ecm]/{})**{})".format(tag, name, Km, n)
                    denomo_string_r = "(1 + (self.{}['{}'][cells.map_mem2ecm])/{})**{})".format(tag, name, Km, n)

                term = "(" + numo_string_r + "/" + denomo_string_r + ")"

                forward_coeff += term

                if i < len(reactant_names) - 1:

                    forward_coeff += "*"

                else:

                    forward_coeff += ")"

            for i, (name, n, Km, tag) in enumerate(zip(product_names, product_coeff, product_Km, prod_transfer_tag)):

                if tag != 'env_concs' or p.sim_ECM is False:

                    numo_string_p = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    denomo_string_p = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)

                else:

                    numo_string_p = "((self.{}['{}'][cells.map_mem2ecm]/{})**{})".format(tag, name, Km, n)
                    denomo_string_p = "(1 + (self.{}['{}'][cells.map_mem2ecm]/{})**{})".format(tag, name, Km, n)


                term = "(" + numo_string_p + "/" + denomo_string_p + ")"

                backward_coeff += term

                if i < len(product_names) - 1:

                    backward_coeff += "*"

                else:

                    backward_coeff += ")"

            # if reaction is reversible deal calculate an equilibrium constant:
            if self.transporters[transp_name].delta_Go is not None:

                # define the reaction equilibrium coefficient expression:
                Keqm = "(np.exp(-(self.transporters['{}'].delta_Go + {})/ (p.R * sim.T)))".format(transp_name,
                                                                                                    echem_string)
            else:

                Q = "0"
                backward_coeff = "0"
                Keqm = "1"

            # Put it all together into a final reaction string (+max rate):

            reversed_term = "(" + Q + "/" + Keqm + ")"

            activator_alpha, inhibitor_alpha = self.get_influencers(a_list, Km_a_list, n_a_list, i_list, Km_i_list,
                n_i_list, reaction_zone='mem')

            vmax = "self.transporters['{}'].vmax".format(transp_name)

            transporter_eval_string = vmax + "*" + activator_alpha + "*" + inhibitor_alpha + "*" + "(" + \
                                   forward_coeff + "-" + "(" + reversed_term + "*" + backward_coeff + ")" + ")"

            # add the composite string describing the reaction math to a new field:
            # Evaluate this expression to assign self.transporters[transp_name].flux field
            self.transporters[transp_name].transporter_eval_string = transporter_eval_string

            # evaluate these expressions to get the correct change in concentration to each reactant/product
            # after assigning .flux as described above
            self.transporters[transp_name].delta_react_eval_strings = delta_strings_reactants
            self.transporters[transp_name].delta_prod_eval_strings = delta_strings_products

            # finally, expressions to do the concentration change after assigning transporters.delta_react
            # and transporters.delta_prod:
            self.transporters[transp_name].react_transport_tag = react_transfer_tag
            self.transporters[transp_name].prod_transport_tag = prod_transfer_tag

    def create_reaction_matrix(self):

        """
        Creates a matrix representation of the system of coupled
        differential equations underlying the reaction network.

        """

        logs.log_info("Writing reaction network matrix...")
        n_reacts = len(self.molecules) + len(self.reactions)
        n_mols = len(self.molecules)

        # initialize the network's reaction matrix:
        self.reaction_matrix = np.zeros((n_mols, n_reacts))

        molecule_keys = list(self.molecules.keys())

        for i, name in enumerate(self.molecules):
            # add in terms referencing the self-growth and decay reaction for each substance
            self.reaction_matrix[i, i] = 1

        for jo, reaction_name in enumerate(self.reactions):

            j = jo + n_mols

            for react_name, coeff in zip(self.reactions[reaction_name].reactants_list,
                self.reactions[reaction_name].reactants_coeff):
                i = molecule_keys.index(react_name)

                self.reaction_matrix[i, j] = -coeff

            for prod_name, coeff in zip(self.reactions[reaction_name].products_list,
                self.reactions[reaction_name].products_coeff):
                i = molecule_keys.index(prod_name)

                self.reaction_matrix[i, j] = coeff

    #------runners------------------------------------------------------------------------------------------------------

    def run_loop(self, t, sim, cells, p):
        """
        Runs the main simulation loop steps for each of the molecules included in the simulation.

        """

        # calculate rates of growth/decay:
        gad_rates = np.asarray(
            [eval(self.molecules[mol].gad_eval_string, self.globals, self.locals) for mol in self.molecules])

        # ... and rates of chemical reactions:
        self.reaction_rates = np.asarray(
            [eval(self.reactions[rn].reaction_eval_string, self.globals, self.locals) for rn in self.reactions])

        # stack into an integrated data structure:
        all_rates = np.vstack((gad_rates, self.reaction_rates))

        # calculate concentration rate of change using linear algebra:
        self.delta_conc = np.dot(self.reaction_matrix, all_rates)

        # Initialize arrays for substance charge contribution:
        net_Q_cell = 0
        net_Q_env = 0

        # get the name of the specific substance:
        for name, deltac in zip(self.molecules, self.delta_conc):

            obj = self.molecules[name]

            # update concentration due to growth/decay and chemical reactions:
            obj.c_cells = obj.c_cells + deltac*p.dt

            # if pumping is enabled:
            if obj.active_pumping:
                obj.pump(sim, cells, p)

            if p.run_sim is True:
                # use the substance as a gating ligand (if desired)
                if obj.ion_channel_gating:
                    obj.gating(sim, self, cells, p)

                # update the global boundary (if desired)
                if obj.change_bounds:
                    obj.update_boundary(t, p)

            # transport the molecule through gap junctions and environment:
            obj.transport(sim, cells, p)

            # update the substance on the inside of the cell:
            obj.updateIntra(sim, self, cells, p)

            # calculate energy charge in the cell:
            self.energy_charge(sim)

            # ensure no negs:
            stb.no_negs(obj.c_mems)

            # calculate the charge density this substance contributes to cell and environment:
            obj_Q_cell = p.F * obj.c_mems * obj.z

            obj_Q_env = p.F * obj.c_env * obj.z
            # add that contribution to the total sum:
            net_Q_cell = net_Q_cell + obj_Q_cell
            net_Q_env = net_Q_env + obj_Q_env

        if p.substances_affect_charge:
            # update charge in the cell and environment of the main bioelectric simulator:
            sim.rho_cells = sim.rho_cells + net_Q_cell
            sim.rho_env = sim.rho_env + net_Q_env

        if self.mit_enabled:  # if enabled, update the mitochondria's voltage and other properties

            self.mit.update(sim, cells, p)

        # manage pH in cells, environment and mitochondria:
        self.pH_handling(sim, cells, p)

    def run_loop_transporters(self, t, sim, sim_metabo, cells, p):

        # call statement to evaluate:
        for name in self.transporters:

            # specific tissue profile regions where the transporter is active:
            targ_mem = self.transporters[name].transporter_targets_mem
            targ_cell = self.transporters[name].transporter_targets_cell
            targ_env = self.transporters[name].transporter_targets_env

            # calculate the flux
            self.transporters[name].flux = eval(self.transporters[name].transporter_eval_string,
                self.globals, self.locals)

            # finally, update the concentrations using the final eval statements:
            for i, (delc, coeff) in enumerate(zip(self.transporters[name].delta_react_eval_strings,
                self.transporters[name].reactants_coeff)):

                # obtain the change for the reactant
                delta_react = coeff*eval(delc,self.globals, self.locals)

                # finally, update the concentrations using the final eval statements:
                if self.transporters[name].react_transport_tag[i] == 'mem_concs':

                    self.mem_concs[self.transporters[name].reactants_list[i]][targ_mem] = \
                        self.mem_concs[self.transporters[name].reactants_list[i]][targ_mem] + \
                        delta_react[targ_mem]*p.dt

                elif self.transporters[name].react_transport_tag[i] == 'env_concs':

                    if p.sim_ECM is True:

                        delta_react_expanded = np.zeros(sim.edl)
                        delta_react_expanded[cells.map_mem2ecm] = delta_react[:]

                        self.env_concs[self.transporters[name].reactants_list[i]][targ_env] = \
                            self.env_concs[self.transporters[name].reactants_list[i]][targ_env] + \
                            delta_react_expanded[targ_env]*p.dt

                    else:

                        self.env_concs[self.transporters[name].reactants_list[i]][targ_mem] = \
                            self.env_concs[self.transporters[name].reactants_list[i]][targ_mem] + \
                            delta_react[targ_mem] * p.dt

                elif self.transporters[name].react_transport_tag[i] == 'cell_concs':

                    self.cell_concs[self.transporters[name].reactants_list[i]][targ_cell] = \
                        self.cell_concs[self.transporters[name].reactants_list[i]][targ_cell] + \
                        delta_react[targ_cell]*p.dt

                elif self.transporters[name].react_transport_tag[i] == 'mit_concs':

                    self.mit_concs[self.transporters[name].reactants_list[i]][targ_cell] = \
                        self.mit_concs[self.transporters[name].reactants_list[i]][targ_cell] + \
                        delta_react[targ_cell] * p.dt

                else:

                    raise BetseParametersException("Internal error: transporter zone not specified correctly!")

            for i, (delc, coeff) in enumerate(zip(self.transporters[name].delta_prod_eval_strings,
                self.transporters[name].products_coeff)):

                # obtain the change for the product
                delta_prod = coeff*eval(delc, self.globals, self.locals)

                # finally, update the concentrations using the final eval statements:
                if self.transporters[name].prod_transport_tag[i] == 'mem_concs':

                    self.mem_concs[self.transporters[name].products_list[i]][targ_mem] = \
                        self.mem_concs[self.transporters[name].products_list[i]][targ_mem] + \
                        delta_prod[targ_mem]*p.dt

                elif self.transporters[name].prod_transport_tag[i] == 'env_concs':

                    if p.sim_ECM is True:

                        delta_prod_expanded = np.zeros(sim.edl)
                        delta_prod_expanded[cells.map_mem2ecm] = delta_prod[:]

                        self.env_concs[self.transporters[name].products_list[i]][targ_env] = \
                            self.env_concs[self.transporters[name].products_list[i]][targ_env] + \
                            delta_prod_expanded[targ_env]*p.dt

                    else:

                        self.env_concs[self.transporters[name].products_list[i]][targ_mem] = \
                            self.env_concs[self.transporters[name].products_list[i]][targ_mem] + \
                            delta_prod[targ_mem] * p.dt

                elif self.transporters[name].prod_transport_tag[i] == 'cell_concs':

                    self.cell_concs[self.transporters[name].products_list[i]][targ_cell] = \
                        self.cell_concs[self.transporters[name].products_list[i]][targ_cell] + \
                        delta_prod[targ_cell]*p.dt

                elif self.transporters[name].prod_transport_tag[i] == 'mit_concs':

                    self.mit_concs[self.transporters[name].products_list[i]][targ_cell] = \
                        self.mit_concs[self.transporters[name].products_list[i]][targ_cell] + \
                        delta_prod[targ_cell] * p.dt

                else:

                    raise BetseParametersException("Internal error: transporter zone not specified correctly!")

    def run_loop_channels(self, sim, sim_metabo, cells, p):

        # get the object corresponding to the specific transporter:
        for i, name in enumerate(self.channels):
            # get the Reaction object:
            obj = self.channels[name]

            # compute the channel activity
            obj.run_channel(sim, sim_metabo, cells, p)

    def run_loop_modulators(self, sim, sim_metabo, cells, p):

        # get the object corresponding to the specific transporter:
        for i, name in enumerate(self.modulators):
            # get the Reaction object:
            obj = self.modulators[name]

            # compute the channel activity
            obj.run_modulator(sim, sim_metabo, cells, p)

    # ------Utility Methods---------------------------------------------------------------------------------------------
    def pH_handling(self, sim, cells, p):
        """
        Molecules may contain dissolved carbon dioxide as a substance,
        and reactions/transporters may act on H+ levels via bicarbonate
        (M-). Therefore, update pH in cells, environment, and
        if enabled, mitochondria.

        """

        if 'CO2' in self.molecules:

            if p.ions_dict['H'] == 1:

                # if the simulation contains sim.cHM_mems, use it and update it!
                sim.cHM_mems = self.CO2.c_cells
                sim.cHM_env = self.CO2.c_env

                # update the cH and pH fields of sim with potentially new value of sim.iM
                sim.cc_cells[sim.iH], sim.pH_cell = stb.bicarbonate_buffer(self.CO2.c_cells,
                    sim.cc_cells[sim.iM])
                sim.cc_env[sim.iH], sim.pH_env = stb.bicarbonate_buffer(self.CO2.c_env, sim.cc_env[sim.iM])

                if self.mit_enabled:
                    # update the cH and pH fields of sim with potentially new value of sim.iM
                    sim.cc_mit[sim.iH], sim.pH_mit = stb.bicarbonate_buffer(self.CO2.c_mit, sim.cc_mit[sim.iM])

            elif p.ions_dict['H'] != 1:

                # update the cH and pH fields of sim with potentially new value of M ion:
                _, sim.pH_cell = stb.bicarbonate_buffer(self.CO2.c_cells, sim.cc_cells[sim.iM])
                _, sim.pH_env = stb.bicarbonate_buffer(self.CO2.c_env, sim.cc_env[sim.iM])

                if self.mit_enabled:
                    # update the cH and pH fields of sim with potentially new value of sim.iM
                    _, sim.pH_mit = stb.bicarbonate_buffer(self.CO2.c_mit, sim.cc_mit[sim.iM])

        else:  # if we're not using CO2 in the simulator, use the default p.CO2*0.03

            CO2 = p.CO2 * 0.03  # get the default concentration of CO2

            if p.ions_dict['H'] == 1:

                # update the cH and pH fields of sim with potentially new value of sim.iM
                sim.cc_cells[sim.iH], sim.pH_cell = stb.bicarbonate_buffer(CO2, sim.cc_cells[sim.iM])
                sim.cc_env[sim.iH], sim.pH_env = stb.bicarbonate_buffer(CO2, sim.cc_env[sim.iM])

                if self.mit_enabled:
                    # update the cH and pH fields of sim with potentially new value of sim.iM
                    sim.cc_mit[sim.iH], sim.pH_mit = stb.bicarbonate_buffer(CO2, sim.cc_mit[sim.iM])

            elif p.ions_dict['H'] != 1:

                # update the cH and pH fields of sim with potentially new value of sim.iM
                _, sim.pH_cell = stb.bicarbonate_buffer(CO2, sim.cc_cells[sim.iM])
                _, sim.pH_env = stb.bicarbonate_buffer(CO2, sim.cc_env[sim.iM])

                if self.mit_enabled:
                    # update the cH and pH fields of sim with potentially new value of sim.iM
                    _, sim.pH_mit = stb.bicarbonate_buffer(CO2, sim.cc_mit[sim.iM])

    def energy_charge(self, sim):

        if 'AMP' in self.molecules:

            numo = (self.ATP.c_cells + 0.5 * self.ADP.c_cells)
            denomo = (self.ATP.c_cells + self.ADP.c_cells + self.AMP.c_cells)

            self.chi = numo / denomo

        else:

            self.chi = np.zeros(sim.cdl)

    def updateInside(self, sim, cells, p):
        """
        Runs the main simulation loop steps for each of the molecules included in the simulation.

        """

        # get the name of the specific substance:
        for name in self.molecules:
            obj = self.molecules[name]
            obj.updateIntra(sim, self, cells, p)

    def mod_after_cut_event(self, target_inds_cell, target_inds_mem, sim, cells, p, met_tag=False):

        # get the name of the specific substance:
        for name in self.molecules:

            obj = self.molecules[name]

            obj.remove_cells(target_inds_cell, target_inds_mem, sim, cells, p)

        if self.mit_enabled:
            self.mit.remove_mits(sim, target_inds_cell)

        if sim.met_concs is not None and met_tag is True:  # update metabolism object if it's being simulated
            sim.met_concs = {'cATP': self.ATP.c_mems,
                'cADP': self.ADP.c_mems,
                'cPi': self.Pi.c_mems}

        for name in self.transporters:
            obj = self.transporters[name]

            obj.update_transporter(sim, cells, p)

        for name in self.channels:
            obj = self.channels[name]

            obj.update_channel(sim, cells, p)

    def clear_cache(self):
        """
        Initializes or clears the time-storage vectors at the beginning of init and sim runs.

        """

        # get the name of the specific substance:
        for name in self.molecules:
            obj = self.molecules[name]

            obj.c_mems_time = []
            obj.c_cells_time = []
            obj.c_env_time = []

            if self.mit_enabled:
                obj.c_mit_time = []

        for name in self.transporters:
            obj = self.transporters[name]

            obj.flux_time = []

        for name in self.reactions:
            obj = self.reactions[name]

            obj.rate_time = []

        if self.mit_enabled:
            self.vmit_time = []
            self.pH_mit_time = []

        self.pH_cells_time = []
        self.pH_env_time = []

        self.chi_time = []

    def write_data(self, sim, p):
        """
        Writes concentration data from a time-step to time-storage vectors.

        """

        # get the name of the specific substance:
        for name in self.molecules:

            obj = self.molecules[name]

            obj.c_mems_time.append(obj.c_mems)
            obj.c_cells_time.append(obj.c_cells)
            obj.c_env_time.append(obj.c_env)

            if self.mit_enabled:
                obj.c_mit_time.append(obj.c_mit)

        # save rates of transporters and reactions
        for name in self.transporters:
            obj = self.transporters[name]

            if obj.flux is not None:
                obj.flux_time.append(obj.flux)

        for i, name in enumerate(self.reactions):
            obj = self.reactions[name]

            if self.reaction_rates is not None:
                obj.rate_time.append(self.reaction_rates[i])

        if self.mit_enabled:
            self.vmit_time.append(self.mit.Vmit[:])
            self.pH_mit_time.append(sim.pH_mit)

        self.pH_cells_time.append(sim.pH_cell)
        self.pH_env_time.append(sim.pH_env)

        self.chi_time.append(self.chi)

    def report(self, sim, p):
        """
        At the end of the simulation, tell user about mean, final concentrations of each molecule.

        """

        for name in self.molecules:

            obj = self.molecules[name]

            if self.mit_enabled:
                logs.log_info('Average concentration of ' + str(name) + ' in the mitochondria: ' +
                              str(np.round(obj.c_mit.mean(), 4)) + ' mmol/L')

            logs.log_info('Average concentration of ' + str(name) + ' in the cell: ' +
                          str(np.round(obj.c_cells.mean(), 4)) + ' mmol/L')

            # logs.log_info('Average concentration of ' + str(name) + ' in the environment: ' +
            #                               str(np.round(obj.c_env.mean(), 4)) + ' mmol/L')

        if self.mit_enabled:
            logs.log_info('Average Vmit: ' + str(np.round(1.0e3 * self.mit.Vmit.mean(), 4)) + ' mV')

            if 'ETC' in self.transporters:

                rate = 0.5 * 3600 * 1e15 * self.mit.mit_vol.mean() * self.ETC.rate.mean()

                if 'ETC_ROS' in self.transporters:
                    rate = rate + 3600 * 1e15 * self.mit.mit_vol.mean() * self.ETC_ROS.rate.mean()

                logs.log_info('Average O2 consumption rate: ' + str(rate) + ' fmol/cell/hr')

        # logs.log_info('Average pH in cell: ' + str(np.round(sim.pH_cell.mean(), 4)))
        # logs.log_info('Average pH in env: ' + str(np.round(sim.pH_env.mean(), 4)))

        # if self.mit_enabled:
        #     logs.log_info('Average pH in mitochondria: ' + str(np.round(sim.pH_mit.mean(), 4)))

        if self.chi.mean() != 0.0:
            logs.log_info('Energy charge of cell ' + str(np.round(self.chi.mean(), 3)))

    def export_all_data(self, sim, cells, p, message='for auxiliary molecules...'):

        """

        Exports concentration data from each molecule to a file for a single cell
        (plot cell defined in params) as a function of time.

        """
        logs.log_info('Exporting raw data for ' + message)
        # get the name of the specific substance:
        for name in self.molecules:
            obj = self.molecules[name]

            obj.export_data(sim, cells, p, self.resultsPath)

        if self.mit_enabled:
            # FIXME we should also save vmit to a file? and pH and vm?
            pass

    def plot(self, sim, cells, p, message='for auxiliary molecules...'):
        """
        Creates plots for each molecule included in the simulation.

        """

        logs.log_info('Plotting 1D and 2D data for ' + message)

        # get the name of the specific substance:
        for name in self.molecules:

            obj = self.molecules[name]

            if p.plot_single_cell_graphs:
                # create line graphs for the substance
                obj.plot_1D(sim, p, self.imagePath)

            # create 2D maps for the substance
            obj.plot_cells(sim, cells, p, self.imagePath)

            # if there's a real environment, plot 2D concentration in the environment
            if p.sim_ECM:
                obj.plot_env(sim, cells, p, self.imagePath)

        # ---------------cell everything plot---------------------------------------------
        data_all1D = []
        fig_all1D = plt.figure()
        ax_all1D = plt.subplot(111)

        # set up the color vector (for plotting complex line graphs)
        maxlen = len(self.molecules)

        lineplots_cm = plt.get_cmap(self.plot_cmap)
        cNorm = colors.Normalize(vmin=0, vmax=maxlen)
        c_names = cm.ScalarMappable(norm=cNorm, cmap=lineplots_cm)

        for i, name in enumerate(self.molecules):
            obj = self.molecules[name]

            c_cells = [arr[p.plot_cell] for arr in obj.c_cells_time]

            ax_all1D.plot(sim.time, c_cells, color=c_names.to_rgba(i), linewidth=2.0, label=name)

        legend = ax_all1D.legend(loc='upper right', shadow=False, frameon=False)

        ax_all1D.set_xlabel('Time [s]')
        ax_all1D.set_ylabel('Concentration [mmol/L]')
        ax_all1D.set_title('Concentration of all substances in cell ' + str(p.plot_cell))

        if p.autosave is True:
            savename = self.imagePath + 'AllCellConcentrations_' + str(p.plot_cell) + '.png'
            plt.savefig(savename, format='png', transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # -------------environment everything plot-------------------------------------------------
        data_all1D = []
        fig_all1D = plt.figure()
        ax_all1D = plt.subplot(111)

        # get a random selection of our chosen colors in the length of our data set:

        for i, name in enumerate(self.molecules):

            obj = self.molecules[name]

            if p.sim_ECM is True:
                c_env = [arr[cells.map_cell2ecm][p.plot_cell] for arr in obj.c_env_time]

            else:

                mem_i = cells.cell_to_mems[p.plot_cell][0]

                c_env = [arr[mem_i] for arr in obj.c_env_time]

            ax_all1D.plot(sim.time, c_env, color=c_names.to_rgba(i), linewidth=2.0, label=name)

        legend = ax_all1D.legend(loc='upper right', shadow=False, frameon=False)

        ax_all1D.set_xlabel('Time [s]')
        ax_all1D.set_ylabel('Concentration [mmol/L]')
        ax_all1D.set_title('Concentration of all substances in environment of cell ' + str(p.plot_cell))

        if p.autosave is True:
            savename = self.imagePath + 'AllEnvConcentrations_' + str(p.plot_cell) + '.png'
            plt.savefig(savename, format='png', transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # ------------------------------------------------------------------------------------------
        if self.mit_enabled:

            # 1 D plot of mitochondrial voltage--------------------------------------------------------
            vmit = [1e3 * arr[p.plot_cell] for arr in self.vmit_time]

            figVmit = plt.figure()
            axVmit = plt.subplot(111)

            axVmit.plot(sim.time, vmit)

            axVmit.set_xlabel('Time [s]')
            axVmit.set_ylabel('Vmit [mV]')
            axVmit.set_title('Mitochondrial transmembrane voltage in cell: ' + str(p.plot_cell))

            if p.autosave is True:
                savename = self.imagePath + 'Vmit_cell_' + str(p.plot_cell) + '.png'
                plt.savefig(savename, format='png', transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

            # 2D plot of mitochondrial voltage ---------------------------------------------------

            fig, ax, cb = viz.plotPolyData(sim, cells, p,
                zdata=self.mit.Vmit * 1e3, number_cells=p.enumerate_cells, clrmap=p.default_cm)

            ax.set_title('Final Mitochondrial Transmembrane Voltage')
            ax.set_xlabel('Spatial distance [um]')
            ax.set_ylabel('Spatial distance [um]')
            cb.set_label('Vmit [mV]')

            if p.autosave is True:
                savename = self.imagePath + '2DVmit.png'
                plt.savefig(savename, format='png', transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

            # plot of all substances in the mitochondria:----------------------------------------------

            data_all1D = []
            fig_all1D = plt.figure()
            ax_all1D = plt.subplot(111)

            # get a random selection of our chosen colors in the length of our data set:

            for i, name in enumerate(self.molecules):
                obj = self.molecules[name]

                c_mit = [arr[p.plot_cell] for arr in obj.c_mit_time]

                ax_all1D.plot(sim.time, c_mit, color=c_names.to_rgba(i), linewidth=2.0, label=name)

            legend = ax_all1D.legend(loc='upper right', shadow=False, frameon=False)

            ax_all1D.set_xlabel('Time [s]')
            ax_all1D.set_ylabel('Concentration [mmol/L]')
            ax_all1D.set_title('Substances in mitochondria of cell ' + str(p.plot_cell))

            if p.autosave is True:
                savename = self.imagePath + 'AllMitConcentrations_' + str(p.plot_cell) + '.png'
                plt.savefig(savename, format='png', transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)
        # ------pH plot------------------------------------------------------------------------------

        # 1 D plot of pH in cell, env and mit ------------------------------------------------------
        pHcell = [arr[p.plot_cell] for arr in self.pH_cells_time]

        if p.sim_ECM:
            pHenv = [arr[cells.map_cell2ecm][p.plot_cell] for arr in self.pH_env_time]

        else:
            avPh = [np.dot(cells.M_sum_mems, arr) / cells.num_mems for arr in self.pH_env_time]
            pHenv = [arr[p.plot_cell] for arr in avPh]

        if self.mit_enabled:
            pHmit = [arr[p.plot_cell] for arr in self.pH_mit_time]

        else:
            pHmit = np.zeros(len(sim.time))

        figpH = plt.figure()
        axpH = plt.subplot(111)

        axpH.plot(sim.time, pHcell, label='cell')
        axpH.plot(sim.time, pHmit, label='mitochondria')
        axpH.plot(sim.time, pHenv, label='env')

        axpH.set_xlabel('Time [s]')
        axpH.set_ylabel('pH')
        axpH.set_title('pH in/near cell : ' + str(p.plot_cell))

        if p.autosave is True:
            savename = self.imagePath + 'pH_' + str(p.plot_cell) + '.png'
            plt.savefig(savename, format='png', transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # -------Reaction rate plot and data export----------------------------------------

        if len(self.reactions):

            # create a suite of single reaction line plots:
            for i, name in enumerate(self.reactions):
                # get the reaction object field
                obj = self.reactions[name]

                # make a 1D plot of this reaction rate:
                obj.plot_1D(sim, cells, p, self.imagePath)

            # now create the "everything" plot and export data to file:

            react_dataM = []
            react_header = 'Time [s], '

            react_dataM.append(sim.time)

            data_all1D = []
            fig_all1D = plt.figure()
            ax_all1D = plt.subplot(111)

            # set up the color vector (for plotting complex line graphs)
            maxlen = len(self.reactions)

            lineplots_cm = plt.get_cmap(self.plot_cmap)
            cNorm = colors.Normalize(vmin=0, vmax=maxlen)
            c_names = cm.ScalarMappable(norm=cNorm, cmap=lineplots_cm)

            for i, name in enumerate(self.reactions):
                # get the reaction object field
                obj = self.reactions[name]

                if len(obj.rate_time) > 0:
                    r_rate = [arr[p.plot_cell] for arr in obj.rate_time]

                    ax_all1D.plot(sim.time, r_rate, color=c_names.to_rgba(i), linewidth=2.0, label=name)

                    react_dataM.append(r_rate)
                    react_header = react_header + name + ' [mM/s]' + ','

            legend = ax_all1D.legend(loc='upper right', shadow=False, frameon=False)

            ax_all1D.set_xlabel('Time [s]')
            ax_all1D.set_ylabel('Rate [mM/s]')
            ax_all1D.set_title('Reaction rates in cell ' + str(p.plot_cell))

            if p.autosave is True:
                savename = self.imagePath + 'AllReactionRates_' + str(p.plot_cell) + '.png'
                plt.savefig(savename, format='png', transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

            react_dataM = np.asarray(react_dataM)

            saveName = 'AllReactionRatesData_' + str(p.plot_cell) + '.csv'

            saveDataReact = os.path.join(self.resultsPath, saveName)

            np.savetxt(saveDataReact, react_dataM.T, delimiter=',', header=react_header)

        # ---Transporter rate plot and data export ------------------------------------------------------

        if len(self.transporters):

            transp_dataM = []
            transp_header = 'Time [s], '

            transp_dataM.append(sim.time)

            # set up the color vector (for plotting complex line graphs)
            maxlen = len(self.transporters)

            lineplots_cm = plt.get_cmap(self.plot_cmap)
            cNorm = colors.Normalize(vmin=0, vmax=maxlen)
            c_names = cm.ScalarMappable(norm=cNorm, cmap=lineplots_cm)

            for i, name in enumerate(self.transporters):
                obj = self.transporters[name]

                # make a 1D plot of this reaction rate:
                obj.plot_1D(sim, cells, p, self.imagePath)

                if len(obj.rate_time) > 0:

                    # check the data structure size for this transporter:
                    if len(obj.rate_time[0]) == sim.cdl:

                        t_rate = [arr[p.plot_cell] for arr in obj.rate_time]

                    elif len(obj.rate_time[0]) == sim.mdl:
                        mem_i = cells.cell_to_mems[p.plot_cell][0]
                        t_rate = [arr[mem_i] for arr in obj.rate_time]

                    else:
                        t_rate = np.zeros(len(sim.time))

                    data_all1D = []
                    fig_all1D = plt.figure()
                    ax_all1D = plt.subplot(111)

                    ax_all1D.plot(sim.time, t_rate, color=c_names.to_rgba(i), linewidth=2.0, label=name)

                    transp_dataM.append(t_rate)
                    transp_header = transp_header + name + ' [mM/s]' + ','

            legend = ax_all1D.legend(loc='upper right', shadow=False, frameon=False)

            ax_all1D.set_xlabel('Time [s]')
            ax_all1D.set_ylabel('Rate [mM/s]')
            ax_all1D.set_title('Transporter rates in cell ' + str(p.plot_cell))

            if p.autosave is True:
                savename = self.imagePath + 'AllTransporterRates_' + str(p.plot_cell) + '.png'
                plt.savefig(savename, format='png', transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

            saveName = 'AllTransporterRatesData_' + str(p.plot_cell) + '.csv'

            saveDataTransp = os.path.join(self.resultsPath, saveName)

            transp_dataM = np.asarray(transp_dataM)

            np.savetxt(saveDataTransp, transp_dataM.T, delimiter=',', header=transp_header)

        # energy charge plots:----------------------------------------------------------
        # 1 D plot of mitochondrial voltage--------------------------------------------------------
        chio = [arr[p.plot_cell] for arr in self.chi_time]

        figChi = plt.figure()
        axChi = plt.subplot(111)

        axChi.plot(sim.time, chio)

        axChi.set_xlabel('Time [s]')
        axChi.set_ylabel('Energy charge')
        axChi.set_title('Energy charge in cell: ' + str(p.plot_cell))

        if p.autosave is True:
            savename = self.imagePath + 'EnergyCharge_cell_' + str(p.plot_cell) + '.png'
            plt.savefig(savename, format='png', transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # ---2D plot--------------------------------------------------------------------



        fig, ax, cb = viz.plotPolyData(sim, cells, p,
            zdata=self.chi, number_cells=p.enumerate_cells, clrmap=p.default_cm)

        ax.set_title('Final Energy Charge of Cell')
        ax.set_xlabel('Spatial distance [um]')
        ax.set_ylabel('Spatial distance [um]')
        cb.set_label('Energy Charge')

        if p.autosave is True:
            savename = self.imagePath + '2DEnergyCharge.png'
            plt.savefig(savename, format='png', transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    def anim(self, sim, cells, p, message='for auxiliary molecules...'):
        """
        Animates 2D data for each molecule in the simulation.
        """

        logs.log_info('Animating data for %s', message)

        # get the name of the specific substance:
        for name in self.molecules:
            obj = self.molecules[name]

            if p.anim.is_after_sim and obj.make_ani is True:
                # create 2D animations for the substance in cells
                obj.anim_cells(sim, cells, p)

                # create 2D animations for the substance in the environment
                if p.sim_ECM:
                    obj.anim_env(sim, cells, p)

    def get_influencers(self, a_list, Km_a_list, n_a_list, i_list, Km_i_list,
        n_i_list, reaction_zone='cell'):

        """
        Given lists of activator and inhibitor names, with associated
        Km and n data, this function constructs string expressions
        representing the net effect of the activators and inhibitors.
        This output can be used as a multiplier for another rate equation.

        a_list:          List of substance names defined in MasterOfNetworks.molecules which serve as activators
        Km_a_list:       List of Hill-function Km's to the activators
        n_a_list:        List of Hill-function exponents (n) to the activators
        i_list:          List of substance names defined in MasterOfNetworks.molecules which serve as inhibitors
        Km_i_list:       List of Hill-function Km's to the inhibitors
        n_i_list:        List of Hill-function exponents (n) to the inhibitors
        reaction_zone:   Where the reaction takes place ('cell', 'mit', 'mem', or 'env')

        Returns
        --------
        activator_alpha, inhibitor_alpha     string expressions representing the activator and inhibitor coeffs

        """

        # initialize string expressions for product of all activators and inhibitors
        activator_alpha = "("
        inhibitor_alpha = "("

        if a_list is not None and a_list != 'None' and len(a_list) > 0:

            # Begin with construction of the activator effects term:
            for i, (name, Km, n) in enumerate(zip(a_list, Km_a_list, n_a_list)):

                # initialize string expressions for activator and inhibitor, respectively
                numo_string_a = ""
                denomo_string_a = ""

                if reaction_zone == 'cell':

                    numo_string_a += "((self.cell_concs['{}']/{})**{})".format(name, Km, n)
                    denomo_string_a += "(1 + (self.cell_concs['{}']/{})**{})".format(name, Km, n)

                elif reaction_zone == 'mem':

                    numo_string_a += "((self.mem_concs['{}']/{})**{})".format(name, Km, n)
                    denomo_string_a += "(1 + (self.mem_concs['{}']/{})**{})".format(name, Km, n)

                elif reaction_zone == 'env':

                    numo_string_a += "((self.env_concs['{}']/{})**{})".format(name, Km, n)
                    denomo_string_a += "(1 + (self.env_concs['{}']/{})**{})".format(name, Km, n)

                elif reaction_zone == 'mit':

                    numo_string_a += "((self.mit_concs['{}']/{})**{})".format(name, Km, n)
                    denomo_string_a += "(1 + (self.mit_concs['{}']/{})**{})".format(name, Km, n)

                term = "(" + numo_string_a + "/" + denomo_string_a + ")"

                activator_alpha += term

                if i < len(a_list) - 1:

                    activator_alpha += "*"

                else:

                    activator_alpha += ")"
        else:

            activator_alpha += "1)"

        if i_list is not None and i_list != 'None' and len(i_list) > 0:

            # Next, construct the inhibitors net effect term:
            for i, (name, Km, n) in enumerate(zip(i_list, Km_i_list, n_i_list)):

                # initialize string expressions for activator and inhibitor, respectively
                numo_string_i = ""
                denomo_string_i = ""

                if reaction_zone == 'cell':

                    numo_string_i += "1"
                    denomo_string_i += "(1 + (self.cell_concs['{}']/{})**{})".format(name, Km, n)

                elif reaction_zone == 'mem':

                    numo_string_i += "1"
                    denomo_string_i += "(1 + (self.mem_concs['{}']/{})**{})".format(name, Km, n)

                elif reaction_zone == 'env':

                    numo_string_i += "1"
                    denomo_string_i += "(1 + (self.env_concs['{}']/{})**{})".format(name, Km, n)

                elif reaction_zone == 'mit':

                    numo_string_i += "1"
                    denomo_string_i += "(1 + (self.mit_concs['{}']/{})**{})".format(name, Km, n)

                term = "(" + numo_string_i + "/" + denomo_string_i + ")"

                inhibitor_alpha += term

                if i < len(i_list) - 1:

                    inhibitor_alpha += "*"

                else:

                    inhibitor_alpha += ")"

        else:

            inhibitor_alpha += "1)"

        return activator_alpha, inhibitor_alpha

class Molecule(object):

    def __init__(self, sim, cells, p):

        self.dummy_dyna = TissueHandler(sim, cells, p)
        self.dummy_dyna.tissueProfiles(sim, cells, p)  # initialize all tissue profiles

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



        self.c_mems, self.c_env, _, _, _, _ = stb.molecule_mover(sim,
                                                                self.c_mems,
                                                                self.c_env,
                                                                cells, p,
                                                                z=self.z,
                                                                Dm = self.Dm,
                                                                Do = self.Do,
                                                                c_bound = self.c_bound,
                                                                ignoreECM = self.ignore_ECM_pump,
                                                                smoothECM = p.smooth_concs,
                                                                ignoreTJ = self.ignoreTJ,
                                                                ignoreGJ = self.ignoreGJ)

    def updateC(self, flux, sim, cells, p):
        """

        General updater for a flux defined on membranes and updating concentrations in
        cells and environment.

        """
        self.c_mems, self.c_env = stb.update_Co(sim, self.c_mems, self.c_cells, flux, cells, p, ignoreECM=True)

    def updateIntra(self, sim, sim_metabo, cells, p):

        self.c_mems, self.c_cells, _ = stb.update_intra(sim, cells, self.c_mems, self.c_cells, self.Do, self.z, p)

        if self.mit_enabled:

            IdCM = np.ones(sim.cdl)

            f_ED = stb.electroflux(self.c_cells, self.c_mit, self.Dm*IdCM, p.tm*IdCM, self.z*IdCM,
                sim_metabo.mit.Vmit, sim.T, p, rho=1)

            # update with flux
            self.c_cells = self.c_cells - f_ED*(sim_metabo.mit.mit_sa/cells.cell_vol)*p.dt
            self.c_mit = self.c_mit + f_ED*(sim_metabo.mit.mit_sa/sim_metabo.mit.mit_vol)*p.dt

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
                                                                     Km_ATP=1.0, met = met_vect, ignoreECM = self.ignore_ECM_pump)
                if p.metabolism_enabled:
                    # update ATP concentrations after pump action:
                    sim.metabo.update_ATP(flux, sim, cells, p)

            else:

                self.c_mems, self.c_env, flux = stb.molecule_transporter(sim, self.c_mems, self.c_env,
                    cells, p, Df=self.Do, z=self.z, pump_into_cell=self.pump_to_cell, alpha_max=self.pump_max_val,
                    Km_X=self.pump_Km, Keq= 1.0, ignoreECM = self.ignore_ECM_pump)

    def gating(self, sim, sim_metabo, cells, p):
        """
        Uses the molecule concentration to open an ion channel in the cell membranes.

        """

        # update membrane permeability if dye targets an ion channel:
        if self.use_gating_ligand:

            # calculate any activators and/or inhibitor effects:
            activator_alpha, inhibitor_alpha = get_influencers(sim, sim_metabo, self.activators_list,
                self.activators_Km, self.activators_n, self.inhibitors_list,
                self.inhibitors_Km, self.inhibitors_n, reaction_zone='mems')


            if self.gating_extracell is False:

                for ion_tag in self.gating_ion:

                    Dm_mod_mol = sim.rho_channel*self.gating_max_val*tb.hill(self.c_mems,
                                                                            self.gating_Hill_K,self.gating_Hill_n)

                    sim.Dm_morpho[ion_tag] = sim.rho_channel*Dm_mod_mol*activator_alpha*inhibitor_alpha

            elif self.gating_extracell is True and p.sim_ECM is True:

                for ion_tag in self.gating_ion:

                    Dm_mod_mol = self.gating_max_val*tb.hill(self.c_env,self.gating_Hill_K,self.gating_Hill_n)

                    sim.Dm_morpho[ion_tag] = (activator_alpha*inhibitor_alpha*sim.rho_channel*
                                              Dm_mod_mol[cells.map_mem2ecm])

    def init_growth(self,cells, p):

        if self.growth_profiles_list is not None and self.growth_profiles_list != 'all':

            self.growth_targets_cell = []
            self.growth_targets_mem = []

            for profile in self.growth_profiles_list:
                targets_cell = self.dummy_dyna.cell_target_inds[profile]
                self.growth_targets_cell.extend(targets_cell)

                targets_mem = self.dummy_dyna.tissue_target_inds[profile]
                self.growth_targets_mem.extend(targets_mem)


        elif self.growth_profiles_list is None or self.growth_profiles_list == 'all':

            self.growth_targets_cell = cells.cell_i
            self.growth_targets_mem = cells.mem_i

    def remove_cells(self, target_inds_cell, target_inds_mem, sim, cells, p):
        """
        During a cutting event, removes the right cells from the simulation network,
        while preserving additional information.

        """

        # remove cells from the cell concentration list:
        ccells2 = np.delete(self.c_cells, target_inds_cell)
        # reassign the new data vector to the object:
        self.c_cells = ccells2[:]

        # remove cells from the mems concentration list:
        cmems2 = np.delete(self.c_mems, target_inds_mem)
        # reassign the new data vector to the object:
        self.c_mems = cmems2[:]

        if self.simple_growth is True and self.growth_mod_function_cells != 1:

            gmfc = np.delete(self.growth_mod_function_cells, target_inds_cell)
            self.growth_mod_function_cells = gmfc[:]

            if len(self.growth_mod_function_cells) == 0:
                self.growth_mod_function_cells = 1

        self.dummy_dyna.tissueProfiles(sim, cells, p)  # re-initialize all tissue profiles
        self.init_growth(cells, p)

        if p.sim_ECM is False:

            cenv2 = np.delete(self.c_env, target_inds_mem)
            self.c_env = cenv2[:]


        if self.mit_enabled:
            # remove cells from the cell concentration list:
            cmit2 = np.delete(self.c_mit, target_inds_cell)
            # reassign the new data vector to the object:
            self.c_mit = cmit2[:]

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

        saveName = 'ExportData_' + self.name + '_' + str(p.plot_cell) + '.csv'

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

        if self.mit_enabled:

            cmit = [arr[ci] for arr in self.c_mit_time]

            headr = headr + 'Mit_Conc_' + self.name + '_mmol/L' + ','

            cmit = np.asarray(cmit)

        time = np.asarray(sim.time)
        ccell = np.asarray(ccell)
        cenv = np.asarray(cenv)

        if self.mit_enabled is False:
            dataM = np.column_stack((time, ccell, cenv))

        else:
            dataM = np.column_stack((time, ccell, cenv, cmit))

        np.savetxt(saveData, dataM, delimiter=',', header=headr)

    def plot_1D(self, sim, p, saveImagePath):
        """
        Create 1D plot of concentration in cell and environment for a single cell (params plot cell)
        as a function of simulation time.

        """

        c_cells = [arr[p.plot_cell] for arr in self.c_cells_time]
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(sim.time, c_cells)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Concentration [mmol/L]')
        ax.set_title('Concentration of ' + self.name + ' in cell ' + str(p.plot_cell))

        if p.autosave is True:
            savename = saveImagePath + 'CellConcentration_' + self.name + '_' + str(p.plot_cell) + '.png'
            plt.savefig(savename, format='png', transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        if self.mit_enabled:

            c_mit = [arr[p.plot_cell] for arr in self.c_mit_time]
            fig = plt.figure()
            ax = plt.subplot(111)
            ax.plot(sim.time, c_mit)
            ax.set_xlabel('Time [s]')
            ax.set_ylabel('Concentration [mmol/L]')
            ax.set_title('Mitochondrial concentration of ' + self.name + ' in cell ' + str(p.plot_cell))

            if p.autosave is True:
                savename = saveImagePath + 'MitConcentration_' + self.name + '_' + str(p.plot_cell) + '.png'
                plt.savefig(savename, format='png', transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

    def plot_cells(self, sim, cells, p, saveImagePath):
        """
        Create 2D plot of cell concentrations.

        """

        fig, ax, cb = viz.plotPrettyPolyData(self.c_mems,
            sim, cells, p,
            number_cells=p.enumerate_cells,
            clrAutoscale=self.plot_autoscale,
            clrMin=self.plot_min,
            clrMax=self.plot_max,
            clrmap=p.default_cm)

        ax.set_title('Final ' + self.name + ' Concentration in Cells')
        ax.set_xlabel('Spatial distance [um]')
        ax.set_ylabel('Spatial distance [um]')
        cb.set_label('Concentration mmol/L')

        if p.autosave is True:
            savename = saveImagePath + '2Dcell_conc_' + self.name + '.png'
            plt.savefig(savename,format='png', transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # mitochondrial plots
        if self.mit_enabled:

            fig, ax, cb = viz.plotPolyData(sim, cells, p, zdata=self.c_mit,
                number_cells=p.enumerate_cells,
                clrAutoscale=self.plot_autoscale,
                clrMin=self.plot_min,
                clrMax=self.plot_max,
                clrmap=p.default_cm)

            ax.set_title('Final ' + self.name + ' Concentration in Mitochondria')
            ax.set_xlabel('Spatial distance [um]')
            ax.set_ylabel('Spatial distance [um]')
            cb.set_label('Concentration mmol/L')

            if p.autosave is True:
                savename = saveImagePath + '2D_mit_conc_' + self.name + '.png'
                plt.savefig(savename, format='png', transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

    def plot_env(self, sim, cells, p, saveImagePath):
        """
        Create 2D plot of environmental concentration.

        """


        fig = plt.figure()
        ax = plt.subplot(111)

        dyeEnv = (self.c_env).reshape(cells.X.shape)

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
        cb.set_label('Concentration mmol/L')

        if p.autosave is True:
            savename = saveImagePath + '2Denv_conc_' + self.name + '.png'
            plt.savefig(savename,format='png',dpi = 300.0, transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    def anim_cells(self, sim, cells, p):
        """
        Create 2D animation of cell concentration.
        """

        AnimCellsTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=[arr for arr in self.c_mems_time],
            label=self.name + '_cells',
            figure_title='Cytosolic ' + self.name,
            colorbar_title='Concentration [mmol/L]',
            is_color_autoscaled=self.plot_autoscale,
            color_min=self.plot_min,
            color_max=self.plot_max)

    def anim_env(self, sim, cells, p):
        """
        Create 2D animation of env concentration.
        """

        env_time_series = [
            env.reshape(cells.X.shape) for env in self.c_env_time]
        AnimEnvTimeSeries(
            sim=sim, cells=cells, p=p,
            time_series=env_time_series,
            label=self.name + '_env',
            figure_title='Environmental ' + self.name,
            colorbar_title='Concentration [mmol/L]',
            is_color_autoscaled=self.plot_autoscale,
            color_min=self.plot_min,
            color_max=self.plot_max)

class Reaction(object):

    def __init__(self, sim, cells, p):

        self.dummy_dyna = TissueHandler(sim, cells, p)
        self.dummy_dyna.tissueProfiles(sim, cells, p)  # initialize all tissue profiles

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

        # activator and inhibitors of the reaction
        self.reaction_activators_list = None
        self.reaction_activators_Km = None
        self.reaction_activators_n = None
        self.reaction_inhibitors_list = None
        self.reaction_inhibitors_Km = None
        self.reaction_inhibitors_n = None

    def plot_1D(self, sim, cells, p, saveImagePath):

        if self.reaction_zone == 'cell':

            r_rate = [arr[p.plot_cell] for arr in self.rate_time]
            fig = plt.figure()
            ax = plt.subplot(111)
            ax.plot(sim.time, r_rate)
            ax.set_xlabel('Time [s]')
            ax.set_ylabel('Rate [mM/s]')
            ax.set_title('Rate of ' + self.name + ' in cell ' + str(p.plot_cell))

        elif self.reaction_zone == 'mit' and self.mit_enabled is True:

            r_rate = [arr[p.plot_cell] for arr in self.rate_time]
            fig = plt.figure()
            ax = plt.subplot(111)
            ax.plot(sim.time, r_rate)
            ax.set_xlabel('Time [s]')
            ax.set_ylabel('Rate [mM/s]')
            ax.set_title('Rate of ' + self.name + ' in mitochondria ' + str(p.plot_cell))

        if p.autosave is True:
            savename = saveImagePath + 'ReactionRate_' + self.name + '.png'
            plt.savefig(savename, format='png', transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

class Transporter(object):

    def __init__(self, sim, cells, p):

        self.dummy_dyna = TissueHandler(sim, cells, p)
        self.dummy_dyna.tissueProfiles(sim, cells, p)  # initialize all tissue profiles

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

        self.flux = None

    def init_reaction(self,cells, p):

        if self.transporter_profiles_list is not None and self.transporter_profiles_list != 'all':

            self.transporter_targets_mem = []
            self.transporter_targets_cell = []
            self.transporter_targets_env = []

            for profile in self.transporter_profiles_list:

                targets_cell = self.dummy_dyna.cell_target_inds[profile]
                self.transporter_targets_cell.extend(targets_cell)

                targets_mem = self.dummy_dyna.tissue_target_inds[profile]
                self.transporter_targets_mem.extend(targets_mem)

                targets_env = self.dummy_dyna.env_target_inds[profile]
                self.transporter_targets_env.extend(targets_env)

        elif self.transporter_profiles_list is None or self.transporter_profiles_list == 'all':

            self.transporter_targets_mem = cells.mem_i
            self.transporter_targets_cell = cells.cell_i
            self.transporter_targets_env = [x for x in range(0, len(cells.xypts))]

    def plot_1D(self, sim, cells, p, saveImagePath):

        if len(self.flux_time) > 0:

            if len(self.flux_time[0]) == sim.cdl:

                r_rate = [arr[p.plot_cell] for arr in self.flux_time]

            else:
                mem_i = cells.cell_to_mems[p.plot_cell][0]
                r_rate = [arr[mem_i] for arr in self.rate_time]

            fig = plt.figure()
            ax = plt.subplot(111)
            ax.plot(sim.time, r_rate)
            ax.set_xlabel('Time [s]')
            ax.set_ylabel('Rate [mM/s]')
            ax.set_title('Rate of ' + self.name + ' in cell ' + str(p.plot_cell))

            if p.autosave is True:
                savename = saveImagePath + 'TransporterRate_' + self.name + '.png'
                plt.savefig(savename, format='png', transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

    def update_transporter(self, sim, cells, p):

        self.dummy_dyna.tissueProfiles(sim, cells, p)  # initialize all tissue profiles
        self.init_reaction(cells, p)

class Channel(object):

    def __init__(self, sim, cells, p):

        self.dummy_dyna = TissueHandler(sim, cells, p)
        self.dummy_dyna.tissueProfiles(sim, cells, p)  # initialize all tissue profiles

    def init_channel(self, ion_string, type_string, max_val, sim, cells, p):

        # get targets for the reaction
        if self.channel_profiles_list is not None and self.channel_profiles_list != 'all':

            self.channel_targets_mem = []

            for profile in self.channel_profiles_list:

                targets_mem = self.dummy_dyna.tissue_target_inds[profile]
                self.channel_targets_mem.extend(targets_mem)

            self.channel_targets_mem = np.asarray(self.channel_targets_mem)

        elif self.channel_profiles_list is None or self.channel_profiles_list == 'all':

            self.channel_targets_mem = np.asarray(cells.mem_i)

        if ion_string == 'Na':

            self.dummy_dyna.maxDmNa = max_val
            self.dummy_dyna.targets_vgNa = self.channel_targets_mem
            class_string = vgna

        elif ion_string == 'NaP':

            self.dummy_dyna.maxDmNaP = max_val
            self.dummy_dyna.targets_vgNaP = self.channel_targets_mem
            class_string = vgnap

        elif ion_string == 'K':

            self.dummy_dyna.maxDmK = max_val
            self.dummy_dyna.targets_vgK = self.channel_targets_mem
            class_string = vgk

        elif ion_string == 'Kir':

            self.dummy_dyna.maxDmKir = max_val
            self.dummy_dyna.targets_vgKir = self.channel_targets_mem
            class_string = vgkir

        elif ion_string == 'Ca':

            self.dummy_dyna.maxDmCa = max_val
            self.dummy_dyna.targets_vgCa = self.channel_targets_mem
            class_string = vgca

        elif ion_string == 'Fun':

            self.dummy_dyna.maxDmFun = max_val
            self.dummy_dyna.targets_vgFun = self.channel_targets_mem
            class_string = vgfun

        else:

            raise BetseParametersException("Substance-modulated ion type not available. "
                                           "Valid choices: Na, K, Ca, NaP, Kir, and Fun")

            # create the desired voltage gated sodium channel instance:

        self.channel_core = getattr(class_string,type_string)()

        if p.run_sim is True:
            # initialize the channel object
            self.channel_core.init(self.dummy_dyna, sim, cells, p)

    def run_channel(self, sim, sim_metabo, cells, p):

        #FIXME do this part up in MoM, then have self.channel_core.modulator = eval(self.channel_core.mod_eval_string)
        # get modulation coefficients by any activating/inhibiting substances:
        activator_alpha, inhibitor_alpha = get_influencers(sim, sim_metabo, self.channel_activators_list,
            self.channel_activators_Km, self.channel_activators_n, self.channel_inhibitors_list,
            self.channel_inhibitors_Km, self.channel_inhibitors_n, reaction_zone='mems')

        # calculate the value of the channel modulation constant:
        self.channel_core.modulator = activator_alpha*inhibitor_alpha

        self.channel_core.run(self.dummy_dyna, sim, cells, p)

    def update_channel(self, sim, cells, p):
        self.dummy_dyna.tissueProfiles(sim, cells, p)  # initialize all tissue profiles
        self.init_channel(self.channel_class, self.channel_type, self.channelMax, sim, cells, p)

class Modulator(object):
    """
    The modulator object allows a substance defined in MasterOfMolecules to
    exert an activating or inhibiting influence over simulation-defined pumps
    and/or gap junctions.

    """
    def __init__(self):

        self.target_label = None
        self.max_val = None

        self.modulator_activators_list = None
        self.modulator_activators_Km = None
        self.modulator_activators_n = None
        self.modulator_inhibitors_list = None
        self.modulator_inhibitors_Km = None
        self.modulator_inhibitors_n = None

    def init_modulator(self, sim, cells, p):

        if self.target_label == 'gj':

            sim.gj_block_o = np.ones(sim.mdl)

        elif self.target_label == 'Na/K-ATPase':

            sim.NaKATP_block_o =  np.ones(sim.mdl)

        elif self.target_label == 'H/K-ATPase':

            sim.HKATP_block_o =  np.ones(sim.mdl)

        elif self.target_label == 'V-ATPase':

            sim.VATP_block_o =  np.ones(sim.mdl)

        else:

            raise BetseParametersException("You have requested a "
                                           "sim modulator that is not "
                                           "available. Available choices "
                                           "are: 'gj', 'Na/K-ATPase', 'H/K-ATPase', "
                                           "and 'V-ATPase' ")

    def run_modulator(self, sim, sim_metabo, cells, p):

        if self.zone == 'env':

            zone_tag = 'env'

        elif self.zone == 'cell':

            zone_tag = 'mems'

        else:

            raise BetseParametersException("You have requested an unavailable modulator zone."
                                           "Available choices are 'env' and 'mems'.")

        # get the coefficients activating and/or inhibiting the sim structure:

        # get modulation coefficients by any activating/inhibiting substances:
        # FIXME do this part up in MoM, then have self.channel_core.modulator = eval(self.channel_core.mod_eval_string)
        activator_alpha, inhibitor_alpha = get_influencers(sim, sim_metabo, self.modulator_activators_list,
            self.modulator_activators_Km, self.modulator_activators_n, self.modulator_inhibitors_list,
            self.modulator_inhibitors_Km, self.modulator_inhibitors_n, reaction_zone=zone_tag)

        # calculate the value of the channel modulation constant:
        modulator = self.max_val*activator_alpha * inhibitor_alpha

        # make size alteration for case of true environment:
        if p.sim_ECM is True and self.zone == 'env':
            modulator = modulator[cells.map_mem2ecm]

        if self.target_label == 'gj':

            sim.gj_block = modulator

        elif self.target_label == 'Na/K-ATPase':

            sim.NaKATP_block = modulator

        elif self.target_label == 'H/K-ATPase':

            sim.HKATP_block = modulator

        elif self.target_label == 'V-ATPase':

            sim.VATP_block = modulator

        elif self.target_label == 'Ca-ATPase':

            sim.CaATP_block = modulator

        elif self.target_label == 'Na/Ca-Exch':

            sim.NaCaExch_block = modulator

        else:

            raise BetseParametersException("You have requested a "
                                           "sim modulator that is not "
                                           "available. Available choices "
                                           "are: 'gj', 'Na/K-ATPase', 'H/K-ATPase', "
                                           "and 'V-ATPase', 'Ca-ATPase', and 'Na/Ca-Exch' ")


