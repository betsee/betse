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

        self.ave_cell_vol = cells.cell_vol.mean()  # average cell volume

        self.reaction_names = []
        self.transporter_names = []
        self.channel_names = []
        self.modulator_names = []

        self.plot_cmap = 'viridis'

    def read_substances(self, sim, cells, config_substances, p):
        """
            Initializes all core data structures and concentration variables
            for all molecules included in the simulation, as well as any ions present in sim.

            config_substances:  dictionary containing BETSE biomolecule template fields

        """

        # Initialize dictionary mapping from molecule names to Molecule objects
        self.molecules = {}
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

    def tissue_init(self, sim, cells, config_substances, p):
        """
        Completes the initialization process of each molecule with additional
        fields, but doesn't touch the concentrations. This is for the case where
        the user changes config file settings after running an init, and ensures
        new parameters are updated.

        """

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

                # factors involving auto-catalytic growth and decay (gad) in the cytoplasm
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
                    mol.growth_activators_k = gad.get('k activators', None)
                    mol.growth_activators_Km = gad.get('Km activators', None)

                    mol.growth_activators_n = gad.get('k activators', None)
                    mol.growth_inhibitors_list = gad.get('inhibitors', None)
                    mol.growth_inhibitors_k = gad.get('k inhibitors', None)
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
            mol = self.molecules['name']

            # assign plotting properties
            pd = mol_dic['plotting']
            mol.make_plots = pd['plot 2D']
            mol.make_ani = pd['animate']

            mol.plot_autoscale = pd['autoscale colorbar']
            mol.plot_max = pd['max val']
            mol.plot_min = pd['min val']

    def read_reactions(self, config_reactions, sim, cells, p):

        """
          Read in and initialize parameters for all user-defined reactions.

          config_options:  dictionary containing BETSE reaction template fields

          """

        # Initialize a list that keeps the Reaction objects:
        self.reactions = {}

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
            obj.reaction_inhibitors_list = react_dic.get('reaction inhibitors', None)
            obj.reaction_inhibitors_Km = react_dic.get('inhibitor Km', None)
            obj.reaction_inhibitors_n = react_dic.get('inhibitor n', None)

            if self.mit_enabled:
                obj.mit_enabled = True
            else:
                obj.mit_enabled = False

            # now we want to load the right concentration data arrays for reactants into the Reaction object:

    def write_reactions(self):
        """
        Reactions are now constructed during the init as strings that are evaluated in eval calls in each time-step.
        This function constructs the evaluation strings for each reaction, given the metadata stored
        in each reaction object (e.g. lists of reactants, products, etc).

        """


        for reaction_name in self.reactions:

            reactant_names = self.reactions[reaction_name].reactants_list
            reactant_coeff = self.reactions[reaction_name].reactants_coeff
            reactant_Km = self.reactions[reaction_name].Km_reactants_list

            product_names = self.reactions[reaction_name].products_list
            product_coeff = self.reactions[reaction_name].products_coeff
            product_Km = self.reactions[reaction_name].Km_products_list

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

            #self.reactions[reaction_name].reaction_eval_string = numo_string + '/' + denomo_string

            # call statement to evaluate:
            # eval(self.reactions['consume_ATP'].reaction_eval_string, globals(), locals())

    # ------Utility Methods---------------------------------------------------------------------------------------------

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

        delta_cells = self.growth_mod_function_cells*self.r_production*inhibitor_alpha*activator_alpha - self.r_decay*cc


        self.c_cells[self.growth_targets_cell] = self.c_cells[self.growth_targets_cell] + \
                                                 delta_cells[self.growth_targets_cell]*p.dt

        self.c_mems[self.growth_targets_mem] = self.c_mems[self.growth_targets_mem] + \
                                               delta_cells[cells.mem_to_cells][self.growth_targets_mem]*p.dt

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

                    raise BetseParametersException('Name of product is not a defined chemical, or is not'
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

                    raise BetseParametersException('Name of product is not a defined chemical, or is not'
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

        elif self.reaction_zone == 'mitochondria' and self.mit_enabled is True:

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

            inds_neg = (Q_deno == 0).nonzero()
            Q_deno[inds_neg] = 1.0e-15  # ensure no division by zero

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

        if self.reaction_zone == 'cell':

            tag = 'cell'

        elif self.reaction_zone == 'mitochondria' and self.mit_enabled is True:

            tag = 'mitochondria'

        else:

            tag = None

        if tag is not None:

            # get net effect of any activators or inhibitors of the reaction:
            activator_alpha, inhibitor_alpha = get_influencers(sim, sim_metabo, self.reaction_activators_list,
                self.reaction_activators_Km, self.reaction_activators_n, self.reaction_inhibitors_list,
                self.reaction_inhibitors_Km, self.reaction_inhibitors_n, reaction_zone=tag)

            # final degree of change, returned for use elsewhere (?):
            flux = activator_alpha*inhibitor_alpha*reaction_rate
            deltaC = activator_alpha*inhibitor_alpha*reaction_rate*p.dt

            if self.reaction_zone == 'cell':

                self.set_reactant_c(deltaC, sim, sim_metabo,'c_cells', 'cc_cells')

                # re-obtain updated concentrations:
                self.get_reactants(sim, sim_metabo, 'c_cells', 'cc_cells')
                self.get_products(sim, sim_metabo, 'c_cells', 'cc_cells')

                self.set_product_c(deltaC, sim, sim_metabo, 'c_cells', 'cc_cells')


            if self.reaction_zone == 'mitochondria' and self.mit_enabled is True:

                self.set_reactant_c(deltaC, sim, sim_metabo,'c_mit', 'cc_mit')

                # obtain updated concentrations:
                self.get_reactants(sim, sim_metabo, 'c_mit', 'cc_mit')
                self.get_products(sim, sim_metabo, 'c_mit', 'cc_mit')

                self.set_product_c(deltaC, sim, sim_metabo, 'c_mit', 'cc_mit')

        else:

            flux = None

        return flux

    def plot_1D(self, sim, cells, p, saveImagePath):

        if self.reaction_zone == 'cell':

            r_rate = [arr[p.plot_cell] for arr in self.rate_time]
            fig = plt.figure()
            ax = plt.subplot(111)
            ax.plot(sim.time, r_rate)
            ax.set_xlabel('Time [s]')
            ax.set_ylabel('Rate [mM/s]')
            ax.set_title('Rate of ' + self.name + ' in cell ' + str(p.plot_cell))

        elif self.reaction_zone == 'mitochondria' and self.mit_enabled is True:

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


