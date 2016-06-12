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

# FIXME we will need another reactions field and read-in case for whether reaction occurs in cell or in mitochondria

class MasterOfMolecules(object):

    def __init__(self, sim, config_substances, p):

        """
         Initializes the MasterOfMolecules object.

        sim:                An instance of simulator
        config_settings     List of dictionaries storing key settings (p.molecules_config)
        p                   An instance of params
        """

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
            gad = mol_dic['growth and decay']

            if gad != 'None':

                obj.simple_growth = True

                obj.r_production = gad['production rate']
                obj.r_decay = gad['decay rate']
                obj.Kgd = gad['Km']

            else:
                obj.simple_growth = False

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

            obj.transfer_membrane = react_dic['transmembrane transfer']

            if obj.transfer_membrane != 'None':
                obj.transfer_V_type = obj.transfer_membrane[0]
                obj.transfer_direction = obj.transfer_membrane[1]

            else:
                obj.transfer_membrane = None

            # now we want to load the right concentration data arrays for reactants into the Reaction object:

            # find out which items (if any) are to be transferred across the membrane -- these are items that occur in
            # both the reactants and the products names list:
            name_set_products = set(obj.products_list)
            name_set_reactants = set(obj.reactants_list)
            # list of substances moved from one side of the membrane to the other
            obj.transfers_list = list(name_set_products & name_set_reactants)

            # store indices of items that are transferred across the membrane from products and reactants lists:--------
            if len(obj.transfers_list) > 0:

                for name in obj.transfers_list:

                    ind_react = obj.reactants_list.index(name)
                    ind_prod = obj.products_list.index(name)

                    obj.transfered_react_inds.append(ind_react)
                    obj.transfered_prod_inds.append(ind_prod)

            # set the reaction states by harvesting the appropriate concentration data:
            self.prime_reactions(obj, sim, cells, p)

    def prime_reactions(self, obj, sim, cells, p):

        # Get into some nasty stuff as we'll define multiple cases for reading in different kinds of concentration
        # arrays: --------------------------------------------------------------------------------------------------

        if obj.reaction_zone == 'cell':

            # Case 1: no transfer across the membrane -- these are in-cytoplasm reactions -------------------------------
            for i, reactant_name in enumerate(obj.reactants_list):

                if reactant_name not in obj.transfers_list: # if we're not dealing with a trans-membrane transfer:

                    self.set_react_concs(obj, sim, cells, p, reactant_name, i, 'c_cells', 'cc_cells')

            # Now load the right concentration data arrays for products into the Reaction object:
            for j, product_name in enumerate(obj.products_list):

                if product_name not in obj.transfers_list:

                    self.set_prod_concs(obj, sim, cells, p, product_name, j, 'c_cells', 'cc_cells')

            # case 2: transfer occurs across a membrane ----------------------------------------------------------------
            if obj.transfer_membrane is not None:

                if obj.transfer_V_type == 'Vmem':  # case 2a: transfer across cell plasma membrane

                    obj.V = sim.vm # link voltage to the membrane voltage

                    if obj.transfer_direction == 'in': # direction into cell

                        for i, reactant_name in enumerate(obj.reactants_list):

                            if reactant_name in obj.transfers_list:

                                self.set_react_concs(obj, sim, cells, p, reactant_name, i, 'c_env', 'cc_env')

                            # Now load the right concentration data arrays for products into the Reaction object:
                        for j, product_name in enumerate(obj.products_list):

                            if product_name in obj.transfers_list:

                                self.set_prod_concs(obj, sim, cells, p, product_name, j, 'c_cells', 'cc_cells')

                    elif obj.transfer_direction == 'out': # directio into environment

                        for i, reactant_name in enumerate(obj.reactants_list):

                            if reactant_name in obj.transfers_list:

                                self.set_react_concs(obj, sim, cells, p, reactant_name, i, 'c_cells', 'cc_cells')

                            # Now load the right concentration data arrays for products into the Reaction object:

                        for j, product_name in enumerate(obj.products_list):

                            if product_name in obj.transfers_list:

                                self.set_prod_concs(obj, sim, cells, p, product_name, j, 'c_env', 'cc_env')

                elif obj.tranfer_V_type == 'Vmit':  # case 2b: transfer across mitochondrial membrane

                    # FIXME need to think about how organelle voltages will be handled
                    obj.V = getattr(sim, 'vmit', None) # try and get a mitochondrial voltage to link

                    if obj.transfer_direction == 'in':  # direction into mitochondria

                        for i, reactant_name in enumerate(obj.reactants_list):

                            if reactant_name in obj.transfers_list:

                                self.set_react_concs(obj, sim, cells, p, reactant_name, i, 'c_cells', 'cc_cells')

                            # Now load the right concentration data arrays for products into the Reaction object:
                        for j, product_name in enumerate(obj.products_list):

                            if product_name in obj.transfers_list:

                                self.set_prod_concs(obj, sim, cells, p, product_name, j, 'c_mit', 'cc_mit')

                    elif obj.transfer_direction == 'out': # direction out to cytoplasm

                        for i, reactant_name in enumerate(obj.reactants_list):

                            if reactant_name in obj.transfers_list:

                                self.set_react_concs(obj, sim, cells, p, reactant_name, i, 'c_mit', 'cc_mit')

                            # Now load the right concentration data arrays for products into the Reaction object:

                        for j, product_name in enumerate(obj.products_list):

                            if product_name in obj.transfers_list:

                                self.set_prod_concs(obj, sim, cells, p, product_name, j, 'c_cells', 'cc_cells')

                elif obj.tranfer_V_type == 'Ver':  # case 2b: transfer across endoplasmic reticulum (ER) membrane

                    # FIXME need to think about how organelle voltages will be handled
                    obj.V = getattr(sim, 'ver', None)  # try and get an er voltage to link-to

                    if obj.transfer_direction == 'in':  # direction into ER

                        for i, reactant_name in enumerate(obj.reactants_list):

                            if reactant_name in obj.transfers_list:

                                self.set_react_concs(obj, sim, cells, p, reactant_name, i, 'c_cells', 'cc_cells')

                            # Now load the right concentration data arrays for products into the Reaction object:
                        for j, product_name in enumerate(obj.products_list):

                            if product_name in obj.transfers_list:

                                self.set_prod_concs(obj, sim, cells, p, product_name, j, 'c_er', 'cc_er')

                    elif obj.transfer_direction == 'out': # direction out to cytoplasm

                        for i, reactant_name in enumerate(obj.reactants_list):

                            if reactant_name in obj.transfers_list:

                                label = 'i' + reactant_name
                                ion_check = getattr(sim, label, None)

                                if ion_check is None:

                                    try:
                                        obj_reactant = getattr(self, reactant_name)
                                        c_react = obj_reactant.c_er
                                        obj.c_reactants.append(c_react)
                                        obj.z_reactants.append(obj_reactant.z)

                                    except:

                                        raise BetseExceptionParameters(
                                            'Name of reactant is not a defined chemical, or is not'
                                            'an ion currently included in the ion profile being used.'
                                            'Please check biomolecule definitions and ion profile'
                                            'settings of your config(s) and try again.')

                                else:
                                    # define the reactant as the ion concentration from the cell concentrations object in sim:
                                    # FIXME we need to make sure this concentration exists...somehow...
                                    obj.c_reactants.append(sim.cc_er[ion_check])
                                    obj.c_reactants.append(sim.zs[ion_check])

                                # add the index of the reaction to the list so we can access modifiers like reaction coefficients
                                # and Km values at a later point:
                                obj.inds_react.append(i)

                            # Now load the right concentration data arrays for products into the Reaction object:

                        for j, product_name in enumerate(obj.products_list):

                            if product_name in obj.transfers_list:

                                # see if the name is an ion defined in sim:
                                label = 'i' + product_name
                                ion_check = getattr(sim, label, None)

                                if ion_check is None:

                                    try:
                                        obj_product = getattr(self, product_name)
                                        c_prod = obj_product.c_cells
                                        obj.c_products.append(c_prod)
                                        obj.z_products.append(obj_product.z)

                                    except:

                                        raise BetseExceptionParameters(
                                            'Name of product is not a defined chemical, or is not'
                                            'an ion currently included in the ion profile being used.'
                                            'Please check biomolecule definitions and ion profile'
                                            'settings of your config(s) and try again.')

                                else:
                                    # define the reactant as the ion concentration from the cell concentrations object in sim:
                                    # FIXME we need to make sure this is a link, so that the sim concentration will be auto-updated here
                                    obj.c_products.append(sim.cc_cells[ion_check])
                                    obj.z_products.append(sim.zs[ion_check])

                                # add the index of the reaction to the list so we can access modifiers like reaction coefficients
                                # and Km values at a later point:
                                obj.inds_prod.append(j)

        elif obj.reaction_zone == 'mitochondria':

            # Case 1: no transfer across the membrane -- these are in-cytoplasm reactions -------------------------------
            for i, reactant_name in enumerate(obj.reactants_list):

                if reactant_name not in obj.transfers_list:  # if we're not dealing with a trans-membrane transfer:

                    self.set_react_concs(obj, sim, cells, p, reactant_name, i, 'c_mit', 'cc_mit')

            # Now load the right concentration data arrays for products into the Reaction object:
            for j, product_name in enumerate(obj.products_list):

                if product_name not in obj.transfers_list:

                    self.set_prod_concs(obj, sim, cells, p, product_name, j, 'c_mit', 'cc_mit')

            # case 2: transfer occurs across a membrane ----------------------------------------------------------------
            if obj.transfer_membrane is not None:

                if obj.tranfer_V_type == 'Vmit':  # case 2b: transfer across mitochondrial membrane

                    # FIXME need to think about how organelle voltages will be handled
                    obj.V = getattr(sim, 'vmit', None)  # try and get a mitochondrial voltage to link

                    if obj.transfer_direction == 'in':  # direction into mitochondria

                        for i, reactant_name in enumerate(obj.reactants_list):

                            if reactant_name in obj.transfers_list:

                                self.set_react_concs(obj, sim, cells, p, reactant_name, i, 'c_cells', 'cc_cells')

                                # Now load the right concentration data arrays for products into the Reaction object:

                        for j, product_name in enumerate(obj.products_list):

                            if product_name in obj.transfers_list:

                                self.set_prod_concs(obj, sim, cells, p, product_name, j, 'c_mit', 'cc_mit')

                    elif obj.transfer_direction == 'out':  # direction out to cytoplasm

                        for i, reactant_name in enumerate(obj.reactants_list):

                            if reactant_name in obj.transfers_list:

                                self.set_react_concs(obj, sim, cells, p, reactant_name, i, 'c_mit', 'cc_mit')

                                # Now load the right concentration data arrays for products into the Reaction object:

                        for j, product_name in enumerate(obj.products_list):

                            if product_name in obj.transfers_list:

                                self.set_prod_concs(obj, sim, cells, p, product_name, j, 'c_cells', 'cc_cells')

        else:
            raise BetseExceptionParameters("You have requested a reaction zone that does not exist."
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

            try:
                obj_reactant = getattr(self, reactant_name)
                c_react = getattr(obj_reactant, reactant_type_self)
                obj.c_reactants.append(c_react)
                obj.z_reactants.append(obj_reactant.z)

                obj.reactant_source_object.append(id(self))
                obj.reactant_source_type.append(reactant_type_self)

            except KeyError:

                raise BetseExceptionParameters('Name of product is not a defined chemical, or is not'
                                               'an ion currently included in the ion profile being used.'
                                               'Please check biomolecule definitions and ion profile'
                                               'settings of your config(s) and try again.')

        else:

            obj.reactant_source_object.append(id(sim))
            obj.reactant_source_type.append(reactant_type_sim)

            # define the reactant as the ion concentration from the cell concentrations object in sim:
            sim_conco = getattr(sim, reactant_type_sim)

            if reactant_type_sim == 'cc_env':
                # define the reactant as the ion concentration from the cell concentrations object in sim:
                if p.sim_ECM:  # get the env concentration from the xy grid to cell centres:
                    sim_conc = sim_conco[ion_check][cells.map_cell2ecm]

                else:  # average the env concentrations to the cell centres
                    sim_conc = np.dot(cells.M_sum_mems, sim_conco[ion_check]) / cells.num_mems

            else:
                sim_conc = sim_conco[ion_check]

            obj.c_reactants.append(sim_conc)
            obj.z_reactants.append(sim.zs[ion_check])

        # add the index of the reaction to the list so we can access modifiers like reaction coefficients
        # and Km values at a later point:
        obj.inds_react.append(i)

    def set_prod_concs(self, obj, sim, cells, p, product_name, j, product_type_self, product_type_sim):

        # FIXME put in a clause for c_mit or cc_mit undefined!

        # Now load the right concentration data arrays for products into the Reaction object:
        # see if the name is an ion defined in sim:
        label = 'i' + product_name
        ion_check = getattr(sim, label, None)

        if ion_check is None:

            try:
                obj_product = getattr(self, product_name)
                c_prod = getattr(obj_product, product_type_self)
                obj.c_products.append(c_prod)
                obj.z_products.append(obj_product.z)

                obj.product_source_object.append(id(self))
                obj.product_source_type.append(product_type_self)

            except KeyError:

                raise BetseExceptionParameters('Name of product is not a defined chemical, or is not'
                                               'an ion currently included in the ion profile being used.'
                                               'Please check biomolecule definitions and ion profile'
                                               'settings of your config(s) and try again.')

        else:

            obj.product_source_object.append(id(sim))
            obj.product_source_type.append(product_type_sim)

            # define the reactant as the ion concentration from the cell concentrations object in sim:
            sim_conco = getattr(sim, product_type_sim)

            if product_type_sim == 'cc_env':
                # define the reactant as the ion concentration from the cell concentrations object in sim:
                if p.sim_ECM:  # get the env concentration from the xy grid to cell centres:
                    sim_conc = sim_conco[ion_check][cells.map_cell2ecm]

                else:  # average the env concentrations to the cell centres
                    sim_conc = np.dot(cells.M_sum_mems, sim_conco[ion_check]) / cells.num_mems

            else:
                sim_conc = sim_conco[ion_check]

            obj.c_products.append(sim_conc)
            obj.z_products.append(sim.zs[ion_check])

        # add the index of the reaction to the list so we can access modifiers like reaction coefficients
        # and Km values at a later point:
        obj.inds_prod.append(j)

    def quick_prime_reactions(self, sim, cells, p):

        for reaction_name in self.reaction_names:

            obj_reaction = getattr(self, reaction_name)

            # clear the reaction and product concentration lists:
            obj_reaction.c_reactants = []
            obj_reaction.c_products = []

            for i, react_name in enumerate(obj_reaction.reactants_list):

                if obj_reaction.reactant_source_object[i] == id(self):

                    molecule = getattr(self, react_name)

                    c_react = getattr(molecule, obj_reaction.reactant_source_type[i])
                    obj_reaction.c_reactants.append(c_react)

                elif obj_reaction.reactant_source_object[i] == id(sim):

                    ion_array = getattr(sim, obj_reaction.reactant_source_type[i])
                    ion_label = 'i' + react_name
                    ion_type = getattr(sim, ion_label)
                    c_reacto = ion_array[ion_type]

                    if obj_reaction.reactant_source_type[i] == 'cc_env':

                        if p.sim_ECM:

                            c_react = c_reacto[cells.map_cell2ecm]

                        else:
                            c_react = np.dot(cells.M_sum_mems, c_reacto)/cells.num_mems

                    else:

                        c_react = c_reacto

                    obj_reaction.c_reactants.append(c_react)


            for j, prod_name in enumerate(obj_reaction.products_list):

                if obj_reaction.product_source_object[j] == id(self):

                    molecule = getattr(self, prod_name)

                    c_prod = getattr(molecule, obj_reaction.product_source_type[j])
                    obj_reaction.c_products.append(c_prod)

                elif obj_reaction.product_source_object[j] == id(sim):

                    ion_array = getattr(sim, obj_reaction.product_source_type[j])
                    ion_label = 'i' + prod_name
                    ion_type = getattr(sim, ion_label)
                    c_prodo = ion_array[ion_type]

                    if obj_reaction.product_source_type[j] == 'cc_env':

                        if p.sim_ECM:

                            c_prod = c_prodo[cells.map_cell2ecm]

                        else:
                            c_prod = np.dot(cells.M_sum_mems, c_prodo) / cells.num_mems

                    else:

                        c_prod = c_prodo

                    obj_reaction.c_products.append(c_prod)

    def assign_new_concentrations(self, obj, delta_c, new_reactants, new_products, sim, cells, p):

        # FIXME need to convert delta_c into a flux to update volumes properly for transmem transfers...!!!
        # if we're not transfering across any membranes; staying in same reaction zone:
        if obj.transfer_membrane is None:

            for i, react_name in enumerate(obj.reactants_list):

                if obj.reactant_source_object[i] == id(self):

                    source_obj = getattr(self, react_name)

                    setattr(source_obj, obj.reactant_source_type[i], new_reactants[i])

                    print('assigning new reactant field to MasterOfMolecules: ', react_name)

                elif obj.reactant_source_object[i] == id(sim):

                    source_obj = getattr(sim, obj.reactant_source_type[i])
                    ion_label = 'i' + react_name
                    ion_type = getattr(sim, ion_label)

                    source_obj[ion_type] = new_reactants[i]

                    print('assigning new reaction field to sim: ', react_name)

            for j, prod_name in enumerate(obj.products_list):

                    if obj.product_source_object[j] == id(self):

                        source_obj = getattr(self, prod_name)

                        setattr(source_obj, obj.product_source_type[j], new_products[j])

                        print('assigning new product field to MasterOfMolecules: ', prod_name)

                    elif obj.product_source_object[j] == id(sim):

                        source_obj = getattr(sim, obj.product_source_type[j])
                        ion_label = 'i' + prod_name
                        ion_type = getattr(sim, ion_label)

                        source_obj[ion_type] = new_products[j]

                        print('assigning new product field to sim: ', prod_name)

        else:

            print("transfer membrane is true, but nothing is worked up yet. Pass on it.")

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

            obj.pump(sim, cells, p)
            obj.transport(sim, cells, p)
            obj.updateIntra(sim, cells, p)

            if obj.simple_growth is True:
                obj.reaction(sim, cells, p)

            if p.run_sim is True:
                obj.gating(sim, cells, p)
                obj.update_boundary(t, p)

    def run_loop_reactions(self, sim, cells, p):

        # get the object corresponding to the specific reaction:
        for i, name in enumerate(self.reaction_names):

            # get the Reaction object:
            obj = getattr(self, name)

            # compute the new reactants and products
            delta_c, new_reactants, new_products = obj.compute_reaction(sim, cells, p)

            # assign new concentrations to the Molecules or sim fields:
            self.assign_new_concentrations(obj, delta_c, new_reactants, new_products, sim, cells, p)

            # re-assign new concentrations to the Reaction object:
            self.quick_prime_reactions(sim, cells, p)

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
        self.pump_to_cell = None
        self.pump_max_val = None
        self.pump_Km = None
        self.use_gating_ligand = None
        self.gating_extracell = None
        self.gating_max_val = None
        self.gating_Hill_K = None
        self.gating_Hill_n = None
        self.gating_ion = None
        self.r_production = None
        self.r_decay = None
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
        self.transfer_membrane = None
        self.transfer_V_type = None
        self.transfer_direction = None

        self.reaction_zone = None
        self.V = None

        self.c_reactants = []  # concentrations of reactions, defined on cell-centres
        self.z_reactants = []  # charge of reactants
        self.inds_react = []  # indices to each reactant, to keep their order
        self.c_products = []  # concentrations of reactions, defined on cell-centres
        self.z_products = []  # charge of products
        self.inds_prod = []  # indices to each reactant, to keep their order
        self.transfered_react_inds = []  # indices of reactants transferred across a membrane
        self.transfered_prod_inds = []

        self.product_source_object = []  # object from which product concentrations come from
        self.product_source_type = []    # type of product concentrations sourced from object
        self.reactant_source_object = []  # object from which reactants concentrations come from
        self.reactant_source_type = []  # type of reactant concentrations sourced from object

    def compute_reaction(self, sim, cells, p):

        # if reaction is reversible, calculate reaction quotient Q, the equilibrium constant, and backwards rate
        if self.delta_Go is not None:

            # define the reaction equilibrium coefficient:
            Keqm = np.exp(-self.delta_Go/(p.R*sim.T))

            # calculate the MM rate coefficient for the backwards reaction direction and Q:
            backwards_term = []
            Q_deno_list = []
            Q_numo_list = []

            # if reaction involves a trans-membrane transfer, init two arrays
            if self.transfer_membrane is not None:
                Qe_prod_list = []
                Qe_react_list = []

            for i, prod in enumerate(self.c_products):
                # calculate a factor in the backwards rate term:
                coeff = self.products_coeff[i]
                Km = self.Km_products_list[i]
                cs = (prod / Km) ** coeff
                term = cs / (1 + cs)
                backwards_term.append(term)

                if self.transfer_membrane is not None and i in self.transfered_prod_inds:

                    if self.transfer_direction is 'in':

                        Qe_prod_list.append(self.z_products[i]*p.F*self.V)

                    elif self.transfer_direction is 'out':

                        Qe_prod_list.append(0)

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

                if self.transfer_membrane is not None and j in self.transfered_react_inds:

                    if self.transfer_direction is 'in':

                        Qe_react_list.append(0)

                    elif self.transfer_direction is 'out':

                        Qe_react_list.append(self.z_reactants[j] * p.F * self.V)

            # get the numerator and denomenator of the reaction quotient:
            Q_deno = np.prod(Q_deno_list, axis=0)
            Q_numo = np.prod(Q_numo_list, axis=0)

            # finally, *the* reaction quotient:
            Q = Q_numo/Q_deno

            # if the reaction involves a transfer across a membrane, alter Keqm appropriately:
            if self.transfer_membrane is not None:
                Qe_react = np.sum(Qe_react_list, axis = 0)
                Qe_prod = np.sum(Qe_prod_list, axis = 0)

                Qe = Qe_prod - Qe_react

                Keqm = np.exp(-(self.delta_Go + Qe)/(p.R*sim.T))  # FIXME check this

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

        # final degree of change, returned for use elsewhere:
        deltaC = reaction_rate*p.dt

        # update concentrations:
        for a, c_react in enumerate(self.c_reactants):

            self.c_reactants[a] = self.c_reactants[a] - deltaC

        for b, c_prod in enumerate(self.c_products):

            self.c_products[b] = self.c_products[b] + deltaC

        return deltaC, self.c_reactants, self.c_products






