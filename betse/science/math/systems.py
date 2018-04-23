#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

#FIXME: This module, which appears to be a partial copy-and-paste of the
#"betse.science.chemistry.networks" module with minor differences, doesn't
#appear to be used anymore. That's great! Would removing this module or perhaps
#shifting its most important functionality into the "networks" module be
#possible? If not, no worries! Translucent lucidity erupts the thunder!

# ....................{ IMPORTS                            }....................
import csv
import matplotlib.pyplot as plt
import numpy as np
from betse.exceptions import BetseSimConfException
from betse.lib import libs
from betse.util.io.log import logs
from betse.util.path import dirs, pathnames
from betse.util.type.mapping.mapcls import DynamicValue, DynamicValueDict
from collections import OrderedDict
from matplotlib import colors

# ....................{ CLASSES                            }....................
class SimMaster(object):

    def __init__(self, config_dic, p, mit_enabled = False):
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
        # Initialize a dict that keeps the Reaction objects specific to environment:
        self.reactions_env = OrderedDict({})
        # Initialize a dict of Transporter objects:
        self.transporters = OrderedDict({})
        # Initialize a dict of Channels:
        self.channels = OrderedDict({})
        # Initialize a dict of modulators:
        self.modulators = OrderedDict({})

        self.reactions_mit = OrderedDict({})

        # Initialize reaction rates array to None (filled in later, if applicable):
        self.reaction_rates = None

        # boolean so that charge will only ever be balanced once:
        self.charge_has_been_balanced = False



        self.vmo = float(p.network_config['optimization']['target Vmem'])

        self.vm = self.vmo
        self.cdl = 1
        self.edl = 1
        self.mdl = 1

        self.div_cell = 5.0e5

        # obtain main dictionaries:
        # obtain specific sub-dictionaries from the config file:
        self.substances_config = config_dic['biomolecules']
        self.reactions_config = config_dic.get('reactions', None)
        self.transporters_config = config_dic.get('transporters', None)
        self.channels_config = config_dic.get('channels', None)

        self.direction_surface_config = config_dic.get('direction surfaces', None)

        # read in substance properties from the config file, and initialize basic properties:
        self.read_substances(self.substances_config, p)
        self.tissue_init(self.substances_config, p)


        if self.reactions_config is not None:
            # initialize the reactions of metabolism:
            self.read_reactions(self.reactions_config, p)
            self.write_reactions()
            self.create_reaction_matrix()
            self.write_reactions_env()
            self.create_reaction_matrix_env()
            self.has_reactions = True

        else:
            self.create_reaction_matrix()
            self.create_reaction_matrix_env()
            self.has_reactions = False

        # initialize transporters, if defined:
        if self.transporters_config is not None:
            self.read_transporters(self.transporters_config, p)
            self.write_transporters(p)

            self.has_transporters = True

        else:
            self.has_transporters = False

        # initialize channels, if desired:
        if self.channels_config is not None:
            self.read_channels(self.channels_config, p)
            self.has_channels = True

        else:
            self.has_channels = False


        # colormap for plotting series of 1D lines:
        self.plot_cmap = 'viridis'

        # default shape identities for network:
        self.reaction_shape = 'rect'
        self.transporter_shape = 'diamond'
        self.channel_shape = 'pentagon'
        self.vmem_shape = 'ellipse'
        self.ed_shape = 'hexagon'
        self.conc_shape = 'oval'

        self.globals = globals()
        self.locals = locals()

        # Build core data structures:

        self.build_sim_matrix(p)

        self.init_saving(p)

        # direction field plots:
        if self.direction_surface_config is not None:

            if len(self.direction_surface_config):

                logs.log_info("Plotting direction surfaces for BIGR networks...")

                self.direction_field_plotter(p)

    def direction_field_plotter(self, p):

        for dsparams in self.direction_surface_config:

            npoints = int(dsparams['number of points'])

            ap = {}
            bp = {}

            ap['name'] = dsparams['substance X name']
            bp['name'] = dsparams['substance Y name']

            ap['A_min'] = float(dsparams['range of X'][0])
            ap['A_max'] = float(dsparams['range of X'][1])

            bp['B_min'] = float(dsparams['range of Y'][0])
            bp['B_max'] = float(dsparams['range of Y'][1])

            alpha = float(dsparams['cmap alpha'])

            self.direction_field(p, ap, bp, npoints=npoints, calpha=alpha)

    def optimize_network(self, config, p):

        pass

    def read_substances(self, config_substances, p):
        """
        Initialize all core data structures and concentration variables for all
        molecules included in the simulation, as well as any ions present in
        sim.

        config_substances:  dictionary containing BETSE biomolecule template fields
        """

        logs.log_info("Reading additional substance data...")

        # Initialize a dictionaries that will eventually hold dynamic values for cell, env and mit concentrations:
        cell_concs_mapping = {}
        env_concs_mapping = {}

        # place to store charge values and membrane diffusivity for all substances AND sim ions:
        self.zmol = {}
        self.Dmem = {}
        self.ED_eval_strings ={}

        # begin by adding all sim ions to the dynamic concentration mapping structures:
        for k, val in p.ions_dict.items():

            if val == 1 and k != 'H': # if the ion is used in the simulation

                # add a field to the MasterOfMolecules corresponding to Molecule object with that name:
                self.molecules[k] = Molecule(p)

                # now set the attributes of that Molecule object with the cornucopia of variables:
                mol = self.molecules[k]

                mol.simple_growth = False

                # assign general properties
                mol.name = k  # let object know who it is

                mol.c_envo = p.env_concs[k]  # initial concentration in the environment [mmol/L]
                mol.c_cello = p.cell_concs[k]  # initial concentration in the cytoplasm [mmol/L]

                # create concentration data arrays:
                mol.c_cells = np.ones(self.cdl) * mol.c_cello
                mol.c_env = np.ones(self.cdl) * mol.c_envo

                name = k

                # create dynamic mappings for the cell and env conc vectors:
                cell_concs_mapping[name] = DynamicValue(
                    lambda name=name: self.molecules[name].c_cells,
                    lambda value, name=name: setattr(self.molecules[name], 'c_cells', value))

                env_concs_mapping[name] = DynamicValue(
                    lambda name=name: self.molecules[name].c_env,
                    lambda value, name=name: setattr(self.molecules[name], 'c_env', value))


                mol.Dm = p.mem_perms[k]
                mol.Do = p.free_diff[k]
                mol.Dgj = p.free_diff[k]*p.gj_surface  # effective diffusion coefficient of substance through GJ
                mol.z = p.ion_charge[k]  # charge (oxidation state)

                self.Dmem[k] = p.mem_perms[k]
                self.zmol[k] = p.ion_charge[k]


        # Now move on to building the data structures for the user-defined substances:
        for q, mol_dic in enumerate(config_substances):
            # get each user-defined name-filed in the dictionary:
            name = str(mol_dic['name'])

            # add a field to the MasterOfMolecules corresponding to Molecule object with that name:
            self.molecules[name] = Molecule(p)

            # now set the attributes of that Molecule object with the cornucopia of variables:
            # define an alias so we don't have to keep typing the long version:
            mol = self.molecules[name]

            # assign general properties
            mol.name = name   # let object know who it is

            mol.c_envo = mol_dic['env conc']  # initial concentration in the environment [mmol/L]
            mol.c_cello = mol_dic['cell conc']  # initial concentration in the cytoplasm [mmol/L]

            # create concentration data arrays:
            mol.c_cells = np.ones(self.cdl)*mol.c_cello

            # create dynamic mappings for the cell and mem conc vectors:
            cell_concs_mapping[name] = DynamicValue(
                lambda name = name: self.molecules[name].c_cells,
                lambda value, name = name: setattr(self.molecules[name], 'c_cells', value))

            # initialize concentration in the environment:
            mol.c_env = np.ones(self.cdl) * mol.c_envo

            env_concs_mapping[name] = DynamicValue(
                lambda name = name: self.molecules[name].c_env,
                lambda value, name = name: setattr(self.molecules[name], 'c_env', value))


        self.cell_concs = DynamicValueDict(cell_concs_mapping)
        self.env_concs = DynamicValueDict(env_concs_mapping)

    def tissue_init(self, config_substances, p):
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

            if name not in p.ions_dict:
                # otherwise, get the name of the Molecule instance we're looking for as an alias:
                mol = self.molecules[name]

                mol.Dm = mol_dic['Dm']  # membrane diffusion coefficient [m2/s]
                mol.Do = mol_dic['Do']  # free diffusion constant in extra and intracellular spaces [m2/s]
                mol.Dgj = mol_dic.get('Dgj', 1.0e-16)  # effective diffusion coefficient of substance through GJ
                mol.z = mol_dic['z']  # charge (oxidation state)

                self.zmol[name] = mol.z
                self.Dmem[name] = mol.Dm

                # factors involving growth and decay (gad) in the cytoplasm
                gad = mol_dic.get('growth and decay', None)

                if gad != 'None' and gad is not None:

                    mol.simple_growth = True

                    mol.r_production = gad['production rate']
                    mol.r_decay = gad['decay rate']

                    mol.growth_profiles_list = gad['apply to']

                    mol.growth_activators_list = gad.get('activators', None)
                    mol.growth_activators_zone = gad.get('zone activators', None)
                    mol.growth_activators_Km = gad.get('Km activators', None)
                    mol.growth_activators_n = gad.get('n activators', None)

                    mol.growth_inhibitors_list = gad.get('inhibitors', None)
                    mol.growth_inhibitors_zone = gad.get('zone inhibitors', None)
                    mol.growth_inhibitors_Km = gad.get('Km inhibitors', None)
                    mol.growth_inhibitors_n = gad.get('n inhibitors', None)

                    # Fill in the blanks if zones aren't supplied:
                    mol.growth_activators_zone, mol.growth_inhibitors_zone = self.default_zones(
                                                                            mol.growth_activators_zone,
                                                                            mol.growth_inhibitors_zone,
                                                                            mol.growth_activators_list,
                                                                            mol.growth_inhibitors_list)

                    mol.growth_mod_function_cells = np.ones(self.cdl)

                else:
                    mol.simple_growth = False
                    mol.growth_profiles_list = None

        # balance charge in cells and env
        if p.substances_affect_charge and self.charge_has_been_balanced is False:
            # calculate the net charge in cells
            Qcells = 0
            Qenv = 0

            for name in self.molecules:

                mol = self.molecules[name]

                Qcells += mol.c_cello * mol.z
                Qenv += mol.c_envo * mol.z

            self.bal_charge(Qcells, 'cell', p)
            self.bal_charge(Qenv, 'env', p)

            self.charge_has_been_balanced = True

        # write substance growth and decay equations:
        self.write_growth_and_decay()

        # write passive electrodiffusion expressions for each ion and molecule:
        logs.log_info("Writing passive electrodiffusion equations...")
        for k, v in self.cell_concs.items():

            name = "'{}'".format(k)

            if self.zmol[k] != 0.0 and self.Dmem[k] != 0.0:

                self.ED_eval_strings[k] = "(self.Dmem[{}]/p.tm)*((p.F*self.vm)/(p.R*p.T))*self.zmol[{}]*" \
                                              "((self.cell_concs[{}] - " \
                                              "self.env_concs[{}]*" \
                                              "np.exp(-(self.zmol[{}]*p.F*self.vm)" \
                                              "/(p.R*p.T)))/(1-np.exp(-(self.zmol[{}]*" \
                                              "p.F*self.vm)/(p.R*p.T))))".format(name, name, name, name, name,
                                                                                     name)

            elif self.zmol[k] == 0.0 and self.Dmem[k] != 0.0:


                self.ED_eval_strings[k] = "(self.Dmem[{}]/p.tm)*(self.cell_concs[{}] - " \
                                          "self.env_concs[{}])".format(name, name, name)

            elif self.Dmem[k] == 0.0:

                self.ED_eval_strings[k] = "np.zeros(self.cdl)"

    def build_indices(self, p):

        # build indices-----------------------------------------------------------------------------------------
        self.molecule_index = {}

        for i, (mol_name, val) in enumerate(self.cell_concs.items()):
            self.molecule_index[mol_name] = i

        self.reaction_index = {}

        for i, rea_name in enumerate(self.reactions):
            self.reaction_index[rea_name] = i

        self.transporter_index = {}

        for i, trans_name in enumerate(self.transporters):
            self.transporter_index[trans_name] = i

        # build a global dictionary of target concentrations:--------------------------------------------------
        self.conc_handler = OrderedDict({})

        for mol_name, val in self.cell_concs.items():

            env_name = mol_name + '_env'

            self.conc_handler[mol_name] = self.cell_concs[mol_name].mean()
            self.conc_handler[env_name] = self.env_concs[mol_name].mean()

        # build a global dictionary of all reactions:---------------------------------------------------------
        self.react_handler = OrderedDict({})

        # define a string term to scale membrane-fluxes to cell-wide concentration change:
        div_string = "*self.div_cell"


        for mol_name, val in self.cell_concs.items():

            # add the transmembrane movement to the mix:
            ed_name = mol_name + '_ed'
            dmem_check = self.Dmem[mol_name]  # check the dmem value

            if dmem_check != 0.0:  # if it's not totally zero, then add passive flux expression to the list
                self.react_handler[ed_name] = self.ED_eval_strings[mol_name] + div_string

            if mol_name not in p.ions_dict:

                mol = self.molecules[mol_name]

                if mol.simple_growth is True:
                    rea_name = mol_name + '_growth'

                    self.react_handler[rea_name] = mol.gad_eval_string_growth

                    rea_name = mol_name + '_decay'

                    self.react_handler[rea_name] = mol.gad_eval_string_decay

        for rea_name in self.reactions:
            self.react_handler[rea_name] = self.reactions[rea_name].reaction_eval_string

        for trans_name in self.transporters:
            self.react_handler[trans_name] = self.transporters[trans_name].transporter_eval_string + div_string

        for chan_name in self.channels:

            self.react_handler[chan_name] = self.channels[chan_name].channel_eval_string + div_string

        self.react_handler_index = {}

        for i, rea_name in enumerate(self.react_handler):
            self.react_handler_index[rea_name] = i

        self.conc_handler_index = {}

        for i, conc_name in enumerate(self.conc_handler):
            self.conc_handler_index[conc_name] = i

    #FIXME: This method has been copy-and-pasted into numerous other submodules,
    #making life substantially harder. Duplicate copies include those in:
    #
    #* "chemistry/networks.py"
    #* "organelles/endo_retic.py".
    #
    #Ideally, a common abstract base class (ABC) shared between these three
    #submodules should handle common code like this. Diaphanous dolphins!
    def init_saving(
        self, p, plot_type='init', nested_folder_name='DynamicSystem'):

        if plot_type == 'sim':
            self.resultsPath = pathnames.join(p.sim_export_dirname, nested_folder_name)
            p.plot_type = 'sim'
        elif plot_type == 'init':
            self.resultsPath = pathnames.join(p.init_export_dirname, nested_folder_name)
            p.plot_type = 'init'

        dirs.make_unless_dir(self.resultsPath)
        self.imagePath = pathnames.join(self.resultsPath, 'fig_')

    def read_reactions(self, config_reactions, p):

        """
          Read in and initialize parameters for all user-defined reactions.

          config_options:  dictionary containing BETSE reaction template fields

          """

        logs.log_info("Reading reaction input data...")

        # Initialize a dict that keeps the Reaction objects for reactions occuring in the 'cell' zone:
        self.reactions = OrderedDict({})
        self.reactions_env = OrderedDict({})

        for q, react_dic in enumerate(config_reactions):

            # get each user-defined name-filed in the dictionary:
            name = str(react_dic['name'])

            zone = react_dic['reaction zone']

            if zone == 'cell':
                # add a Reaction object to MasterOfReactions reaction dictionary:
                self.reactions[name] = Reaction(p)
                obj = self.reactions[name]

            elif zone == 'env':

                self.reactions_env[name] = Reaction(p)
                obj = self.reactions_env[name]

            else:
                logs.log_warning("---------------------------------------------------------------")
                logs.log_warning("Reaction {} defined on inappropriate zone (or mit not enabled)".format(name))
                logs.log_warning("Reverting this reaction to the 'cell' zone!")
                logs.log_warning("---------------------------------------------------------------")
                # add a Reaction object to MasterOfReactions reaction dictionary:
                self.reactions[name] = Reaction(p)
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

            # Fill in the blanks if zones aren't supplied:
            obj.reaction_activators_zone, obj.reaction_inhibitors_zone = self.default_zones(
                obj.reaction_activators_zone,
                obj.reaction_inhibitors_zone,
                obj.reaction_activators_list,
                obj.reaction_inhibitors_list)

        for name in self.reactions:

            msg = "Including the cell-zone reaction: {}".format(name)
            logs.log_info(msg)

    def read_transporters(self, config_transporters, p):

        """
            Read in and initialize parameters for all user-defined transporters.

            config_options:  dictionary containing BETSE transporter template fields

            """

        logs.log_info("Reading transporter input data...")

        # Initialize a dict of Transporter objects:
        self.transporters = OrderedDict({})

        for q, trans_dic in enumerate(config_transporters):
            # get each user-defined name-filed in the dictionary:
            name = str(trans_dic['name'])

            # add a field to the MasterOfReactions corresponding to Transporter object with that name:

            self.transporters[name] = Transporter(p)

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

            # Fill in the blanks if zones aren't supplied:
            obj.transporter_activators_zone, obj.transporter_inhibitors_zone = self.default_zones(
                obj.transporter_activators_zone,
                obj.transporter_inhibitors_zone,
                obj.transporter_activators_list,
                obj.transporter_inhibitors_list)

        for name in self.transporters:

            msg = "Including the network transporter: {}".format(name)
            logs.log_info(msg)

    def read_channels(self, config_channels, p):

        logs.log_info("Reading channel input data...")

        # Initialize a dict of Channels:
        self.channels = OrderedDict({})

        for q, chan_dic in enumerate(config_channels):
            # get each user-defined name-filed in the dictionary:
            name = str(chan_dic['name'])

            # add a field to the MasterOfReactions corresponding to Channel object with that name:
            self.channels[name] = Channel(p)

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

            # Fill in the blanks if zones aren't supplied:
            obj.channel_activators_zone, obj.channel_inhibitors_zone = self.default_zones(
                obj.channel_activators_zone,
                obj.channel_inhibitors_zone,
                obj.channel_activators_list,
                obj.channel_inhibitors_list)

            # activator/inhibitor lists and associated data:
            a_list = obj.channel_activators_list
            Km_a_list = obj.channel_activators_Km
            n_a_list = obj.channel_activators_n
            zone_a = obj.channel_activators_zone

            i_list = obj.channel_inhibitors_list
            Km_i_list = obj.channel_inhibitors_Km
            n_i_list = obj.channel_inhibitors_n
            zone_i = obj.channel_inhibitors_zone

            all_alpha, _, _ = self.get_influencers(a_list, Km_a_list, n_a_list, i_list, Km_i_list, n_i_list,
                                                    reaction_zone='cell', zone_tags_a=zone_a,
                                                    zone_tags_i=zone_i, in_mem_tag=True)


            obj.alpha_eval_string = "(" + all_alpha + ")"

            if obj.channel_class == 'Na':
                obj.ion_strings = ['Na']

            elif obj.channel_class == 'K':
                obj.ion_strings = ['K']

            elif obj.channel_class == 'Ca':
                obj.ion_strings = ['Ca']

            elif obj.channel_class == 'Cl':
                obj.ion_strings = ['Cl']

            elif obj.channel_class == 'Cat':
                obj.ion_strings = ['Na']

            elif obj.channel_class == 'Fun':
                obj.ion_strings = ['Na']

            obj.init_channel(obj.ion_strings, obj.channel_type, obj.channelMax, obj.alpha_eval_string, p)

    def write_growth_and_decay(self):

        logs.log_info("Writing substance growth/decay equations...")

        for mol_name in self.molecules:


            if self.molecules[mol_name].simple_growth is True:


                cc = "(self.cell_concs['{}'])".format(mol_name)

                # get activators and inhibitors and associated data:
                a_list = self.molecules[mol_name].growth_activators_list
                Km_a_list = self.molecules[mol_name].growth_activators_Km
                n_a_list = self.molecules[mol_name].growth_activators_n
                zone_a = self.molecules[mol_name].growth_activators_zone

                i_list = self.molecules[mol_name].growth_inhibitors_list
                Km_i_list = self.molecules[mol_name].growth_inhibitors_Km
                n_i_list = self.molecules[mol_name].growth_inhibitors_n
                zone_i = self.molecules[mol_name].growth_inhibitors_zone

                all_alpha, alpha_tex, gad_tex_var_list = \
                             self.get_influencers(a_list, Km_a_list, n_a_list,
                             i_list, Km_i_list, n_i_list,reaction_zone='cell',
                             zone_tags_a=zone_a, zone_tags_i=zone_i, in_mem_tag=False)

                # define remaining portion of substance's growth and decay expression:

                mod_funk = "(self.molecules['{}'].growth_mod_function_cells)".format(mol_name)
                r_prod =  "(self.molecules['{}'].r_production)".format(mol_name)
                r_decay = "(self.molecules['{}'].r_decay)".format(mol_name)

                gad_eval_string = mod_funk + "*" + r_prod + "*" + all_alpha + "-" + \
                                  r_decay + "*" + cc


                #------------------------------------------------------------------------------

                gad_eval_string_growth = all_alpha
                gad_eval_string_decay = r_decay + "*" + cc

                self.molecules[mol_name].gad_eval_string = gad_eval_string

                # add in the separate growth decay eval strings for case of optimization:
                self.molecules[mol_name].gad_eval_string_growth = gad_eval_string_growth
                self.molecules[mol_name].gad_eval_string_decay = gad_eval_string_decay

            else:

                self.molecules[mol_name].gad_eval_string = "np.zeros(self.cdl)"
                self.molecules[mol_name].gad_eval_string_growth = "np.zeros(self.cdl)"
                self.molecules[mol_name].gad_eval_string_decay = "np.zeros(self.cdl)"

                self.molecules[mol_name].gad_tex_string = ""
                self.molecules[mol_name].gad_tex_vars = ""

    def write_reactions(self):
        """
        Reactions are now constructed during the init as strings that are evaluated in eval calls in each time-step.
        This function constructs the evaluation strings for each reaction, given the metadata stored
        in each reaction object (e.g. lists of reactants, products, etc).

        """

        logs.log_info("Writing reaction equations for cell zone...")

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
            zone_a = self.reactions[reaction_name].reaction_activators_zone

            i_list = self.reactions[reaction_name].reaction_inhibitors_list
            Km_i_list = self.reactions[reaction_name].reaction_inhibitors_Km
            n_i_list = self.reactions[reaction_name].reaction_inhibitors_n
            zone_i = self.reactions[reaction_name].reaction_inhibitors_zone

            r_zone = self.reactions[reaction_name].reaction_zone

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

                # tex_name = name

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
                Keqm = "(np.exp(-self.reactions['{}'].delta_Go / (p.R * p.T)))".format(reaction_name)


            else:

                Q = "0"
                backward_coeff = "0"
                Keqm = "1"


            # Put it all together into a final reaction string (+max rate):

            reversed_term = "(" + Q + "/" + Keqm + ")"


            all_alpha, alpha_tex, rea_tex_var_list = self.get_influencers(a_list, Km_a_list,
                                                                    n_a_list,
                                                                    i_list, Km_i_list, n_i_list,
                                                                    reaction_zone= r_zone, zone_tags_a=zone_a,
                                                                    zone_tags_i=zone_i, in_mem_tag= False)

            vmax = "self.reactions['{}'].vmax".format(reaction_name)

            reaction_eval_string = vmax + "*" + all_alpha + "*" + "(" + \
                                   forward_coeff + "-" + "(" + reversed_term + "*" + backward_coeff + ")" + ")"


            # add the composite string describing the reaction math to a new field:
            self.reactions[reaction_name].reaction_eval_string = reaction_eval_string

    def write_reactions_mit(self):
        """
        Reactions are now constructed during the init as strings that are evaluated in eval calls in each time-step.
        This function constructs the evaluation strings for each reaction, given the metadata stored
        in each reaction object (e.g. lists of reactants, products, etc).

        """

        logs.log_info("Writing reaction equations for mit zone...")

        for reaction_name in self.reactions_mit:

        # define aliases for convenience:

            reactant_names = self.reactions_mit[reaction_name].reactants_list
            reactant_coeff = self.reactions_mit[reaction_name].reactants_coeff
            reactant_Km = self.reactions_mit[reaction_name].Km_reactants_list

            product_names = self.reactions_mit[reaction_name].products_list
            product_coeff = self.reactions_mit[reaction_name].products_coeff
            product_Km = self.reactions_mit[reaction_name].Km_products_list

            # activator/inhibitor lists and associated data:
            a_list = self.reactions_mit[reaction_name].reaction_activators_list
            Km_a_list = self.reactions_mit[reaction_name].reaction_activators_Km
            n_a_list = self.reactions_mit[reaction_name].reaction_activators_n
            zone_a = self.reactions_mit[reaction_name].reaction_activators_zone

            i_list = self.reactions_mit[reaction_name].reaction_inhibitors_list
            Km_i_list = self.reactions_mit[reaction_name].reaction_inhibitors_Km
            n_i_list = self.reactions_mit[reaction_name].reaction_inhibitors_n
            zone_i = self.reactions_mit[reaction_name].reaction_inhibitors_zone

            r_zone = self.reactions_mit[reaction_name].reaction_zone

            # first calculate a reaction coefficient Q, as a string expression
            numo_string_Q = "("
            denomo_string_Q = "("


            for i, (name, coeff) in enumerate(zip(reactant_names, reactant_coeff)):


                denomo_string_Q += "(self.mit_concs['{}']".format(name)

                denomo_string_Q += "**{})".format(coeff)


                if i < len(reactant_names) - 1:

                    denomo_string_Q += "*"

                else:

                    denomo_string_Q += ")"


            for i, (name, coeff) in enumerate(zip(product_names, product_coeff)):

                # tex_name = name + "_{mit}"
                numo_string_Q += "(self.mit_concs['{}']".format(name)

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

                numo_string_r = "((self.mit_concs['{}']/{})**{})".format(name, Km, n)
                denomo_string_r = "(1 + (self.mit_concs['{}']/{})**{})".format(name, Km, n)

                term = "(" + numo_string_r + "/" + denomo_string_r + ")"

                forward_coeff += term


                if i < len(reactant_names) - 1:

                    forward_coeff += "*"
                else:

                    forward_coeff += ")"

            for i, (name, n, Km) in enumerate(zip(product_names, product_coeff, product_Km)):

                # tex_name = name + '_{mit}'

                numo_string_p = "((self.mit_concs['{}']/{})**{})".format(name, Km, n)
                denomo_string_p = "(1 + (self.mit_concs['{}']/{})**{})".format(name, Km, n)

                term = "(" + numo_string_p + "/" + denomo_string_p + ")"


                if i < len(product_names) - 1:

                    backward_coeff += "*"

                else:

                    backward_coeff += ")"

            # if reaction is reversible deal calculate an equilibrium constant:
            if self.reactions_mit[reaction_name].delta_Go is not None:

                # define the reaction equilibrium coefficient expression:
                Keqm = "(np.exp(-self.reactions_mit['{}'].delta_Go / (p.R * p.T)))".format(reaction_name)

            else:

                Q = "0"
                backward_coeff = "0"
                Keqm = "1"

            # Put it all together into a final reaction string (+max rate):

            reversed_term = "(" + Q + "/" + Keqm + ")"

            all_alpha, alpha_tex, rea_tex_var_list = self.get_influencers(a_list, Km_a_list,
                                                                    n_a_list,
                                                                    i_list, Km_i_list, n_i_list,
                                                                    reaction_zone= r_zone, zone_tags_a=zone_a,
                                                                    zone_tags_i=zone_i, in_mem_tag=False)

            vmax = "self.reactions_mit['{}'].vmax".format(reaction_name)

            reaction_eval_string = vmax + "*" + all_alpha + "*" + "(" + \
                                   forward_coeff + "-" + "(" + reversed_term + "*" + backward_coeff + ")" + ")"


            # add the composite string describing the reaction math to a new field:
            self.reactions_mit[reaction_name].reaction_eval_string = reaction_eval_string

    def write_reactions_env(self):

        """
        Reactions are now constructed during the init as strings that are evaluated in eval calls in each time-step.
        This function constructs the evaluation strings for each reaction, given the metadata stored
        in each reaction object (e.g. lists of reactants, products, etc).

        """

        logs.log_info("Writing reaction equations for env zone...")

        for reaction_name in self.reactions_env:


            # define aliases for convenience:

            reactant_names = self.reactions_env[reaction_name].reactants_list
            reactant_coeff = self.reactions_env[reaction_name].reactants_coeff
            reactant_Km = self.reactions_env[reaction_name].Km_reactants_list

            product_names = self.reactions_env[reaction_name].products_list
            product_coeff = self.reactions_env[reaction_name].products_coeff
            product_Km = self.reactions_env[reaction_name].Km_products_list

            # activator/inhibitor lists and associated data:
            a_list = self.reactions_env[reaction_name].reaction_activators_list
            Km_a_list = self.reactions_env[reaction_name].reaction_activators_Km
            n_a_list = self.reactions_env[reaction_name].reaction_activators_n
            zone_a = self.reactions_env[reaction_name].reaction_activators_zone

            i_list = self.reactions_env[reaction_name].reaction_inhibitors_list
            Km_i_list = self.reactions_env[reaction_name].reaction_inhibitors_Km
            n_i_list = self.reactions_env[reaction_name].reaction_inhibitors_n
            zone_i = self.reactions_env[reaction_name].reaction_inhibitors_zone

            r_zone = self.reactions_env[reaction_name].reaction_zone

            # first calculate a reaction coefficient Q, as a string expression
            numo_string_Q = "("
            denomo_string_Q = "("

            for i, (name, coeff) in enumerate(zip(reactant_names, reactant_coeff)):

                denomo_string_Q += "(self.env_concs['{}']".format(name)

                denomo_string_Q += "**{})".format(coeff)


                if i < len(reactant_names) - 1:

                    denomo_string_Q += "*"

                else:

                    denomo_string_Q += ")"

            for i, (name, coeff) in enumerate(zip(product_names, product_coeff)):

                numo_string_Q += "(self.env_concs['{}']".format(name)

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

                numo_string_r = "((self.env_concs['{}']/{})**{})".format(name, Km, n)
                denomo_string_r = "(1 + (self.env_concs['{}']/{})**{})".format(name, Km, n)

                term = "(" + numo_string_r + "/" + denomo_string_r + ")"

                forward_coeff += term


                if i < len(reactant_names) - 1:

                    forward_coeff += "*"

                else:

                    forward_coeff += ")"

            for i, (name, n, Km) in enumerate(zip(product_names, product_coeff, product_Km)):


                numo_string_p = "((self.env_concs['{}']/{})**{})".format(name, Km, n)
                denomo_string_p = "(1 + (self.env_concs['{}']/{})**{})".format(name, Km, n)

                term = "(" + numo_string_p + "/" + denomo_string_p + ")"


                backward_coeff += term


                if i < len(product_names) - 1:

                    backward_coeff += "*"

                else:

                    backward_coeff += ")"

            # if reaction is reversible deal calculate an equilibrium constant:
            if self.reactions_env[reaction_name].delta_Go is not None:

                # define the reaction equilibrium coefficient expression:
                Keqm = "(np.exp(-self.reactions_env['{}'].delta_Go / (p.R * p.T)))".format(reaction_name)

            else:

                Q = "0"
                backward_coeff = "0"
                Keqm = "1"


            # Put it all together into a final reaction string (+max rate):

            reversed_term = "(" + Q + "/" + Keqm + ")"


            all_alpha, alpha_tex, rea_tex_var_list = self.get_influencers(a_list, Km_a_list,
                                                                          n_a_list,
                                                                          i_list, Km_i_list, n_i_list,
                                                                          reaction_zone=r_zone, zone_tags_a=zone_a,
                                                                          zone_tags_i=zone_i, in_mem_tag=False)

            vmax = "self.reactions_env['{}'].vmax".format(reaction_name)

            reaction_eval_string = vmax + "*" + all_alpha + "*" + "(" + \
                                   forward_coeff + "-" + "(" + reversed_term + "*" + backward_coeff + ")" + ")"


            # add the composite string describing the reaction math to a new field:
            self.reactions_env[reaction_name].reaction_eval_string = reaction_eval_string

    def write_transporters(self, p):
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

                vmem = "self.vm"   # get the transmembrane voltage for this category

                div_term = "self.div_cell"

                in_delta_term_react = "-self.transporters['{}'].flux*{}".format(transp_name, div_term)
                in_delta_term_prod = "self.transporters['{}'].flux*{}".format(transp_name, div_term)


                out_delta_term_react = "-self.transporters['{}'].flux*{}".format(transp_name, div_term)

                out_delta_term_prod = "self.transporters['{}'].flux*{}".format(transp_name, div_term)

                all_alpha, _, _ = self.get_influencers(a_list, Km_a_list,
                                                                n_a_list, i_list,
                                                                Km_i_list, n_i_list,
                                                                reaction_zone='cell', in_mem_tag=False)


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

                if tag == 'mem_concs':

                    tag2 = 'cell_concs'

                    denomo_string_Q += "(self.{}['{}']".format(tag2, name)

                elif tag == 'mit_concs':
                    denomo_string_Q += "(self.{}['{}']".format(tag, name)

                elif tag == 'cell_concs':
                    denomo_string_Q += "(self.{}['{}']".format(tag, name)

                # get the concentration from the environment, mapped to respective membranes:
                elif tag == 'env_concs':

                    denomo_string_Q += "(self.{}['{}']".format(tag, name)

                else:

                    raise BetseSimConfException("Transporter tag not properly defined!")

                denomo_string_Q += "**{})".format(coeff)

                if i < len(reactant_names) - 1:

                    denomo_string_Q += "*"

                else:

                    denomo_string_Q += ")"

            for i, (name, coeff, tag) in enumerate(zip(product_names, product_coeff, prod_transfer_tag)):

                if tag == 'mem_concs':

                    tag2 = 'cell_concs'

                    numo_string_Q += "(self.{}['{}']".format(tag2, name)

                elif tag == 'mit_concs':
                    numo_string_Q += "(self.{}['{}']".format(tag, name)

                elif tag == 'cell_concs':
                    numo_string_Q += "(self.{}['{}']".format(tag, name)


                # get the concentration from the environment mapped to the respective membranes:
                elif tag == 'env_concs':

                    numo_string_Q += "(self.{}['{}']".format(tag, name)

                else:

                    raise BetseSimConfException("Transporter tag not properly defined!")

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

                if tag == 'mem_concs':

                    tag2 = 'cell_concs'

                    numo_string_r = "((self.{}['{}']/{})**{})".format(tag2, name, Km, n)
                    denomo_string_r = "(1 + (self.{}['{}']/{})**{})".format(tag2, name, Km, n)

                elif tag == 'mit_concs':

                    numo_string_r = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    denomo_string_r = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)

                elif tag == 'cell_concs':

                    numo_string_r = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    denomo_string_r = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)

                elif tag == 'env_concs':

                    numo_string_r = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    denomo_string_r = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)

                else:

                    raise BetseSimConfException("Transporter tag not properly defined!")

                term = "(" + numo_string_r + "/" + denomo_string_r + ")"

                forward_coeff += term

                if i < len(reactant_names) - 1:

                    forward_coeff += "*"

                else:

                    forward_coeff += ")"

            for i, (name, n, Km, tag) in enumerate(zip(product_names, product_coeff, product_Km, prod_transfer_tag)):

                if tag == 'mem_concs':

                    tag2 = 'cell_concs'

                    numo_string_p = "((self.{}['{}']/{})**{})".format(tag2, name, Km, n)
                    denomo_string_p = "(1 + (self.{}['{}']/{})**{})".format(tag2, name, Km, n)


                elif tag == 'mit_concs':

                    numo_string_p = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    denomo_string_p = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)

                elif tag == 'cell_concs':

                    numo_string_p = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    denomo_string_p = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)

                elif tag == 'env_concs':

                    numo_string_p = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    denomo_string_p = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)


                else:

                    raise BetseSimConfException("Transporter tag not properly defined!")


                term = "(" + numo_string_p + "/" + denomo_string_p + ")"

                backward_coeff += term

                if i < len(product_names) - 1:

                    backward_coeff += "*"

                else:

                    backward_coeff += ")"

            # if reaction is reversible, calculate an equilibrium constant:
            if self.transporters[transp_name].delta_Go is not None:

                # define the reaction equilibrium coefficient expression:
                Keqm = "(np.exp(-(self.transporters['{}'].delta_Go + {})/ (p.R * p.T)))".format(transp_name,
                                                                                                    echem_string)

            else:

                Q = "0"
                backward_coeff = "0"
                Keqm = "1"

            # Put it all together into a final reaction string (+max rate):
            reversed_term = "(" + Q + "/" + Keqm + ")"

            vmax = "self.transporters['{}'].vmax".format(transp_name)


            # write the final LaTeX expressions:
            if self.transporters[transp_name].delta_Go is not None:

                # calculate the evaluation string expression for the transporter:
                transporter_eval_string = vmax + "*" + all_alpha + "*" + "(" + \
                                          forward_coeff + "-"  + reversed_term + "*" + backward_coeff  + ")"


            else:

                # calculate the evaluation string expression for the transporter:
                transporter_eval_string = vmax + "*" + all_alpha + "*" + forward_coeff


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

        logs.log_info("Writing reaction network matrix for cell zone...")
        n_reacts = len(self.molecules) + len(self.reactions)  # (this is + len(molecules) due to added gad reactions
        n_mols = len(self.cell_concs)

        n_subs = len(self.molecules)

        # initialize the network's reaction matrix:
        self.reaction_matrix = np.zeros((n_mols, n_reacts))

        molecule_keys = list(self.cell_concs.keys())

        for i, name in enumerate(self.molecules):

            jj = molecule_keys.index(name)
            # add in terms referencing the self-growth and decay reaction for each substance
            self.reaction_matrix[jj, i] += 1

        for jo, reaction_name in enumerate(self.reactions):

            j = jo + n_subs

            for react_name, coeff in zip(self.reactions[reaction_name].reactants_list,
                self.reactions[reaction_name].reactants_coeff):
                i = molecule_keys.index(react_name)

                self.reaction_matrix[i, j] += -coeff

            for prod_name, coeff in zip(self.reactions[reaction_name].products_list,
                self.reactions[reaction_name].products_coeff):
                i = molecule_keys.index(prod_name)

                self.reaction_matrix[i, j] += coeff

    def create_reaction_matrix_env(self):

        """
        Creates a matrix representation of the system of coupled
        differential equations underlying the reaction network for environment.

        """

        logs.log_info("Writing reaction network matrix for env zone...")
        n_reacts = len(self.reactions_env)
        n_mols = len(self.env_concs)

        # initialize the network's reaction matrix:
        self.reaction_matrix_env = np.zeros((n_mols, n_reacts))

        molecule_keys = list(self.env_concs.keys())

        for j, reaction_name in enumerate(self.reactions_env):

            for react_name, coeff in zip(self.reactions_env[reaction_name].reactants_list,
                                         self.reactions_env[reaction_name].reactants_coeff):
                i = molecule_keys.index(react_name)

                self.reaction_matrix_env[i, j] += -coeff

            for prod_name, coeff in zip(self.reactions_env[reaction_name].products_list,
                                        self.reactions_env[reaction_name].products_coeff):
                i = molecule_keys.index(prod_name)

                self.reaction_matrix_env[i, j] += coeff

    def export_eval_strings(self, p):

        # Absolute path of the text file to write this solution to:
        saveData = pathnames.join(
            self.resultsPath, 'NetworkModelEvalStrings.csv')

        with open(saveData, 'w', newline='') as csvfile:
            eqwriter = csv.writer(csvfile, delimiter='\t',
                quotechar='|', quoting=csv.QUOTE_NONE)

            for mol in self.molecules:

                if self.molecules[mol].simple_growth is True:
                    eq_var = self.molecules[mol].gad_eval_string

                    eqwriter.writerow([mol, eq_var])

            for rea in self.reactions:
                eq_var = self.reactions[rea].reaction_eval_string

                eqwriter.writerow([rea, eq_var])

            for rea in self.reactions_mit:
                eq_var = self.reactions_mit[rea].reaction_eval_string
                eqwriter.writerow([rea, eq_var])

            for trans in self.transporters:
                eq_var = self.transporters[trans].transporter_eval_string

                eqwriter.writerow([trans, eq_var])

    def build_reaction_network(self, p):
        """
        Uses pydot to create and return a directed graph representing only growth/decay and chemical reactions,
        omitting any channels and activation/inhibition relationships considered in this reaction network object.

        """

        # If PyDot is unimportable, raise an exception.
        libs.die_unless_runtime_optional('pydot')

        # Delay import of pydot in case the user doesn't have it and needs to
        # turn this functionality off.
        import pydot

        # alpha value to decrease saturation of graph node colors
        # alpha_val = 0.5

        # define some basic colormap scaling properties for the dataset:
        vals = np.asarray([v.c_cells.mean() for (c, v) in self.molecules.items()])
        minc = vals.min()
        maxc = vals.max()
        colors.Normalize(vmin=minc, vmax=maxc)

        # create a graph object
        graphicus_maximus = pydot.Dot(
            graph_type='digraph', concentrate=True, nodesep=0.1, ranksep=0.3,
            overlap='compress', splines=True)

        # add each substance as a node in the graph:
        for i, (name, val) in enumerate(self.cell_concs.items()):

            if name not in p.ions_dict:

                mol = self.molecules[name]

                nde = pydot.Node(name, shape = self.conc_shape)
                graphicus_maximus.add_node(nde)

                if mol.simple_growth:

                    # add node & edge for growth reaction component:
                    rea_name = name + '_growth'
                    rea_node = pydot.Node(rea_name, style = 'filled', shape = self.reaction_shape)
                    graphicus_maximus.add_node(rea_node)

                    # if the substance has autocatalytic growth capacity add the edge in:
                    graphicus_maximus.add_edge(pydot.Edge(rea_name, name, arrowhead='normal', coeff = 1.0,
                                                          ))

                    # add node & edge for decay reaction component:
                    rea_name = name + '_decay'
                    rea_node = pydot.Node(rea_name, style = 'filled', shape = self.reaction_shape)
                    graphicus_maximus.add_node(rea_node)

                    # if the substance has autocatalytic growth capacity add the edge in:
                    graphicus_maximus.add_edge(pydot.Edge(name, rea_name, arrowhead='normal', coeff =1.0,
                                                          ))

            else:  # add the ion as a node

                nde = pydot.Node(name, shape = self.conc_shape)
                graphicus_maximus.add_node(nde)


            # add expression nodes for electrodiffusion/diffusion:
            # check the Dmem value:
            dmem_check = self.Dmem[name]

            if dmem_check != 0.0:
                react_label = name + '_ed'
                nde = pydot.Node(react_label, style='filled', shape=self.ed_shape)
                graphicus_maximus.add_node(nde)

                name_env = name + '_env'
                nde = pydot.Node(name_env, shape=self.conc_shape)
                graphicus_maximus.add_node(nde)

                graphicus_maximus.add_edge(pydot.Edge(name, react_label, arrowhead='normal', coeff=1.0,
                                                      ))
                graphicus_maximus.add_edge(pydot.Edge(react_label, name_env, arrowhead='normal', coeff=1.0,
                                                      ))

        # if there are any reactions in the cytosol, add them to the graph
        if len(self.reactions) > 0:

            for i, name in enumerate(self.reactions):
                nde = pydot.Node(name, style='filled', shape=self.reaction_shape)
                graphicus_maximus.add_node(nde)

        # if there are any reactions, plot their edges on the graph--------------------------------------------------
        if len(self.reactions) > 0:

            for i, name in enumerate(self.reactions):

                rea = self.reactions[name]

                for i, react_name in enumerate(rea.reactants_list):

                    rea_coeff = rea.reactants_coeff[i]
                    graphicus_maximus.add_edge(pydot.Edge(react_name, name, arrowhead='normal', coeff = rea_coeff,
                                                          ))

                for j, prod_name in enumerate(rea.products_list):

                    prod_coeff = rea.products_coeff[j]
                    graphicus_maximus.add_edge(pydot.Edge(name, prod_name, arrowhead='normal', coeff = prod_coeff,
                                                          ))


        # if there are any transporters, plot them on the graph:
        if len(self.transporters) > 0:

            for name in self.transporters:
                nde = pydot.Node(name, style='filled', shape=self.transporter_shape)
                graphicus_maximus.add_node(nde)

            for name in self.transporters:

                trans = self.transporters[name]

                for i, (react_name, tag) in enumerate(zip(trans.reactants_list, trans.react_transport_tag)):


                    rea_coeff = trans.reactants_coeff[i]

                    if tag == 'cell_concs' or tag == 'mem_concs':


                        graphicus_maximus.add_edge(pydot.Edge(react_name, name, arrowhead='normal',coeff = rea_coeff,
                                                              ))

                    else:

                        if tag == 'env_concs':

                            react_name += '_env'

                        elif tag == 'mit_concs':

                            react_name += '_mit'

                        nde = pydot.Node(react_name, shape = self.conc_shape)
                        graphicus_maximus.add_node(nde)

                        graphicus_maximus.add_edge(pydot.Edge(react_name, name, arrowhead='normal',coeff=rea_coeff,
                                                              ))

                for j, (prod_name, tag) in enumerate(zip(trans.products_list, trans.prod_transport_tag)):

                    prod_coeff = trans.products_coeff[j]

                    if tag == 'cell_concs' or tag == 'mem_concs':

                        graphicus_maximus.add_edge(pydot.Edge(name, prod_name, arrowhead='normal', coeff= prod_coeff,
                                                              ))

                    else:

                        if tag == 'env_concs':

                            prod_name += '_env'

                        elif tag == 'mit_concs':

                            prod_name += '_mit'

                        nde = pydot.Node(prod_name, shape = self.conc_shape)
                        graphicus_maximus.add_node(nde)

                        graphicus_maximus.add_edge(pydot.Edge(name, prod_name, arrowhead='normal', coeff = prod_coeff,
                                                              ))

        # if there are any channels, plot them on the graph:
        if len(self.channels) > 0:

            for ch_k, ch_v in self.channels.items():

                nde = pydot.Node(ch_k, style='filled', shape=self.channel_shape)
                graphicus_maximus.add_node(nde)

                ch_type = ch_v.channel_class

                if ch_type == 'K':

                    graphicus_maximus.add_edge(pydot.Edge('K', ch_k, arrowhead='normal', coeff=1.0,
                                                          ))
                    graphicus_maximus.add_edge(pydot.Edge(ch_k, 'K_env', arrowhead='normal', coeff=1.0,
                                                          ))

                elif ch_type == 'Na':

                    graphicus_maximus.add_edge(pydot.Edge('Na', ch_k, arrowhead='normal', coeff=1.0,
                                                          ))
                    graphicus_maximus.add_edge(pydot.Edge(ch_k, 'Na_env', arrowhead='normal', coeff=1.0,
                                                          ))

                elif ch_type == 'Cl':

                    graphicus_maximus.add_edge(pydot.Edge('Cl', ch_k, arrowhead='normal', coeff=1.0,
                                                          ))
                    graphicus_maximus.add_edge(pydot.Edge(ch_k, 'Cl_env', arrowhead='normal', coeff=1.0,
                                                          ))

                elif ch_type == 'Fun':

                    graphicus_maximus.add_edge(pydot.Edge('K', ch_k, arrowhead='normal', coeff=1.0,
                                                          ))
                    graphicus_maximus.add_edge(pydot.Edge(ch_k, 'K_env', arrowhead='normal', coeff=1.0,
                                                          ))

                    graphicus_maximus.add_edge(pydot.Edge('Na', ch_k, arrowhead='normal', coeff=0.2,
                                                          ))
                    graphicus_maximus.add_edge(pydot.Edge(ch_k, 'Na_env', arrowhead='normal', coeff=0.2,
                                                          ))

                elif ch_type == 'Cat':

                    graphicus_maximus.add_edge(pydot.Edge('K', ch_k, arrowhead='normal', coeff=1.0,
                                                          ))
                    graphicus_maximus.add_edge(pydot.Edge(ch_k, 'K_env', arrowhead='normal', coeff=1.0,
                                                          ))

                    graphicus_maximus.add_edge(pydot.Edge('Na', ch_k, arrowhead='normal', coeff=1.0,
                                                          ))
                    graphicus_maximus.add_edge(pydot.Edge(ch_k, 'Na_env', arrowhead='normal', coeff=1.0,
                                                          ))

                elif ch_type == 'Ca':

                    graphicus_maximus.add_edge(pydot.Edge('Ca', ch_k, arrowhead='normal', coeff=1.0,
                                                          ))
                    graphicus_maximus.add_edge(pydot.Edge(ch_k, 'Ca_env', arrowhead='normal', coeff=1.0,
                                                          ))


        return graphicus_maximus

    def get_influencers(self, a_list, Km_a_list, n_a_list, i_list, Km_i_list,
        n_i_list, reaction_zone='cell', zone_tags_a = None, zone_tags_i = None, in_mem_tag = False):

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
        reaction_zone:   Where the reaction/transporter/growth takes place ('cell', 'mit', 'mem')
        zone_tags_a:     Place where the activator concentration works from ('cell' or 'env'). See *.
        zone_tags_i:     Place where the inhibitor concentration works from ('cell' or 'env'). See *.
        in_mem_tag:      Flag to specify that the activator/inhibitors are acting on membrane-bound channels or
                         transporters (in which case they're voltage sensitive due to interactions with Vmem).

        * zone tags only apply to reactions occurring in 'cell' or 'mem' regions.

        Returns
        --------
        activator_alpha, inhibitor_alpha     string expressions representing the activator and inhibitor coeffs

        """

        # if zone tag lists aren't supplied or are invalid entries, create defaults:
        zone_tags_a, zone_tags_i = self.default_zones(zone_tags_a, zone_tags_i, a_list, i_list)

        # initialize string expressions for product of all activators and inhibitors
        activator_alpha = "("
        inhibitor_alpha = "("

        # initialize lists for special "direct terms" that are the concentrations unmodified:

        direct_terms_list = []

        direct_term_alpha = "("

        # set the appropriate value for the absence of an activator or inhibitor expression:
        if reaction_zone == 'cell':
            dl = "np.ones(self.cdl))"

        elif reaction_zone == 'mem':
            dl = "np.ones(self.cdl))"

        elif reaction_zone == 'mit':
            dl = 'np.ones(self.cdl))'

        elif reaction_zone == 'env':
            dl = "np.ones(self.cdl))"

        else:
            raise BetseSimConfException("Reaction zone requested not an existing option!")


        if a_list is not None and a_list != 'None' and len(a_list) > 0:

            # initialize an empty list that will contain independently acting activator terms:
            independent_terms = []

            # Begin with construction of the activator effects term:
            for i, (name, Kmo, n, zone_tag) in enumerate(zip(a_list, Km_a_list, n_a_list, zone_tags_a)):

                # check and see if the user requested voltage sensitivity using a '*' prefix to the name:
                if name.endswith('*'):
                    # capture the original name that will work in dictionaries:
                    name = name[:-1]
                    vs_request = True
                else:
                    vs_request = False

                independence_tag = False  # force independence and direct tags to be false by default
                direct_tag = False

                # check and see if name ends with a bang (meaning it's an independently acting activator that
                # is added to, rather than multiplied, in the main expression)
                if name.endswith('!'):
                    # capture the original name that will work in dictionaries:
                    name = name[:-1]
                    # set the independence tag to convey bang-effect
                    independence_tag = True

                elif name.endswith('&'):
                    name = name[:-1]
                    direct_tag = True

                # First deal with the fact that Km may be voltage sensitive if the substance is charged:
                # get the charge value for the substance:
                za = self.molecules[name].z

                if in_mem_tag is True and za != 0.0 and vs_request is True:

                    msg = "Vmem-sensitive gating used for activator: "+ name

                    logs.log_info(msg)

                    Kmc = '({}*np.exp(-(self.vm*{}*p.F)/(p.R*p.T)))'.format(Kmo, za)
                    Kme = '({}*np.exp((self.vm*{}*p.F)/(p.R*p.T)))'.format(Kmo, za)

                else:
                    Kmc = '{}'.format(Kmo)
                    Kme = '{}'.format(Kmo)

                # initialize string expressions for activator and inhibitor, respectively
                numo_string_a = ""
                denomo_string_a = ""

                direct_string_a = ""

                if reaction_zone == 'cell':

                    if zone_tag == 'cell':

                        numo_string_a += "((self.cell_concs['{}']/{})**{})".format(name, Kmc, n)
                        denomo_string_a += "(1 + (self.cell_concs['{}']/{})**{})".format(name, Kmc, n)

                        direct_string_a += "self.cell_concs['{}']".format(name)


                    elif zone_tag == 'env':

                        numo_string_a += "((self.env_concs['{}']/{})**{})".format(name, Kme, n)
                        denomo_string_a += "(1 + (self.env_concs['{}']/{})**{})".format(name, Kme, n)

                        direct_string_a += "self.env_concs['{}']".format(name)

                elif reaction_zone == 'mem':

                    if zone_tag == 'cell':

                        numo_string_a += "((self.cell_concs['{}']/{})**{})".format(name, Kmc, n)
                        denomo_string_a += "(1 + (self.cell_concs['{}']/{})**{})".format(name, Kmc, n)

                        direct_string_a += "self.cell_concs['{}']".format(name)

                    elif zone_tag == 'env':


                        numo_string_a += "((self.env_concs['{}']/{})**{})".format(name, Kme, n)
                        denomo_string_a += "(1 + (self.env_concs['{}']/{})**{})".format(name, Kme, n)

                        direct_string_a += "self.env_concs['{}']".format(name)

                elif reaction_zone == 'mit':

                    numo_string_a += "((self.mit_concs['{}']/{})**{})".format(name, Kmo, n)
                    denomo_string_a += "(1 + (self.mit_concs['{}']/{})**{})".format(name, Kmo, n)

                    direct_string_a += "self.mit_concs['{}']".format(name)

                elif reaction_zone == 'env':

                    if zone_tag == 'cell':

                        numo_string_a += "1"

                        denomo_string_a += "1"

                        direct_string_a += "1"

                    elif zone_tag == 'env':

                        numo_string_a += "((self.env_concs['{}']/{})**{})".format(name, Kme, n)
                        denomo_string_a += "(1 + (self.env_concs['{}']/{})**{})".format(name, Kme, n)

                        direct_string_a += "self.env_concs['{}']".format(name)

                else:
                    raise BetseSimConfException("You've asked for a reaction zone"
                                                   "that doesn't exist. Ensure all "
                                                   "reaction and transporter zones are 'cell' or 'env'.")


                term = "(" + numo_string_a + "/" + denomo_string_a + ")"


                #---Update the final activator alpha term.
                if independence_tag is False and direct_tag is False:
                    # include the new activator term as a multiplier:
                    activator_alpha += term


                elif independence_tag is True:
                    # include the term in the list; we will add them to the final coefficient at the end:
                    independent_terms.append(term)


                elif direct_tag is True:
                    # include the direct term (no modifications) to the string as an activator
                    direct_terms_list.append(direct_string_a)



                if i < len(a_list) - 1 and direct_tag is False and independence_tag is False:

                    activator_alpha += "*"



                elif i >= len(a_list) -1: # if we've reached the end of the activators list:

                    # add any independently acting terms from the list:
                    for indt in independent_terms:

                        if activator_alpha == "(": # if nothing has been added to the activator alpha string
                            activator_alpha = '(1'


                        elif activator_alpha.endswith('*'): # otherwise, if it's a running list, remove the * operator
                            activator_alpha = activator_alpha[:-1]

                        activator_alpha += "+" + indt


                    if activator_alpha.endswith("*"):
                        activator_alpha = activator_alpha[:-1]

                    # cap things off with a final parens:
                    activator_alpha += ")"

        else:

            activator_alpha += dl

        if i_list is not None and i_list != 'None' and len(i_list) > 0:

            # Next, construct the inhibitors net effect term:
            for i, (name, Kmo, n, zone_tag) in enumerate(zip(i_list, Km_i_list, n_i_list, zone_tags_i)):


                # Check for and deal with accessory name tag flags:
                if name.endswith('*'):
                    # capture the original name that will work in dictionaries:
                    name = name[:-1]
                    vs_request = True
                else:
                    vs_request = False

                direct_tag = False  # Force the "direct" tag to be false by default

                if name.endswith('!'):  # remove any bangs a user might specify; inhibitors always act multiplicatively

                    name = name[:-1]

                elif name.endswith('&'):
                    name = name[:-1]
                    direct_tag = True

                zi =self.molecules[name].z

                if in_mem_tag is True and zi != 0.0 and vs_request is True:

                    msg = "Vmem-sensitive gating used for inhibitor: " + name

                    logs.log_info(msg)

                    Kmc = '({}*np.exp(-(self.vm*{}*p.F)/(p.R*p.T)))'.format(Kmo, zi)
                    Kme = '({}*np.exp((self.vm*{}*p.F)/(p.R*p.T)))'.format(Kmo, zi)

                else:
                    Kmc = '{}'.format(Kmo)
                    Kme = '{}'.format(Kmo)



                # initialize string expressions for activator and inhibitor, respectively
                numo_string_i = ""
                denomo_string_i = ""
                direct_string_i = ""

                if reaction_zone == 'cell':

                    if zone_tag == 'cell':


                        numo_string_i += "1"
                        denomo_string_i += "(1 + (self.cell_concs['{}']/{})**{})".format(name, Kmc, n)

                        direct_string_i += "-self.cell_concs['{}']".format(name)

                    elif zone_tag == 'env':

                        numo_string_i += "1"
                        denomo_string_i += "(1 + (self.env_concs['{}']/{})**{})".format(name, Kme, n)
                        direct_string_i += "-self.env_concs['{}']".format(name)

                elif reaction_zone == 'mem':

                    if zone_tag == 'cell':

                        numo_string_i += "1"
                        denomo_string_i += "(1 + (self.cell_concs['{}']/{})**{})".format(name, Kmc, n)

                        direct_string_i += "-self.cell_concs['{}']".format(name)

                    elif zone_tag == 'env':

                        numo_string_i += "1"
                        denomo_string_i += "(1 + (self.env_concs['{}']/{})**{})".format(name, Kme, n)

                        direct_string_i += "-self.env_concs['{}']".format(name)

                elif reaction_zone == 'mit':

                    numo_string_i += "1"
                    denomo_string_i += "(1 + (self.mit_concs['{}']/{})**{})".format(name, Kmo, n)

                    direct_string_i += "-self.mit_concs['{}']".format(name)

                elif reaction_zone == 'env':

                    if zone_tag == 'cell':

                        numo_string_i += "1"

                        denomo_string_i += "1"

                        direct_string_i += "1"

                    elif zone_tag == 'env':

                        numo_string_i += "1"
                        denomo_string_i += "(1 + (self.env_concs['{}']/{})**{})".format(name, Kme, n)
                        direct_string_i += "-self.env_concs['{}']".format(name)


                term = "(" + numo_string_i + "/" + denomo_string_i + ")"


                #-------------------------------------------------

                if direct_tag is False:
                    inhibitor_alpha += term


                elif direct_tag is True:

                    direct_terms_list.append(direct_string_i)


                if i < len(i_list) - 1:

                    inhibitor_alpha += "*"


                elif i >= len(i_list) -1:   # else if we've reached the end:

                    if inhibitor_alpha.endswith("*"):

                        inhibitor_alpha = inhibitor_alpha[:-1]

                    inhibitor_alpha += ")"


        else:

            inhibitor_alpha += dl

        # finalize the multiplier:

        if inhibitor_alpha == "()":
            inhibitor_alpha = ""
            alpha = activator_alpha


        elif activator_alpha == "()":
            activator_alpha = ""
            alpha = inhibitor_alpha

        else:

            alpha = activator_alpha + '*' + inhibitor_alpha


        if len(direct_terms_list) > 0: # if there's anything in the direct terms list

            for q, dit in enumerate(direct_terms_list):

                direct_term_alpha += dit

                if q < len(direct_terms_list) - 1:
                    direct_term_alpha += '+'

                else:
                    direct_term_alpha += ")"  # finish it off with a parens or space:

        alpha_tex = ""
        tex_list = []


        return alpha, alpha_tex, tex_list

    def default_zones(self, zone_tags_a, zone_tags_i, a_list, i_list):

        # if zone tag lists aren't supplied or are invalid entries, create defaults:
        if zone_tags_a is None or zone_tags_a == 'None' or len(zone_tags_a) == 0:
            if a_list is not None and a_list != 'None' and len(a_list)>0:
                zone_tags_a = ['cell' for x in a_list]

        if zone_tags_i is None or zone_tags_i == 'None' or len(zone_tags_i) == 0:

            if i_list is not None and i_list != 'None' and len(i_list) > 0:
                zone_tags_i = ['cell' for x in i_list]

        return zone_tags_a, zone_tags_i

    def bal_charge(self, Q, tag, p):

        if tag == 'cell':

            if Q < 0 and np.abs(Q) <= self.cell_concs['P'].mean():  # if net charge is anionic

                self.cell_concs['P'] = self.cell_concs['P'] - np.abs(Q)

            elif Q > 0 and np.abs(Q) <= self.cell_concs['K'].mean():

                self.cell_concs['K'] = self.cell_concs['K'] - np.abs(Q)

            elif Q < 0 and np.abs(Q) > self.cell_concs['P'].mean():  # if net charge is anionic
                raise BetseSimConfException("You've defined way more anionic charge in"
                                               "the extra substances than we can "
                                               "compensate for. Either turn 'substances "
                                               "affect Vmem' off, or try again.")
            elif Q > 0 and np.abs(Q) > self.cell_concs['K'].mean():
                raise BetseSimConfException("You've defined way more cationic charge in"
                                               "the extra substances than we can "
                                               "compensate for. Either turn 'substances "
                                               "affect Vmem' off, or try again.")

        elif tag == 'env':

            if Q < 0 and np.abs(Q) <= self.cell_concs['P'].mean():  # if net charge is anionic
                self.cell_concs['P'] = self.cell_concs['P'] - np.abs(Q)
            elif Q > 0 and np.abs(Q) <= self.cell_concs['K'].mean():
                self.cell_concs['K']= self.cell_concs['K'] - np.abs(Q)
            elif Q < 0 and np.abs(Q) > self.cell_concs['P'].mean():  # if net charge is anionic
                raise BetseSimConfException("You've defined way more anionic charge in"
                                               "the extra substances than we can "
                                               "compensate for. Either turn 'substances "
                                               "affect Vmem' off, or try again.")
            elif Q > 0 and np.abs(Q) > self.cell_concs['K'].mean():
                raise BetseSimConfException("You've defined way more cationic charge in"
                                               "the extra substances than we can "
                                               "compensate for. Either turn 'substances "
                                               "affect Vmem' off, or try again.")

    def build_sim_matrix(self, p):

        # report
        logs.log_info("Building simulation matrix for whole BETSE network!")

        # If either NetworkX or PyDot are unimportable, raise an exception.
        libs.die_unless_runtime_optional('networkx', 'pydot')

        # Import NetworkX and NetworkX's PyDot interface:
        from networkx import nx_pydot

        # set the Vmem to target value requested by user:
        self.vm = self.vmo

        self.init_saving(p, plot_type='init', nested_folder_name='DynamicSystem')

        # create the reaction and concentration handlers, which organize different kinds of reactions and
        # concentrations:

        # create a complete graph using pydot and the network plotting method:
        grapha = self.build_reaction_network(p)

        # if saving is enabled, export the graph of the network used in the optimization:
        if p.autosave is True and p.plot_network is True:
            savename = self.imagePath + 'NetworkGraph' + '.svg'
            grapha.write_svg(savename, prog='dot')

        # convert the graph into a networkx format so it's easy to manipulate:
        network = nx_pydot.from_pydot(grapha)

        self.net_graph = network

        self.build_indices(p)

        # Optimization Matrix for Concentrations:----------------------------------------------------------------------

        # build a network matrix in order to easily organize reaction relationships needed for the optimization:
        self.network_opt_M = np.zeros((len(self.conc_handler), len(self.react_handler)))

        # build the reaction matrix based on the network reaction graph:
        for node_a, node_b in network.edges():

            # get the coefficient of stoichiometry used in the reaction relationship:
            edge_coeff = network.edge[node_a][node_b][0]['coeff']

            # when building the graph, the node shape was used to signify the item:
            node_type_a = network.node[node_a].get('shape', None)
            node_type_b = network.node[node_b].get('shape', None)

            # if node a is a concentration and b is a reaction, we must be dealing with a reactant:
            if node_type_a is self.conc_shape and node_type_b == self.reaction_shape:

                # find the numerical index of the reactant (node_a) and reaction (node_b) in the handlers:
                row_i = self.conc_handler_index[node_a]
                col_j = self.react_handler_index[node_b]

                # add them to the optimization matrix:
                self.network_opt_M[row_i, col_j] = -1 * edge_coeff

            # perform a similar type of reasoning for the opposite relationship (indicating a product):
            elif node_type_a == self.reaction_shape and node_type_b is self.conc_shape:

                row_i = self.conc_handler_index[node_b]
                col_j = self.react_handler_index[node_a]

                self.network_opt_M[row_i, col_j] = 1 * edge_coeff

            # the next two blocks to do exactly the same thing as above, except with transporters, which
            # have a diamond shape as they're treated differently:
            elif node_type_a is self.conc_shape and node_type_b == self.transporter_shape:

                row_i = self.conc_handler_index[node_a]
                col_j = self.react_handler_index[node_b]

                self.network_opt_M[row_i, col_j] = -1.0 * edge_coeff

            elif node_type_a == self.transporter_shape and node_type_b is self.conc_shape:

                row_i = self.conc_handler_index[node_b]
                col_j = self.react_handler_index[node_a]

                self.network_opt_M[row_i, col_j] = 1.0 * edge_coeff

            elif node_type_a == self.conc_shape and node_type_b is self.ed_shape:

                row_i = self.conc_handler_index[node_a]
                col_j = self.react_handler_index[node_b]

                self.network_opt_M[row_i, col_j] = -1.0 * edge_coeff

            elif node_type_a == self.ed_shape and node_type_b is self.conc_shape:

                row_i = self.conc_handler_index[node_b]
                col_j = self.react_handler_index[node_a]

                self.network_opt_M[row_i, col_j] = 1.0 * edge_coeff

            elif node_type_a == self.conc_shape and node_type_b is self.channel_shape:

                row_i = self.conc_handler_index[node_a]
                col_j = self.react_handler_index[node_b]

                self.network_opt_M[row_i, col_j] = -1.0 * edge_coeff

            elif node_type_a == self.channel_shape and node_type_b is self.conc_shape:

                row_i = self.conc_handler_index[node_b]
                col_j = self.react_handler_index[node_a]

                self.network_opt_M[row_i, col_j] = 1.0 * edge_coeff

        self.network_opt_M_full = self.network_opt_M * 1

        # --------------------------------------------------
        # This is idiotic, but refine the above matrix, assuming environmental concentrations are constant, and
        # adding in a row to take care of zero transmembrane current at steady-state

        # restructure the optimization matrix, first by removing the environmental values as we'll assume they're constant:
        MM = []
        self.output_handler = OrderedDict({})

        for i, (k, v) in enumerate(self.conc_handler.items()):

            # Don't include environmental concentrations in the optimization (assume they're constant)

            if k.endswith("_env"):

                pass

            else:
                MM.append(self.network_opt_M[i, :].tolist())

                name_tag = 'd/dt ' + str(k)
                self.output_handler[name_tag] = v

        # next define a current expression as the final row:
        Jrow = [0 for nn in self.react_handler]

        # create another vector that will give total sum of pump/transporter currents only
        self.Jpumps = [0 for nn in self.react_handler]

        # current scaling factor requires conversion back to flux units:
        ff = -1.0 * (1/self.div_cell) * p.F

        for i, (k, v) in enumerate(self.react_handler.items()):

            if k.endswith("_ed"):  # if we're dealing with an electrodiffusing item...

                kk = k[:-3]  # get the proper keyname...
                z = self.zmol[kk]  # get the charge value...
                Jrow[i] = z * ff

            if k in self.transporters:

                for li in self.transporters[k].transport_in_list:
                    z = self.zmol[li]
                    coeff = self.net_graph.edge[k][li][0]['coeff']
                    Jrow[i] += -coeff * z * ff
                    self.Jpumps[i] += -coeff * z * ff

                for lo in self.transporters[k].transport_out_list:
                    z = self.zmol[lo]
                    coeff = self.net_graph.edge[lo][k][0]['coeff']
                    Jrow[i] += coeff * z * ff
                    self.Jpumps[i] += coeff * z * ff

            if k in self.channels:

                chan_type = self.channels[k].channel_class

                if chan_type == 'K' or chan_type == 'Na':

                    Jrow[i] = ff

                elif chan_type == 'Cl':

                    Jrow[i] = -ff

                elif chan_type == 'Fun' or chan_type == 'Cat':

                    Jrow[i] = ff

        MM.append(Jrow)
        self.output_handler['Jmem'] = 0

        MM = np.array(MM)

        self.network_opt_M = MM

        origin_o = np.ones(len(self.react_handler))

        for i, rea_name in enumerate(self.react_handler):

            if rea_name.endswith('_growth'):

                name = rea_name[:-7]

                origin_o[i] = self.molecules[name].r_production

            elif rea_name.endswith('_decay'):

                name = rea_name[:-6]

                origin_o[i] = self.molecules[name].r_decay

            elif rea_name.endswith('_ed'):

                name = rea_name[:-3]

                origin_o[i] = self.Dmem[name]

            elif rea_name in self.reactions:

                origin_o[i] = self.reactions[rea_name].vmax

            elif rea_name in self.transporters:

                origin_o[i] = self.transporters[rea_name].vmax

            elif rea_name in self.channels:

                origin_o[i] = self.channels[rea_name].channelMax

        self.origin_o = origin_o

    def optimize(self, sim, cells, p):
        pass

    def direction_field(self, p, aa_params, bb_params, npoints = 25, calpha = 0.5):
        """
        Plots a direction field for 2 user-specified parameters of the network.
        :param sim:
        :param cells:
        :param p:
        :return:
        """


        Amin = aa_params['A_min']
        Amax = aa_params['A_max']
        Bmin = bb_params['B_min']
        Bmax = bb_params['B_max']

        aa = np.linspace(Amin, Amax, npoints)
        bb = np.linspace(Bmin, Bmax, npoints)

        AA, BB = np.meshgrid(aa, bb)

        outputdA = np.zeros(AA.shape)
        outputdB = np.zeros(BB.shape)

        name_A = aa_params['name']
        name_B = bb_params['name']

        self.vm = self.vmo

        # for rea in self.react_handler:
        #
        #     print(self.react_handler[rea])
        #     print('-----')

        for i in range(AA.shape[0]):

            for j in range(BB.shape[1]):

                if name_B == 'Vmem' and name_A in self.cell_concs:

                    self.vm = BB[i, j]
                    self.cell_concs[name_A] = AA[i, j]

                    tagA = 'd/dt ' + name_A
                    ind_A = list(self.output_handler).index(tagA)

                    ind_B = list(self.output_handler).index('Jmem')

                elif name_A == 'Vmem' and name_B in self.cell_concs:

                    self.cell_concs[name_B] = BB[i, j]
                    self.vm = AA[i, j]

                    tagB = 'd/dt ' + name_B
                    ind_B = list(self.output_handler).index(tagB)

                    ind_A = list(self.output_handler).index('Jmem')

                elif name_A in self.cell_concs and name_B in self.cell_concs:


                    self.cell_concs[name_A] = AA[i,j]
                    self.cell_concs[name_B] = BB[i,j]

                    tagA = 'd/dt ' + name_A
                    ind_A = list(self.output_handler).index(tagA)


                    tagB = 'd/dt ' + name_B
                    ind_B = list(self.output_handler).index(tagB)

                else:

                    raise BetseSimConfException("Something's not right with the way direction surface "
                                                  "entities have been specified. Please check the config "
                                                  "settings and try again.")

                r_base = [eval(self.react_handler[rea], self.globals, self.locals).mean() for rea in
                          self.react_handler]

                outputs = np.dot(self.network_opt_M, r_base)

                outputdB[i, j] = outputs[ind_B]
                outputdA[i, j] = outputs[ind_A]

        MagM = (np.hypot(outputdA, outputdB))
        # avoid zero division errors
        MagM[MagM == 0] = 1.
        # normalize each arrow
        outputdA2 = outputdA / MagM
        outputdB2 = outputdB / MagM

        UU = np.sqrt(outputdA ** 2 + (outputdB) ** 2)

        fig = plt.figure()
        ax = plt.subplot(111)
        plt.pcolormesh(AA, BB, UU, cmap = p.background_cm, alpha = calpha)
        ax.quiver(AA, BB, outputdA2, outputdB2)
        plt.axis([Amin, Amax, Bmin, Bmax])

        ax.set_title('Direction field for ' + name_A + ' versus ' + name_B)
        ax.set_xlabel('Value of ' + name_A)
        ax.set_ylabel('Value of ' + name_B)

        if p.autosave:
            savename = self.imagePath + 'DirectionField_' + name_A + '_' + name_B + '.png'
            plt.savefig(savename, format='png', transparent=True)

        plt.close()

        #-----

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        surf = ax.plot_surface(AA, BB, UU, linewidth=0.0, rstride=1, cstride=1,
                               shade=True, cmap=p.background_cm, antialiased=True)

        # ax.set_title('Direction Surface for ' + name_A + ' versus ' + name_B)
        ax.set_xlabel('Value of ' + name_A)
        ax.set_ylabel('Value of ' + name_B)
        # ax.view_init(elev=10, azim=-70)

        if p.autosave is True:
            savename = self.imagePath + 'DirectionSurface3D_' + name_A + '_' + name_B + '.png'
            plt.savefig(savename, format='png', transparent=True)

        plt.close()

        # if p.turn_all_plots_off is False:
        #     plt.show(block=False)

class Molecule(object):

    def __init__(self, p):

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
        self.n_decay = None

        self.growth_activators_list = None
        self.growth_activators_Km = None
        self.growth_activators_n = None
        self.growth_inhibitors_list = None
        self.growth_inhibitors_Km = None
        self.growth_inhibitors_n = None

        self.ion_activators_list = None
        self.ion_activators_Km = None
        self.ion_activators_n = None
        self.ion_inhibitors_list = None
        self.ion_inhibitors_Km = None
        self.ion_inhibitors_n = None

        self.growth_targets_cell = None

        self.gating_mod = 1.0

class Reaction(object):

    def __init__(self, p):

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

class Transporter(object):

    def __init__(self, p):

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

class Channel(object):

    def __init__(self, p):

        self.channel_eval_string = ""

    def init_channel(self, ion_string, type_string, max_val, alpha_mod, p):

        for ion in ion_string:

            c_in = "self.cell_concs['{}']".format(ion)
            c_out = "self.env_concs['{}']".format(ion)

            zz = "self.zmol['{}']".format(ion)

            front_term = "((({}*1e-9)/p.tm)*((p.F*self.vm)/(p.R*p.T))*{})".format(max_val, zz)

            alpha_v = "(({}*self.vm*p.F)/(p.R*p.T))".format(zz)

            numo_term = "({} - {}*np.exp(-{}))".format(c_in, c_out, alpha_v)
            deno_term = "(1 - np.exp(-{}))".format(alpha_v)

            # write the eval string for the channel:
            self.channel_eval_string = "{}*{}*({}/{})".format(alpha_mod, front_term, numo_term, deno_term)


