#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# FIXME this is the old version of what's currently networks.py. This is kept around for reference, just in case!

"""

Create and electrodiffuses a suite of customizable general molecule in the BETSE ecosystem,
including functionality to pump the molecule, use it as a gating ligand, produce and consume it,
and use it an enzyme to facilitate another reaction. The molecule is assumed to be at low
concentrations and to not have a significant effect on system voltages or currents. This
module creates a structure containing all user-defined molecules, along with the facilities
to initialize, define the core computations for a simulation loop, save and report on data,
and plot.

"""

# import os
# import os.path
# import numpy as np
# from betse.science import toolbox as tb
# from betse.science import sim_toolbox as stb
# from betse.science.tissue.handler import TissueHandler
# from betse.science.event import modulators as mods
# from betse.util.io.log import logs
# import matplotlib.pyplot as plt
# from betse.exceptions import BetseSimConfigException
# from betse.science.visual import plot as viz
# from betse.science.visual.anim.anim import AnimCellsMembranesData, AnimEnvTimeSeries
# from betse.science.organelles.mitochondria import Mito
# from matplotlib import colors
# from matplotlib import cm
#
# from betse.science.tissue.channels import vg_na as vgna
# from betse.science.tissue.channels import vg_nap as vgnap
# from betse.science.tissue.channels import vg_k as vgk
# from betse.science.tissue.channels import vg_kir as vgkir
# from betse.science.tissue.channels import vg_funny as vgfun
# from betse.science.tissue.channels import vg_ca as vgca
#
#
# class MasterOfMolecules(object):
#
#     def __init__(self, sim, cells, config_substances, p, mit_enabled = False):
#
#         """
#          Initializes the MasterOfMolecules object.
#
#         sim:                An instance of simulator
#         config_settings     List of dictionaries storing key settings (p.molecules_config)
#         p                   An instance of params
#         """
#         # all metabolic simulations require ATP, ADP and Pi. Initialize these fields to None so that we can test
#         # for their pressence in metabolic sims:
#         self.ATP = None
#         self.ADP = None
#         self.Pi = None
#
#         # set the key controlling presence of mitochondria:
#         self.mit_enabled = mit_enabled
#
#         self.read_substances(sim, cells, config_substances, p)
#         self.tissue_init(sim, cells, config_substances, p)
#
#         self.ave_cell_vol = cells.cell_vol.mean()  # average cell volume
#
#         self.reaction_names = []
#         self.transporter_names = []
#         self.channel_names = []
#         self.modulator_names = []
#
#         self.plot_cmap = 'viridis'
#
#     def read_substances(self, sim, cells, config_substances, p):
#         """
#             Initializes all of the main variables for all molecules included in the simulation.
#
#             config_substances:  dictionary containing BETSE biomolecule template fields
#
#         """
#         # Initialize a list that will keep track of molecule names in the simulation
#         self.molecule_names = []
#         # Initialize a list that keeps track of molecule index
#         self.molecule_index = []
#
#         for q, mol_dic in enumerate(config_substances):
#             # get each user-defined name-filed in the dictionary:
#             name = str(mol_dic['name'])
#
#             # add the name to the name catalogue:
#             self.molecule_names.append(name)
#             self.molecule_index.append(q)
#
#             # add a field to the MasterOfMolecules corresponding to Molecule object with that name
#             setattr(self, name, Molecule(sim, cells, p))
#
#             # now set the attributes of that Molecule object with the cornucopia of variables:
#
#             # get MasterOfMolecules.name
#             obj = getattr(self, name)
#
#             # assign general properties
#             obj.name = name   # let object know who it is
#
#             obj.c_envo = mol_dic['env conc']  # initial concentration in the environment [mmol/L]
#             obj.c_cello = mol_dic['cell conc']  # initial concentration in the cytoplasm [mmol/L]
#
#             obj.c_mito = mol_dic.get('mit conc',None)  # initialized to None if optional fields not present
#             obj.c_ero = mol_dic.get('er conc', None)
#
#             # create data structures to use with sim --------
#             # initialize concentrations in cells:
#             obj.c_cells = np.ones(sim.cdl)*obj.c_cello
#             obj.c_mems = np.ones(sim.mdl) * obj.c_cello
#
#             # if there is an initial concentration for mitochondria, define a conc vector for it:
#             if obj.c_mito is not None and self.mit_enabled:
#
#                 obj.c_mit = np.ones(sim.cdl)*obj.c_mito
#
#             elif self.mit_enabled and obj.c_mito is None:
#                 obj.c_mit = np.zeros(sim.cdl)
#
#             # if there is an initial concentration for endo retic, define a conc vector for it:
#             if obj.c_ero is not None:
#
#                 obj.c_er = np.ones(sim.cdl)*obj.c_ero
#
#             # initialize concentration in the environment:
#             if p.sim_ECM is False:
#                 obj.c_env = np.ones(sim.mdl) * obj.c_envo
#             else:
#                 obj.c_env = np.ones(sim.edl) * obj.c_envo
#
#             # initialize concentration at the boundary
#             obj.c_bound = obj.c_envo
#
#             if self.mit_enabled:
#                 obj.mit_enabled = True
#             else:
#                 obj.mit_enabled = False
#
#         # if mitochondria are enabled:
#         if self.mit_enabled:
#
#             self.mit = Mito(sim, cells, p)
#         else:
#             self.mit = None
#
#     def tissue_init(self, sim, cells, config_substances, p):
#
#         for q, mol_dic in enumerate(config_substances):
#             # get each user-defined name-filed in the dictionary:
#             name = str(mol_dic['name'])
#
#             if name not in self.molecule_names:
#
#                 # add the name to the name catalogue:
#                 self.molecule_names.append(name)
#                 self.molecule_index.append(q)
#
#                 # add a field to the MasterOfMolecules corresponding to Molecule object with that name
#                 setattr(self, name, Molecule(sim, cells, p))
#
#                 # now set the attributes of that Molecule object with the cornucopia of variables:
#
#                 # get MasterOfMolecules.name
#                 obj = getattr(self, name)
#
#                 # assign general properties
#                 obj.name = name  # let object know who it is
#
#                 obj.c_envo = mol_dic['env conc']  # initial concentration in the environment [mmol/L]
#                 obj.c_cello = mol_dic['cell conc']  # initial concentration in the cytoplasm [mmol/L]
#
#                 obj.c_mito = mol_dic.get('mit conc', None)  # initialized to None if optional fields not present
#                 obj.c_ero = mol_dic.get('er conc', None)
#
#                 # create data structures to use with sim --------
#                 # initialize concentrations in cells:
#                 obj.c_cells = np.ones(sim.cdl) * obj.c_cello
#                 obj.c_mems = np.ones(sim.mdl) * obj.c_cello
#
#                 # if there is an initial concentration for mitochondria, define a conc vector for it:
#                 if obj.c_mito is not None and self.mit_enabled:
#
#                     obj.c_mit = np.ones(sim.cdl) * obj.c_mito
#
#                 elif self.mit_enabled and obj.c_mito is None:
#                     obj.c_mit = np.zeros(sim.cdl)
#
#                 # if there is an initial concentration for endo retic, define a conc vector for it:
#                 if obj.c_ero is not None:
#                     obj.c_er = np.ones(sim.cdl) * obj.c_ero
#
#                 # initialize concentration in the environment:
#                 if p.sim_ECM is False:
#                     obj.c_env = np.ones(sim.mdl) * obj.c_envo
#                 else:
#                     obj.c_env = np.ones(sim.edl) * obj.c_envo
#
#                 # initialize concentration at the boundary
#                 obj.c_bound = obj.c_envo
#
#                 if self.mit_enabled:
#                     obj.mit_enabled = True
#                 else:
#                     obj.mit_enabled = False
#
#             # get MasterOfMolecules.name
#             obj = getattr(self, name)
#
#             obj.Dm = mol_dic['Dm']  # membrane diffusion coefficient [m2/s]
#             obj.Do = mol_dic['Do']  # free diffusion constant in extra and intracellular spaces [m2/s]
#             obj.z = mol_dic['z']  # charge (oxidation state)
#
#             obj.ignore_ECM_pump = mol_dic['ignore ECM']
#
#             obj.ignoreTJ = mol_dic['TJ permeable']  # ignore TJ?
#             obj.ignoreGJ = mol_dic['GJ impermeable'] # ignore GJ?
#
#             obj.TJ_factor = float(mol_dic['TJ factor'])
#
#             # factors involving auto-catalytic growth and decay in the cytoplasm
#             gad = mol_dic.get('growth and decay', None)
#
#             if gad != 'None' and gad is not None:
#
#                 obj.simple_growth = True
#
#                 obj.r_production = gad['production rate']
#                 obj.r_decay = gad['decay rate']
#                 obj.Kgd = gad['Km']
#                 obj.n_production = gad['n']
#
#                 obj.growth_profiles_list = gad['apply to']
#
#                 modulator_function_name = gad.get('modulator function', None)
#
#                 obj.growth_activators_list = gad.get('activators', None)
#                 obj.growth_activators_k = gad.get('k activators', None)
#                 obj.growth_activators_Km = gad.get('Km activators', None)
#
#                 obj.growth_activators_n = gad.get('k activators', None)
#                 obj.growth_inhibitors_list = gad.get('inhibitors', None)
#                 obj.growth_inhibitors_k = gad.get('k inhibitors', None)
#                 obj.growth_inhibitors_Km = gad.get('Km inhibitors', None)
#                 obj.growth_inhibitors_n = gad.get('k inhibitors', None)
#
#                 obj.init_growth(cells, p)
#
#                 if modulator_function_name != 'None' and modulator_function_name is not None:
#                     obj.growth_mod_function_mems, _ = getattr(mods, modulator_function_name)(obj.growth_targets_mem,
#                                                                                               cells, p)
#                     obj.growth_mod_function_cells, _ = getattr(mods, modulator_function_name)(obj.growth_targets_cell,
#                                                                                                cells, p)
#
#                 else:
#                     obj.growth_mod_function_mems = 1
#                     obj.growth_mod_function_cells = 1
#
#             else:
#                 obj.simple_growth = False
#                 obj.growth_profiles_list = None
#
#             # assign ion channel gating properties
#             icg = mol_dic.get('ion channel gating', None)
#
#             if icg is not None:
#
#                 obj.ion_channel_gating = True
#
#                 gating_ion_o = icg['ion channel target']  # get a target ion label to gate membrane to (or 'None')
#
#                 if gating_ion_o != 'None':
#                     obj.use_gating_ligand = True
#
#                     obj.gating_ion = []
#
#                     for ion_o in gating_ion_o:
#                         obj.gating_ion.append(sim.get_ion(ion_o))
#
#                 else:
#                     obj.use_gating_ligand = False
#                     obj.gating_ion = []
#
#                 obj.gating_Hill_K = float(icg['target Hill coefficient'])
#                 obj.gating_Hill_n = float(icg['target Hill exponent'])
#                 obj.gating_max_val = float(icg['peak channel opening'])
#                 obj.gating_extracell = icg['acts extracellularly']
#
#                 # get any optional activators and inhibitors for the channel:
#                 obj.activators_list = icg.get('activators', None)
#                 obj.activators_Km = icg.get('Km activators', None)
#                 obj.activators_n = icg.get('n activators', None)
#
#                 obj.inhibitors_list = icg.get('inhibitors', None)
#                 obj.inhibitors_Km = icg.get('Km inhibitors', None)
#                 obj.inhibitors_n = icg.get('n inhibitors', None)
#
#             else:
#
#                 obj.ion_channel_gating = False
#
#             # assign active pumping properties
#             ap = mol_dic.get('active pumping', None)
#
#             if ap is not None:
#
#                 obj.active_pumping = True
#                 obj.use_pumping = ap['turn on']
#                 obj.pump_to_cell = ap['pump to cell']
#                 obj.pump_max_val = ap['maximum rate']
#                 obj.pump_Km = ap['pump Km']
#                 obj.pumps_use_ATP = ap['uses ATP']
#
#
#             else:
#                 obj.active_pumping = False
#
#             # assign boundary change event properties
#             cab = mol_dic.get('change at bounds', None)
#
#             if cab is not None:
#                 obj.change_bounds = True
#                 obj.change_at_bounds = cab['event happens']
#                 obj.change_bounds_start = cab['change start']
#                 obj.change_bounds_end = cab['change finish']
#                 obj.change_bounds_rate = cab['change rate']
#                 obj.change_bounds_target = cab['concentration']
#
#             else:
#                 obj.change_bounds = False
#
#             # assign plotting properties
#             pd = mol_dic['plotting']
#             obj.make_plots = pd['plot 2D']
#             obj.make_ani = pd['animate']
#
#             obj.plot_autoscale = pd['autoscale colorbar']
#             obj.plot_max = pd['max val']
#             obj.plot_min = pd['min val']
#
#     def plot_init(self, config_substances):
#
#         for q, mol_dic in enumerate(config_substances):
#             # get each user-defined name-filed in the dictionary:
#             name = str(mol_dic['name'])
#
#             # get MasterOfMolecules.name
#             obj = getattr(self, name)
#
#             # assign plotting properties
#             pd = mol_dic['plotting']
#             obj.make_plots = pd['plot 2D']
#             obj.make_ani = pd['animate']
#
#             obj.plot_autoscale = pd['autoscale colorbar']
#             obj.plot_max = pd['max val']
#             obj.plot_min = pd['min val']
#
#     def read_reactions(self, config_reactions, sim, cells, p):
#
#         """
#           Read in and initialize parameters for all user-defined reactions.
#
#           config_options:  dictionary containing BETSE reaction template fields
#
#           """
#
#         # Initialize a list that keeps the index of the reaction:
#         self.reaction_index = []
#         self.reaction_names = []
#
#         for q, react_dic in enumerate(config_reactions):
#             # get each user-defined name-filed in the dictionary:
#             name = str(react_dic['name'])
#
#             # add the name to the name catalogue:
#             self.reaction_names.append(name)
#             self.reaction_index.append(q)
#
#             # add a field to the MasterOfReactions corresponding to Reaction object with that name:
#             #self._reactions[name] = Reaction()
#             setattr(self, name, Reaction(sim, cells, p))
#
#             # now set the attributes of that Reaction object with the cornucopia of variables:
#
#             # get MasterOfMolecules.name
#             obj = getattr(self, name)
#
#             # assign general properties
#             obj.name = name  # let object know who it is
#
#             # list where the reaction takes place; if field not specified default to 'cell':
#             obj.reaction_zone = react_dic['reaction zone']
#
#             # set the main fields of the reaction object
#             obj.reactants_list = react_dic['reactants']
#             obj.reactants_coeff = react_dic['reactant multipliers']
#
#             obj.products_list = react_dic['products']
#             obj.products_coeff = react_dic['product multipliers']
#
#             obj.Km_reactants_list = react_dic['Km reactants']
#             obj.Km_products_list = react_dic['Km products']
#             obj.vmax = float(react_dic['max rate'])
#
#             obj.delta_Go = react_dic['standard free energy']
#
#             if obj.delta_Go == 'None':
#                 obj.delta_Go = None     # make the field a proper None variable
#
#             else:
#                 obj.delta_Go = float(obj.delta_Go)
#
#             obj.reaction_activators_list = react_dic.get('reaction activators', None)
#             obj.reaction_activators_Km = react_dic.get('activator Km', None)
#             obj.reaction_activators_n = react_dic.get('activator n', None)
#             obj.reaction_inhibitors_list = react_dic.get('reaction inhibitors', None)
#             obj.reaction_inhibitors_Km = react_dic.get('inhibitor Km', None)
#             obj.reaction_inhibitors_n = react_dic.get('inhibitor n', None)
#
#             if self.mit_enabled:
#                 obj.mit_enabled = True
#             else:
#                 obj.mit_enabled = False
#
#             # now we want to load the right concentration data arrays for reactants into the Reaction object:
#
#             # check the reaction states to make sure all reagents/products are named biomolecules or ions:
#             # self.check_reactions(obj, sim, cells, p)
#
#     def read_transporters(self, config_transporters, sim, cells, p):
#
#         """
#             Read in and initialize parameters for all user-defined transporters.
#
#             config_options:  dictionary containing BETSE transporter template fields
#
#             """
#
#         # Initialize a list that keeps the index of the reaction:
#         self.transporter_index = []
#         self.transporter_names = []
#
#         for q, trans_dic in enumerate(config_transporters):
#             # get each user-defined name-filed in the dictionary:
#             name = str(trans_dic['name'])
#
#             # add the name to the name catalogue:
#             self.transporter_names.append(name)
#             self.transporter_index.append(q)
#
#             # add a field to the MasterOfReactions corresponding to Transporter object with that name:
#
#             setattr(self, name, Transporter(sim, cells, p))
#
#             # now set the attributes of that Transporter object with the cornucopia of variables:
#
#             # get MasterOfMolecules.name
#             obj = getattr(self, name)
#
#             # assign general properties
#             obj.name = name  # let object know who it is
#
#             # list where the reaction takes place; if field not specified default to 'cell':
#             obj.reaction_zone = trans_dic['reaction zone']
#
#             # set the main fields of the reaction object
#             obj.reactants_list = trans_dic['reactants']
#             obj.reactants_coeff = trans_dic['reactant multipliers']
#
#             obj.products_list = trans_dic['products']
#             obj.products_coeff = trans_dic['product multipliers']
#
#             obj.Km_reactants_list = trans_dic['Km reactants']
#             obj.Km_products_list = trans_dic['Km products']
#
#             obj.transport_out_list = trans_dic['transfered out of cell']
#             obj.transport_in_list = trans_dic['transfered into cell']
#
#             obj.vmax = float(trans_dic['max rate'])
#
#             obj.delta_Go = trans_dic['standard free energy']
#             obj.ignore_ECM_transporter = trans_dic['ignore ECM']
#
#             obj.transporter_profiles_list = trans_dic['apply to']
#             obj.init_reaction(cells, p)
#
#             if obj.delta_Go == 'None':
#                 obj.delta_Go = None  # make the field a proper None variable
#
#             else:
#                 obj.delta_Go = float(obj.delta_Go)
#
#             obj.transporter_activators_list = trans_dic.get('transporter activators', None)
#             obj.transporter_activators_Km = trans_dic.get('activator Km', None)
#             obj.transporter_activators_n = trans_dic.get('activator n', None)
#
#             obj.transporter_inhibitors_list = trans_dic.get('transporter inhibitors', None)
#             obj.transporter_inhibitors_Km = trans_dic.get('inhibitor Km', None)
#             obj.transporter_inhibitors_n = trans_dic.get('inhibitor n', None)
#
#             if self.mit_enabled:
#                 obj.mit_enabled = True
#             else:
#                 obj.mit_enabled = False
#
#     def read_channels(self, config_channels, sim, cells, p):
#
#         # Initialize a list that keeps the index of the channel:
#         self.channel_index = []
#         self.channel_names = []
#
#         for q, chan_dic in enumerate(config_channels):
#             # get each user-defined name-filed in the dictionary:
#             name = str(chan_dic['name'])
#
#             # add the name to the name catalogue:
#             self.channel_names.append(name)
#             self.channel_index.append(q)
#
#             # add a field to the MasterOfReactions corresponding to Channel object with that name:
#             setattr(self, name, Channel(sim, cells, p))
#
#             # now set the attributes of that channel:
#
#             # get MasterOfMolecules.name
#             obj = getattr(self, name)
#
#             # assign general properties
#             obj.name = name  # let object know who it is
#
#             # list where the reaction takes place; if field not specified default to 'cell':
#             obj.reaction_zone = 'cell'
#
#             obj.channel_class = chan_dic['channel class']
#             obj.channel_type = chan_dic['channel type']
#             obj.channelMax = chan_dic['max conductivity']
#             obj.channel_profiles_list = chan_dic['apply to']
#
#             obj.channel_activators_list = chan_dic.get('channel activators', None)
#             obj.channel_activators_Km = chan_dic.get('activator Km', None)
#             obj.channel_activators_n = chan_dic.get('activator n', None)
#             obj.channel_inhibitors_list = chan_dic.get('channel inhibitors', None)
#             obj.channel_inhibitors_Km = chan_dic.get('inhibitor Km', None)
#             obj.channel_inhibitors_n = chan_dic.get('inhibitor n', None)
#
#
#             obj.init_channel(obj.channel_class, obj.channel_type, obj.channelMax, sim, cells, p)
#
#     def read_modulators(self, config_modulators, sim, cells, p):
#
#         self.modulator_index = []
#         self.modulator_names = []
#
#         for q, mod_dic in enumerate(config_modulators):
#
#             name = str(mod_dic['name'])
#
#             self.modulator_names.append(name)
#             self.modulator_index.append(q)
#
#             # add a field to the MasterOfReactions corresponding to Channel object with that name:
#             setattr(self, name, Modulator())
#
#             # now set the attributes of that channel:
#
#             # get MasterOfMolecules.name
#             obj = getattr(self, name)
#
#             obj.target_label = str(mod_dic['target'])
#             obj.zone = str(mod_dic.get('zone', 'cell'))
#             obj.max_val = float(mod_dic['max effect'])
#             obj.modulator_activators_list = mod_dic.get('activators', None)
#             obj.modulator_activators_Km = mod_dic.get('activator Km', None)
#             obj.modulator_activators_n = mod_dic.get('activator n', None)
#             obj.modulator_inhibitors_list = mod_dic.get('inhibitors', None)
#             obj.modulator_inhibitors_Km = mod_dic.get('inhibitor Km', None)
#             obj.modulator_inhibitors_n = mod_dic.get('inhibitor n', None)
#
#             obj.init_modulator(sim, cells, p)
#
#     def set_react_sources(self, obj, sim, cells, p, reactant_name, i, reactant_type_self, reactant_type_sim):
#
#         """
#
#         obj:                    Molecule field of MasterOfMolecules
#         sim:                    Instance of BETSE sim
#         reactant_type_self:     Concentration field of Molecule: 'c_cells', 'c_env', 'c_mit' or 'c_er'
#         reactant_type_sim:      Concentration field of sim: 'cc_cells', 'cc_env',  'cc_mit' or 'cc_er'
#
#         """
#
#         label = 'i' + reactant_name
#         ion_check = getattr(sim, label, None)
#
#         if ion_check is None:
#
#             obj.reactant_source_object.append(id(self))
#             obj.reactant_source_type.append(reactant_type_self)
#
#         else:
#
#             obj.reactant_source_object.append(id(sim))
#             obj.reactant_source_type.append(reactant_type_sim)
#
#         # add the index of the reaction to the list so we can access modifiers like reaction coefficients
#         # and Km values at a later point:
#         obj.inds_react.append(i)
#
#     def set_prod_sources(self, obj, sim, cells, p, product_name, j, product_type_self, product_type_sim):
#
#         # Now load the right concentration data arrays for products into the Reaction object:
#         # see if the name is an ion defined in sim:
#         label = 'i' + product_name
#         ion_check = getattr(sim, label, None)
#
#         if ion_check is None:
#
#             obj.product_source_object.append(id(self))
#             obj.product_source_type.append(product_type_self)
#
#         else:
#
#             obj.product_source_object.append(id(sim))
#             obj.product_source_type.append(product_type_sim)
#
#         # add the index of the reaction to the list so we can access modifiers like reaction coefficients
#         # and Km values at a later point:
#         obj.inds_prod.append(j)
#
#     def init_saving(self, cells, p, plot_type = 'init', nested_folder_name = 'Molecules'):
#
#             # init files
#             if p.autosave is True:
#
#                 if plot_type == 'sim':
#                     results_path = os.path.join(p.sim_results, nested_folder_name)
#                     p.plot_type = 'sim'
#
#                 elif plot_type == 'init':
#                     results_path = os.path.join(p.init_results, nested_folder_name)
#                     p.plot_type = 'init'
#
#                 self.resultsPath = os.path.expanduser(results_path)
#                 os.makedirs(self.resultsPath, exist_ok=True)
#
#                 self.imagePath = os.path.join(self.resultsPath, 'fig_')
#
#             # check that the plot cell is in range of the available cell indices:
#             if p.plot_cell not in cells.cell_i:
#                 raise BetseSimConfigException(
#                     'The "plot cell" defined in the "results" section of your '
#                     'configuration file does not exist in your cluster. '
#                     'Choose a plot cell number smaller than the maximum cell number.')
#
#     def run_loop(self, t, sim, cells, p):
#         """
#         Runs the main simulation loop steps for each of the molecules included in the simulation.
#
#         """
#
#         # Initialize arrays for substance charge contribution:
#         net_Q_cell = 0
#         net_Q_env = 0
#
#         # get the name of the specific substance:
#         for name in self.molecule_names:
#
#             obj = getattr(self,name)
#
#             # if pumping is enabled:
#             if obj.active_pumping:
#                 obj.pump(sim, cells, p)
#
#             # update the production-decay regulatory network (if defined):
#             if obj.simple_growth is True:
#                 obj.growth_and_decay(self, sim, cells, p)
#
#             if p.run_sim is True:
#                 # use the substance as a gating ligand (if desired)
#                 if obj.ion_channel_gating:
#                     obj.gating(sim, self, cells, p)
#
#                 # update the global boundary (if desired)
#                 if obj.change_bounds:
#                     obj.update_boundary(t, p)
#
#             # transport the molecule through gap junctions and environment:
#             obj.transport(sim, cells, p)
#
#             # update the substance on the inside of the cell:
#             obj.updateIntra(sim, self, cells, p)
#
#             # calculate energy charge in the cell:
#             self.energy_charge(sim)
#
#             # ensure no negs:
#             stb.no_negs(obj.c_mems)
#
#             # calculate the charge density this substance contributes to cell and environment:
#             obj_Q_cell = p.F * obj.c_mems * obj.z
#
#             obj_Q_env = p.F * obj.c_env * obj.z
#             # add that contribution to the total sum:
#             net_Q_cell = net_Q_cell + obj_Q_cell
#             net_Q_env = net_Q_env + obj_Q_env
#
#         if p.substances_affect_charge:
#             # update charge in the cell and environment of the main bioelectric simulator:
#             sim.rho_cells = sim.rho_cells + net_Q_cell
#             sim.rho_env = sim.rho_env + net_Q_env
#
#         if self.mit_enabled:  # if enabled, update the mitochondria's voltage and other properties
#
#             self.mit.update(sim, cells, p)
#
#         # manage pH in cells, environment and mitochondria:
#         self.pH_handling(sim, cells, p)
#
#     def pH_handling(self, sim, cells, p):
#         """
#         Molecules may contain dissolved carbon dioxide as a substance,
#         and reactions/transporters may act on H+ levels via bicarbonate
#         (M-). Therefore, update pH in cells, environment, and
#         if enabled, mitochondria.
#
#         """
#
#
#         if 'CO2' in self.molecule_names:
#
#             if p.ions_dict['H'] == 1:
#
#                 # if the simulation contains sim.cHM_mems, use it and update it!
#                 sim.cHM_mems = self.CO2.c_cells
#                 sim.cHM_env = self.CO2.c_env
#
#                 # update the cH and pH fields of sim with potentially new value of sim.iM
#                 sim.cc_cells[sim.iH], sim.pH_cell = stb.bicarbonate_buffer(self.CO2.c_cells, sim.cc_cells[sim.iM])
#                 sim.cc_env[sim.iH], sim.pH_env = stb.bicarbonate_buffer(self.CO2.c_env, sim.cc_env[sim.iM])
#
#                 if self.mit_enabled:
#                     # update the cH and pH fields of sim with potentially new value of sim.iM
#                     sim.cc_mit[sim.iH], sim.pH_mit = stb.bicarbonate_buffer(self.CO2.c_mit, sim.cc_mit[sim.iM])
#
#             elif p.ions_dict['H'] != 1:
#
#                 # update the cH and pH fields of sim with potentially new value of M ion:
#                 _, sim.pH_cell = stb.bicarbonate_buffer(self.CO2.c_cells, sim.cc_cells[sim.iM])
#                 _, sim.pH_env = stb.bicarbonate_buffer(self.CO2.c_env, sim.cc_env[sim.iM])
#
#                 if self.mit_enabled:
#                     # update the cH and pH fields of sim with potentially new value of sim.iM
#                     _, sim.pH_mit = stb.bicarbonate_buffer(self.CO2.c_mit, sim.cc_mit[sim.iM])
#
#         else: # if we're not using CO2 in the simulator, use the default p.CO2*0.03
#
#             CO2 = p.CO2*0.03 # get the default concentration of CO2
#
#             if p.ions_dict['H'] == 1:
#
#                 # update the cH and pH fields of sim with potentially new value of sim.iM
#                 sim.cc_cells[sim.iH], sim.pH_cell = stb.bicarbonate_buffer(CO2, sim.cc_cells[sim.iM])
#                 sim.cc_env[sim.iH], sim.pH_env = stb.bicarbonate_buffer(CO2, sim.cc_env[sim.iM])
#
#                 if self.mit_enabled:
#                     # update the cH and pH fields of sim with potentially new value of sim.iM
#                     sim.cc_mit[sim.iH], sim.pH_mit = stb.bicarbonate_buffer(CO2, sim.cc_mit[sim.iM])
#
#             elif p.ions_dict['H'] != 1:
#
#                 # update the cH and pH fields of sim with potentially new value of sim.iM
#                 _, sim.pH_cell = stb.bicarbonate_buffer(CO2, sim.cc_cells[sim.iM])
#                 _, sim.pH_env = stb.bicarbonate_buffer(CO2, sim.cc_env[sim.iM])
#
#                 if self.mit_enabled:
#                     # update the cH and pH fields of sim with potentially new value of sim.iM
#                     _, sim.pH_mit = stb.bicarbonate_buffer(CO2, sim.cc_mit[sim.iM])
#
#     def energy_charge(self, sim):
#
#         if 'AMP' in self.molecule_names:
#
#             numo = (self.ATP.c_cells + 0.5*self.ADP.c_cells)
#             denomo = (self.ATP.c_cells + self.ADP.c_cells + self.AMP.c_cells)
#
#             self.chi = numo/denomo
#
#         else:
#
#             self.chi = np.zeros(sim.cdl)
#
#     def updateInside(self, sim, cells, p):
#         """
#         Runs the main simulation loop steps for each of the molecules included in the simulation.
#
#         """
#
#         # get the name of the specific substance:
#         for name in self.molecule_names:
#
#             obj = getattr(self, name)
#             obj.updateIntra(sim, self, cells, p)
#
#     def run_loop_reactions(self, t, sim, sim_metabo, cells, p):
#
#         # get the object corresponding to the specific reaction:
#         for i, name in enumerate(self.reaction_names):
#
#             # get the Reaction object:
#             obj = getattr(self, name)
#
#             # compute the new reactants and products
#             obj.rate = obj.compute_reaction(sim, sim_metabo, cells, p)
#
#     def run_loop_transporters(self, t, sim, sim_metabo, cells, p):
#
#         # get the object corresponding to the specific transporter:
#         for i, name in enumerate(self.transporter_names):
#
#             # get the Reaction object:
#             obj = getattr(self, name)
#
#             # compute the new reactants and products
#             obj.rate = obj.compute_reaction(sim, sim_metabo, cells, p)
#
#     def run_loop_channels(self, sim, sim_metabo, cells, p):
#
#         # get the object corresponding to the specific transporter:
#         for i, name in enumerate(self.channel_names):
#
#             # get the Reaction object:
#             obj = getattr(self, name)
#
#             # compute the channel activity
#             obj.run_channel(sim, sim_metabo, cells, p)
#
#     def run_loop_modulators(self, sim, sim_metabo, cells, p):
#
#         # get the object corresponding to the specific transporter:
#         for i, name in enumerate(self.modulator_names):
#             # get the Reaction object:
#             obj = getattr(self, name)
#
#             # compute the channel activity
#             obj.run_modulator(sim, sim_metabo, cells, p)
#
#     def mod_after_cut_event(self,target_inds_cell, target_inds_mem, sim, cells, p, met_tag = False):
#
#         # get the name of the specific substance:
#         for name in self.molecule_names:
#             obj = getattr(self, name)
#
#             obj.remove_cells(target_inds_cell, target_inds_mem, sim, cells, p)
#
#         if self.mit_enabled:
#
#             self.mit.remove_mits(sim, target_inds_cell)
#
#         if sim.met_concs is not None and met_tag is True: #update metabolism object if it's being simulated
#             sim.met_concs = {'cATP': self.ATP.c_mems,
#                 'cADP': self.ADP.c_mems,
#                 'cPi': self.Pi.c_mems}
#
#         for name in self.transporter_names:
#             obj = getattr(self, name)
#
#             obj.update_transporter(sim, cells, p)
#
#         for name in self.channel_names:
#             obj = getattr(self, name)
#
#             obj.update_channel(sim, cells, p)
#
#     def clear_cache(self):
#         """
#         Initializes or clears the time-storage vectors at the beginning of init and sim runs.
#
#         """
#
#         # get the name of the specific substance:
#         for name in self.molecule_names:
#             obj = getattr(self, name)
#
#             obj.c_mems_time = []
#             obj.c_cells_time = []
#             obj.c_env_time = []
#
#             if self.mit_enabled:
#                 obj.c_mit_time = []
#
#         for name in self.transporter_names:
#             obj = getattr(self, name)
#
#             obj.rate_time = []
#
#
#         for name in self.reaction_names:
#             obj = getattr(self, name)
#
#             obj.rate_time = []
#
#         if self.mit_enabled:
#             self.vmit_time = []
#             self.pH_mit_time = []
#
#         self.pH_cells_time = []
#         self.pH_env_time = []
#
#         self.chi_time = []
#
#     def write_data(self, sim, p):
#         """
#         Writes concentration data from a time-step to time-storage vectors.
#
#         """
#
#         # get the name of the specific substance:
#         for name in self.molecule_names:
#             obj = getattr(self, name)
#
#             obj.c_mems_time.append(obj.c_mems)
#             obj.c_cells_time.append(obj.c_cells)
#             obj.c_env_time.append(obj.c_env)
#
#             if self.mit_enabled:
#                 obj.c_mit_time.append(obj.c_mit)
#
#
#         # save rates of transporters and reactions
#         for name in self.transporter_names:
#             obj = getattr(self, name)
#
#             if obj.rate is not None:
#
#                 obj.rate_time.append(obj.rate)
#
#         for name in self.reaction_names:
#             obj = getattr(self, name)
#
#             if obj.rate is not None:
#
#                 obj.rate_time.append(obj.rate)
#
#
#         if self.mit_enabled:
#             self.vmit_time.append(self.mit.Vmit[:])
#             self.pH_mit_time.append(sim.pH_mit)
#
#         self.pH_cells_time.append(sim.pH_cell)
#         self.pH_env_time.append(sim.pH_env)
#
#         self.chi_time.append(self.chi)
#
#     def report(self, sim, p):
#         """
#         At the end of the simulation, tell user about mean, final concentrations of each molecule.
#
#         """
#
#         for name in self.molecule_names:
#
#             obj = getattr(self, name)
#
#             if self.mit_enabled:
#
#                 logs.log_info('Average concentration of ' + str(name) + ' in the mitochondria: ' +
#                           str(np.round(obj.c_mit.mean(), 4)) + ' mmol/L')
#
#             logs.log_info('Average concentration of ' + str(name) + ' in the cell: ' +
#                                            str(np.round(obj.c_cells.mean(), 4)) + ' mmol/L')
#
#             # logs.log_info('Average concentration of ' + str(name) + ' in the environment: ' +
#             #                               str(np.round(obj.c_env.mean(), 4)) + ' mmol/L')
#
#         if self.mit_enabled:
#             logs.log_info('Average Vmit: ' + str(np.round(1.0e3*self.mit.Vmit.mean(), 4)) + ' mV')
#
#             if 'ETC' in self.transporter_names:
#
#                 rate =  0.5*3600*1e15*self.mit.mit_vol.mean()*self.ETC.rate.mean()
#
#                 if 'ETC_ROS' in self.transporter_names:
#                     rate = rate + 3600*1e15*self.mit.mit_vol.mean()*self.ETC_ROS.rate.mean()
#
#                 logs.log_info('Average O2 consumption rate: ' + str(rate) + ' fmol/cell/hr')
#
#         # logs.log_info('Average pH in cell: ' + str(np.round(sim.pH_cell.mean(), 4)))
#         # logs.log_info('Average pH in env: ' + str(np.round(sim.pH_env.mean(), 4)))
#
#         # if self.mit_enabled:
#         #     logs.log_info('Average pH in mitochondria: ' + str(np.round(sim.pH_mit.mean(), 4)))
#
#         if self.chi.mean() != 0.0:
#
#             logs.log_info('Energy charge of cell ' + str(np.round(self.chi.mean(), 3)))
#
#     def export_all_data(self, sim, cells, p, message = 'for auxiliary molecules...'):
#
#         """
#
#         Exports concentration data from each molecule to a file for a single cell
#         (plot cell defined in params) as a function of time.
#
#         """
#         logs.log_info('Exporting raw data for ' + message)
#         # get the name of the specific substance:
#         for name in self.molecule_names:
#             obj = getattr(self, name)
#
#             obj.export_data(sim, cells, p, self.resultsPath)
#
#         if self.mit_enabled:
#             # FIXME we should also save vmit to a file? and pH and vm?
#             pass
#
#     def plot(self, sim, cells, p, message = 'for auxiliary molecules...'):
#         """
#         Creates plots for each molecule included in the simulation.
#
#         """
#
#         logs.log_info('Plotting 1D and 2D data for ' + message)
#
#         # get the name of the specific substance:
#         for name in self.molecule_names:
#             obj = getattr(self, name)
#
#             if p.plot_single_cell_graphs:
#             # create line graphs for the substance
#                 obj.plot_1D(sim, p, self.imagePath)
#
#             # create 2D maps for the substance
#             obj.plot_cells(sim, cells, p, self.imagePath)
#
#             # if there's a real environment, plot 2D concentration in the environment
#             if p.sim_ECM:
#
#                 obj.plot_env(sim, cells, p, self.imagePath)
#
#         #---------------cell everything plot---------------------------------------------
#         data_all1D = []
#         fig_all1D = plt.figure()
#         ax_all1D = plt.subplot(111)
#
#
#         # set up the color vector (for plotting complex line graphs)
#         maxlen = len(self.molecule_names)
#
#         lineplots_cm = plt.get_cmap(self.plot_cmap)
#         cNorm = colors.Normalize(vmin=0, vmax=maxlen)
#         c_names = cm.ScalarMappable(norm=cNorm, cmap=lineplots_cm)
#
#
#         for i, name in enumerate(self.molecule_names):
#             obj = getattr(self, name)
#
#             c_cells = [arr[p.plot_cell] for arr in obj.c_cells_time]
#
#             ax_all1D.plot(sim.time, c_cells, color = c_names.to_rgba(i), linewidth = 2.0, label=name)
#
#         legend = ax_all1D.legend(loc = 'upper right', shadow = False, frameon = False)
#
#         ax_all1D.set_xlabel('Time [s]')
#         ax_all1D.set_ylabel('Concentration [mmol/L]')
#         ax_all1D.set_title('Concentration of all substances in cell '  + str(p.plot_cell))
#
#         if p.autosave is True:
#             savename = self.imagePath + 'AllCellConcentrations_' + str(p.plot_cell) + '.png'
#             plt.savefig(savename, format='png', transparent=True)
#
#         if p.turn_all_plots_off is False:
#             plt.show(block=False)
#
#         #-------------environment everything plot-------------------------------------------------
#         data_all1D = []
#         fig_all1D = plt.figure()
#         ax_all1D = plt.subplot(111)
#
#         # get a random selection of our chosen colors in the length of our data set:
#         # c_names = np.random.choice(self.c_string, len(self.molecule_names))
#
#         for i, name in enumerate(self.molecule_names):
#             obj = getattr(self, name)
#
#             if p.sim_ECM is True:
#                 c_env = [arr[cells.map_cell2ecm][p.plot_cell] for arr in obj.c_env_time]
#
#             else:
#
#                 mem_i = cells.cell_to_mems[p.plot_cell][0]
#
#                 c_env = [arr[mem_i] for arr in obj.c_env_time]
#
#
#             ax_all1D.plot(sim.time, c_env, color = c_names.to_rgba(i), linewidth = 2.0, label=name)
#
#         legend = ax_all1D.legend(loc = 'upper right', shadow = False, frameon = False)
#
#         ax_all1D.set_xlabel('Time [s]')
#         ax_all1D.set_ylabel('Concentration [mmol/L]')
#         ax_all1D.set_title('Concentration of all substances in environment of cell ' + str(p.plot_cell))
#
#         if p.autosave is True:
#             savename = self.imagePath + 'AllEnvConcentrations_' + str(p.plot_cell) + '.png'
#             plt.savefig(savename, format='png', transparent=True)
#
#         if p.turn_all_plots_off is False:
#             plt.show(block=False)
#
#         #------------------------------------------------------------------------------------------
#         if self.mit_enabled:
#
#             # 1 D plot of mitochondrial voltage--------------------------------------------------------
#             vmit = [1e3*arr[p.plot_cell] for arr in self.vmit_time]
#
#             figVmit = plt.figure()
#             axVmit = plt.subplot(111)
#
#             axVmit.plot(sim.time, vmit)
#
#             axVmit.set_xlabel('Time [s]')
#             axVmit.set_ylabel('Vmit [mV]')
#             axVmit.set_title('Mitochondrial transmembrane voltage in cell: ' + str(p.plot_cell))
#
#             if p.autosave is True:
#                 savename = self.imagePath + 'Vmit_cell_' + str(p.plot_cell) + '.png'
#                 plt.savefig(savename, format='png', transparent=True)
#
#             if p.turn_all_plots_off is False:
#                 plt.show(block=False)
#
#             # 2D plot of mitochondrial voltage ---------------------------------------------------
#
#             fig, ax, cb = viz.plotPolyData(sim, cells, p,
#                 zdata=self.mit.Vmit*1e3, number_cells=p.enumerate_cells, clrmap=p.default_cm)
#
#             ax.set_title('Final Mitochondrial Transmembrane Voltage')
#             ax.set_xlabel('Spatial distance [um]')
#             ax.set_ylabel('Spatial distance [um]')
#             cb.set_label('Vmit [mV]')
#
#             if p.autosave is True:
#                 savename = self.imagePath + '2DVmit.png'
#                 plt.savefig(savename, format='png', transparent=True)
#
#             if p.turn_all_plots_off is False:
#                 plt.show(block=False)
#
#             # plot of all substances in the mitochondria:----------------------------------------------
#
#             data_all1D = []
#             fig_all1D = plt.figure()
#             ax_all1D = plt.subplot(111)
#
#             # get a random selection of our chosen colors in the length of our data set:
#             # c_names = np.random.choice(self.c_string, len(self.molecule_names))
#
#             for i, name in enumerate(self.molecule_names):
#                 obj = getattr(self, name)
#
#                 c_mit = [arr[p.plot_cell] for arr in obj.c_mit_time]
#
#                 ax_all1D.plot(sim.time, c_mit, color = c_names.to_rgba(i), linewidth=2.0, label=name)
#
#             legend = ax_all1D.legend(loc='upper right', shadow=False, frameon=False)
#
#             ax_all1D.set_xlabel('Time [s]')
#             ax_all1D.set_ylabel('Concentration [mmol/L]')
#             ax_all1D.set_title('Substances in mitochondria of cell ' + str(p.plot_cell))
#
#             if p.autosave is True:
#                 savename = self.imagePath + 'AllMitConcentrations_' + str(p.plot_cell) + '.png'
#                 plt.savefig(savename, format='png', transparent=True)
#
#             if p.turn_all_plots_off is False:
#                 plt.show(block=False)
#         #------pH plot------------------------------------------------------------------------------
#
#         # 1 D plot of pH in cell, env and mit ------------------------------------------------------
#         pHcell = [arr[p.plot_cell] for arr in self.pH_cells_time]
#
#         if p.sim_ECM:
#             pHenv = [arr[cells.map_cell2ecm][p.plot_cell] for arr in self.pH_env_time]
#
#         else:
#             avPh = [np.dot(cells.M_sum_mems, arr)/cells.num_mems for arr in self.pH_env_time]
#             pHenv = [arr[p.plot_cell] for arr in avPh]
#
#         if self.mit_enabled:
#             pHmit = [arr[p.plot_cell] for arr in self.pH_mit_time]
#
#         else:
#             pHmit = np.zeros(len(sim.time))
#
#         figpH = plt.figure()
#         axpH = plt.subplot(111)
#
#         axpH.plot(sim.time, pHcell, label = 'cell')
#         axpH.plot(sim.time, pHmit, label='mitochondria')
#         axpH.plot(sim.time, pHenv, label='env')
#
#         axpH.set_xlabel('Time [s]')
#         axpH.set_ylabel('pH')
#         axpH.set_title('pH in/near cell : ' + str(p.plot_cell))
#
#         if p.autosave is True:
#             savename = self.imagePath + 'pH_' + str(p.plot_cell) + '.png'
#             plt.savefig(savename, format='png', transparent=True)
#
#         if p.turn_all_plots_off is False:
#             plt.show(block=False)
#
#         #-------Reaction rate plot and data export----------------------------------------
#
#         if len(self.reaction_names):
#
#             # create a suite of single reaction line plots:
#             for i, name in enumerate(self.reaction_names):
#                 # get the reaction object field
#                 obj = getattr(self, name)
#
#                 # make a 1D plot of this reaction rate:
#                 obj.plot_1D(sim, cells, p, self.imagePath)
#
#             # now create the "everything" plot and export data to file:
#
#             react_dataM = []
#             react_header = 'Time [s], '
#
#             react_dataM.append(sim.time)
#
#             data_all1D = []
#             fig_all1D = plt.figure()
#             ax_all1D = plt.subplot(111)
#
#             # set up the color vector (for plotting complex line graphs)
#             maxlen = len(self.reaction_names)
#
#             lineplots_cm = plt.get_cmap(self.plot_cmap)
#             cNorm = colors.Normalize(vmin=0, vmax=maxlen)
#             c_names = cm.ScalarMappable(norm=cNorm, cmap=lineplots_cm)
#
#             for i, name in enumerate(self.reaction_names):
#                 # get the reaction object field
#                 obj = getattr(self, name)
#
#                 if len(obj.rate_time) > 0:
#
#                     r_rate = [arr[p.plot_cell] for arr in obj.rate_time]
#
#                     ax_all1D.plot(sim.time, r_rate, color = c_names.to_rgba(i), linewidth=2.0, label=name)
#
#                     react_dataM.append(r_rate)
#                     react_header = react_header + name + ' [mM/s]'+ ','
#
#             legend = ax_all1D.legend(loc='upper right', shadow=False, frameon=False)
#
#             ax_all1D.set_xlabel('Time [s]')
#             ax_all1D.set_ylabel('Rate [mM/s]')
#             ax_all1D.set_title('Reaction rates in cell ' + str(p.plot_cell))
#
#             if p.autosave is True:
#                 savename = self.imagePath + 'AllReactionRates_' + str(p.plot_cell) + '.png'
#                 plt.savefig(savename, format='png', transparent=True)
#
#             if p.turn_all_plots_off is False:
#                 plt.show(block=False)
#
#             react_dataM = np.asarray(react_dataM)
#
#             saveName = 'AllReactionRatesData_' + str(p.plot_cell) + '.csv'
#
#             saveDataReact = os.path.join(self.resultsPath, saveName)
#
#             np.savetxt(saveDataReact, react_dataM.T, delimiter=',', header=react_header)
#
#         #---Transporter rate plot and data export ------------------------------------------------------
#
#         if len(self.transporter_names):
#
#             transp_dataM = []
#             transp_header = 'Time [s], '
#
#             transp_dataM.append(sim.time)
#
#             # set up the color vector (for plotting complex line graphs)
#             maxlen = len(self.transporter_names)
#
#             lineplots_cm = plt.get_cmap(self.plot_cmap)
#             cNorm = colors.Normalize(vmin=0, vmax=maxlen)
#             c_names = cm.ScalarMappable(norm=cNorm, cmap=lineplots_cm)
#
#             for i, name in enumerate(self.transporter_names):
#                 obj = getattr(self, name)
#
#                 # make a 1D plot of this reaction rate:
#                 obj.plot_1D(sim, cells, p, self.imagePath)
#
#                 if len(obj.rate_time) > 0:
#
#                     # check the data structure size for this transporter:
#                     if len(obj.rate_time[0]) == sim.cdl:
#
#                         t_rate = [arr[p.plot_cell] for arr in obj.rate_time]
#
#                     elif len(obj.rate_time[0]) == sim.mdl:
#                         mem_i = cells.cell_to_mems[p.plot_cell][0]
#                         t_rate = [arr[mem_i] for arr in obj.rate_time]
#
#                     else:
#                         t_rate = np.zeros(len(sim.time))
#
#                     data_all1D = []
#                     fig_all1D = plt.figure()
#                     ax_all1D = plt.subplot(111)
#
#                     ax_all1D.plot(sim.time, t_rate, color = c_names.to_rgba(i), linewidth=2.0, label=name)
#
#                     transp_dataM.append(t_rate)
#                     transp_header = transp_header + name + ' [mM/s]' + ','
#
#             legend = ax_all1D.legend(loc='upper right', shadow=False, frameon=False)
#
#             ax_all1D.set_xlabel('Time [s]')
#             ax_all1D.set_ylabel('Rate [mM/s]')
#             ax_all1D.set_title('Transporter rates in cell ' + str(p.plot_cell))
#
#             if p.autosave is True:
#                 savename = self.imagePath + 'AllTransporterRates_' + str(p.plot_cell) + '.png'
#                 plt.savefig(savename, format='png', transparent=True)
#
#             if p.turn_all_plots_off is False:
#                 plt.show(block=False)
#
#             saveName = 'AllTransporterRatesData_' + str(p.plot_cell) + '.csv'
#
#             saveDataTransp = os.path.join(self.resultsPath, saveName)
#
#             transp_dataM = np.asarray(transp_dataM)
#
#             np.savetxt(saveDataTransp, transp_dataM.T, delimiter=',', header=transp_header)
#
#         # energy charge plots:----------------------------------------------------------
#         # 1 D plot of mitochondrial voltage--------------------------------------------------------
#         chio = [arr[p.plot_cell] for arr in self.chi_time]
#
#         figChi = plt.figure()
#         axChi = plt.subplot(111)
#
#         axChi.plot(sim.time, chio)
#
#         axChi.set_xlabel('Time [s]')
#         axChi.set_ylabel('Energy charge')
#         axChi.set_title('Energy charge in cell: ' + str(p.plot_cell))
#
#         if p.autosave is True:
#             savename = self.imagePath + 'EnergyCharge_cell_' + str(p.plot_cell) + '.png'
#             plt.savefig(savename, format='png', transparent=True)
#
#         if p.turn_all_plots_off is False:
#             plt.show(block=False)
#
#         #---2D plot--------------------------------------------------------------------
#
#
#
#         fig, ax, cb = viz.plotPolyData(sim, cells, p,
#             zdata=self.chi, number_cells=p.enumerate_cells, clrmap=p.default_cm)
#
#         ax.set_title('Final Energy Charge of Cell')
#         ax.set_xlabel('Spatial distance [um]')
#         ax.set_ylabel('Spatial distance [um]')
#         cb.set_label('Energy Charge')
#
#         if p.autosave is True:
#             savename = self.imagePath + '2DEnergyCharge.png'
#             plt.savefig(savename, format='png', transparent=True)
#
#         if p.turn_all_plots_off is False:
#             plt.show(block=False)
#
#     def anim(self, sim, cells, p, message = 'for auxiliary molecules...'):
#         """
#         Animates 2D data for each molecule in the simulation.
#
#         """
#
#         logs.log_info('Animating data for ' + message)
#         # get the name of the specific substance:
#         for name in self.molecule_names:
#
#             obj = getattr(self, name)
#
#             if p.anim.is_after_sim and obj.make_ani is True:
#
#                 # create 2D animations for the substance in cells
#                 obj.anim_cells(sim, cells, p)
#
#                 # create 2D animations for the substance in the environment
#                 if p.sim_ECM:
#
#                     obj.anim_env(sim, cells, p)
#
# class Molecule(object):
#
#     def __init__(self, sim, cells, p):
#
#         self.dummy_dyna = TissueHandler(sim, cells, p)
#         self.dummy_dyna.tissueProfiles(sim, cells, p)  # initialize all tissue profiles
#
#
#
#         # Set all fields to None -- these will be dynamically set by MasterOfMolecules
#
#         self.c_cello = None
#         self.c_memo = None
#         self.c_env = None
#         self.z = None
#         self.Dm = None
#         self.Do = None
#         self.c_bound = None
#         self.Kgd = None
#         self.use_pumping = None
#         self.pumps_use_ATP = None
#         self.pump_to_cell = None
#         self.pump_max_val = None
#         self.pump_Km = None
#         self.use_gating_ligand = None
#         self.gating_extracell = None
#         self.gating_max_val = None
#         self.gating_Hill_K = None
#         self.gating_Hill_n = None
#         self.gating_ion = None
#
#         self.change_at_bounds = None
#         self.change_bounds_start = None
#         self.change_bounds_end = None
#         self.change_bounds_rate = None
#         self.change_bounds_target = None
#         self.make_plots = None
#         self.make_ani = None
#         self.plot_autoscale = None
#         self.plot_max = None
#         self.plot_min = None
#
#         self.r_production = None
#         self.r_decay = None
#         self.n_production = None
#
#         self.growth_activators_list = None
#         self.growth_activators_k = None
#         self.growth_activators_Km = None
#         self.growth_activators_n = None
#         self.growth_inhibitors_list = None
#         self.growth_inhibitors_k = None
#         self.growth_inhibitors_Km = None
#         self.growth_inhibitors_n = None
#
#     def transport(self, sim, cells, p):
#         """
#         Transports the molecule across the membrane,
#         through gap junctions, and if p.sim_ECM is true,
#         through extracellular spaces and the environment.
#
#         """
#
#
#
#         self.c_mems, self.c_env, _, _, _, _ = stb.molecule_mover(sim,
#                                                                 self.c_mems,
#                                                                 self.c_env,
#                                                                 cells, p,
#                                                                 z=self.z,
#                                                                 Dm = self.Dm,
#                                                                 Do = self.Do,
#                                                                 c_bound = self.c_bound,
#                                                                 ignoreECM = self.ignore_ECM_pump,
#                                                                 smoothECM = p.smooth_concs,
#                                                                 ignoreTJ = self.ignoreTJ,
#                                                                 ignoreGJ = self.ignoreGJ)
#
#     def updateC(self, flux, sim, cells, p):
#         """
#
#         General updater for a flux defined on membranes and updating concentrations in
#         cells and environment.
#
#         """
#         self.c_mems, self.c_env = stb.update_Co(sim, self.c_mems, self.c_cells, flux, cells, p, ignoreECM=True)
#
#     def updateIntra(self, sim, sim_metabo, cells, p):
#
#         self.c_mems, self.c_cells, _ = stb.update_intra(sim, cells, self.c_mems, self.c_cells, self.Do, self.z, p)
#
#         if self.mit_enabled:
#
#             IdCM = np.ones(sim.cdl)
#
#             f_ED = stb.electroflux(self.c_cells, self.c_mit, self.Dm*IdCM, p.tm*IdCM, self.z*IdCM,
#                 sim_metabo.mit.Vmit, sim.T, p, rho=1)
#
#             # update with flux
#             self.c_cells = self.c_cells - f_ED*(sim_metabo.mit.mit_sa/cells.cell_vol)*p.dt
#             self.c_mit = self.c_mit + f_ED*(sim_metabo.mit.mit_sa/sim_metabo.mit.mit_vol)*p.dt
#
#     def pump(self, sim, cells, p):
#
#         """
#         Defines a generic active transport pump that can be used to move
#         a general molecule (such as serotonin or glutamate)
#         into or out of the cell by active transport.
#
#         Works on the basic premise of enzymatic pumps defined elsewhere:
#
#         pump_out is True:
#
#         cX_cell + cATP  -------> cX_env + cADP + cPi
#
#         pump_out is False:
#
#         cX_env + cATP  <------- cX_cell + cADP + cPi
#
#         """
#
#         if self.use_pumping:
#
#             if self.pumps_use_ATP:
#
#                 if p.metabolism_enabled:
#                     met_vect = sim.met_concs
#                 else:
#                     met_vect = None
#
#                 self.c_mems, self.c_env, flux = stb.molecule_pump(sim, self.c_mems, self.c_env,
#                                                                      cells, p, Df=self.Do, z=self.z,
#                                                                      pump_into_cell=self.pump_to_cell,
#                                                                      alpha_max=self.pump_max_val, Km_X=self.pump_Km,
#                                                                      Km_ATP=1.0, met = met_vect, ignoreECM = self.ignore_ECM_pump)
#                 if p.metabolism_enabled:
#                     # update ATP concentrations after pump action:
#                     sim.metabo.update_ATP(flux, sim, cells, p)
#
#             else:
#
#                 self.c_mems, self.c_env, flux = stb.molecule_transporter(sim, self.c_mems, self.c_env,
#                     cells, p, Df=self.Do, z=self.z, pump_into_cell=self.pump_to_cell, alpha_max=self.pump_max_val,
#                     Km_X=self.pump_Km, Keq= 1.0, ignoreECM = self.ignore_ECM_pump)
#
#     def gating(self, sim, sim_metabo, cells, p):
#         """
#         Uses the molecule concentration to open an ion channel in the cell membranes.
#
#         """
#
#         # update membrane permeability if dye targets an ion channel:
#         if self.use_gating_ligand:
#
#             # calculate any activators and/or inhibitor effects:
#             activator_alpha, inhibitor_alpha = get_influencers(sim, sim_metabo, self.activators_list,
#                 self.activators_Km, self.activators_n, self.inhibitors_list,
#                 self.inhibitors_Km, self.inhibitors_n, reaction_zone='mems')
#
#
#             if self.gating_extracell is False:
#
#                 for ion_tag in self.gating_ion:
#
#                     Dm_mod_mol = sim.rho_channel*self.gating_max_val*tb.hill(self.c_mems,
#                                                                             self.gating_Hill_K,self.gating_Hill_n)
#
#                     sim.Dm_morpho[ion_tag] = sim.rho_channel*Dm_mod_mol*activator_alpha*inhibitor_alpha
#
#             elif self.gating_extracell is True and p.sim_ECM is True:
#
#                 for ion_tag in self.gating_ion:
#
#                     Dm_mod_mol = self.gating_max_val*tb.hill(self.c_env,self.gating_Hill_K,self.gating_Hill_n)
#
#                     sim.Dm_morpho[ion_tag] = (activator_alpha*inhibitor_alpha*sim.rho_channel*
#                                               Dm_mod_mol[cells.map_mem2ecm])
#
#     def init_growth(self,cells, p):
#
#         if self.growth_profiles_list is not None and self.growth_profiles_list != 'all':
#
#             self.growth_targets_cell = []
#             self.growth_targets_mem = []
#
#             for profile in self.growth_profiles_list:
#                 targets_cell = self.dummy_dyna.cell_target_inds[profile]
#                 self.growth_targets_cell.extend(targets_cell)
#
#                 targets_mem = self.dummy_dyna.tissue_target_inds[profile]
#                 self.growth_targets_mem.extend(targets_mem)
#
#
#         elif self.growth_profiles_list is None or self.growth_profiles_list == 'all':
#
#             self.growth_targets_cell = cells.cell_i
#             self.growth_targets_mem = cells.mem_i
#
#     def growth_and_decay(self, super_self, sim, cells, p):
#         """
#         Grows and/or decays the molecule concentration in the cell cytosol using a simple rate equation
#         representing saturating autocatalytic production/decay via Michaelis-Menten kinetics.
#
#         """
#
#         cc = self.c_cells/self.Kgd
#
#         activator_alpha, inhibitor_alpha = get_influencers_grn(sim, super_self, self.growth_activators_list,
#             self.growth_activators_k, self.growth_activators_Km, self.growth_activators_n,
#             self.growth_inhibitors_list, self.growth_inhibitors_k, self.growth_inhibitors_Km,
#             self.growth_inhibitors_n, reaction_zone='cell')
#
#         delta_cells = self.growth_mod_function_cells*self.r_production*inhibitor_alpha*activator_alpha - self.r_decay*cc
#
#
#         self.c_cells[self.growth_targets_cell] = self.c_cells[self.growth_targets_cell] + \
#                                                  delta_cells[self.growth_targets_cell]*p.dt
#
#         self.c_mems[self.growth_targets_mem] = self.c_mems[self.growth_targets_mem] + \
#                                                delta_cells[cells.mem_to_cells][self.growth_targets_mem]*p.dt
#
#         # make sure the concs inside the cell are evenly mixed after production/decay:
#         # self.updateIntra(sim, cells, p)
#
#     def remove_cells(self, target_inds_cell, target_inds_mem, sim, cells, p):
#         """
#         During a cutting event, removes the right cells from the simulation network,
#         while preserving additional information.
#
#         """
#
#         # remove cells from the cell concentration list:
#         ccells2 = np.delete(self.c_cells, target_inds_cell)
#         # reassign the new data vector to the object:
#         self.c_cells = ccells2[:]
#
#         # remove cells from the mems concentration list:
#         cmems2 = np.delete(self.c_mems, target_inds_mem)
#         # reassign the new data vector to the object:
#         self.c_mems = cmems2[:]
#
#         if self.simple_growth is True and self.growth_mod_function_cells != 1:
#
#             gmfc = np.delete(self.growth_mod_function_cells, target_inds_cell)
#             self.growth_mod_function_cells = gmfc[:]
#
#             if len(self.growth_mod_function_cells) == 0:
#                 self.growth_mod_function_cells = 1
#
#         self.dummy_dyna.tissueProfiles(sim, cells, p)  # re-initialize all tissue profiles
#         self.init_growth(cells, p)
#
#         if p.sim_ECM is False:
#
#             cenv2 = np.delete(self.c_env, target_inds_mem)
#             self.c_env = cenv2[:]
#
#
#         if self.mit_enabled:
#             # remove cells from the cell concentration list:
#             cmit2 = np.delete(self.c_mit, target_inds_cell)
#             # reassign the new data vector to the object:
#             self.c_mit = cmit2[:]
#
#     def update_boundary(self, t, p):
#         """
#         Run a dynamic event in the sim, which alters concentration at the global boundary.
#
#         t:          simulation world time
#         p:          parameters instance
#
#         """
#
#         if self.change_at_bounds:
#
#             effector_MorphEnv = tb.pulse(t,self.change_bounds_start,self.change_bounds_end,self.change_bounds_rate)
#
#             if p.sim_ECM is False:
#                 self.c_env[:] = self.conc_MorphEnv*effector_MorphEnv + self.c_envo*(1-effector_MorphEnv)
#
#             elif p.sim_ECM is True: # simulate addition of counter salt to maintain charge neutrality:
#                 self.c_bound = self.conc_MorphEnv*effector_MorphEnv + self.c_envo*(1-effector_MorphEnv)
#
#     def export_data(self, sim, cells, p, savePath):
#
#         saveName = 'ExportData_' + self.name + '_' + str(p.plot_cell) + '.csv'
#
#         saveData = os.path.join(savePath, saveName)
#
#         ci = p.plot_cell  # index of cell to get time-dependent data for
#
#         # create the header, first entry will be time:
#         headr = 'time_s' + ','
#
#         ccell = [arr[ci] for arr in self.c_cells_time]
#
#         headr = headr + 'Cell_Conc_' + self.name + '_mmol/L' + ','
#
#         if p.sim_ECM is True:
#
#             cenv = [obj_cenv[cells.map_cell2ecm][ci] for obj_cenv in self.c_env_time]
#
#         else:
#
#             cenv = [np.dot(cells.M_sum_mems, obj_cenv) / cells.num_mems for obj_cenv in self.c_env_time]
#
#         headr = headr + 'Env_Conc_' + self.name + '_mmol/L' + ','
#
#         if self.mit_enabled:
#
#             cmit = [arr[ci] for arr in self.c_mit_time]
#
#             headr = headr + 'Mit_Conc_' + self.name + '_mmol/L' + ','
#
#             cmit = np.asarray(cmit)
#
#         time = np.asarray(sim.time)
#         ccell = np.asarray(ccell)
#         cenv = np.asarray(cenv)
#
#         if self.mit_enabled is False:
#             dataM = np.column_stack((time, ccell, cenv))
#
#         else:
#             dataM = np.column_stack((time, ccell, cenv, cmit))
#
#         np.savetxt(saveData, dataM, delimiter=',', header=headr)
#
#     def plot_1D(self, sim, p, saveImagePath):
#         """
#         Create 1D plot of concentration in cell and environment for a single cell (params plot cell)
#         as a function of simulation time.
#
#         """
#
#         c_cells = [arr[p.plot_cell] for arr in self.c_cells_time]
#         fig = plt.figure()
#         ax = plt.subplot(111)
#         ax.plot(sim.time, c_cells)
#         ax.set_xlabel('Time [s]')
#         ax.set_ylabel('Concentration [mmol/L]')
#         ax.set_title('Concentration of ' + self.name + ' in cell ' + str(p.plot_cell))
#
#         if p.autosave is True:
#             savename = saveImagePath + 'CellConcentration_' + self.name + '_' + str(p.plot_cell) + '.png'
#             plt.savefig(savename, format='png', transparent=True)
#
#         if p.turn_all_plots_off is False:
#             plt.show(block=False)
#
#         if self.mit_enabled:
#
#             c_mit = [arr[p.plot_cell] for arr in self.c_mit_time]
#             fig = plt.figure()
#             ax = plt.subplot(111)
#             ax.plot(sim.time, c_mit)
#             ax.set_xlabel('Time [s]')
#             ax.set_ylabel('Concentration [mmol/L]')
#             ax.set_title('Mitochondrial concentration of ' + self.name + ' in cell ' + str(p.plot_cell))
#
#             if p.autosave is True:
#                 savename = saveImagePath + 'MitConcentration_' + self.name + '_' + str(p.plot_cell) + '.png'
#                 plt.savefig(savename, format='png', transparent=True)
#
#             if p.turn_all_plots_off is False:
#                 plt.show(block=False)
#
#     def plot_cells(self, sim, cells, p, saveImagePath):
#         """
#         Create 2D plot of cell concentrations.
#
#         """
#
#         fig, ax, cb = viz.plotPrettyPolyData(self.c_mems,
#             sim, cells, p,
#             number_cells=p.enumerate_cells,
#             clrAutoscale=self.plot_autoscale,
#             clrMin=self.plot_min,
#             clrMax=self.plot_max,
#             clrmap=p.default_cm)
#
#         ax.set_title('Final ' + self.name + ' Concentration in Cells')
#         ax.set_xlabel('Spatial distance [um]')
#         ax.set_ylabel('Spatial distance [um]')
#         cb.set_label('Concentration mmol/L')
#
#         if p.autosave is True:
#             savename = saveImagePath + '2Dcell_conc_' + self.name + '.png'
#             plt.savefig(savename,format='png', transparent=True)
#
#         if p.turn_all_plots_off is False:
#             plt.show(block=False)
#
#         # mitochondrial plots
#         if self.mit_enabled:
#
#             fig, ax, cb = viz.plotPolyData(sim, cells, p, zdata=self.c_mit,
#                 number_cells=p.enumerate_cells,
#                 clrAutoscale=self.plot_autoscale,
#                 clrMin=self.plot_min,
#                 clrMax=self.plot_max,
#                 clrmap=p.default_cm)
#
#             ax.set_title('Final ' + self.name + ' Concentration in Mitochondria')
#             ax.set_xlabel('Spatial distance [um]')
#             ax.set_ylabel('Spatial distance [um]')
#             cb.set_label('Concentration mmol/L')
#
#             if p.autosave is True:
#                 savename = saveImagePath + '2D_mit_conc_' + self.name + '.png'
#                 plt.savefig(savename, format='png', transparent=True)
#
#             if p.turn_all_plots_off is False:
#                 plt.show(block=False)
#
#     def plot_env(self, sim, cells, p, saveImagePath):
#         """
#         Create 2D plot of environmental concentration.
#
#         """
#
#
#         fig = plt.figure()
#         ax = plt.subplot(111)
#
#         dyeEnv = (self.c_env).reshape(cells.X.shape)
#
#         xmin = cells.xmin*p.um
#         xmax = cells.xmax*p.um
#         ymin = cells.ymin*p.um
#         ymax = cells.ymax*p.um
#
#         bkgPlot = ax.imshow(dyeEnv,origin='lower',extent=[xmin,xmax,ymin,ymax], cmap=p.default_cm)
#
#
#         if self.plot_autoscale is False:
#             bkgPlot.set_clim(self.plot_min, self.plot_max)
#
#         cb = fig.colorbar(bkgPlot)
#
#         ax.axis('equal')
#
#         ax.axis([xmin,xmax,ymin,ymax])
#
#         ax.set_title('Final ' + self.name + ' Concentration in Environment')
#         ax.set_xlabel('Spatial distance [um]')
#         ax.set_ylabel('Spatial distance [um]')
#         cb.set_label('Concentration mmol/L')
#
#         if p.autosave is True:
#             savename = saveImagePath + '2Denv_conc_' + self.name + '.png'
#             plt.savefig(savename,format='png',dpi = 300.0, transparent=True)
#
#         if p.turn_all_plots_off is False:
#             plt.show(block=False)
#
#     def anim_cells(self, sim, cells, p):
#         """
#         Create 2D animation of cell concentration.
#         """
#
#         AnimCellsMembranesData(
#             sim=sim, cells=cells, p=p,
#             times_membranes_midpoint_data=[arr for arr in self.c_mems_time],
#             label=self.name + 'cells',
#             figure_title='Cytosolic ' + self.name,
#             colorbar_title='Concentration [mmol/L]',
#             is_color_autoscaled=self.plot_autoscale,
#             color_min=self.plot_min,
#             color_max=self.plot_max)
#
#     def anim_env(self, sim, cells, p):
#         """
#         Create 2D animation of env concentration.
#         """
#
#         env_time_series = [
#             env.reshape(cells.X.shape) for env in self.c_env_time]
#         AnimEnvTimeSeries(
#             sim=sim, cells=cells, p=p,
#             time_series=env_time_series,
#             label=self.name + '_env',
#             figure_title='Environmental ' + self.name,
#             colorbar_title='Concentration [mmol/L]',
#             is_color_autoscaled=self.plot_autoscale,
#             color_min=self.plot_min,
#             color_max=self.plot_max)
#
# class Reaction(object):
#
#     def __init__(self, sim, cells, p):
#
#         self.dummy_dyna = TissueHandler(sim, cells, p)
#         self.dummy_dyna.tissueProfiles(sim, cells, p)  # initialize all tissue profiles
#
#         # pre-populate the object with fields that will be assigned by MasterOfMolecules
#
#         self.name = None
#         self.reactants_list = None
#         self.reactants_coeff = None
#         self.products_list = None
#         self.products_coeff = None
#         self.Km_reactants_list = None
#         self.Km_products_list = None
#         self.vmax = None
#         self.delta_Go = None
#
#         self.reaction_zone = None
#
#         self.c_reactants = []  # concentrations of reactions, defined on cell-centres
#         self.z_reactants = []  # charge of reactants
#         self.inds_react = []  # indices to each reactant, to keep their order
#         self.c_products = []  # concentrations of reactions, defined on cell-centres
#         self.z_products = []  # charge of products
#         self.inds_prod = []  # indices to each reactant, to keep their order
#
#         self.product_source_object = []  # object from which product concentrations come from
#         self.product_source_type = []    # type of product concentrations sourced from object
#         self.reactant_source_object = []  # object from which reactants concentrations come from
#         self.reactant_source_type = []  # type of reactant concentrations sourced from object
#
#         # activator and inhibitors of the reaction
#         self.reaction_activators_list = None
#         self.reaction_activators_Km = None
#         self.reaction_activators_n = None
#         self.reaction_inhibitors_list = None
#         self.reaction_inhibitors_Km = None
#         self.reaction_inhibitors_n = None
#
#     def get_reactants(self, sim, sim_metabo, reactant_type_self, reactant_type_sim):
#
#         """
#         Get updated concentrations direct from sources for chemical reaction's reactants.
#
#         sim:                    An instance of simulator
#         sim_metabo:             Sim instance of reaction block (sim.metabo or sim.reacto for general reactions)
#         reactant_type_self:     Data type to retrieve: c_cells, c_mems, c_mit
#         reactant_type_sim:      Data typpe to retrieve: cc_cells, cc_mems, cc_mit
#
#         """
#
#         self.c_reactants = []
#
#         for reactant_name in self.reactants_list:
#
#             label = 'i' + reactant_name
#             ion_check = getattr(sim, label, None)
#
#             if ion_check is None:
#
#                 try:
#                     obj_reactant = getattr(sim_metabo, reactant_name)
#                     c_react = getattr(obj_reactant, reactant_type_self)
#                     self.c_reactants.append(c_react)
#
#                 except KeyError:
#
#                     raise BetseSimConfigException('Name of product is not a defined chemical, or is not'
#                                                    'an ion currently included in the ion profile being used.'
#                                                    'Please check biomolecule definitions and ion profile'
#                                                    'settings of your config(s) and try again.')
#
#             else:
#
#                 # define the reactant as the ion concentration from the cell concentrations object in sim:
#                 sim_conco = getattr(sim, reactant_type_sim)
#
#                 sim_conc = sim_conco[ion_check]
#
#                 self.c_reactants.append(sim_conc)
#
#     def get_products(self, sim, sim_metabo, product_type_self, product_type_sim):
#
#         """
#         Get updated concentrations direct from sources for chemical reaction's products.
#
#         sim:                    An instance of simulator
#         sim_metabo:             Sim instance of reaction block (sim.metabo or sim.reacto for general reactions)
#         product_type_self:     Data type to retrieve: c_cells, c_mems, c_mit
#         product_type_sim:      Data typpe to retrieve: cc_cells, cc_mems, cc_mit
#
#         """
#
#         self.c_products = []
#
#         for product_name in self.products_list:
#
#             label = 'i' + product_name
#             ion_check = getattr(sim, label, None)
#
#             if ion_check is None:
#
#                 try:
#                     obj_prod = getattr(sim_metabo, product_name)
#                     c_prod = getattr(obj_prod, product_type_self)
#                     self.c_products.append(c_prod)
#
#                 except KeyError:
#
#                     raise BetseSimConfigException('Name of product is not a defined chemical, or is not'
#                                                    'an ion currently included in the ion profile being used.'
#                                                    'Please check biomolecule definitions and ion profile'
#                                                    'settings of your config(s) and try again.')
#
#             else:
#
#                 # define the reactant as the ion concentration from the cell concentrations object in sim:
#                 sim_conco = getattr(sim, product_type_sim)
#
#                 sim_conc = sim_conco[ion_check]
#
#                 self.c_products.append(sim_conc)
#
#     def set_reactant_c(self, deltac, sim, sim_metabo, reactant_type_self, reactant_type_sim):
#
#         for i, reactant_name in enumerate(self.reactants_list):
#
#             conc = self.c_reactants[i] - deltac*self.reactants_coeff[i]
#
#             label = 'i' + reactant_name
#             ion_check = getattr(sim, label, None)
#
#             if ion_check is None:
#
#                 obj_reactant = getattr(sim_metabo, reactant_name)
#                 setattr(obj_reactant, reactant_type_self, conc)
#
#             else:
#                 # define the reactant as the ion concentration from the cell concentrations object in sim:
#                 sim_conc = getattr(sim, reactant_type_sim)
#
#                 sim_conc[ion_check] = conc
#
#     def set_product_c(self, deltac, sim, sim_metabo, product_type_self, product_type_sim):
#
#         for i, product_name in enumerate(self.products_list):
#
#             conc = self.c_products[i] + deltac*self.products_coeff[i]
#
#             label = 'i' + product_name
#             ion_check = getattr(sim, label, None)
#
#             if ion_check is None:
#
#                 obj_reactant = getattr(sim_metabo, product_name)
#                 setattr(obj_reactant, product_type_self, conc)
#
#             else:
#                 # define the reactant as the ion concentration from the cell concentrations object in sim:
#                 sim_conc = getattr(sim, product_type_sim)
#
#                 sim_conc[ion_check] = conc
#
#     def compute_reaction(self, sim, sim_metabo, cells, p):
#
#
#         # get up-to-date concentration data for the reaction:
#         if self.reaction_zone == 'cell':
#
#             self.get_reactants(sim, sim_metabo, 'c_cells', 'cc_cells')
#             self.get_products(sim, sim_metabo, 'c_cells', 'cc_cells')
#
#         elif self.reaction_zone == 'mitochondria' and self.mit_enabled is True:
#
#             self.get_reactants(sim, sim_metabo, 'c_mit', 'cc_mit')
#             self.get_products(sim, sim_metabo, 'c_mit', 'cc_mit')
#
#         # if reaction is reversible, calculate reaction quotient Q, the equilibrium constant, and backwards rate
#         if self.delta_Go is not None:
#
#             # define the reaction equilibrium coefficient:
#             Keqm = np.exp(-self.delta_Go/(p.R*sim.T))
#
#             # calculate the MM rate coefficient for the backwards reaction direction and Q:
#             backwards_term = []
#             Q_deno_list = []
#             Q_numo_list = []
#
#             for i, prod in enumerate(self.c_products):
#                 # calculate a factor in the backwards rate term:
#                 coeff = self.products_coeff[i]
#                 Km = self.Km_products_list[i]
#                 cs = (prod / Km) ** coeff
#                 term = cs / (1 + cs)
#                 backwards_term.append(term)
#
#                 # calculate a factor in the reaction quotient numerator term:
#                 ci = prod**coeff
#                 Q_numo_list.append(ci)
#
#             backwards_term = np.asarray(backwards_term)
#
#             backwards_rate = self.vmax * np.prod(backwards_term, axis=0)
#
#             for j, react in enumerate(self.c_reactants):
#                 # calculate a factor in the reaction quotient denominator term:
#                 coeff = self.reactants_coeff[j]
#                 cj = react**coeff
#                 Q_deno_list.append(cj)
#
#             # get the numerator and denomenator of the reaction quotient:
#             Q_deno = np.prod(Q_deno_list, axis=0)
#             Q_numo = np.prod(Q_numo_list, axis=0)
#
#             inds_neg = (Q_deno == 0).nonzero()
#             Q_deno[inds_neg] = 1.0e-15  # ensure no division by zero
#
#             # finally, *the* reaction quotient:
#             Q = Q_numo/Q_deno
#
#         else:
#
#             Q  = 0
#             backwards_rate = 0
#             Keqm = 1
#
#         # Next, calculate the MM rate coefficient for the forward reaction direction:
#         forward_term = []
#
#         for i, react in enumerate(self.c_reactants):
#
#             coeff = self.reactants_coeff[i]
#             Km = self.Km_reactants_list[i]
#             cs = (react/Km)**coeff
#
#             term = cs/(1 + cs)
#
#             forward_term.append(term)
#
#         forward_term = np.asarray(forward_term)
#
#         forward_rate = self.vmax*np.prod(forward_term, axis=0)
#
#         reaction_rate = forward_rate - (Q/Keqm)*backwards_rate
#
#         if self.reaction_zone == 'cell':
#
#             tag = 'cell'
#
#         elif self.reaction_zone == 'mitochondria' and self.mit_enabled is True:
#
#             tag = 'mitochondria'
#
#         else:
#
#             tag = None
#
#         if tag is not None:
#
#             # get net effect of any activators or inhibitors of the reaction:
#             activator_alpha, inhibitor_alpha = get_influencers(sim, sim_metabo, self.reaction_activators_list,
#                 self.reaction_activators_Km, self.reaction_activators_n, self.reaction_inhibitors_list,
#                 self.reaction_inhibitors_Km, self.reaction_inhibitors_n, reaction_zone=tag)
#
#             # final degree of change, returned for use elsewhere (?):
#             flux = activator_alpha*inhibitor_alpha*reaction_rate
#             deltaC = activator_alpha*inhibitor_alpha*reaction_rate*p.dt
#
#             if self.reaction_zone == 'cell':
#
#                 self.set_reactant_c(deltaC, sim, sim_metabo,'c_cells', 'cc_cells')
#
#                 # re-obtain updated concentrations:
#                 self.get_reactants(sim, sim_metabo, 'c_cells', 'cc_cells')
#                 self.get_products(sim, sim_metabo, 'c_cells', 'cc_cells')
#
#                 self.set_product_c(deltaC, sim, sim_metabo, 'c_cells', 'cc_cells')
#
#
#             if self.reaction_zone == 'mitochondria' and self.mit_enabled is True:
#
#                 self.set_reactant_c(deltaC, sim, sim_metabo,'c_mit', 'cc_mit')
#
#                 # obtain updated concentrations:
#                 self.get_reactants(sim, sim_metabo, 'c_mit', 'cc_mit')
#                 self.get_products(sim, sim_metabo, 'c_mit', 'cc_mit')
#
#                 self.set_product_c(deltaC, sim, sim_metabo, 'c_mit', 'cc_mit')
#
#         else:
#
#             flux = None
#
#         return flux
#
#     def plot_1D(self, sim, cells, p, saveImagePath):
#
#         if self.reaction_zone == 'cell':
#
#             r_rate = [arr[p.plot_cell] for arr in self.rate_time]
#             fig = plt.figure()
#             ax = plt.subplot(111)
#             ax.plot(sim.time, r_rate)
#             ax.set_xlabel('Time [s]')
#             ax.set_ylabel('Rate [mM/s]')
#             ax.set_title('Rate of ' + self.name + ' in cell ' + str(p.plot_cell))
#
#         elif self.reaction_zone == 'mitochondria' and self.mit_enabled is True:
#
#             r_rate = [arr[p.plot_cell] for arr in self.rate_time]
#             fig = plt.figure()
#             ax = plt.subplot(111)
#             ax.plot(sim.time, r_rate)
#             ax.set_xlabel('Time [s]')
#             ax.set_ylabel('Rate [mM/s]')
#             ax.set_title('Rate of ' + self.name + ' in mitochondria ' + str(p.plot_cell))
#
#         if p.autosave is True:
#             savename = saveImagePath + 'ReactionRate_' + self.name + '.png'
#             plt.savefig(savename, format='png', transparent=True)
#
#         if p.turn_all_plots_off is False:
#             plt.show(block=False)
#
# class Transporter(object):
#
#     def __init__(self, sim, cells, p):
#
#         self.dummy_dyna = TissueHandler(sim, cells, p)
#         self.dummy_dyna.tissueProfiles(sim, cells, p)  # initialize all tissue profiles
#
#         # pre-populate the object with fields that will be assigned by MasterOfMolecules
#
#         self.name = None
#         self.reactants_list = None
#         self.reactants_coeff = None
#         self.products_list = None
#         self.products_coeff = None
#         self.Km_reactants_list = None
#         self.Km_products_list = None
#         self.vmax = None
#         self.delta_Go = None
#
#         self.reaction_zone = None
#
#         self.c_reactants = []  # concentrations of reactions, defined on cell-centres
#         self.z_reactants = []  # charge of reactants
#         self.inds_react = []  # indices to each reactant, to keep their order
#         self.c_products = []  # concentrations of reactions, defined on cell-centres
#         self.z_products = []  # charge of products
#         self.inds_prod = []  # indices to each reactant, to keep their order
#
#         self.product_source_object = []  # object from which product concentrations come from
#         self.product_source_type = []    # type of product concentrations sourced from object
#         self.reactant_source_object = []  # object from which reactants concentrations come from
#         self.reactant_source_type = []  # type of reactant concentrations sourced from object
#
#         # activator and inhibitors of the reaction
#         self.transporter_activators_list = None
#         self.transporter_activators_Km = None
#         self.transporter_activators_n = None
#         self.transporter_inhibitors_list = None
#         self.transporter_inhibitors_Km = None
#         self.transporter_inhibitors_n = None
#
#         self.transport_out_list = None
#         self.transport_in_list = None
#
#     def get_reactants(self, sim, sim_metabo, reactant_type_self, reactant_type_sim):
#
#         """
#         Get updated concentrations direct from sources for chemical reaction's reactants.
#
#         sim:                    An instance of simulator
#         sim_metabo:             Sim instance of reaction block (sim.metabo or sim.reacto for general reactions)
#         reactant_type_self:     Data type to retrieve: c_cells, c_mems, c_mit
#         reactant_type_sim:      Data typpe to retrieve: cc_cells, cc_mems, cc_mit
#
#         """
#
#         self.c_reactants = []
#
#         for reactant_name in self.reactants_list:
#
#             label = 'i' + reactant_name
#             ion_check = getattr(sim, label, None)
#
#             if ion_check is None:
#
#                 try:
#                     obj_reactant = getattr(sim_metabo, reactant_name)
#                     c_react = getattr(obj_reactant, reactant_type_self)
#                     self.c_reactants.append(c_react)
#
#                 except KeyError:
#
#                     raise BetseSimConfigException('Name of product is not a defined chemical, or is not'
#                                                    'an ion currently included in the ion profile being used.'
#                                                    'Please check biomolecule definitions and ion profile'
#                                                    'settings of your config(s) and try again.')
#
#             else:
#
#                 # define the reactant as the ion concentration from the cell concentrations object in sim:
#                 sim_conco = getattr(sim, reactant_type_sim)
#
#                 sim_conc = sim_conco[ion_check]
#
#                 self.c_reactants.append(sim_conc)
#
#     def get_products(self, sim, sim_metabo, product_type_self, product_type_sim):
#
#         """
#         Get updated concentrations direct from sources for chemical reaction's products.
#
#         sim:                    An instance of simulator
#         sim_metabo:             Sim instance of reaction block (sim.metabo or sim.reacto for general reactions)
#         product_type_self:     Data type to retrieve: c_cells, c_mems, c_mit
#         product_type_sim:      Data typpe to retrieve: cc_cells, cc_mems, cc_mit
#
#         """
#
#         self.c_products = []
#
#         for product_name in self.products_list:
#
#             label = 'i' + product_name
#             ion_check = getattr(sim, label, None)
#
#             if ion_check is None:
#
#                 try:
#                     obj_prod = getattr(sim_metabo, product_name)
#                     c_prod = getattr(obj_prod, product_type_self)
#                     self.c_products.append(c_prod)
#
#                 except KeyError:
#
#                     raise BetseSimConfigException('Name of product is not a defined chemical, or is not'
#                                                    'an ion currently included in the ion profile being used.'
#                                                    'Please check biomolecule definitions and ion profile'
#                                                    'settings of your config(s) and try again.')
#
#             else:
#
#                 # define the reactant as the ion concentration from the cell concentrations object in sim:
#                 sim_conco = getattr(sim, product_type_sim)
#
#                 sim_conc = sim_conco[ion_check]
#
#                 self.c_products.append(sim_conc)
#
#     def set_reactant_c(self, deltaMoles, sim, sim_metabo, reactant_tags, cells, p, ignoreECM = True):
#
#         targ_cells = self.transporter_targets_cell
#         targ_mems = self.transporter_targets_mem
#         targ_env = self.transporter_targets_env
#
#         for i, reactant_name in enumerate(self.reactants_list):
#
#             type_tag = reactant_tags[i]
#
#             if type_tag != 'c_env':
#
#                 if type_tag == 'c_mems':
#
#                     deltaC = deltaMoles/cells.mem_vol
#
#                     conc = np.zeros(sim.mdl)
#
#                     conc[targ_mems] = self.c_reactants[i][targ_mems] - deltaC[targ_mems]*self.reactants_coeff[i]
#
#                     label = 'i' + reactant_name
#                     ion_check = getattr(sim, label, None)
#
#                     if ion_check is None:
#
#                         obj_reactant = getattr(sim_metabo, reactant_name)
#                         setattr(obj_reactant, 'c_mems', conc)
#
#                     else:
#                         # define the reactant as the ion concentration from the cell concentrations object in sim:
#                         sim_conc = getattr(sim, 'cc_mems')
#                         sim_conc[ion_check] = conc
#
#                 elif type_tag == 'c_cells':
#
#                     deltaC = deltaMoles / cells.cell_vol
#
#                     conc = np.zeros(sim.cdl)
#
#                     conc[targ_cells] = self.c_reactants[i][targ_cells] - deltaC[targ_cells]*self.reactants_coeff[i]
#
#                     label = 'i' + reactant_name
#                     ion_check = getattr(sim, label, None)
#
#                     if ion_check is None:
#
#                         obj_reactant = getattr(sim_metabo, reactant_name)
#                         setattr(obj_reactant, 'c_cells', conc)
#
#                     else:
#                         # define the reactant as the ion concentration from the cell concentrations object in sim:
#                         sim_conc = getattr(sim, 'cc_cells')
#
#                         sim_conc[ion_check] = conc
#
#                 elif type_tag == 'c_mit':
#
#                     deltaC = deltaMoles / sim_metabo.mit.mit_vol
#
#                     conc = np.zeros(sim.cdl)
#
#                     conc[targ_cells] = self.c_reactants[i][targ_cells] - deltaC[targ_cells] * self.reactants_coeff[i]
#
#                     label = 'i' + reactant_name
#                     ion_check = getattr(sim, label, None)
#
#                     if ion_check is None:
#
#                         # print('setting mitochondria reactant')
#                         obj_reactant = getattr(sim_metabo, reactant_name)
#                         setattr(obj_reactant, 'c_mit', conc)
#                         # print('obj_reactant')
#                         # print('-----------')
#
#                     else:
#                         # define the reactant as the ion concentration from the cell concentrations object in sim:
#                         sim_conc = getattr(sim, 'cc_mit')
#                         sim_conc[ion_check] = conc
#
#             elif type_tag == 'c_env':
#
#                 if p.sim_ECM is False:
#
#                     deltaC = deltaMoles/cells.mem_sa
#
#                     conc = np.zeros(sim.mdl)
#
#                     conc[targ_mems] = self.c_reactants[i][targ_mems] - deltaC[targ_mems] * self.reactants_coeff[i]
#
#                     mean_conc = conc.mean()
#                     conc[:] = mean_conc
#
#                     label = 'i' + reactant_name
#                     ion_check = getattr(sim, label, None)
#
#                     if ion_check is None:
#
#                         obj_reactant = getattr(sim_metabo, reactant_name)
#                         setattr(obj_reactant, 'c_env', conc)
#
#                     else:
#
#                         # define the reactant as the ion concentration from the cell concentrations object in sim:
#                         sim_conc = getattr(sim, 'cc_env')
#                         sim_conc[ion_check] = conc
#
#
#                 elif p.sim_ECM is True:
#
#                     flux = deltaMoles/cells.mem_sa
#
#                     flux_env = np.zeros(sim.edl)
#                     flux_env[cells.map_mem2ecm] = flux
#
#                     # save values at the cluster boundary:
#                     bound_vals = flux_env[cells.ecm_bound_k]
#
#                     # set the values of the global environment to zero:
#                     flux_env[cells.inds_env] = 0
#
#                     # finally, ensure that the boundary values are restored:
#                     flux_env[cells.ecm_bound_k] = bound_vals
#
#                     if ignoreECM is True:
#
#                         delta_env = (flux_env * cells.memSa_per_envSquare) / cells.ecm_vol
#
#                     else:
#
#                         delta_env = (flux_env * cells.memSa_per_envSquare) / cells.true_ecm_vol
#
#                     label = 'i' + reactant_name
#                     ion_check = getattr(sim, label, None)
#
#                     if ion_check is None:
#
#                         obj_reactant = getattr(sim_metabo, reactant_name)
#
#                         conco = np.zeros(sim.edl)
#
#                         conco[targ_env] = obj_reactant.c_env[targ_env] - delta_env[targ_env]*self.reactants_coeff[i]
#
#                         setattr(obj_reactant, 'c_env', conco)
#
#                     else:
#                         # define the reactant as the ion concentration from the cell concentrations object in sim:
#                         sim_conc = getattr(sim, 'cc_env')
#
#                         conco = np.zeros(sim.edl)
#
#                         conco[targ_env] = sim_conc[ion_check][targ_env] - delta_env[targ_env]*self.reactants_coeff[i]
#
#                         sim_conc[ion_check] = conco[:]
#
#     def set_product_c(self, deltaMoles, sim, sim_metabo, product_tags, cells, p, ignoreECM = True):
#
#         targ_cells = self.transporter_targets_cell
#         targ_mems = self.transporter_targets_mem
#         targ_env = self.transporter_targets_env
#
#         for i, product_name in enumerate(self.products_list):
#
#             type_tag = product_tags[i]
#
#             if type_tag != 'c_env':
#
#                 if type_tag == 'c_mems':
#
#                     deltaC = deltaMoles / cells.mem_vol
#
#                     conc = np.zeros(sim.mdl)
#
#                     conc[targ_mems] = self.c_products[i][targ_mems] + deltaC[targ_mems] * self.products_coeff[i]
#
#                     label = 'i' + product_name
#                     ion_check = getattr(sim, label, None)
#
#                     if ion_check is None:
#
#                         obj_product = getattr(sim_metabo, product_name)
#                         setattr(obj_product, 'c_mems', conc)
#
#                     else:
#                         # define the reactant as the ion concentration from the cell concentrations object in sim:
#                         sim_conc = getattr(sim, 'cc_mems')
#
#                         sim_conc[ion_check] = conc
#
#                 if type_tag == 'c_cells':
#
#                     deltaC = deltaMoles / cells.cell_vol
#
#                     conc = np.zeros(sim.cdl)
#
#                     conc[targ_cells] = self.c_products[i][targ_cells] + deltaC[targ_cells] * self.products_coeff[i]
#
#                     label = 'i' + product_name
#                     ion_check = getattr(sim, label, None)
#
#                     if ion_check is None:
#
#                         obj_product = getattr(sim_metabo, product_name)
#                         setattr(obj_product, 'c_cells', conc)
#
#                     else:
#                         # define the reactant as the ion concentration from the cell concentrations object in sim:
#                         sim_conc = getattr(sim, 'cc_cells')
#
#                         sim_conc[ion_check] = conc
#
#                 elif type_tag == 'c_mit':
#
#                     deltaC = deltaMoles / sim_metabo.mit.mit_vol
#
#                     conc = np.zeros(sim.cdl)
#
#                     conc[targ_cells] = self.c_products[i][targ_cells] + deltaC * self.products_coeff[i]
#
#                     label = 'i' + product_name
#                     ion_check = getattr(sim, label, None)
#
#                     if ion_check is None:
#
#                         obj_product = getattr(sim_metabo, product_name)
#                         setattr(obj_product, 'c_mit', conc)
#
#                     else:
#                         # define the reactant as the ion concentration from the cell concentrations object in sim:
#                         sim_conc = getattr(sim, 'cc_mit')
#
#                         sim_conc[ion_check] = conc
#
#             elif type_tag == 'c_env':
#
#                 if p.sim_ECM is False:
#
#                     deltaC = deltaMoles / cells.mem_sa
#
#                     conc = np.zeros(sim.mdl)
#
#                     conc[targ_mems] = self.c_products[i][targ_mems] + deltaC[targ_mems] * self.products_coeff[i]
#
#                     # mean_conc = conc.mean()
#                     # conc[:] = mean_conc
#
#                     label = 'i' + product_name
#                     ion_check = getattr(sim, label, None)
#
#                     if ion_check is None:
#
#                         obj_product = getattr(sim_metabo, product_name)
#                         setattr(obj_product, 'c_env', conc)
#
#                     else:
#
#
#                         # define the reactant as the ion concentration from the cell concentrations object in sim:
#                         sim_conc = getattr(sim, 'cc_env')
#                         sim_conc[ion_check] = conc
#
#
#
#                 elif p.sim_ECM is True:
#
#                     flux = deltaMoles / cells.mem_sa
#
#                     flux_env = np.zeros(sim.edl)
#                     flux_env[cells.map_mem2ecm] = flux
#
#                     # save values at the cluster boundary:
#                     bound_vals = flux_env[cells.ecm_bound_k]
#
#                     # set the values of the global environment to zero:
#                     flux_env[cells.inds_env] = 0
#
#                     # finally, ensure that the boundary values are restored:
#                     flux_env[cells.ecm_bound_k] = bound_vals
#
#                     if ignoreECM is True:
#
#                         delta_env = (flux_env * cells.memSa_per_envSquare) / cells.ecm_vol
#
#                     else:
#
#                         delta_env = (flux_env * cells.memSa_per_envSquare) / cells.true_ecm_vol
#
#                     label = 'i' + product_name
#                     ion_check = getattr(sim, label, None)
#
#                     if ion_check is None:
#
#                         obj_product = getattr(sim_metabo, product_name)
#
#                         conco = np.zeros(sim.edl)
#
#                         conco[targ_env] = obj_product.c_env[targ_env] + delta_env[targ_env]*self.products_coeff[i]
#
#                         setattr(obj_product, 'c_env', conco)
#
#                     else:
#                         # define the reactant as the ion concentration from the cell concentrations object in sim:
#                         sim_conc = getattr(sim, 'cc_env')
#
#                         conco = np.zeros(sim.edl)
#
#                         conco[targ_env] = sim_conc[ion_check][targ_env] + delta_env[targ_env] * self.products_coeff[i]
#
#                         sim_conc[ion_check] = conco[:]
#
#     def init_reaction(self,cells, p):
#
#         if self.transporter_profiles_list is not None and self.transporter_profiles_list != 'all':
#
#             self.transporter_targets_mem = []
#             self.transporter_targets_cell = []
#             self.transporter_targets_env = []
#
#             for profile in self.transporter_profiles_list:
#
#                 targets_cell = self.dummy_dyna.cell_target_inds[profile]
#                 self.transporter_targets_cell.extend(targets_cell)
#
#                 targets_mem = self.dummy_dyna.tissue_target_inds[profile]
#                 self.transporter_targets_mem.extend(targets_mem)
#
#                 targets_env = self.dummy_dyna.env_target_inds[profile]
#                 self.transporter_targets_env.extend(targets_env)
#
#         elif self.transporter_profiles_list is None or self.transporter_profiles_list == 'all':
#
#             self.transporter_targets_mem = cells.mem_i
#             self.transporter_targets_cell = cells.cell_i
#             self.transporter_targets_env = [x for x in range(0, len(cells.xypts))]
#
#     def compute_reaction(self, sim, sim_metabo, cells, p):
#
#         if self.reaction_zone == 'cell':
#
#
#             type_self_out = 'c_env'
#             type_sim_out = 'cc_env'
#
#             type_self_in = 'c_mems'
#             type_sim_in = 'cc_mems'
#
#             type_self = 'c_mems'
#             type_sim = 'cc_mems'
#
#             vmem = sim.vm   # get the transmembrane voltage for this category
#
#         elif self.reaction_zone == 'mitochondria' and self.mit_enabled is True:
#             type_self_out = 'c_cells'
#             type_sim_out = 'cc_cells'
#
#             type_self_in = 'c_mit'
#             type_sim_in = 'cc_mit'
#
#             type_self = 'c_mit'
#             type_sim = 'cc_mit'
#
#             vmem = sim_metabo.mit.Vmit  # get the transmembrane voltage for this category
#
#         else:
#             return  # don't proceed any further, get out of this function.
#
#         echem_terms = []
#         c_reactants_trans = []   # substances transferred across membrane -- start state concs
#         c_products_trans = []    # substances transferred across membrane -- end state concs
#         trans_react_index = []
#         trans_prod_index = []
#         react_transfer_tag = []
#         prod_transfer_tag = []
#
#         if self.transport_out_list != None:
#
#             for out_name in self.transport_out_list:
#
#                 _, c_in = get_conc(sim, sim_metabo, out_name, type_self_in, type_sim_in, cells, p)
#                 z_out, c_out = get_conc(sim, sim_metabo, out_name, type_self_out, type_sim_out, cells, p)
#
#                 c_reactants_trans.append(c_in)
#                 react_transfer_tag.append(type_self_in)
#
#                 c_products_trans.append(c_out)
#                 prod_transfer_tag.append(type_self_out)
#
#                 ind_r = self.reactants_list.index(out_name)
#                 coeff = self.reactants_coeff[ind_r]
#
#                 ind_p = self.products_list.index(out_name)
#
#                 trans_react_index.append(ind_r)
#                 trans_prod_index.append(ind_p)
#
#                 # get the electrochemical potential term for this reagent
#                 # it's negative because it starts inside the cell
#                 out_term = -z_out*coeff*vmem*p.F
#
#                 echem_terms.append(out_term)
#
#         if self.transport_in_list != None:
#
#             for in_name in self.transport_in_list:
#
#                 z_in, c_in = get_conc(sim, sim_metabo, in_name, type_self_in, type_sim_in, cells, p)
#                 _, c_out = get_conc(sim, sim_metabo, in_name, type_self_out, type_sim_out, cells, p)
#
#                 c_reactants_trans.append(c_out)
#                 react_transfer_tag.append(type_self_out)
#
#                 c_products_trans.append(c_in)
#                 prod_transfer_tag.append(type_self_in)
#
#                 ind_r = self.reactants_list.index(in_name)
#                 coeff = self.reactants_coeff[ind_r]
#
#                 ind_p = self.products_list.index(in_name)
#
#                 trans_react_index.append(ind_r)
#                 trans_prod_index.append(ind_p)
#
#                 # get the electrochemical potential term for this reagent
#                 # it's positive because it ends inside the cell
#                 in_term = z_in*coeff*vmem*p.F
#
#                 echem_terms.append(in_term)
#
#         echem_terms = np.asarray(echem_terms)
#         vmem_term = np.sum(echem_terms, axis = 0)
#
#         # modification factor for the equilibrium constant due to transfer of charged item across
#         #  transmembrane voltage:
#
#         # Kmod = np.exp(-vmem_term/(p.R*sim.T))
#
#         deltaGi = -vmem_term/(p.R * sim.T)
#
#         # get up-to-date concentration data for the reaction:
#         self.get_reactants(sim, sim_metabo, type_self, type_sim)
#         self.get_products(sim, sim_metabo, type_self, type_sim)
#
#         if self.reaction_zone == 'cell':
#
#             self.reactant_transfer_tag = ['c_mems' for x in range(0, len(self.c_reactants))]
#             self.product_transfer_tag = ['c_mems' for x in range(0, len(self.c_products))]
#
#         elif self.reaction_zone == 'mitochondria' and self.mit_enabled is True:
#
#             self.reactant_transfer_tag = ['c_mit' for x in range(0, len(self.c_reactants))]
#             self.product_transfer_tag = ['c_mit' for x in range(0, len(self.c_products))]
#
#         for i, ind_r in enumerate(trans_react_index):
#
#             self.c_reactants[ind_r] = c_reactants_trans[i]
#             self.reactant_transfer_tag[ind_r] = react_transfer_tag[i]
#
#         for j, ind_p in enumerate(trans_prod_index):
#
#             self.c_products[ind_p] = c_products_trans[j]
#             self.product_transfer_tag[ind_p] = prod_transfer_tag[j]
#
#         # define the reaction equilibrium coefficient:
#         # Ko = np.exp(-self.delta_Go/(p.R*sim.T))
#
#         Keqm = np.exp(-self.delta_Go/(p.R*sim.T) + deltaGi)
#
#         # if self.name == 'ETC':
#         #
#         #     print(Keqm.min(), Keqm.max(), Keqm.mean())
#             # print(deltaGi.mean())
#         #
#         #     deltaGo = -self.delta_Go / (p.R * sim.T)
#         #     deltaGi = - vmem_term / (p.R * sim.T)
#         #     deltaGG = -self.delta_Go / (p.R * sim.T) - vmem_term / (p.R * sim.T)
#         #
#         #     print(deltaGo)
#
#             # print(deltaGi.min(), deltaGi.max(), deltaGi.mean())
#
#             # print(Keqm.min(), Keqm.max(), Keqm.mean())
#             # print(Ko.min(), Ko.max(), Ko.mean())
#             # print(Kmod.min(), Kmod.max(), Kmod.mean())
#
#         # calculate the MM rate coefficient for the backwards reaction direction and Q:
#         backwards_term = []
#         Q_deno_list = []
#         Q_numo_list = []
#
#         for i, prod in enumerate(self.c_products):
#             # calculate a factor in the backwards rate term:
#             coeff = self.products_coeff[i]
#             Km = self.Km_products_list[i]
#             cs = (prod / Km) ** coeff
#             term = cs / (1 + cs)
#             backwards_term.append(term)
#
#             # calculate a factor in the reaction quotient numerator term:
#             ci = prod**coeff
#             Q_numo_list.append(ci)
#
#         backwards_term = np.asarray(backwards_term)
#
#         backwards_rate = self.vmax * np.prod(backwards_term, axis=0)
#
#         for j, react in enumerate(self.c_reactants):
#             # calculate a factor in the reaction quotient denominator term:
#             coeff = self.reactants_coeff[j]
#             cj = react**coeff
#             Q_deno_list.append(cj)
#
#         # get the numerator and denomenator of the reaction quotient:
#         Q_deno = np.prod(Q_deno_list, axis=0)
#         Q_numo = np.prod(Q_numo_list, axis=0)
#
#         inds_neg = (Q_deno == 0).nonzero()
#         Q_deno[inds_neg] = 1.0e-15   # ensure no division by zero
#
#         # finally, *the* reaction quotient:
#         Q = Q_numo/Q_deno
#
#         # Next, calculate the MM rate coefficient for the forward reaction direction:
#         forward_term = []
#
#         for i, react in enumerate(self.c_reactants):
#
#             coeff = self.reactants_coeff[i]
#             Km = self.Km_reactants_list[i]
#             cs = (react/Km)**coeff
#
#             term = cs/(1 + cs)
#
#             forward_term.append(term)
#
#         forward_term = np.asarray(forward_term)
#
#         forward_rate = self.vmax*np.prod(forward_term, axis=0)
#
#         reaction_rate = forward_rate - (Q/Keqm)*backwards_rate
#
#         # get net effect of any activators or inhibitors of the reaction:
#
#         if self.reaction_zone == 'cell':
#
#             tag = 'mems'
#
#         elif self.reaction_zone == 'mitochondria' and self.mit_enabled is True:
#
#             tag = 'mitochondria'
#
#
#         activator_alpha, inhibitor_alpha = get_influencers(sim, sim_metabo, self.transporter_activators_list,
#             self.transporter_activators_Km, self.transporter_activators_n, self.transporter_inhibitors_list,
#             self.transporter_inhibitors_Km, self.transporter_inhibitors_n, reaction_zone=tag)
#
#
#         deltaC = activator_alpha*inhibitor_alpha*reaction_rate*p.dt
#
#         # multiply by compartment volume to get a mole transfer:
#         if self.reaction_zone == 'cell':
#             deltaMoles = deltaC*cells.mem_vol
#
#         elif self.reaction_zone == 'mitochondria' and self.mit_enabled is True:
#             deltaMoles = deltaC*sim_metabo.mit.mit_vol
#
#         else:
#             deltaMoles = None
#
#
#         if deltaMoles is not None:
#
#             self.set_reactant_c(deltaMoles, sim, sim_metabo,self.reactant_transfer_tag, cells, p,
#                                 ignoreECM = self.ignore_ECM_transporter)
#
#             # get up-to-date concentration data for the reaction:
#             self.get_reactants(sim, sim_metabo, type_self, type_sim)
#             self.get_products(sim, sim_metabo, type_self, type_sim)
#
#             self.set_product_c(deltaMoles, sim, sim_metabo, self.product_transfer_tag, cells, p,
#                                ignoreECM = self.ignore_ECM_transporter)
#
#         else:
#             deltaC = None
#
#         return deltaC
#
#     def plot_1D(self, sim, cells, p, saveImagePath):
#
#         if len(self.rate_time) > 0:
#
#             if len(self.rate_time[0]) == sim.cdl:
#
#                 r_rate = [arr[p.plot_cell] for arr in self.rate_time]
#
#             else:
#                 mem_i = cells.cell_to_mems[p.plot_cell][0]
#                 r_rate = [arr[mem_i] for arr in self.rate_time]
#
#             fig = plt.figure()
#             ax = plt.subplot(111)
#             ax.plot(sim.time, r_rate)
#             ax.set_xlabel('Time [s]')
#             ax.set_ylabel('Rate [mM/s]')
#             ax.set_title('Rate of ' + self.name + ' in cell ' + str(p.plot_cell))
#
#             if p.autosave is True:
#                 savename = saveImagePath + 'TransporterRate_' + self.name + '.png'
#                 plt.savefig(savename, format='png', transparent=True)
#
#             if p.turn_all_plots_off is False:
#                 plt.show(block=False)
#
#     def update_transporter(self, sim, cells, p):
#
#         self.dummy_dyna.tissueProfiles(sim, cells, p)  # initialize all tissue profiles
#         self.init_reaction(cells, p)
#
# class Channel(object):
#
#     def __init__(self, sim, cells, p):
#
#         self.dummy_dyna = TissueHandler(sim, cells, p)
#         self.dummy_dyna.tissueProfiles(sim, cells, p)  # initialize all tissue profiles
#
#     def init_channel(self, ion_string, type_string, max_val, sim, cells, p):
#
#         # get targets for the reaction
#         if self.channel_profiles_list is not None and self.channel_profiles_list != 'all':
#
#             self.channel_targets_mem = []
#
#             for profile in self.channel_profiles_list:
#
#                 targets_mem = self.dummy_dyna.tissue_target_inds[profile]
#                 self.channel_targets_mem.extend(targets_mem)
#
#             self.channel_targets_mem = np.asarray(self.channel_targets_mem)
#
#         elif self.channel_profiles_list is None or self.channel_profiles_list == 'all':
#
#             self.channel_targets_mem = np.asarray(cells.mem_i)
#
#
#         if ion_string == 'Na':
#
#             self.dummy_dyna.maxDmNa = max_val
#             self.dummy_dyna.targets_vgNa = self.channel_targets_mem
#             class_string = vgna
#
#         elif ion_string == 'NaP':
#
#             self.dummy_dyna.maxDmNaP = max_val
#             self.dummy_dyna.targets_vgNaP = self.channel_targets_mem
#             class_string = vgnap
#
#         elif ion_string == 'K':
#
#             self.dummy_dyna.maxDmK = max_val
#             self.dummy_dyna.targets_vgK = self.channel_targets_mem
#             class_string = vgk
#
#         elif ion_string == 'Kir':
#
#             self.dummy_dyna.maxDmKir = max_val
#             self.dummy_dyna.targets_vgKir = self.channel_targets_mem
#             class_string = vgkir
#
#         elif ion_string == 'Ca':
#
#             self.dummy_dyna.maxDmCa = max_val
#             self.dummy_dyna.targets_vgCa = self.channel_targets_mem
#             class_string = vgca
#
#         elif ion_string == 'Fun':
#
#             self.dummy_dyna.maxDmFun = max_val
#             self.dummy_dyna.targets_vgFun = self.channel_targets_mem
#             class_string = vgfun
#
#         else:
#
#             raise BetseSimConfigException("Substance-modulated ion type not available. "
#                                            "Valid choices: Na, K, Ca, NaP, Kir, and Fun")
#
#             # create the desired voltage gated sodium channel instance:
#
#         self.channel_core = getattr(class_string,type_string)()
#
#
#         if p.run_sim is True:
#             # initialize the channel object
#             self.channel_core.init(self.dummy_dyna, sim, cells, p)
#
#     def run_channel(self, sim, sim_metabo, cells, p):
#
#         # get modulation coefficients by any activating/inhibiting substances:
#         activator_alpha, inhibitor_alpha = get_influencers(sim, sim_metabo, self.channel_activators_list,
#             self.channel_activators_Km, self.channel_activators_n, self.channel_inhibitors_list,
#             self.channel_inhibitors_Km, self.channel_inhibitors_n, reaction_zone='mems')
#
#         # calculate the value of the channel modulation constant:
#         self.channel_core.modulator = activator_alpha*inhibitor_alpha
#
#         self.channel_core.run(self.dummy_dyna, sim, cells, p)
#
#     def update_channel(self, sim, cells, p):
#         self.dummy_dyna.tissueProfiles(sim, cells, p)  # initialize all tissue profiles
#         self.init_channel(self.channel_class, self.channel_type, self.channelMax, sim, cells, p)
#
# class Modulator(object):
#     """
#     The modulator object allows a substance defined in MasterOfMolecules to
#     exert an activating or inhibiting influence over simulation-defined pumps
#     and/or gap junctions.
#
#     """
#     def __init__(self):
#
#         self.target_label = None
#         self.max_val = None
#
#         self.modulator_activators_list = None
#         self.modulator_activators_Km = None
#         self.modulator_activators_n = None
#         self.modulator_inhibitors_list = None
#         self.modulator_inhibitors_Km = None
#         self.modulator_inhibitors_n = None
#
#     def init_modulator(self, sim, cells, p):
#
#         if self.target_label == 'gj':
#
#             sim.gj_block_o = np.ones(sim.mdl)
#
#         elif self.target_label == 'Na/K-ATPase':
#
#             sim.NaKATP_block_o =  np.ones(sim.mdl)
#
#         elif self.target_label == 'H/K-ATPase':
#
#             sim.HKATP_block_o =  np.ones(sim.mdl)
#
#         elif self.target_label == 'V-ATPase':
#
#             sim.VATP_block_o =  np.ones(sim.mdl)
#
#         else:
#
#             raise BetseSimConfigException("You have requested a "
#                                            "sim modulator that is not "
#                                            "available. Available choices "
#                                            "are: 'gj', 'Na/K-ATPase', 'H/K-ATPase', "
#                                            "and 'V-ATPase' ")
#
#     def run_modulator(self, sim, sim_metabo, cells, p):
#
#         if self.zone == 'env':
#
#             zone_tag = 'env'
#
#         elif self.zone == 'cell':
#
#             zone_tag = 'mems'
#
#         else:
#
#             raise BetseSimConfigException("You have requested an unavailable modulator zone."
#                                            "Available choices are 'env' and 'mems'.")
#
#         # get the coefficients activating and/or inhibiting the sim structure:
#
#         # get modulation coefficients by any activating/inhibiting substances:
#         activator_alpha, inhibitor_alpha = get_influencers(sim, sim_metabo, self.modulator_activators_list,
#             self.modulator_activators_Km, self.modulator_activators_n, self.modulator_inhibitors_list,
#             self.modulator_inhibitors_Km, self.modulator_inhibitors_n, reaction_zone=zone_tag)
#
#         # calculate the value of the channel modulation constant:
#         modulator = self.max_val*activator_alpha * inhibitor_alpha
#
#         # make size alteration for case of true environment:
#         if p.sim_ECM is True and self.zone == 'env':
#             modulator = modulator[cells.map_mem2ecm]
#
#         if self.target_label == 'gj':
#
#             sim.gj_block = modulator
#
#         elif self.target_label == 'Na/K-ATPase':
#
#             sim.NaKATP_block = modulator
#
#         elif self.target_label == 'H/K-ATPase':
#
#             sim.HKATP_block = modulator
#
#         elif self.target_label == 'V-ATPase':
#
#             sim.VATP_block = modulator
#
#         elif self.target_label == 'Ca-ATPase':
#
#             sim.CaATP_block = modulator
#
#         elif self.target_label == 'Na/Ca-Exch':
#
#             sim.NaCaExch_block = modulator
#
#         else:
#
#             raise BetseSimConfigException("You have requested a "
#                                            "sim modulator that is not "
#                                            "available. Available choices "
#                                            "are: 'gj', 'Na/K-ATPase', 'H/K-ATPase', "
#                                            "and 'V-ATPase', 'Ca-ATPase', and 'Na/Ca-Exch' ")
#
# def get_influencers(sim, sim_metabo, a_list, Km_a_list, n_a_list, i_list, Km_i_list,
#                     n_i_list, reaction_zone='cell'):
#     """
#     Get coefficients representing the net effect of all activators and inhibitors on a particular reaction.
#
#     Parameters
#     ------------
#     sim                 Instance of BETSE simulator
#     sim_metabo:         Instance of MasterOfMetabolism
#     a_list:             activator names list
#     Km_a_list:          activator half-max constants
#     n_a_list:           activator Hill exponents
#     i_list:             inhibitor names list
#     Km_i_list:          inhibitor half-max constants
#     n_i_list:           inhibitor Hill exponents
#     reaction_zone:      Reaction occurring in 'cell' or 'mitochondria'
#
#     Returns
#     ------------
#     activator_alpha         Coefficient of net effect of activators
#     inhibitor_alpha         Coefficient of net effect of inhibitors
#     """
#
#     if reaction_zone == 'cell':
#         type_self = 'c_cells'
#         type_sim = 'cc_cells'
#
#     elif reaction_zone == 'mems':
#         type_self = 'c_mems'
#         type_sim = 'cc_mems'
#
#     elif reaction_zone == 'mitochondria':
#
#         type_self = 'c_mit'
#         type_sim = 'cc_mit'
#
#     elif reaction_zone == 'env':
#
#         type_self = 'c_env'
#         type_sim = 'cc_env'
#
#     # initialize a blank list
#     activator_terms = []
#
#     if a_list is not None and a_list != 'None' and len(a_list) > 0:  # if user specified activators for growth/decay
#
#         # get reaction zone for data type:
#
#         # get the activator concentration for the substance, and
#         # create a term based on Hill form:
#         for i, activator_name in enumerate(a_list):
#
#             label = 'i' + activator_name
#             ion_check = getattr(sim, label, None)
#
#             if ion_check is None:
#
#                 try:
#                     obj_activator = getattr(sim_metabo, activator_name)
#                     c_act = getattr(obj_activator, type_self)
#
#                 except KeyError:
#
#                     raise BetseSimConfigException('Name of reaction activator is not a defined chemical, '
#                                                    'or is not an ion currently included in the ion profile '
#                                                    'being used.'
#                                                    'Please check biomolecule definitions and ion profile'
#                                                    'settings of your config(s) and try again.')
#
#             else:
#                 # define the reactant as the ion concentration from the cell concentrations object in sim:
#                 sim_conco = getattr(sim, type_sim)
#                 c_act = sim_conco[ion_check]
#
#             Km_act = Km_a_list[i]
#             n_act = n_a_list[i]
#
#             cs = (c_act / Km_act)**n_act
#
#             act_term = cs / (1 + cs)
#
#             activator_terms.append(act_term)
#
#         activator_terms = np.asarray(activator_terms)
#
#         # calculate the net effect of all activator terms:
#         activator_alpha = np.prod(activator_terms, axis=0)
#
#
#     else:
#
#         activator_alpha = 1
#
#
#     # initialize a blank list
#     inhibitor_terms = []
#
#     if i_list is not None and i_list != 'None' and len(i_list) > 0:  # if user specified inhibitors for growth/decay
#
#         # get the inhibitor concentration for the substance, and
#         # create a term based on Hill form:
#         for j, inhibitor_name in enumerate(i_list):
#
#             label = 'i' + inhibitor_name
#             ion_check = getattr(sim, label, None)
#
#             if ion_check is None:
#
#                 try:
#                     obj_inhibitor = getattr(sim_metabo, inhibitor_name)
#                     c_inh = getattr(obj_inhibitor, type_self)
#
#                 except KeyError:
#
#                     raise BetseSimConfigException('Name of substance is not a defined chemical, '
#                                                    'or is not an ion currently included in the ion profile '
#                                                    'being used.'
#                                                    'Please check biomolecule definitions and ion profile'
#                                                    'settings of your config(s) and try again.')
#
#             else:
#                 # define the reactant as the ion concentration from the cell concentrations object in sim:
#                 sim_conco = getattr(sim, type_sim)
#                 c_inh = sim_conco[ion_check]
#
#             # print(i_list, Km_i_list, n_i_list)
#
#             Km_inh = Km_i_list[j]
#             n_inh = n_i_list[j]
#
#             cs = (c_inh / Km_inh) ** n_inh
#
#             inh_term = 1 / (1 + cs)
#
#             inhibitor_terms.append(inh_term)
#
#         inhibitor_terms = np.asarray(inhibitor_terms)
#
#         # calculate the net effect of all activator terms:
#         inhibitor_alpha = np.prod(inhibitor_terms, axis = 0)
#
#     else:
#         inhibitor_alpha = 1
#
#
#     return activator_alpha, inhibitor_alpha
#
# def get_influencers_grn(sim, sim_metabo, a_list, k_a_list, Km_a_list, n_a_list, i_list, k_i_list, Km_i_list,
#         n_i_list, reaction_zone='cell'):
#
#     """
#     Get coefficients representing the net effect of all activators and inhibitors on a particular reaction.
#
#     Parameters
#     ------------
#     sim                 Instance of BETSE simulator
#     sim_metabo:         Instance of MasterOfMetabolism
#     a_list:             activator names list
#     Km_a_list:          activator half-max constants
#     n_a_list:           activator Hill exponents
#     i_list:             inhibitor names list
#     Km_i_list:          inhibitor half-max constants
#     n_i_list:           inhibitor Hill exponents
#     reaction_zone:      Reaction occurring in 'cell' or 'mitochondria'
#
#     Returns
#     ------------
#     activator_alpha         Coefficient of net effect of activators
#     inhibitor_alpha         Coefficient of net effect of inhibitors
#     """
#
#     if reaction_zone == 'cell':
#         type_self = 'c_cells'
#         type_sim = 'cc_cells'
#
#     if reaction_zone == 'mems':
#         type_self = 'c_mems'
#         type_sim = 'cc_mems'
#
#     elif reaction_zone == 'mitochondria':
#
#         type_self = 'c_mit'
#         type_sim = 'cc_mit'
#
#     # initialize a blank list
#     activator_terms = []
#
#     if a_list is not None and a_list != 'None' and len(a_list) > 0:  # if user specified activators for growth/decay
#
#         # get reaction zone for data type:
#
#         # get the activator concentration for the substance, and
#         # create a term based on Hill form:
#         for i, activator_name in enumerate(a_list):
#
#             label = 'i' + activator_name
#             ion_check = getattr(sim, label, None)
#
#             if ion_check is None:
#
#                 try:
#                     obj_activator = getattr(sim_metabo, activator_name)
#                     c_act = getattr(obj_activator, type_self)
#
#                 except KeyError:
#
#                     raise BetseSimConfigException('Name of reaction activator is not a defined chemical, '
#                                                    'or is not an ion currently included in the ion profile '
#                                                    'being used.'
#                                                    'Please check biomolecule definitions and ion profile'
#                                                    'settings of your config(s) and try again.')
#
#             else:
#                 # define the reactant as the ion concentration from the cell concentrations object in sim:
#                 sim_conco = getattr(sim, type_sim)
#                 c_act = sim_conco[ion_check]
#
#             k_act = k_a_list[i]
#             Km_act = Km_a_list[i]
#             n_act = n_a_list[i]
#
#             cs = ((c_act*k_act)/Km_act) ** n_act
#
#             act_term = cs / (1 + cs)
#
#             activator_terms.append(act_term)
#
#         activator_terms = np.asarray(activator_terms)
#
#         # calculate the net effect of all activator terms:
#         activator_alpha = np.prod(activator_terms, axis=0)
#
#     else:
#
#         activator_alpha = 1
#
#     # initialize a blank list
#     inhibitor_terms = []
#
#     if i_list is not None and i_list != 'None' and len(i_list) > 0:  # if user specified inhibitors for growth/decay
#
#         # get the inhibitor concentration for the substance, and
#         # create a term based on Hill form:
#         for j, inhibitor_name in enumerate(i_list):
#
#             label = 'i' + inhibitor_name
#             ion_check = getattr(sim, label, None)
#
#             if ion_check is None:
#
#                 try:
#                     obj_inhibitor = getattr(sim_metabo, inhibitor_name)
#                     c_inh = getattr(obj_inhibitor, type_self)
#
#                 except KeyError:
#
#                     raise BetseSimConfigException('Name of substance is not a defined chemical, '
#                                                    'or is not an ion currently included in the ion profile '
#                                                    'being used.'
#                                                    'Please check biomolecule definitions and ion profile'
#                                                    'settings of your config(s) and try again.')
#
#             else:
#                 # define the reactant as the ion concentration from the cell concentrations object in sim:
#                 sim_conco = getattr(sim, type_sim)
#                 c_inh = sim_conco[ion_check]
#
#             k_inh = k_i_list[j]
#             Km_inh = Km_i_list[j]
#             n_inh = n_i_list[j]
#
#             cs = ((c_inh*k_inh)/Km_inh) ** n_inh
#
#             inh_term = 1 / (1 + cs)
#
#             inhibitor_terms.append(inh_term)
#
#         inhibitor_terms = np.asarray(inhibitor_terms)
#
#         # calculate the net effect of all activator terms:
#         inhibitor_alpha = np.prod(inhibitor_terms, axis=0)
#
#     else:
#         inhibitor_alpha = 1
#
#     return activator_alpha, inhibitor_alpha
#
# def get_conc(sim, sim_metabo, name, type_self, type_sim, cells, p):
#
#         label = 'i' + name
#         ion_check = getattr(sim, label, None)
#
#         if ion_check is None:
#
#             try:
#                 obj = getattr(sim_metabo, name)
#                 c = getattr(obj, type_self)
#                 z = obj.z
#
#             except KeyError:
#
#                 raise BetseSimConfigException('Name of substance is not a defined chemical, '
#                                                'or is not an ion currently included in the ion profile '
#                                                'being used.'
#                                                'Please check biomolecule definitions and ion profile'
#                                                'settings of your config(s) and try again.')
#
#         else:
#             # define the reactant as the ion concentration from the cell concentrations object in sim:
#             sim_conco = getattr(sim, type_sim)
#             c = sim_conco[ion_check]
#             z = sim.zs[ion_check]
#
#         if type_self == 'c_env' or type_sim == 'cc_env':
#
#             if p.sim_ECM:
#                 c = c[cells.map_mem2ecm]
#
#         return z, c
#

#------------------------------------------------------------------------------------------------------------------
# WASTELANDS
#------------------------------------------------------------------------------------------------------------------
# def assign_new_concentrations(self, obj, delta_c, new_reactants, new_products, sim, cells, p):
#
#     for i, react_name in enumerate(obj.reactants_list):
#
#         if obj.reactant_source_object[i] == id(self):
#
#             source_obj = getattr(self, react_name)
#
#             setattr(source_obj, obj.reactant_source_type[i], new_reactants[i])
#
#         elif obj.reactant_source_object[i] == id(sim):
#
#             source_obj = getattr(sim, obj.reactant_source_type[i])
#             ion_label = 'i' + react_name
#             ion_type = getattr(sim, ion_label)
#
#             source_obj[ion_type] = new_reactants[i]
#
#     for j, prod_name in enumerate(obj.products_list):
#
#             if obj.product_source_object[j] == id(self):
#
#                 source_obj = getattr(self, prod_name)
#
#                 setattr(source_obj, obj.product_source_type[j], new_products[j])
#
#             elif obj.product_source_object[j] == id(sim):
#
#                 source_obj = getattr(sim, obj.product_source_type[j])
#                 ion_label = 'i' + prod_name
#                 ion_type = getattr(sim, ion_label)
#
#                 source_obj[ion_type] = new_products[j]

# def check_reactions(self, obj, sim, cells, p):

#
#     if obj.reaction_zone == 'cell':
#
#         # Case 1: no transfer across the membrane -- these are in-cytoplasm reactions -------------------------------
#         for i, reactant_name in enumerate(obj.reactants_list):
#
#             self.set_react_concs(obj, sim, cells, p, reactant_name, i, 'c_cells', 'cc_cells')
#
#         # Now load the right concentration data arrays for products into the Reaction object:
#         for j, product_name in enumerate(obj.products_list):
#
#             self.set_prod_concs(obj, sim, cells, p, product_name, j, 'c_cells', 'cc_cells')
#
#
#     elif obj.reaction_zone == 'mitochondria':
#
#         # Case 1: no transfer across the membrane -- these are in-cytoplasm reactions -------------------------------
#         for i, reactant_name in enumerate(obj.reactants_list):
#
#             self.set_react_concs(obj, sim, cells, p, reactant_name, i, 'c_mit', 'cc_mit')
#
#         # Now load the right concentration data arrays for products into the Reaction object:
#         for j, product_name in enumerate(obj.products_list):
#
#             self.set_prod_concs(obj, sim, cells, p, product_name, j, 'c_mit', 'cc_mit')
#
#     else:
#         raise BetseSimConfigException("You have requested a reaction zone that does not exist."
#                                        "Valid options include: 'cell' and 'mitochondria'."
#                                        "Please check your config(s) settings defining reactions"
#                                        " and try again. ")
