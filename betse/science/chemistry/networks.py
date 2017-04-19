#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS                            }....................
import csv
import math
import os
import os.path
import matplotlib.pyplot as plt
import numpy as np
from collections import OrderedDict
from matplotlib import cm
from matplotlib import colors
from scipy.ndimage.filters import gaussian_filter
# from scipy.optimize import basinhopping
from betse.exceptions import BetseSimConfigException, BetseSimInstabilityException
from betse.lib import libs
from betse.science import sim_toolbox as stb
from betse.science.channels import cation as vgcat
from betse.science.channels import vg_ca as vgca
from betse.science.channels import vg_cl as vgcl
from betse.science.channels import vg_funny as vgfun
from betse.science.channels import vg_k as vgk
from betse.science.channels import vg_na as vgna
from betse.science.chemistry.netplot import plot_master_network, set_net_opts
from betse.science.config.export.confvis import SimConfVisualCellsNonYAML
from betse.science.math import modulate as mods
from betse.science.math import toolbox as tb
from betse.science.organelles.mitochondria import Mito
from betse.science.simulate.simphase import SimPhase
from betse.science.tissue.handler import TissueHandler
from betse.science.visual.anim.anim import AnimFlatCellsTimeSeries, AnimEnvTimeSeries
from betse.science.visual.plot import plotutil as viz
from betse.util.io.log import logs
from betse.util.path import paths
from betse.util.type.mappings import DynamicValue, DynamicValueDict
from betse.util.type.types import type_check

# FIXME: if moving to have unpacked membrane concs, update transporters...

# ....................{ CLASSES                            }....................
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

        # set the key controlling presence of mitochondria:
        self.mit_enabled = mit_enabled

        # read in substance properties from the config file, and initialize basic properties:
        self.read_substances(sim, cells, config_substances, p)
        self.tissue_init(sim, cells, config_substances, p)

        # # write substance growth and decay equations:
        # self.write_growth_and_decay()

        self.ave_cell_vol = cells.cell_vol.mean()  # average cell volume

        # colormap for plotting series of 1D lines:
        self.plot_cmap = 'viridis'

        # default shape identities for network:
        self.reaction_shape = 'rect'
        self.transporter_shape = 'diamond'
        self.channel_shape = 'pentagon'
        self.vmem_shape = 'ellipse'
        self.ed_shape = 'hexagon'
        self.conc_shape = 'oval'

        self.extra_J_mem = np.zeros(sim.mdl)
        self.extra_J_env = np.zeros(sim.edl)


    #------------Initializers-------------------------------------------------------------------------------------------
    def read_substances(self, sim, cells, config_substances, p):
        """
        Initialize all core data structures and concentration variables for all
        molecules included in the simulation, as well as any ions present in
        sim.

        config_substances:  dictionary containing BETSE biomolecule template fields
        """

        logs.log_info("Reading additional substance data...")

        # Initialize a dictionaries that will eventually hold dynamic values for cell, env and mit concentrations:
        cell_concs_mapping = {}
        mem_concs_mapping = {}
        env_concs_mapping = {}
        bound_concs_mapping = {}

        # place to store charge values and membrane diffusivity for all substances AND sim ions:
        self.zmol = {}
        self.Dmem = {}
        self.ED_eval_strings ={}


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
                    lambda value, ion_index=ion_index: sim.cc_cells.__setitem__(ion_index, value))

                mem_concs_mapping[k] = DynamicValue(
                    lambda ion_index=ion_index: sim.cc_at_mem[ion_index],
                    lambda value, ion_index=ion_index: sim.cc_at_mem.__setitem__(ion_index, value))

                env_concs_mapping[k] = DynamicValue(
                    lambda ion_index=ion_index: sim.cc_env[ion_index],
                    lambda value, ion_index=ion_index: sim.cc_env.__setitem__(ion_index, value))

                bound_concs_mapping[k] = DynamicValue(
                    lambda ion_index=ion_index: sim.c_env_bound[ion_index],
                    lambda value, ion_index=ion_index: sim.c_env_bound.__setitem__(ion_index, value))

                self.zmol[k] = sim.zs[ion_index]
                self.Dmem[k] = sim.Dm_cells[ion_index].mean()

                if self.mit_enabled:
                    mit_concs_mapping[k] = DynamicValue(
                        lambda ion_index=ion_index: sim.cc_mit[ion_index],
                        lambda value, ion_index=ion_index: sim.cc_mit.__setitem__(ion_index, value))



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

            # determine if the user has requested an asymmetric initial condition:
            init_asym = mol_dic.get('initial asymmetry', None)

            if init_asym != 'None' and init_asym is not None:

                mol.init_asym, _ = getattr(mods, init_asym)(cells.cell_i, cells, p)

            else:

                mol.init_asym = np.ones(sim.cdl)


            # create concentration data arrays:
            mol.c_cells = np.ones(sim.cdl) * mol.c_cello*mol.init_asym
            # mol.c_mems = np.ones(sim.mdl) * mol.c_cello

            #FIXME: Add explicit checks to avoid name conflicts between
            #user-defined substances defined below and ions defined above.
            #Na, K, Cl, Ca, H, M, P == reserved names for ions!

            # create dynamic mappings for the cell conc vectors:
            cell_concs_mapping[name] = DynamicValue(
                lambda name = name: self.molecules[name].c_cells,
                lambda value, name = name: setattr(self.molecules[name], 'c_cells', value))

            # initialize array for data stored at membranes:
            mol.cc_at_mem = mol.c_cells[cells.mem_to_cells]

            # create dynamic mappings for the mem conc vectors:
            mem_concs_mapping[name] = DynamicValue(
                lambda name = name: self.molecules[name].cc_at_mem,
                lambda value, name = name: setattr(self.molecules[name], 'cc_at_mem', value))

            # storage array for intracellular flux:
            mol.flux_intra = np.zeros(sim.mdl)

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
                mol.Dgj = mol_dic.get('Dgj', 1.0e-16)  # effective diffusion coefficient of substance through GJ
                mol.z = mol_dic['z']  # charge (oxidation state)


                # factors describing potential transport along current-aligned microtubules
                mol.u_mt = float(mol_dic.get('u_mtube', 0.0))

                self.zmol[name] = mol.z
                self.Dmem[name] = mol.Dm

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
                    # mol.Kgd = gad['Km']
                    # mol.n_decay = gad['n']

                    mol.growth_profiles_list = gad['apply to']

                    modulator_function_name = gad.get('modulator function', None)

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

                    mol.init_growth(cells, p)

                    # the modulator function, if requested, makes gad occur modulated by a function.
                    # Make this happen, if it's requested:

                    if modulator_function_name != 'None' and modulator_function_name is not None:
                        # mol.growth_mod_function_mems, _ = getattr(mods, modulator_function_name)(cells.mem_i,
                        #                                                                           cells, p)
                        mol.growth_mod_function_cells, _ = getattr(mods, modulator_function_name)(cells.cell_i,
                                                                                                   cells, p)

                    else:
                        # mol.growth_mod_function_mems = np.ones(sim.mdl)
                        mol.growth_mod_function_cells = np.ones(sim.cdl)

                else:
                    mol.simple_growth = False
                    mol.growth_profiles_list = None

                # assign ion channel gating properties
                icg = mol_dic.get('ion channel gating', None)

                if icg is not None:

                    mol.ion_channel_gating = True

                    mol.gating_channel_name = icg.get('channel name', 'Gated channel')

                    mol.gating_ion_name = icg['ion channel target']  # get a target ion label to gate membrane to (or 'None')

                    if mol.gating_ion_name != 'None':
                        mol.use_gating_ligand = True

                        mol.gating_ion = []

                        for ion_o in mol.gating_ion_name:
                            mol.gating_ion.append(sim.get_ion(ion_o))

                    else:
                        mol.use_gating_ligand = False
                        mol.gating_ion = []

                    mol.gating_Hill_K = float(icg['target Hill coefficient'])
                    mol.gating_Hill_n = float(icg['target Hill exponent'])
                    mol.gating_max_val = float(icg['peak channel opening'])
                    mol.gating_extracell = icg['acts extracellularly']

                    # get any optional activators and inhibitors for the channel:
                    mol.ion_activators_list = icg.get('activators', None)
                    mol.ion_activators_Km = icg.get('Km activators', None)
                    mol.ion_activators_n = icg.get('n activators', None)
                    mol.ion_activators_zone = icg.get('zone activators', None)

                    mol.ion_inhibitors_list = icg.get('inhibitors', None)
                    mol.ion_inhibitors_Km = icg.get('Km inhibitors', None)
                    mol.ion_inhibitors_n = icg.get('n inhibitors', None)
                    mol.ion_inhibitors_zone = icg.get('zone inhibitors', None)

                    # Fill in the blanks if zones aren't supplied:
                    mol.ion_activators_zone, mol.ion_inhibitors_zone = self.default_zones(
                                                                            mol.ion_activators_zone,
                                                                            mol.ion_inhibitors_zone,
                                                                            mol.ion_activators_list,
                                                                            mol.ion_inhibitors_list)

                    tex_vars = []

                    # define the modulating coefficients:
                    alpha_ion, alpha_tex, tex_vars = \
                        self.get_influencers(mol.ion_activators_list,
                        mol.ion_activators_Km, mol.ion_activators_n,
                        mol.ion_inhibitors_list, mol.ion_inhibitors_Km,
                        mol.ion_inhibitors_n, reaction_zone='mem', tex_list = tex_vars,
                        zone_tags_a=mol.ion_activators_zone,
                        zone_tags_i=mol.ion_inhibitors_zone,
                        in_mem_tag=True)

                    mol.gating_mod_eval_string = "(" + alpha_ion + ")"


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

        # balance charge in cells, env, and mits
        if p.substances_affect_charge and self.charge_has_been_balanced is False:
            # calculate the net charge in cells
            Qcells = 0
            Qenv = 0
            Qmit = 0

            for name in self.molecules:

                mol = self.molecules[name]

                Qcells += mol.c_cello * mol.z
                Qenv += mol.c_envo * mol.z

                if self.mit_enabled:
                    Qmit += mol.c_mito * mol.z

            self.bal_charge(Qcells, sim, 'cell', p)
            self.bal_charge(Qenv, sim, 'env', p)

            if self.mit_enabled:
                self.bal_charge(Qmit, sim, 'mit', p)

            self.charge_has_been_balanced = True

        # write substance growth and decay equations:
        self.write_growth_and_decay()

        # write passive electrodiffusion expressions for each ion and molecule:
        logs.log_info("Writing passive electrodiffusion equations...")
        for k, v in self.cell_concs.items():

            name = "'{}'".format(k)

            if self.zmol[k] != 0.0 and self.Dmem[k] != 0.0:

                if p.sim_ECM is True:

                    self.ED_eval_strings[k] = "(self.Dmem[{}]/p.tm)*((p.F*sim.vm)/(p.R*sim.T))*self.zmol[{}]*" \
                                      "((self.mem_concs[{}] - " \
                                      "self.env_concs[{}][cells.map_mem2ecm]*" \
                                      "np.exp(-(self.zmol[{}]*p.F*sim.vm)" \
                                      "/(p.R*sim.T)))/(1-np.exp(-(self.zmol[{}]*" \
                                      "p.F*sim.vm)/(p.R*sim.T))))".format(name, name, name, name, name, name)

                else:

                    self.ED_eval_strings[k] = "(self.Dmem[{}]/p.tm)*((p.F*sim.vm)/(p.R*sim.T))*self.zmol[{}]*" \
                                              "((self.mem_concs[{}] - " \
                                              "self.env_concs[{}][cells.mem_to_cells]*" \
                                              "np.exp(-(self.zmol[{}]*p.F*sim.vm)" \
                                              "/(p.R*sim.T)))/(1-np.exp(-(self.zmol[{}]*" \
                                              "p.F*sim.vm)/(p.R*sim.T))))".format(name, name, name, name, name, name)



            elif self.zmol[k] == 0.0 and self.Dmem[k] != 0.0:

                if p.sim_ECM is True:

                    self.ED_eval_strings[k] = "(self.Dmem[{}]/p.tm)*(self.mem_concs[{}] - " \
                                          "self.env_concs[{}][cells.map_mem2ecm])".format(name, name, name)

                else:

                    self.ED_eval_strings[k] = "(self.Dmem[{}]/p.tm)*(self.mem_concs[{}] - " \
                                          "self.env_concs[{}][cells.mem_to_cells])".format(name, name, name)

            elif self.Dmem[k] == 0.0:

                self.ED_eval_strings[k] = "np.zeros(sim.mdl)"

    def init_network(self, sim, cells, p):

        # initialize network optimization, etc:
        if p.network_config is not None:

            opti = p.network_config['optimize network']

            self.opti_method = p.network_config['optimization method']
            self.opti_N = p.network_config['optimization steps']
            # after primary initialization, check and see if optimization required:

            if opti is True:
                logs.log_info("The Metabolic Network is being analyzed for optimal rates...")
                self.optimizer(sim, cells, p)

    def build_indices(self, sim, cells, p):

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

            if self.mit_enabled:

                mit_name = mol_name + '_mit'

                if self.mit_concs[mol_name] is None:

                    self.conc_handler[mit_name] = 0

                else:
                    self.conc_handler[mit_name] = self.mit_concs[mol_name]

        # build a global dictionary of all reactions:---------------------------------------------------------
        self.react_handler = OrderedDict({})

        # define a string term to scale membrane-fluxes to cell-wide concentration change:
        div_string = "*(cells.cell_sa.mean()/cells.cell_vol.mean())"
        # div_string = '*1.0'

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

        if self.mit_enabled:
            for rea_name_mit in self.reactions_mit:
                self.react_handler[rea_name_mit] = self.reactions_mit[rea_name_mit].reaction_eval_string

        for trans_name in self.transporters:
            self.react_handler[trans_name] = self.transporters[trans_name].transporter_eval_string + div_string

        for chan_name in self.channels:

            self.react_handler[chan_name] = "self.channels['{}'].channel_core.chan_flux".format(chan_name) + div_string

        self.react_handler_index = {}

        for i, rea_name in enumerate(self.react_handler):
            self.react_handler_index[rea_name] = i

        self.conc_handler_index = {}

        for i, conc_name in enumerate(self.conc_handler):
            self.conc_handler_index[conc_name] = i

    def plot_init(self, config_dic, p):

        config_substances = config_dic['biomolecules']

        # read in network plotting options:
        self.net_plot_opts = config_dic.get('network plotting', None)

        # set plotting options for the network:
        set_net_opts(self, self.net_plot_opts, p)

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
            raise BetseSimConfigException(
                'The "plot cell" defined in the "results" section of your '
                'configuration file does not exist in your cluster. '
                'Choose a plot cell number smaller than the maximum cell number.')

    def read_reactions(self, config_reactions, sim, cells, p):

        """
          Read in and initialize parameters for all user-defined reactions.

          config_options:  dictionary containing BETSE reaction template fields

          """

        logs.log_info("Reading reaction input data...")

        # Initialize a dict that keeps the Reaction objects for reactions occuring in the 'cell' zone:
        self.reactions = OrderedDict({})
        self.reactions_env = OrderedDict({})

        if self.mit_enabled is True:
            # Initialize a dict that keeps the Reaction objects for reactions occuring in the 'mit' zone:
            self.reactions_mit = OrderedDict({})

        for q, react_dic in enumerate(config_reactions):

            # get each user-defined name-filed in the dictionary:
            name = str(react_dic['name'])

            zone = react_dic['reaction zone']

            if zone == 'cell':
                # add a Reaction object to MasterOfReactions reaction dictionary:
                self.reactions[name] = Reaction(sim, cells, p)
                obj = self.reactions[name]

            elif zone == 'mit' and self.mit_enabled is True:
                # add a Reaction object to MasterOfReactions reaction dictionary:
                self.reactions_mit[name] = Reaction(sim, cells, p)
                obj = self.reactions_mit[name]

            elif zone == 'env':

                self.reactions_env[name] = Reaction(sim, cells, p)
                obj = self.reactions_env[name]

            else:
                logs.log_warning("---------------------------------------------------------------")
                logs.log_warning("Reaction {} defined on inappropriate zone (or mit not enabled)".format(name))
                logs.log_warning("Reverting this reaction to the 'cell' zone!")
                logs.log_warning("---------------------------------------------------------------")
                # add a Reaction object to MasterOfReactions reaction dictionary:
                self.reactions[name] = Reaction(sim, cells, p)
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


            if self.mit_enabled:
                obj.mit_enabled = True
            else:
                obj.mit_enabled = False

        for name in self.reactions:

            msg = "Including the cell-zone reaction: {}".format(name)
            logs.log_info(msg)

        if self.mit_enabled is True:

            for name in self.reactions_mit:

                msg = "Including the mit-zone reaction: {}".format(name)
                logs.log_info(msg)

    def read_transporters(self, config_transporters, sim, cells, p):

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

            # Fill in the blanks if zones aren't supplied:
            obj.transporter_activators_zone, obj.transporter_inhibitors_zone = self.default_zones(
                obj.transporter_activators_zone,
                obj.transporter_inhibitors_zone,
                obj.transporter_activators_list,
                obj.transporter_inhibitors_list)

            if self.mit_enabled:
                obj.mit_enabled = True
            else:
                obj.mit_enabled = False

        for name in self.transporters:

            msg = "Including the network transporter: {}".format(name)
            logs.log_info(msg)

    def read_channels(self, config_channels, sim, cells, p):

        logs.log_info("Reading channel input data...")

        # Initialize a dict of Channels:
        self.channels = OrderedDict({})

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

            tex_vars = []

            all_alpha, alpha_tex, tex_vars = self.get_influencers(
                                                                    a_list, Km_a_list, n_a_list, i_list,
                                                                    Km_i_list, n_i_list, tex_list = tex_vars,
                                                                    reaction_zone='mem', zone_tags_a=zone_a,
                                                                    zone_tags_i=zone_i, in_mem_tag=True)


            obj.alpha_eval_string = "(" + all_alpha + ")"

        # for name in self.channels:
        #
        #     msg = "Including the network channel: {}".format(name)
        #     logs.log_info(msg)

    def read_modulators(self, config_modulators, sim, cells, p):

        logs.log_info("Reading modulator input data...")

        # Initialize a dict of modulators:
        self.modulators = OrderedDict({})

        for q, mod_dic in enumerate(config_modulators):

            name = str(mod_dic['name'])

            # add a field to the MasterOfReactions corresponding to Channel object with that name:
            self.modulators[name] = Modulator()

            # now set the attributes of that channel:

            # get MasterOfMolecules.name
            obj = self.modulators[name]

            obj.target_label = str(mod_dic['target'])
            # obj.zone = str(mod_dic['zone'])
            obj.max_val = float(mod_dic['max effect'])
            obj.modulator_activators_list = mod_dic.get('activators', None)
            obj.modulator_activators_Km = mod_dic.get('activator Km', None)
            obj.modulator_activators_n = mod_dic.get('activator n', None)
            obj.modulator_activators_zone = mod_dic.get('activator zone', None)

            obj.modulator_inhibitors_list = mod_dic.get('inhibitors', None)
            obj.modulator_inhibitors_Km = mod_dic.get('inhibitor Km', None)
            obj.modulator_inhibitors_n = mod_dic.get('inhibitor n', None)
            obj.modulator_inhibitors_zone = mod_dic.get('inhibitor zone', None)

            # Fill in the blanks if zones aren't supplied:
            obj.modulator_activators_zone, obj.modulator_inhibitors_zone = self.default_zones(
                obj.modulator_activators_zone,
                obj.modulator_inhibitors_zone,
                obj.modulator_activators_list,
                obj.modulator_inhibitors_list)

            obj.init_modulator(sim, cells, p)

            # activator/inhibitor lists and associated data:
            a_list = obj.modulator_activators_list
            Km_a_list = obj.modulator_activators_Km
            n_a_list = obj.modulator_activators_n
            zone_a = obj.modulator_activators_zone

            i_list = obj.modulator_inhibitors_list
            Km_i_list = obj.modulator_inhibitors_Km
            n_i_list = obj.modulator_inhibitors_n
            zone_i = obj.modulator_inhibitors_zone

            tex_vars = []

            all_alpha, alpha_tex, tex_vars = self.get_influencers(a_list, Km_a_list,
                                                                    n_a_list, i_list,
                                                                    Km_i_list, n_i_list, tex_list=tex_vars,
                                                                    reaction_zone='mem', zone_tags_a = zone_a,
                                                                    zone_tags_i=zone_i, in_mem_tag=False)

            obj.alpha_eval_string = "(" + all_alpha + ")"

        for name in self.modulators:
            msg = "Including the modulation: {}".format(name)
            logs.log_info(msg)

    def write_growth_and_decay(self):

        logs.log_info("Writing substance growth/decay equations...")

        for mol_name in self.molecules:


            if self.molecules[mol_name].simple_growth is True:

                # initialize an empty list that will hold strings defining fixed parameter values as LaTeX math string
                gad_tex_var_list = []


                cc = "(self.molecules['{}'].c_cells)".format(mol_name)

                cc_tex = r"[%s]" % (mol_name)

                # Kgd = self.molecules[mol_name].Kgd
                # n_d = self.molecules[mol_name].n_decay

                # write fixed parameter values to LaTeX:
                # kval = tex_val(Kgd)
                # Kgd_tex = "K_{%s}^{decay} & =" % (mol_name)
                # Kgd_tex += kval

                # nval = tex_val(n_d)
                # n_d_tex = "n_{%s}^{decay} & =" % (mol_name)
                # n_d_tex += nval

                # gad_tex_var_list.append(Kgd_tex)
                # gad_tex_var_list.append(n_d_tex)

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
                             i_list, Km_i_list, n_i_list, tex_list = gad_tex_var_list,reaction_zone='cell',
                             zone_tags_a=zone_a, zone_tags_i=zone_i, in_mem_tag=False)

                # define remaining portion of substance's growth and decay expression:

                mod_funk = "(self.molecules['{}'].growth_mod_function_cells)".format(mol_name)
                r_prod =  "(self.molecules['{}'].r_production)".format(mol_name)
                r_decay = "(self.molecules['{}'].r_decay)".format(mol_name)

                gad_eval_string = mod_funk + "*" + r_prod + "*" + all_alpha + "-" + \
                                  r_decay + "*" + cc

                # Formatting for LaTeX equation export ---------------------------------------

                r_prod_tex = r"r_{%s}^{max}\," % mol_name
                r_dec_tex = r" - \delta_{%s}\," % mol_name

                gad_tex_initial = r"r_{%s}=" % mol_name

                gad_tex_string = gad_tex_initial + r_prod_tex + alpha_tex + r_dec_tex + cc_tex

                # write fixed parameter names and values to the LaTeX storage list:
                rpval = tex_val(self.molecules[mol_name].r_production)
                r_p_tex = "r_{%s}^{max} & =" % (mol_name)
                r_p_tex += rpval

                rdval = tex_val(self.molecules[mol_name].r_decay)
                r_d_tex = r"\delta_{%s} & =" % (mol_name)
                r_d_tex += rdval

                gad_tex_var_list.append(r_p_tex)
                gad_tex_var_list.append(r_d_tex)

                #------------------------------------------------------------------------------

                gad_eval_string_growth = all_alpha
                gad_eval_string_decay = r_decay + "*" + cc

                self.molecules[mol_name].gad_eval_string = gad_eval_string

                # add in the separate growth decay eval strings for case of optimization:
                self.molecules[mol_name].gad_eval_string_growth = gad_eval_string_growth
                self.molecules[mol_name].gad_eval_string_decay = gad_eval_string_decay

                # add in the tex string describing the growth and decay expression:
                self.molecules[mol_name].gad_tex_string = gad_tex_string

                # finalize the fixed parameters LaTeX list:
                tex_params = r"\begin{aligned}"

                for i, tex_str in enumerate(gad_tex_var_list):

                    tex_params += tex_str

                    if i < len(gad_tex_var_list) -1:

                        tex_params += r"\\"

                    else:

                        tex_params += r"\end{aligned}"

                self.molecules[mol_name].gad_tex_vars = tex_params

            else:

                self.molecules[mol_name].gad_eval_string = "np.zeros(sim.cdl)"
                self.molecules[mol_name].gad_eval_string_growth = "np.zeros(sim.cdl)"
                self.molecules[mol_name].gad_eval_string_decay = "np.zeros(sim.cdl)"

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

            # initialize an empty list that will hold strings defining fixed parameter values as LaTeX math string
            rea_tex_var_list = []

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

            numo_tex_Q = ""
            deno_tex_Q = ""

            for i, (name, coeff) in enumerate(zip(reactant_names, reactant_coeff)):

                tex_name = name

                denomo_string_Q += "(self.cell_concs['{}']".format(name)

                denomo_string_Q += "**{})".format(coeff)

                deno_tex_Q += "[%s]^{%.1f}" % (tex_name, coeff)

                if i < len(reactant_names) - 1:

                    denomo_string_Q += "*"
                    deno_tex_Q += r"\,"

                else:

                    denomo_string_Q += ")"


            for i, (name, coeff) in enumerate(zip(product_names, product_coeff)):

                tex_name = name
                numo_string_Q += "(self.cell_concs['{}']".format(name)

                numo_string_Q += "**{})".format(coeff)

                numo_tex_Q += "[%s]^{%.1f}" % (tex_name, coeff)

                if i < len(product_names) - 1:

                    numo_string_Q += "*"
                    numo_tex_Q += r"\,"

                else:

                    numo_string_Q += ")"

            # define the final reaction quotient string, Q:
            Q = "(" + numo_string_Q + '/' + denomo_string_Q + ")"

            Q_tex = r"\frac{%s}{%s}" % (numo_tex_Q, deno_tex_Q)  # LaTeX version

            # next calculate the forward and backward reaction rate coefficients:---------------------------------------

            forward_coeff = "("
            backward_coeff = "("

            fwd_coeff_tex = ""
            bwd_coeff_tex = ""

            for i, (name, n, Km) in enumerate(zip(reactant_names, reactant_coeff, reactant_Km)):

                tex_name = name
                numo_string_r = "((self.cell_concs['{}']/{})**{})".format(name, Km, n)
                denomo_string_r = "(1 + (self.cell_concs['{}']/{})**{})".format(name, Km, n)

                cc_tex = r"\left(\frac{[%s]}{K_{%s_f}^{%s}}\right)" % (tex_name, reaction_name, tex_name)
                tex_term = r"\left(\frac{%s}{1 + %s}\right)" % (cc_tex, cc_tex)

                term = "(" + numo_string_r + "/" + denomo_string_r + ")"

                forward_coeff += term

                # tex equation:

                fwd_coeff_tex += tex_term

                # write fixed parameter names and values to the LaTeX storage list:
                kmval = tex_val(Km)
                km_tex = "K_{%s_f}^{%s} & =" % (reaction_name, name)
                km_tex += kmval

                nval = tex_val(n)
                n_tex = "n_{%s_f}^{%s} & =" % (reaction_name, name)
                n_tex += nval

                rea_tex_var_list.append(km_tex)
                rea_tex_var_list.append(n_tex)

                if i < len(reactant_names) - 1:

                    forward_coeff += "*"

                    fwd_coeff_tex += r"\,"

                else:

                    forward_coeff += ")"

            for i, (name, n, Km) in enumerate(zip(product_names, product_coeff, product_Km)):

                tex_name = name

                numo_string_p = "((self.cell_concs['{}']/{})**{})".format(name, Km, n)
                denomo_string_p = "(1 + (self.cell_concs['{}']/{})**{})".format(name, Km, n)

                term = "(" + numo_string_p + "/" + denomo_string_p + ")"

                cc_tex = r"\left(\frac{[%s]}{K_{%s_r}^{%s}}\right)" % (tex_name, reaction_name, tex_name)
                tex_term = r"\left(\frac{%s}{1 + %s}\right)" % (cc_tex, cc_tex)

                backward_coeff += term

                bwd_coeff_tex += tex_term

                if self.reactions[reaction_name].delta_Go is not None:
                    # write fixed parameter names and values to the LaTeX storage list:
                    kmval = tex_val(Km)
                    km_tex = "K_{%s_r}^{%s} & =" % (reaction_name, name)
                    km_tex += kmval

                    nval = tex_val(n)
                    n_tex = "n_{%s_r}^{%s} & =" % (reaction_name, name)
                    n_tex += nval

                    rea_tex_var_list.append(km_tex)
                    rea_tex_var_list.append(n_tex)

                if i < len(product_names) - 1:

                    backward_coeff += "*"

                    bwd_coeff_tex += r"\,"

                else:

                    backward_coeff += ")"

            # if reaction is reversible deal calculate an equilibrium constant:
            if self.reactions[reaction_name].delta_Go is not None:

                # define the reaction equilibrium coefficient expression:
                Keqm = "(np.exp(-self.reactions['{}'].delta_Go / (p.R * sim.T)))".format(reaction_name)

                Keqm_tex = r"exp\left(-\frac{deltaG_{%s}^{o}}{R\,T}\right)" % reaction_name

                # write fixed parameter names and values to the LaTeX storage list:
                gval = tex_val(self.reactions[reaction_name].delta_Go)
                g_tex = "deltaG_{%s}^{o} & =" % reaction_name
                g_tex += gval
                rea_tex_var_list.append(g_tex)

            else:

                Q = "0"
                backward_coeff = "0"
                Keqm = "1"

                Q_tex = ""
                Keqm_tex = ""
                bwd_coeff_tex = ""


            # Put it all together into a final reaction string (+max rate):

            reversed_term = "(" + Q + "/" + Keqm + ")"

            reversed_term_tex = r"\frac{%s}{%s}" % (Q_tex, Keqm_tex)

            all_alpha, alpha_tex, rea_tex_var_list = self.get_influencers(a_list, Km_a_list,
                                                                    n_a_list,
                                                                    i_list, Km_i_list, n_i_list,
                                                                    tex_list = rea_tex_var_list,
                                                                    reaction_zone= r_zone, zone_tags_a=zone_a,
                                                                    zone_tags_i=zone_i, in_mem_tag= False)

            vmax = "self.reactions['{}'].vmax".format(reaction_name)

            reaction_eval_string = vmax + "*" + all_alpha + "*" + "(" + \
                                   forward_coeff + "-" + "(" + reversed_term + "*" + backward_coeff + ")" + ")"

            r_tex = "r_{%s}^{max}" % reaction_name
            r_main = "r_{%s} = " % reaction_name

            if self.reactions[reaction_name].delta_Go is not None:

                reaction_tex_string = r_main + r_tex + r"\," + alpha_tex + r"\,\left(" + fwd_coeff_tex + "-" +  \
                                      reversed_term_tex + r"\," + bwd_coeff_tex + r"\right)"

            else:

                reaction_tex_string = r_main + r_tex + r"\," + alpha_tex + r"\," + fwd_coeff_tex

            # write fixed parameter names and values to the LaTeX storage list:
            rval = tex_val(self.reactions[reaction_name].vmax)
            r_tex = "r_{%s}^{max} & =" % (reaction_name)
            r_tex += rval
            rea_tex_var_list.append(r_tex)

            # finalize the fixed parameters LaTeX list:
            tex_params = r"\begin{aligned}"
            for i, tex_str in enumerate(rea_tex_var_list):
                tex_params += tex_str
                if i < len(rea_tex_var_list) - 1:
                    tex_params += r"\\"
                else:
                    tex_params += r"\end{aligned}"

            self.reactions[reaction_name].reaction_tex_vars = tex_params


            # add the composite string describing the reaction math to a new field:
            self.reactions[reaction_name].reaction_eval_string = reaction_eval_string

            self.reactions[reaction_name].reaction_tex_string = reaction_tex_string

    def write_reactions_mit(self):
        """
        Reactions are now constructed during the init as strings that are evaluated in eval calls in each time-step.
        This function constructs the evaluation strings for each reaction, given the metadata stored
        in each reaction object (e.g. lists of reactants, products, etc).

        """

        logs.log_info("Writing reaction equations for mit zone...")

        for reaction_name in self.reactions_mit:

            # initialize an empty list that will hold strings defining fixed parameter values as LaTeX math string
            rea_tex_var_list = []

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

            numo_tex_Q = ""
            deno_tex_Q = ""

            for i, (name, coeff) in enumerate(zip(reactant_names, reactant_coeff)):

                tex_name = name + '_{mit}'

                denomo_string_Q += "(self.mit_concs['{}']".format(name)

                denomo_string_Q += "**{})".format(coeff)

                deno_tex_Q += "[%s]^{%.1f}" % (tex_name, coeff)

                if i < len(reactant_names) - 1:

                    denomo_string_Q += "*"
                    deno_tex_Q += r"\,"

                else:

                    denomo_string_Q += ")"


            for i, (name, coeff) in enumerate(zip(product_names, product_coeff)):

                tex_name = name + "_{mit}"
                numo_string_Q += "(self.mit_concs['{}']".format(name)

                numo_string_Q += "**{})".format(coeff)

                numo_tex_Q += "[%s]^{%.1f}" % (tex_name, coeff)

                if i < len(product_names) - 1:

                    numo_string_Q += "*"
                    numo_tex_Q += r"\,"

                else:

                    numo_string_Q += ")"

            # define the final reaction quotient string, Q:
            Q = "(" + numo_string_Q + '/' + denomo_string_Q + ")"

            Q_tex = r"\frac{%s}{%s}" % (numo_tex_Q, deno_tex_Q)  # LaTeX version

            # next calculate the forward and backward reaction rate coefficients:---------------------------------------

            forward_coeff = "("
            backward_coeff = "("

            fwd_coeff_tex = ""
            bwd_coeff_tex = ""

            for i, (name, n, Km) in enumerate(zip(reactant_names, reactant_coeff, reactant_Km)):

                tex_name = name + "_{mit}"

                numo_string_r = "((self.mit_concs['{}']/{})**{})".format(name, Km, n)
                denomo_string_r = "(1 + (self.mit_concs['{}']/{})**{})".format(name, Km, n)

                cc_tex = r"\left(\frac{[%s]}{K_{%s_f}^{%s}}\right)" % (tex_name, reaction_name, tex_name)
                tex_term = r"\left(\frac{%s}{1 + %s}\right)" % (cc_tex, cc_tex)

                term = "(" + numo_string_r + "/" + denomo_string_r + ")"

                forward_coeff += term

                # tex equation:

                fwd_coeff_tex += tex_term

                # write fixed parameter names and values to the LaTeX storage list:
                kmval = tex_val(Km)
                km_tex = "K_{%s_f}^{%s} & =" % (reaction_name, name)
                km_tex += kmval

                nval = tex_val(n)
                n_tex = "n_{%s_f}^{%s} & =" % (reaction_name, name)
                n_tex += nval

                rea_tex_var_list.append(km_tex)
                rea_tex_var_list.append(n_tex)

                if i < len(reactant_names) - 1:

                    forward_coeff += "*"

                    fwd_coeff_tex += r"\,"

                else:

                    forward_coeff += ")"

            for i, (name, n, Km) in enumerate(zip(product_names, product_coeff, product_Km)):

                tex_name = name + '_{mit}'

                numo_string_p = "((self.mit_concs['{}']/{})**{})".format(name, Km, n)
                denomo_string_p = "(1 + (self.mit_concs['{}']/{})**{})".format(name, Km, n)

                term = "(" + numo_string_p + "/" + denomo_string_p + ")"

                cc_tex = r"\left(\frac{[%s]}{K_{%s_r}^{%s}}\right)" % (tex_name, reaction_name, tex_name)
                tex_term = r"\left(\frac{%s}{1 + %s}\right)" % (cc_tex, cc_tex)

                backward_coeff += term

                bwd_coeff_tex += tex_term

                if self.reactions_mit[reaction_name].delta_Go is not None:
                    # write fixed parameter names and values to the LaTeX storage list:
                    kmval = tex_val(Km)
                    km_tex = "K_{%s_r}^{%s} & =" % (reaction_name, name)
                    km_tex += kmval

                    nval = tex_val(n)
                    n_tex = "n_{%s_r}^{%s} & =" % (reaction_name, name)
                    n_tex += nval

                    rea_tex_var_list.append(km_tex)
                    rea_tex_var_list.append(n_tex)

                if i < len(product_names) - 1:

                    backward_coeff += "*"

                    bwd_coeff_tex += r"\,"

                else:

                    backward_coeff += ")"

            # if reaction is reversible deal calculate an equilibrium constant:
            if self.reactions_mit[reaction_name].delta_Go is not None:

                # define the reaction equilibrium coefficient expression:
                Keqm = "(np.exp(-self.reactions_mit['{}'].delta_Go / (p.R * sim.T)))".format(reaction_name)

                Keqm_tex = r"exp\left(-\frac{deltaG_{%s}^{o}}{R\,T}\right)" % reaction_name

                # write fixed parameter names and values to the LaTeX storage list:
                gval = tex_val(self.reactions_mit[reaction_name].delta_Go)
                g_tex = "deltaG_{%s}^{o} & =" % reaction_name
                g_tex += gval
                rea_tex_var_list.append(g_tex)

            else:

                Q = "0"
                backward_coeff = "0"
                Keqm = "1"

                Q_tex = ""
                Keqm_tex = ""
                bwd_coeff_tex = ""


            # Put it all together into a final reaction string (+max rate):

            reversed_term = "(" + Q + "/" + Keqm + ")"

            reversed_term_tex = r"\frac{%s}{%s}" % (Q_tex, Keqm_tex)

            all_alpha, alpha_tex, rea_tex_var_list = self.get_influencers(a_list, Km_a_list,
                                                                    n_a_list,
                                                                    i_list, Km_i_list, n_i_list,
                                                                    tex_list = rea_tex_var_list,
                                                                    reaction_zone= r_zone, zone_tags_a=zone_a,
                                                                    zone_tags_i=zone_i, in_mem_tag=False)

            vmax = "self.reactions_mit['{}'].vmax".format(reaction_name)

            reaction_eval_string = vmax + "*" + all_alpha + "*" + "(" + \
                                   forward_coeff + "-" + "(" + reversed_term + "*" + backward_coeff + ")" + ")"

            r_tex = "r_{%s}^{max}" % reaction_name
            r_main = "r_{%s} = " % reaction_name

            if self.reactions_mit[reaction_name].delta_Go is not None:

                reaction_tex_string = r_main + r_tex + r"\," + alpha_tex + r"\,\left(" + fwd_coeff_tex + "-" +  \
                                      reversed_term_tex + r"\," + bwd_coeff_tex + r"\right)"

            else:

                reaction_tex_string = r_main + r_tex + r"\," + alpha_tex + r"\," + fwd_coeff_tex

            # write fixed parameter names and values to the LaTeX storage list:
            rval = tex_val(self.reactions_mit[reaction_name].vmax)
            r_tex = "r_{%s}^{max} & =" % (reaction_name)
            r_tex += rval
            rea_tex_var_list.append(r_tex)

            # finalize the fixed parameters LaTeX list:
            tex_params = r"\begin{aligned}"
            for i, tex_str in enumerate(rea_tex_var_list):
                tex_params += tex_str
                if i < len(rea_tex_var_list) - 1:
                    tex_params += r"\\"
                else:
                    tex_params += r"\end{aligned}"

            self.reactions_mit[reaction_name].reaction_tex_vars = tex_params


            # add the composite string describing the reaction math to a new field:
            self.reactions_mit[reaction_name].reaction_eval_string = reaction_eval_string

            self.reactions_mit[reaction_name].reaction_tex_string = reaction_tex_string

    def write_reactions_env(self):

        """
        Reactions are now constructed during the init as strings that are evaluated in eval calls in each time-step.
        This function constructs the evaluation strings for each reaction, given the metadata stored
        in each reaction object (e.g. lists of reactants, products, etc).

        """

        logs.log_info("Writing reaction equations for env zone...")

        for reaction_name in self.reactions_env:

            # initialize an empty list that will hold strings defining fixed parameter values as LaTeX math string
            rea_tex_var_list = []

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

            numo_tex_Q = ""
            deno_tex_Q = ""

            for i, (name, coeff) in enumerate(zip(reactant_names, reactant_coeff)):

                tex_name = name

                denomo_string_Q += "(self.env_concs['{}']".format(name)

                denomo_string_Q += "**{})".format(coeff)

                deno_tex_Q += "[%s]^{%.1f}" % (tex_name, coeff)

                if i < len(reactant_names) - 1:

                    denomo_string_Q += "*"
                    deno_tex_Q += r"\,"

                else:

                    denomo_string_Q += ")"

            for i, (name, coeff) in enumerate(zip(product_names, product_coeff)):

                tex_name = name
                numo_string_Q += "(self.env_concs['{}']".format(name)

                numo_string_Q += "**{})".format(coeff)

                numo_tex_Q += "[%s]^{%.1f}" % (tex_name, coeff)

                if i < len(product_names) - 1:

                    numo_string_Q += "*"
                    numo_tex_Q += r"\,"

                else:

                    numo_string_Q += ")"

            # define the final reaction quotient string, Q:
            Q = "(" + numo_string_Q + '/' + denomo_string_Q + ")"

            Q_tex = r"\frac{%s}{%s}" % (numo_tex_Q, deno_tex_Q)  # LaTeX version

            # next calculate the forward and backward reaction rate coefficients:---------------------------------------

            forward_coeff = "("
            backward_coeff = "("

            fwd_coeff_tex = ""
            bwd_coeff_tex = ""

            for i, (name, n, Km) in enumerate(zip(reactant_names, reactant_coeff, reactant_Km)):

                tex_name = name
                numo_string_r = "((self.env_concs['{}']/{})**{})".format(name, Km, n)
                denomo_string_r = "(1 + (self.env_concs['{}']/{})**{})".format(name, Km, n)

                cc_tex = r"\left(\frac{[%s]}{K_{%s_f}^{%s}}\right)" % (tex_name, reaction_name, tex_name)
                tex_term = r"\left(\frac{%s}{1 + %s}\right)" % (cc_tex, cc_tex)

                term = "(" + numo_string_r + "/" + denomo_string_r + ")"

                forward_coeff += term

                # tex equation:

                fwd_coeff_tex += tex_term

                # write fixed parameter names and values to the LaTeX storage list:
                kmval = tex_val(Km)
                km_tex = "K_{%s_f}^{%s} & =" % (reaction_name, name)
                km_tex += kmval

                nval = tex_val(n)
                n_tex = "n_{%s_f}^{%s} & =" % (reaction_name, name)
                n_tex += nval

                rea_tex_var_list.append(km_tex)
                rea_tex_var_list.append(n_tex)

                if i < len(reactant_names) - 1:

                    forward_coeff += "*"

                    fwd_coeff_tex += r"\,"

                else:

                    forward_coeff += ")"

            for i, (name, n, Km) in enumerate(zip(product_names, product_coeff, product_Km)):

                tex_name = name

                numo_string_p = "((self.env_concs['{}']/{})**{})".format(name, Km, n)
                denomo_string_p = "(1 + (self.env_concs['{}']/{})**{})".format(name, Km, n)

                term = "(" + numo_string_p + "/" + denomo_string_p + ")"

                cc_tex = r"\left(\frac{[%s]}{K_{%s_r}^{%s}}\right)" % (tex_name, reaction_name, tex_name)
                tex_term = r"\left(\frac{%s}{1 + %s}\right)" % (cc_tex, cc_tex)

                backward_coeff += term

                bwd_coeff_tex += tex_term

                if self.reactions_env[reaction_name].delta_Go is not None:
                    # write fixed parameter names and values to the LaTeX storage list:
                    kmval = tex_val(Km)
                    km_tex = "K_{%s_r}^{%s} & =" % (reaction_name, name)
                    km_tex += kmval

                    nval = tex_val(n)
                    n_tex = "n_{%s_r}^{%s} & =" % (reaction_name, name)
                    n_tex += nval

                    rea_tex_var_list.append(km_tex)
                    rea_tex_var_list.append(n_tex)

                if i < len(product_names) - 1:

                    backward_coeff += "*"

                    bwd_coeff_tex += r"\,"

                else:

                    backward_coeff += ")"

            # if reaction is reversible deal calculate an equilibrium constant:
            if self.reactions_env[reaction_name].delta_Go is not None:

                # define the reaction equilibrium coefficient expression:
                Keqm = "(np.exp(-self.reactions_env['{}'].delta_Go / (p.R * sim.T)))".format(reaction_name)

                Keqm_tex = r"exp\left(-\frac{deltaG_{%s}^{o}}{R\,T}\right)" % reaction_name

                # write fixed parameter names and values to the LaTeX storage list:
                gval = tex_val(self.reactions_env[reaction_name].delta_Go)
                g_tex = "deltaG_{%s}^{o} & =" % reaction_name
                g_tex += gval
                rea_tex_var_list.append(g_tex)

            else:

                Q = "0"
                backward_coeff = "0"
                Keqm = "1"

                Q_tex = ""
                Keqm_tex = ""
                bwd_coeff_tex = ""

            # Put it all together into a final reaction string (+max rate):

            reversed_term = "(" + Q + "/" + Keqm + ")"

            reversed_term_tex = r"\frac{%s}{%s}" % (Q_tex, Keqm_tex)

            all_alpha, alpha_tex, rea_tex_var_list = self.get_influencers(a_list, Km_a_list,
                                                                          n_a_list,
                                                                          i_list, Km_i_list, n_i_list,
                                                                          tex_list=rea_tex_var_list,
                                                                          reaction_zone=r_zone, zone_tags_a=zone_a,
                                                                          zone_tags_i=zone_i, in_mem_tag=False)

            vmax = "self.reactions_env['{}'].vmax".format(reaction_name)

            reaction_eval_string = vmax + "*" + all_alpha + "*" + "(" + \
                                   forward_coeff + "-" + "(" + reversed_term + "*" + backward_coeff + ")" + ")"

            r_tex = "r_{%s}^{max}" % reaction_name
            r_main = "r_{%s} = " % reaction_name

            if self.reactions_env[reaction_name].delta_Go is not None:

                reaction_tex_string = r_main + r_tex + r"\," + alpha_tex + r"\,\left(" + fwd_coeff_tex + "-" + \
                                      reversed_term_tex + r"\," + bwd_coeff_tex + r"\right)"

            else:

                reaction_tex_string = r_main + r_tex + r"\," + alpha_tex + r"\," + fwd_coeff_tex

            # write fixed parameter names and values to the LaTeX storage list:
            rval = tex_val(self.reactions_env[reaction_name].vmax)
            r_tex = "r_{%s}^{max} & =" % (reaction_name)
            r_tex += rval
            rea_tex_var_list.append(r_tex)

            # finalize the fixed parameters LaTeX list:
            tex_params = r"\begin{aligned}"
            for i, tex_str in enumerate(rea_tex_var_list):
                tex_params += tex_str
                if i < len(rea_tex_var_list) - 1:
                    tex_params += r"\\"
                else:
                    tex_params += r"\end{aligned}"

            self.reactions_env[reaction_name].reaction_tex_vars = tex_params

            # add the composite string describing the reaction math to a new field:
            self.reactions_env[reaction_name].reaction_eval_string = reaction_eval_string

            self.reactions_env[reaction_name].reaction_tex_string = reaction_tex_string

    def write_transporters(self, sim, cells, p):
        """
        Reactions are now constructed during the init as strings that are evaluated in eval calls in each time-step.
        This function constructs the evaluation strings for each reaction, given the metadata stored
        in each reaction object (e.g. lists of reactants, products, etc).

        """

        logs.log_info("Writing transporter equations...")

        for transp_name in self.transporters:

            # initialize an empty list that will hold strings defining fixed parameter values as LaTeX math string
            trans_tex_var_list = []

            # initialize a field that will report net charge transported
            self.transporters[transp_name].net_z = 0

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
            echem_terms_list_tex = []

            if reaction_zone == 'cell':

                type_out = 'env_concs'

                vmem = "sim.vm"   # get the transmembrane voltage for this category

                vmem_tex = "V_{mem}"

                in_delta_term_react = "(np.dot(cells.M_sum_mems, -self.transporters['{}'].flux*cells.mem_sa)/cells.cell_vol)".format(transp_name)
                in_delta_term_prod = "(np.dot(cells.M_sum_mems, self.transporters['{}'].flux*cells.mem_sa)/cells.cell_vol)".format(transp_name)

                if p.sim_ECM is True:

                    # if self.transporters[transp_name].ignore_ECM_transporter is True:

                    out_delta_term_react = "-self.transporters['{}'].flux*(cells.memSa_per_envSquare" \
                                           "[cells.map_mem2ecm]/cells.ecm_vol)".format(transp_name)

                    out_delta_term_prod = "self.transporters['{}'].flux*(cells.memSa_per_envSquare" \
                                          "[cells.map_mem2ecm]/cells.ecm_vol)".format(transp_name)

                else:
                    out_delta_term_react = "(np.dot(cells.M_sum_mems, -self.transporters['{}'].flux*cells.mem_sa)/cells.cell_vol)".format(transp_name)

                    out_delta_term_prod = "(np.dot(cells.M_sum_mems, self.transporters['{}'].flux*cells.mem_sa)/cells.cell_vol)".format(transp_name)

                all_alpha, alpha_tex, trans_tex_var_list = self.get_influencers(a_list, Km_a_list,
                                                                n_a_list, i_list,
                                                                Km_i_list, n_i_list, tex_list=trans_tex_var_list,
                                                                reaction_zone='mem', in_mem_tag=False)

            elif reaction_zone == 'mit' and self.mit_enabled is True:

                vmem_tex = "V_{mit}"

                # initialize lists to hold the reactants and product transfer tags (initialized to zone):
                react_transfer_tag = ['mit_concs' for x in reactant_names]
                prod_transfer_tag = ['mit_concs' for x in product_names]

                type_out = 'cell_concs'

                vmem = "self.mit.Vmit"  # get the transmembrane voltage for this category

                in_delta_term_react = "-self.transporters['{}'].flux*(self.mit.mit_sa/self.mit.mit_vol)".\
                                                                        format(transp_name)

                in_delta_term_prod = "self.transporters['{}'].flux*(self.mit.mit_sa/self.mit.mit_vol)".\
                                                                format(transp_name)

                out_delta_term_react = "-self.transporters['{}'].flux*(self.mit.mit_sa/cells.cell_vol)". \
                                                                        format(transp_name)

                out_delta_term_prod = "self.transporters['{}'].flux*(self.mit.mit_sa/cells.cell_vol)". \
                                                                     format(transp_name)

                all_alpha, alpha_tex, trans_tex_var_list = self.get_influencers(a_list, Km_a_list,
                                                                n_a_list, i_list,
                                                                Km_i_list, n_i_list, reaction_zone='mit',
                                                                tex_list=trans_tex_var_list, in_mem_tag=False)

        # initialize list that will hold expression for calculating net concentration change

            delta_strings_reactants = [in_delta_term_react for x in reactant_names]
            delta_strings_products = [in_delta_term_prod for x in product_names]

            # calculate a reactants zone tag list and terms affecting transporter free energy via Vmem:

            if transport_out_list != 'None':

                for out_name in transport_out_list:

                    if out_name in p.ions_dict:

                        if p.ions_dict[out_name] == 1:

                            ion_ind = sim.get_ion(out_name)

                            # get the index for the substance in the molecules database
                            prod_i = product_names.index(out_name)

                            prod_transfer_tag[prod_i] = type_out

                            coeff = product_coeff[prod_i]

                            # calculate the effect of transfer on transporter's free energy:
                            eterm = "-{}*sim.zs[{}]*p.F*{}".format(coeff, ion_ind, vmem)

                            # add term to net charge list:
                            self.transporters[transp_name].net_z += coeff*sim.zs[ion_ind]

                            # LaTeX version:
                            eterm_tex = r"-%s\,[%s]\,z_{%s}\,F\,%s" % (coeff, out_name, out_name, vmem_tex)

                            # write fixed parameter names and values to the LaTeX storage list:
                            zval = tex_val(sim.zs[ion_ind])
                            z_tex = "z_{%s} & =" % (out_name)
                            z_tex += zval

                    else:

                        # get the index for the substance in the molecules database
                        prod_i = product_names.index(out_name)

                        prod_transfer_tag[prod_i] = type_out

                        coeff = product_coeff[prod_i]

                        # calculate the effect of transfer on transporter's free energy:
                        eterm = "-{}*self.molecules['{}'].z*p.F*{}".format(coeff, out_name, vmem)

                        self.transporters[transp_name].net_z += coeff*self.molecules[out_name].z

                        # LaTeX version:
                        eterm_tex = r"-%s\,[%s]\,z_{%s}\,F\,%s" % (coeff, out_name, out_name, vmem_tex)

                        # write fixed parameter names and values to the LaTeX storage list:
                        zval = tex_val(self.molecules[out_name].z)
                        z_tex = "z_{%s} & =" % (out_name)
                        z_tex += zval

                    trans_tex_var_list.append(z_tex)

                    echem_terms_list.append(eterm)
                    echem_terms_list_tex.append(eterm_tex)

                    # update delta string for correct transfer:
                    delta_strings_products[prod_i] = out_delta_term_prod


            if transport_in_list != 'None':

                for in_name in transport_in_list:

                    if in_name in p.ions_dict:

                        if p.ions_dict[in_name] == 1:

                            ion_ind = sim.get_ion(in_name)

                            # get the index for the substance in the molecules database
                            react_i = reactant_names.index(in_name)
                            react_transfer_tag[react_i] = type_out

                            coeff = reactant_coeff[react_i]

                            # calculate the effect of transfer on transporter's free energy:
                            eterm = "{}*sim.zs[{}]*p.F*{}".format(coeff, ion_ind, vmem)

                            self.transporters[transp_name].net_z += -coeff * sim.zs[ion_ind]

                            # LaTeX version:
                            eterm_tex = r"%s\,[%s]\,z_{%s}\,F\,%s" % (coeff, in_name, in_name, vmem_tex)

                            # write fixed parameter names and values to the LaTeX storage list:
                            zval = tex_val(sim.zs[ion_ind])
                            z_tex = "z_{%s} & =" % (in_name)
                            z_tex += zval

                    else:

                        # get the index for the substance in the molecules database
                        react_i = reactant_names.index(in_name)
                        react_transfer_tag[react_i] = type_out

                        coeff = reactant_coeff[react_i]

                        # calculate the effect of transfer on transporter's free energy:
                        eterm = "{}*self.molecules['{}'].z*p.F*{}".format(coeff, in_name, vmem)

                        self.transporters[transp_name].net_z += -coeff * self.molecules[in_name].z

                        # LaTeX version:
                        eterm_tex = r"%s\,[%s]\,z_{%s}\,F\,%s" % (coeff, in_name, in_name, vmem_tex)

                        # write fixed parameter names and values to the LaTeX storage list:
                        zval = tex_val(self.molecules[in_name].z)
                        z_tex = "z_{%s} & =" % (in_name)
                        z_tex += zval

                    trans_tex_var_list.append(z_tex)

                    echem_terms_list.append(eterm)
                    echem_terms_list_tex.append(eterm_tex)

                    # update delta string for correct transfer:
                    delta_strings_reactants[react_i] = out_delta_term_react

            # create the eterms string expression describing net effect of trans-membrane fluxes on free energy:
            echem_string = "("
            echem_tex_string = ""

            for i, (et, etx) in enumerate(zip(echem_terms_list, echem_terms_list_tex)):

                echem_string += et
                echem_tex_string += etx

                if i < len(echem_terms_list) -1:

                    echem_string += "+"
                    echem_tex_string += "+"


                else:
                    echem_string += ")"

            # first calculate a reaction coefficient Q, as a string expression
            numo_string_Q = "("
            denomo_string_Q = "("

            numo_tex = ""
            deno_tex = ""

            for i, (name, coeff, tag) in enumerate(zip(reactant_names, reactant_coeff, react_transfer_tag)):

                if tag == 'mem_concs':

                    denomo_string_Q += "(self.{}['{}']".format(tag, name)
                    tex_name = name

                elif tag == 'mit_concs':
                    denomo_string_Q += "(self.{}['{}']".format(tag, name)
                    tex_name = name + '_{mit}'

                elif tag == 'cell_concs':
                    denomo_string_Q += "(self.{}['{}']".format(tag, name)
                    tex_name = name

                # get the concentration from the environment, mapped to respective membranes:
                elif tag == 'env_concs':

                    if p.sim_ECM is True:

                        denomo_string_Q += "(self.{}['{}'][cells.map_mem2ecm]".format(tag, name)

                    else:
                        denomo_string_Q += "(self.{}['{}']".format(tag, name)

                    tex_name = name + '_{env}'

                else:

                    raise BetseSimConfigException("Transporter tag not properly defined!")

                denomo_string_Q += "**{})".format(coeff)

                deno_tex += "[%s]^{%s}" % (tex_name, coeff)

                if i < len(reactant_names) - 1:

                    denomo_string_Q += "*"
                    deno_tex += r"\,"

                else:

                    denomo_string_Q += ")"

            for i, (name, coeff, tag) in enumerate(zip(product_names, product_coeff, prod_transfer_tag)):

                if tag == 'mem_concs':

                    numo_string_Q += "(self.{}['{}']".format(tag, name)
                    tex_name = name

                elif tag == 'mit_concs':

                    numo_string_Q += "(self.{}['{}']".format(tag, name)
                    tex_name = name + "_{mit}"

                elif tag == 'cell_concs':

                    numo_string_Q += "(self.{}['{}']".format(tag, name)
                    tex_name = name

                # get the concentration from the environment mapped to the respective membranes:
                elif tag == 'env_concs':

                    if p.sim_ECM is True:
                        numo_string_Q += "(self.{}['{}'][cells.map_mem2ecm]".format(tag, name)

                    else:
                        numo_string_Q += "(self.{}['{}']".format(tag, name)

                    tex_name = name + "_{env}"

                else:

                    raise BetseSimConfigException("Transporter tag not properly defined!")

                numo_string_Q += "**{})".format(coeff)

                numo_tex += "[%s]^{%s}" % (tex_name, coeff)

                if i < len(product_names) - 1:

                    numo_string_Q += "*"

                    numo_tex += r"\,"

                else:

                    numo_string_Q += ")"

            # define the final reaction quotient string, Q:
            Q = "(" + numo_string_Q + '/' + denomo_string_Q + ")"

            Q_tex = r"\frac{%s}{%s}" % (numo_tex, deno_tex)

            # next calculate the forward and backward reaction rate coefficients:---------------------------------------
            forward_coeff = "("
            backward_coeff = "("

            fwd_tex_coeff = ""
            bwd_tex_coeff = ""

            for i, (name, n, Km, tag) in enumerate(zip(reactant_names, reactant_coeff, reactant_Km, react_transfer_tag)):

                if tag == 'mem_concs':

                    numo_string_r = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    denomo_string_r = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    tex_name = name

                elif tag == 'mit_concs':

                    numo_string_r = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    denomo_string_r = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)

                    tex_name = name + "_{mit}"

                elif tag == 'cell_concs':

                    numo_string_r = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    denomo_string_r = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    tex_name = name

                elif tag == 'env_concs':

                    if p.sim_ECM is True:

                        numo_string_r = "((self.{}['{}'][cells.map_mem2ecm]/{})**{})".format(tag, name, Km, n)
                        denomo_string_r = "(1 + (self.{}['{}'][cells.map_mem2ecm]/{})**{})".format(tag, name, Km, n)

                    else:

                        numo_string_r = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                        denomo_string_r = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)

                    tex_name = name + "_{env}"

                else:

                    raise BetseSimConfigException("Transporter tag not properly defined!")

                term = "(" + numo_string_r + "/" + denomo_string_r + ")"

                forward_coeff += term

                cc_tex = r"\left(\frac{[%s]}{K_{%s_f}^{%s}}\right)^{n_{%s_f}^{%s}}" % (tex_name, transp_name, tex_name,
                                                                                    transp_name, tex_name)
                fwd_tex_coeff += r"\left(\frac{%s}{1+%s}\right)" % (cc_tex, cc_tex)

                if i < len(reactant_names) - 1:

                    forward_coeff += "*"
                    fwd_tex_coeff += r"\,"

                else:

                    forward_coeff += ")"

                # write fixed parameter names and values to the LaTeX storage list:
                kval = tex_val(Km)
                k_tex = "K_{%s}^{%s} & =" % (transp_name, tex_name)
                k_tex += kval
                trans_tex_var_list.append(k_tex)

                nval = tex_val(n)
                n_tex = "n_{%s}^{%s} & =" % (transp_name, tex_name)
                n_tex += nval
                trans_tex_var_list.append(n_tex)

            for i, (name, n, Km, tag) in enumerate(zip(product_names, product_coeff, product_Km, prod_transfer_tag)):

                if tag == 'mem_concs':

                    # tag2 = 'cell_concs'

                    numo_string_p = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    denomo_string_p = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)

                    tex_name = name

                elif tag == 'mit_concs':

                    numo_string_p = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    denomo_string_p = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)

                    tex_name = name + "_{mit}"

                elif tag == 'cell_concs':

                    numo_string_p = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                    denomo_string_p = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)

                    tex_name = name

                elif tag == 'env_concs':

                    if p.sim_ECM is True:

                        numo_string_p = "((self.{}['{}'][cells.map_mem2ecm]/{})**{})".format(tag, name, Km, n)
                        denomo_string_p = "(1 + (self.{}['{}'][cells.map_mem2ecm]/{})**{})".format(tag, name, Km, n)

                    else:

                        numo_string_p = "((self.{}['{}']/{})**{})".format(tag, name, Km, n)
                        denomo_string_p = "(1 + (self.{}['{}']/{})**{})".format(tag, name, Km, n)

                    tex_name = name + "_{env}"

                else:

                    raise BetseSimConfigException("Transporter tag not properly defined!")


                term = "(" + numo_string_p + "/" + denomo_string_p + ")"

                backward_coeff += term

                # LaTeX version:
                cc_tex = r"\left(\frac{[%s]}{K_{%s_r}^{%s}}\right)^{n_{%s_r}^{%s}}" % (tex_name, transp_name, tex_name,
                                                                                    transp_name, tex_name)
                bwd_tex_coeff += r"\left(\frac{%s}{1+%s}\right)" % (cc_tex, cc_tex)

                if i < len(product_names) - 1:

                    backward_coeff += "*"

                    bwd_tex_coeff += r"\,"

                else:

                    backward_coeff += ")"

                if self.transporters[transp_name].delta_Go is not None:

                    # write fixed parameter names and values to the LaTeX storage list:
                    kval = tex_val(Km)
                    k_tex = "K_{%s}^{%s} & =" % (transp_name, tex_name)
                    k_tex += kval
                    trans_tex_var_list.append(k_tex)

                    nval = tex_val(n)
                    n_tex = "n_{%s}^{%s} & =" % (transp_name, tex_name)
                    n_tex += nval
                    trans_tex_var_list.append(n_tex)

            # if reaction is reversible, calculate an equilibrium constant:
            if self.transporters[transp_name].delta_Go is not None:

                # define the reaction equilibrium coefficient expression:
                Keqm = "(np.exp(-(self.transporters['{}'].delta_Go + {})/ (p.R * sim.T)))".format(transp_name,
                                                                                                    echem_string)

                Keqm_tex = r"\frac{exp(-deltaG_{%s}^{o} + %s)}{R\,T}" % (transp_name, echem_tex_string)

                # write fixed parameter names and values to the LaTeX storage list:
                gval = tex_val(self.transporters[transp_name].delta_Go)
                g_tex = "deltaG_{%s}^{o} & =" % (transp_name)
                g_tex += gval
                trans_tex_var_list.append(g_tex)

            else:

                Q = "0"
                backward_coeff = "0"
                Keqm = "1"

                Keqm_tex = ""

            # Put it all together into a final reaction string (+max rate):
            reversed_term = "(" + Q + "/" + Keqm + ")"

            rev_term_tex = r"\left(\frac{%s}{%s}\right)" % (Q_tex, Keqm_tex)

            vmax = "self.transporters['{}'].vmax".format(transp_name)

            # Write the LaTeX version and append fixed parameter name to tex list:
            v_texo = "r_{%s}^{max}" % (transp_name)
            vm_tex = "r_{%s}" % (transp_name)
            vval = tex_val(self.transporters[transp_name].vmax)
            v_tex = "v_{%s}^{max} & =" % (transp_name)
            v_tex += vval
            trans_tex_var_list.append(v_tex)

            # write the final LaTeX expressions:
            if self.transporters[transp_name].delta_Go is not None:

                transporter_tex_string = vm_tex + " = " + v_texo + r"\," + alpha_tex + r"\,\left(" + fwd_tex_coeff + \
                                         "-" + rev_term_tex + r"\," + bwd_tex_coeff + r"\right)"

                # calculate the evaluation string expression for the transporter:
                transporter_eval_string = vmax + "*" + all_alpha + "*" + "(" + \
                                          forward_coeff + "-"  + reversed_term + "*" + backward_coeff  + ")"


            else:
                transporter_tex_string = vm_tex + " = " + v_tex + r"\," + alpha_tex + r"\," + fwd_tex_coeff

                # calculate the evaluation string expression for the transporter:
                transporter_eval_string = vmax + "*" + all_alpha + "*" + forward_coeff


            # finalize the fixed parameters LaTeX list:
            tex_params = r"\begin{aligned}"

            for i, tex_str in enumerate(trans_tex_var_list):

                tex_params += tex_str

                if i < len(trans_tex_var_list)-1:
                    tex_params += r"\\"
                elif i == len(trans_tex_var_list)-1:
                    tex_params += r"\end{aligned}"

            # add the composite string describing the reaction math to a new field:
            # Evaluate this expression to assign self.transporters[transp_name].flux field
            self.transporters[transp_name].transporter_eval_string = transporter_eval_string

            # write LaTeX expression for transporter:
            self.transporters[transp_name].transporter_tex_string = transporter_tex_string
            self.transporters[transp_name].transporter_tex_vars = tex_params

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

    def create_reaction_matrix_mit(self):

        """
        Creates a matrix representation of the system of coupled
        differential equations underlying the reaction network for mitochondria.

        """

        logs.log_info("Writing reaction network matrix for mit zone...")
        n_reacts = len(self.reactions_mit)
        n_mols = len(self.mit_concs)

        # initialize the network's reaction matrix:
        self.reaction_matrix_mit = np.zeros((n_mols, n_reacts))

        molecule_keys = list(self.mit_concs.keys())

        for j, reaction_name in enumerate(self.reactions_mit):

            for react_name, coeff in zip(self.reactions_mit[reaction_name].reactants_list,
                self.reactions_mit[reaction_name].reactants_coeff):
                i = molecule_keys.index(react_name)

                self.reaction_matrix_mit[i, j] += -coeff

            for prod_name, coeff in zip(self.reactions_mit[reaction_name].products_list,
                self.reactions_mit[reaction_name].products_coeff):
                i = molecule_keys.index(prod_name)

                self.reaction_matrix_mit[i, j] += coeff

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

    #------runners------------------------------------------------------------------------------------------------------
    def run_loop(self, t, sim, cells, p):
        """
        Runs the main simulation loop steps for each of the molecules included in the simulation.

        """

        gad_rates_o = []

        gad_targs = []

        init_rates = []

        gad_rates = []

        self.extra_rho_cells = np.zeros(sim.cdl)
        self.extra_rho_mems = np.zeros(sim.mdl)
        self.extra_rho_env = np.zeros(sim.edl)
        self.extra_rho_mit = np.zeros(sim.cdl)
        self.extra_J_mem = np.zeros(sim.mdl)
        self.extra_J_env = np.zeros(sim.edl)
        
        globalo = globals()
        localo = locals()

        for mol in self.molecules:

            # calculate concentrations at membranes:
            obj = self.molecules[mol]
            obj.update_intra(sim, cells, p)

            # calculate rates of growth/decay:
            gad_rates_o.append(eval(self.molecules[mol].gad_eval_string, globalo, localo))

            gad_targs.append(self.molecules[mol].growth_targets_cell)

            init_rates.append(np.zeros(sim.cdl))

        for mat, trgs, rts in zip(init_rates, gad_targs, gad_rates_o):

            mat[trgs] = rts[trgs]

            gad_rates.append(mat)

        gad_rates = np.asarray(gad_rates)

        # ... and rates of chemical reactions in cell:
        self.reaction_rates = np.asarray(
            [eval(self.reactions[rn].reaction_eval_string, globalo, localo) for rn in self.reactions])

        # stack into an integrated data structure:
        if len(self.reaction_rates) > 0:
            all_rates = np.vstack((gad_rates, self.reaction_rates))

        else:
            all_rates = gad_rates

        # calculate concentration rate of change using linear algebra:
        self.delta_conc = np.dot(self.reaction_matrix, all_rates)

        if self.mit_enabled and len(self.reactions_mit)>0:
            # ... rates of chemical reactions in mitochondria:
            self.reaction_rates_mit = np.asarray(
                [eval(self.reactions_mit[rn].reaction_eval_string, globalo, localo) for
                    rn in self.reactions_mit])

            # calculate concentration rate of change using linear algebra:
            self.delta_conc_mit = np.dot(self.reaction_matrix_mit, self.reaction_rates_mit)

        # update environmental concentrations:
        if len(self.reactions_env)>0:
            # ... rates of chemical reactions in env:
            self.reaction_rates_env = np.asarray(
                [eval(self.reactions_env[rn].reaction_eval_string, globalo, localo) for
                    rn in self.reactions_env])

            # calculate concentration rate of change using linear algebra:
            self.delta_conc_env = np.dot(self.reaction_matrix_env, self.reaction_rates_env)

        else:

            self.delta_conc_env = None

        if self.delta_conc_env is not None:

            for name_env, deltae in zip(self.env_concs, self.delta_conc_env):

                conce = self.env_concs[name_env]
                self.env_concs[name_env] = conce + deltae*p.dt

        for ii, (name, deltac) in enumerate(zip(self.cell_concs, self.delta_conc)):

            conco = self.cell_concs[name]

            if name not in p.ions_dict:

                obj = self.molecules[name]

            else:
                obj = None

            # update concentration due to growth/decay and chemical reactions:
            self.cell_concs[name] = conco + deltac * p.dt


            if obj is not None:

                # if pumping is enabled:
                if obj.active_pumping:
                    obj.pump(sim, cells, p)

            # if p.run_sim is True:
            # use the substance as a gating ligand (if desired)

                if obj.ion_channel_gating:
                    obj.gating_mod = eval(obj.gating_mod_eval_string, globalo, localo)
                    obj.gating(sim, cells, p)

                if p.run_sim is True:
                    # update the global boundary (if desired)
                    if obj.change_bounds:
                        obj.update_boundary(t, p)

                # transport the molecule through gap junctions and environment:
                obj.transport(sim, cells, p)

                # print(obj.name, obj.f_mem.mean())

                # ensure no negs:
                stb.no_negs(obj.c_cells)
                stb.no_negs(obj.c_env)

                if p.substances_affect_charge:
                    # calculate the charge density this substance contributes to cell and environment:
                    self.extra_rho_cells[:] += p.F*obj.c_cells*obj.z
                    self.extra_rho_env[:] += p.F*obj.c_env*obj.z
                    self.extra_J_mem[:] +=  -obj.z*obj.f_mem*p.F + obj.z*obj.f_gj*p.F
                    self.extra_rho_mems[:] += p.F*obj.cc_at_mem*obj.z


            if self.mit_enabled and len(self.reactions_mit)>0:

                concm = self.mit_concs[name]
                # update concentration due to chemical reactions in mitochondria:
                self.mit_concs[name] = concm + self.delta_conc_mit[ii]*p.dt

                # ensure no negative values:
                stb.no_negs(self.mit_concs[name])


            if self.mit_enabled and obj is not None:
                # calculate the charge density this substance contributes to mit:
                self.extra_rho_mit[:] += p.F*obj.c_mit*obj.z

        # calculate energy charge in the cell:
        self.energy_charge(sim)

        if p.substances_affect_charge:
            sim.extra_rho_cells = self.extra_rho_cells
            sim.extra_rho_mems = self.extra_rho_mems*cells.diviterm[cells.mem_to_cells]
            sim.extra_rho_env = self.extra_rho_env
            sim.extra_J_mem = self.extra_J_mem


        if self.mit_enabled:  # if enabled, update the mitochondria's voltage and other properties
            self.mit.extra_rho = self.extra_rho_mit[:]
            self.mit.update(sim, cells, p)

    def run_loop_transporters(self, t, sim, cells, p):

        globalo = globals()
        localo = locals()

        # call statement to evaluate:
        for name in self.transporters:

            # specific tissue profile regions where the transporter is active:
            targ_mem = self.transporters[name].transporter_targets_mem
            targ_cell = self.transporters[name].transporter_targets_cell
            targ_env = self.transporters[name].transporter_targets_env

            # calculate the flux
            self.transporters[name].flux = sim.rho_pump*eval(self.transporters[name].transporter_eval_string,
                globalo, localo)


            self.extra_J_mem += self.transporters[name].net_z*self.transporters[name].flux*p.F
            sim.extra_J_mem = self.extra_J_mem

            # finally, update the concentrations using the final eval statements:
            for i, (delc, coeff) in enumerate(zip(self.transporters[name].delta_react_eval_strings,
                self.transporters[name].reactants_coeff)):

                # obtain the change for the reactant

                delta_react = coeff*eval(delc, globalo, localo)

                # finally, update the concentrations using the final eval statements:
                if self.transporters[name].react_transport_tag[i] == 'mem_concs':

                    # print("--> ", self.transporters[name].reactants_list[i], delta_react.mean())

                    self.cell_concs[self.transporters[name].reactants_list[i]][targ_cell] = \
                        self.cell_concs[self.transporters[name].reactants_list[i]][targ_cell] + \
                        delta_react[targ_cell]*p.dt

                    self.mem_concs[self.transporters[name].reactants_list[i]][targ_mem] = \
                        self.mem_concs[self.transporters[name].reactants_list[i]][targ_mem] - \
                        self.transporters[name].flux[targ_mem]*(cells.mem_sa[targ_mem]/cells.mem_vol[targ_mem])*p.dt

                elif self.transporters[name].react_transport_tag[i] == 'env_concs':

                    if p.sim_ECM is True:

                        delta_react_expanded = np.zeros(sim.edl)
                        delta_react_expanded[cells.map_mem2ecm] = delta_react[:]

                        self.env_concs[self.transporters[name].reactants_list[i]][targ_env] = (
                            self.env_concs[self.transporters[name].reactants_list[i]][targ_env] +
                            delta_react_expanded[targ_env]*p.dt)

                    else:

                        self.env_concs[self.transporters[name].reactants_list[i]][targ_mem] = \
                            self.env_concs[self.transporters[name].reactants_list[i]][targ_mem] + \
                            delta_react[targ_cell][cells.mem_to_cells] * p.dt

                elif self.transporters[name].react_transport_tag[i] == 'cell_concs':

                    self.cell_concs[self.transporters[name].reactants_list[i]][targ_cell] = \
                        self.cell_concs[self.transporters[name].reactants_list[i]][targ_cell] + \
                        delta_react[targ_cell]*p.dt

                elif self.transporters[name].react_transport_tag[i] == 'mit_concs':

                    self.mit_concs[self.transporters[name].reactants_list[i]][targ_cell] = \
                        self.mit_concs[self.transporters[name].reactants_list[i]][targ_cell] + \
                        delta_react[targ_cell] * p.dt

                else:

                    raise BetseSimConfigException("Internal error: transporter zone not specified correctly!")

            for i, (delc, coeff) in enumerate(zip(self.transporters[name].delta_prod_eval_strings,
                self.transporters[name].products_coeff)):

                # obtain the change for the product
                delta_prod = coeff*eval(delc, globalo, localo)

                # finally, update the concentrations using the final eval statements:
                if self.transporters[name].prod_transport_tag[i] == 'mem_concs':

                    self.cell_concs[self.transporters[name].products_list[i]][targ_cell] = \
                        self.cell_concs[self.transporters[name].products_list[i]][targ_cell] + \
                        delta_prod[targ_cell]*p.dt

                    self.mem_concs[self.transporters[name].products_list[i]][targ_mem] = \
                        self.mem_concs[self.transporters[name].products_list[i]][targ_mem] + \
                        self.transporters[name].flux[targ_mem]*(cells.mem_sa[targ_mem]/cells.mem_vol[targ_mem])*p.dt



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
                            delta_prod[targ_cell][cells.mem_to_cells] * p.dt

                elif self.transporters[name].prod_transport_tag[i] == 'cell_concs':

                    self.cell_concs[self.transporters[name].products_list[i]][targ_cell] = \
                        self.cell_concs[self.transporters[name].products_list[i]][targ_cell] + \
                        delta_prod[targ_cell]*p.dt

                elif self.transporters[name].prod_transport_tag[i] == 'mit_concs':

                    self.mit_concs[self.transporters[name].products_list[i]][targ_cell] = \
                        self.mit_concs[self.transporters[name].products_list[i]][targ_cell] + \
                        delta_prod[targ_cell] * p.dt

                else:

                    raise BetseSimConfigException("Internal error: transporter zone not specified correctly!")

            # print("---------------------------------")

    def run_loop_channels(self, sim, cells, p):

        globalo = globals()
        localo = locals()

        # get the object corresponding to the specific transporter:
        for i, name in enumerate(self.channels):

            # compute the channel activity
            # calculate the value of the channel modulation constant:
            moddy = eval(self.channels[name].alpha_eval_string, globalo, localo)

            self.channels[name].channel_core.modulator = moddy

            self.channels[name].channel_core.run(self.channels[name].dummy_dyna, sim, cells, p)

    def run_loop_modulators(self, sim, cells, p):

        globalo = globals()
        localo = locals()

        # get the object corresponding to the specific transporter:
        for i, name in enumerate(self.modulators):

            obj = self.modulators[name]

            # calculate the value of the channel modulation constant:
            modulator = obj.max_val*eval(obj.alpha_eval_string, globalo, localo)

            # # make size alteration for case of true environment:
            # if p.sim_ECM is True and obj.zone == 'env':
            #     modulator = modulator[cells.map_mem2ecm]

            # print(modulator.min())
            # print("-----------")

            if obj.target_label == 'gj':

                sim.gj_block = modulator

            elif obj.target_label == 'Na/K-ATPase':

                sim.NaKATP_block = modulator

            elif obj.target_label == 'H/K-ATPase':

                sim.HKATP_block = modulator

            elif obj.target_label == 'V-ATPase':

                sim.VATP_block = modulator

            elif obj.target_label == 'Ca-ATPase':

                sim.CaATP_block = modulator

            elif obj.target_label == 'Na/Ca-Exch':

                sim.NaCaExch_block = modulator

            else:

                raise BetseSimConfigException("You have requested a "
                                               "sim modulator that is not "
                                               "available. Available choices "
                                               "are: 'gj', 'Na/K-ATPase', 'H/K-ATPase', "
                                               "and 'V-ATPase', 'Ca-ATPase', and 'Na/Ca-Exch' ")

    # ------Utility Methods--------------------------------------------------------------------------------------------
    def bal_charge(self, Q, sim, tag, p):

        if tag == 'cell':

            if Q < 0 and np.abs(Q) <= sim.cc_cells[sim.iP].mean():  # if net charge is anionic

                sim.cc_cells[sim.iP] = sim.cc_cells[sim.iP] - np.abs(Q)

            elif Q > 0 and np.abs(Q) <= sim.cc_cells[sim.iK].mean():

                sim.cc_cells[sim.iK] = sim.cc_cells[sim.iK] - np.abs(Q)

            elif Q < 0 and np.abs(Q) > sim.cc_mems[sim.iP].mean():  # if net charge is anionic
                raise BetseSimConfigException("You've defined way more anionic charge in"
                                               "the extra substances than we can "
                                               "compensate for. Either turn 'substances "
                                               "affect Vmem' off, or try again.")
            elif Q > 0 and np.abs(Q) > sim.cc_cells[sim.iK].mean():
                raise BetseSimConfigException("You've defined way more cationic charge in"
                                               "the extra substances than we can "
                                               "compensate for. Either turn 'substances "
                                               "affect Vmem' off, or try again.")

            sim.extra_rho_cells = p.F*Q*np.ones(sim.cdl)

        elif tag == 'env':

            if Q < 0 and np.abs(Q) <= sim.cc_env[sim.iP].mean():  # if net charge is anionic
                sim.cc_env[sim.iP] = sim.cc_env[sim.iP] - np.abs(Q)
            elif Q > 0 and np.abs(Q) <= sim.cc_env[sim.iK].mean():
                sim.cc_env[sim.iK] = sim.cc_env[sim.iK] - np.abs(Q)
            elif Q < 0 and np.abs(Q) > sim.cc_env[sim.iP].mean():  # if net charge is anionic
                raise BetseSimConfigException("You've defined way more anionic charge in"
                                               "the extra substances than we can "
                                               "compensate for. Either turn 'substances "
                                               "affect Vmem' off, or try again.")
            elif Q > 0 and np.abs(Q) > sim.cc_env[sim.iK].mean():
                raise BetseSimConfigException("You've defined way more cationic charge in"
                                               "the extra substances than we can "
                                               "compensate for. Either turn 'substances "
                                               "affect Vmem' off, or try again.")

            sim.extra_rho_env = p.F*Q*np.ones(sim.edl)

        elif tag == 'mit':

            if Q < 0 and np.abs(Q) <= sim.cc_mit[sim.iP].mean():  # if net charge is anionic
                sim.cc_mit[sim.iP] = sim.cc_mit[sim.iP] - np.abs(Q)

            elif Q > 0 and np.abs(Q) <= sim.cc_mit[sim.iK].mean():
                sim.cc_mit[sim.iK] = sim.cc_mit[sim.iK] - np.abs(Q)

            elif Q < 0 and np.abs(Q) > sim.cc_mit[sim.iP].mean():  # if net charge is anionic
                raise BetseSimConfigException("You've defined way more anionic charge in"
                                               "the extra substances than we can "
                                               "compensate for. Either turn 'substances "
                                               "affect Vmem' off, or try again.")
            elif Q > 0 and np.abs(Q) > sim.cc_mit[sim.iK].mean():
                raise BetseSimConfigException("You've defined way more cationic charge in"
                                               "the extra substances than we can "
                                               "compensate for. Either turn 'substances "
                                               "affect Vmem' off, or try again.")

            self.mit.extra_rho = p.F*Q*np.ones(sim.cdl)

    def get_rho_mem(self, cells, p):

        # FIXME add this to if substances affect charge loop!!

        self.rho_at_mem = np.zeros(len(cells.mem_mids_flat))

        for ind, mol in self.molecules.items():

            self.rho_at_mem += mol.z*p.F*mol.cc_at_mem*cells.diviterm[cells.mem_to_cells]

        self.rho_cells = np.dot(cells.M_sum_mems, self.rho_at_mem)/cells.num_mems

    def energy_charge(self, sim):

        if 'AMP' in self.molecules:

            numo = (self.cell_concs['ATP'] + 0.5 * self.cell_concs['ADP'])
            denomo = (self.cell_concs['ATP'] + self.cell_concs['ADP'] + self.cell_concs['AMP'])

            self.chi = numo / denomo

        else:

            self.chi = np.zeros(sim.cdl)

    def mod_after_cut_event(self, target_inds_cell, target_inds_mem, sim, cells, p, met_tag=False):

        self.extra_J_mem = np.zeros(sim.mdl)

        # get the name of the specific substance:
        for name in self.molecules:

            obj = self.molecules[name]

            obj.remove_cells(target_inds_cell, target_inds_mem, sim, cells, p)

        if self.mit_enabled:
            self.mit.remove_mits(sim, target_inds_cell)

        if sim.met_concs is not None and met_tag is True:  # update metabolism object if it's being simulated
            sim.met_concs = {'cATP': self.cell_concs['ATP'][cells.mem_to_cells],
                'cADP': self.cell_concs['ADP'][cells.mem_to_cells],
                'cPi': self.cell_concs['Pi'][cells.mem_to_cells]}

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

        self.chi_time = []

    def write_data(self, sim, cells, p):
        """
        Writes concentration data from a time-step to time-storage vectors.

        """

        # get the name of the specific substance:
        for name in self.molecules:

            obj = self.molecules[name]

            obj.c_mems_time.append(obj.cc_at_mem)
            obj.c_cells_time.append(obj.c_cells)

            # smooth env data if necessary:
            if p.smooth_concs is False and p.sim_ECM is True:
                cc_env = gaussian_filter(obj.c_env.reshape(cells.X.shape), 1.0).ravel()

            else:
                cc_env = np.copy(obj.c_env)

            obj.c_env_time.append(cc_env)

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

        self.chi_time.append(self.chi)

    def report(self, sim, p):
        """
        At the end of the simulation, tell user about mean, final concentrations of each molecule.

        """



        logs.log_info('time: '+ str(np.round(sim.time[-1], 2)) +
                      ' s of ' + str(np.round(p.total_time, 2)) + ' s')

        for name in self.molecules:

            obj = self.molecules[name]

            if self.mit_enabled:
                logs.log_info('Average ' + str(name) + ' in the mitochondria: ' +
                              str(np.round(obj.c_mit.mean(), 4)) + ' mmol/L')

            logs.log_info('Average ' + str(name) + ' in the cell: ' +
                          str(np.round(obj.c_cells.mean(), 4)) + ' mmol/L')

            #
            # logs.log_info('Average ' + str(name) + ' at mems: ' +
            #               str(np.round(obj.c_mems.mean(), 4)) + ' mmol/L')

            # logs.log_info('Average concentration of ' + str(name) + ' in the environment: ' +
            #                               str(np.round(obj.c_env.mean(), 4)) + ' mmol/L')

        if self.mit_enabled:
            logs.log_info('Average Vmit: ' + str(np.round(1.0e3 * self.mit.Vmit.mean(), 4)) + ' mV')

            # if 'ETC' in self.transporters:
            #
            #     rate = 0.5 * 3600 * 1e15 * self.mit.mit_vol.mean() * self.transporters['ETC'].rate.mean()
            #
            #     if 'ETC_ROS' in self.transporters:
            #         rate = rate + 3600 * 1e15 * self.mit.mit_vol.mean() * self.transporters['ETC_ROS'].rate.mean()
            #
            #     logs.log_info('Average O2 consumption rate: ' + str(rate) + ' fmol/cell/hr')


        if self.chi.mean() != 0.0:
            logs.log_info('Energy charge: ' + str(np.round(self.chi.mean(), 3)))

        if 'H+' in self.molecules:

            obj = self.molecules['H+']

            ph = -np.log10(1.0e-3*obj.c_cells.mean())

            logs.log_info('Average pH in the cell: ' +
                          str(np.round(ph, 2)))





        logs.log_info("-------------------------------------------------------------------")

    @type_check
    def export_all_data(self, sim, cells, p, message: str):
        """
        Exports concentration data from each molecule to a file for a single cell
        (plot cell defined in params) as a function of time.
        """

        logs.log_info('Exporting raw data for %s...', message)

        # get the name of the specific substance:
        for name in self.molecules:
            obj = self.molecules[name]
            obj.export_data(sim, cells, p, self.resultsPath)

        # FIXME we should also save vmit to a file? and pH and vm?
        if self.mit_enabled:
            pass

        # write all network model LaTeX equations to a text file:
        self.export_equations(p)
        self.export_eval_strings(p)

    @type_check
    def plot(self, sim, cells, p, message: str):
        """
        Creates plots for each molecule included in the simulation.
        """

        logs.log_info('Plotting 1D and 2D data for %s...', message)

        # network plot
        if p.plot_network is True:
            # whole_graph = self.plot_network(p)
            whole_graph = plot_master_network(self, p)
            savename = self.imagePath + 'NetworkGraph_Cell_' + str(p.plot_cell) + '.svg'
            whole_graph.write_svg(savename)

        # get the name of the specific substance:
        for name in self.molecules:
            obj = self.molecules[name]

            if p.plot_networks_single_cell:
                # create line graphs for the substance
                obj.plot_1D(sim, p, self.imagePath)

            # create 2D maps for the substance
            obj.plot_cells(sim, cells, p, self.imagePath)

            # if there's a real environment, plot 2D concentration in the environment
            if p.sim_ECM:
                obj.plot_env(sim, cells, p, self.imagePath)

        # ---------------cell everything plot---------------------------------------------
        # data_all1D = []
        plt.figure()
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

        ax_all1D.legend(loc='upper right', shadow=False, frameon=False)

        ax_all1D.set_xlabel('Time [s]')
        ax_all1D.set_ylabel('Concentration [mmol/L]')
        ax_all1D.set_title('Concentration of all substances in cell ' + str(p.plot_cell))

        if p.autosave is True:
            savename = self.imagePath + 'AllCellConcentrations_' + str(p.plot_cell) + '.png'
            plt.savefig(savename, format='png', transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # -------------environment everything plot-------------------------------------------------
        # data_all1D = []
        plt.figure()
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

        ax_all1D.legend(loc='upper right', shadow=False, frameon=False)

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

            plt.figure()
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
            # data_all1D = []
            plt.figure()
            ax_all1D = plt.subplot(111)

            # get a random selection of our chosen colors in the length of our data set:

            for i, name in enumerate(self.molecules):
                obj = self.molecules[name]

                c_mit = [arr[p.plot_cell] for arr in obj.c_mit_time]

                ax_all1D.plot(sim.time, c_mit, color=c_names.to_rgba(i), linewidth=2.0, label=name)

            ax_all1D.legend(loc='upper right', shadow=False, frameon=False)

            ax_all1D.set_xlabel('Time [s]')
            ax_all1D.set_ylabel('Concentration [mmol/L]')
            ax_all1D.set_title('Substances in mitochondria of cell ' + str(p.plot_cell))

            if p.autosave is True:
                savename = self.imagePath + 'AllMitConcentrations_' + str(p.plot_cell) + '.png'
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

            # data_all1D = []
            plt.figure()
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

            ax_all1D.legend(loc='upper right', shadow=False, frameon=False)

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

                if len(obj.flux_time) > 0:

                    # check the data structure size for this transporter:
                    if len(obj.flux_time[0]) == sim.cdl:

                        t_rate = [arr[p.plot_cell] for arr in obj.flux_time]

                    elif len(obj.flux_time[0]) == sim.mdl:
                        mem_i = cells.cell_to_mems[p.plot_cell][0]
                        t_rate = [arr[mem_i] for arr in obj.flux_time]

                    else:
                        t_rate = np.zeros(len(sim.time))

                    # data_all1D = []
                    plt.figure()
                    ax_all1D = plt.subplot(111)

                    ax_all1D.plot(sim.time, t_rate, color=c_names.to_rgba(i), linewidth=2.0, label=name)

                    transp_dataM.append(t_rate)
                    transp_header = transp_header + name + ' [mM/s]' + ','

            ax_all1D.legend(loc='upper right', shadow=False, frameon=False)

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


        #----Channels flux plot-------------------------------------------------------------------
        if len(self.channels):

            for kch, vch in self.channels.items():

                vch.plot_2D(sim, cells, p, self.imagePath)

        # energy charge plots:----------------------------------------------------------
        # 1 D plot of mitochondrial voltage--------------------------------------------------------
        chio = [arr[p.plot_cell] for arr in self.chi_time]

        plt.figure()
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

    @type_check
    def anim(self, phase: SimPhase, message: str) -> None:
        """
        Animates 2D data for each molecule in the simulation.
        """

        logs.log_info('Animating data for %s...', message)

        # get the name of the specific substance:
        for name in self.molecules:
            obj = self.molecules[name]

            if phase.p.anim.is_after_sim and obj.make_ani:
                # create 2D animations for the substance in cells
                obj.anim_cells(phase)

                # create 2D animations for the substance in the environment
                if phase.p.sim_ECM:
                    obj.anim_env(phase)

    def default_zones(self, zone_tags_a, zone_tags_i, a_list, i_list):

        # if zone tag lists aren't supplied or are invalid entries, create defaults:
        if zone_tags_a is None or zone_tags_a == 'None' or len(zone_tags_a) == 0:
            if a_list is not None and a_list != 'None' and len(a_list)>0:
                zone_tags_a = ['cell' for x in a_list]

        if zone_tags_i is None or zone_tags_i == 'None' or len(zone_tags_i) == 0:

            if i_list is not None and i_list != 'None' and len(i_list) > 0:
                zone_tags_i = ['cell' for x in i_list]

        return zone_tags_a, zone_tags_i

    def export_equations(self, p):

        # Absolute path of the text file to write this solution to:
        saveData = paths.join(self.resultsPath, 'NetworkModelLatexEquations.csv')

        with open(saveData, 'w', newline='') as csvfile:
            eqwriter = csv.writer(csvfile, delimiter='\t',
                quotechar='|', quoting=csv.QUOTE_NONE)

            for mol in self.molecules:

                if self.molecules[mol].simple_growth is True:
                    text_var = self.molecules[mol].gad_tex_vars
                    eq_var = self.molecules[mol].gad_tex_string

                    eqwriter.writerow([mol, eq_var, text_var])

            for rea in self.reactions:
                text_var = self.reactions[rea].reaction_tex_vars
                eq_var = self.reactions[rea].reaction_tex_string

                eqwriter.writerow([rea, eq_var, text_var])

            for rea in self.reactions_mit:
                text_var = self.reactions_mit[rea].reaction_tex_vars
                eq_var = self.reactions_mit[rea].reaction_tex_string

                eqwriter.writerow([rea, eq_var, text_var])

            for trans in self.transporters:

                text_var = self.transporters[trans].transporter_tex_vars
                eq_var = self.transporters[trans].transporter_tex_string

                eqwriter.writerow([trans, eq_var, text_var])

    def export_eval_strings(self, p):

        # results_path = os.path.join(p.init_results, 'model_eval_strings')
        #
        # self.resultsPath = os.path.expanduser(results_path)
        # os.makedirs(self.resultsPath, exist_ok=True)

        # Absolute path of the text file to write this solution to:
        saveData = paths.join(self.resultsPath, 'NetworkModelEvalStrings.csv')

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

        # # add Vmem to the graph
        # nde = pydot.Node('Vmem', style='filled', shape=self.vmem_shape)
        # graphicus_maximus.add_node(nde)

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

        # if there are any reactions in the mitochondria, add them to the graph
        if len(self.reactions_mit) > 0:

            for i, name in enumerate(self.reactions_mit):
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

        # if there are any mitochondria zone reactions, plot their edges on the graph (and react/prod nodes):
        if len(self.reactions_mit) > 0:

            for i, name in enumerate(self.reactions_mit):

                rea = self.reactions_mit[name]

                for i, react_name in enumerate(rea.reactants_list):

                    react_name += '_mit'

                    nde = pydot.Node(react_name, shape = self.conc_shape)
                    graphicus_maximus.add_node(nde)

                    rea_coeff = rea.reactants_coeff[i]
                    graphicus_maximus.add_edge(pydot.Edge(react_name, name, arrowhead='normal', coeff = rea_coeff,
                                                          ))

                for j, prod_name in enumerate(rea.products_list):

                    prod_name += '_mit'

                    nde = pydot.Node(prod_name, shape = self.conc_shape)
                    graphicus_maximus.add_node(nde)

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
        n_i_list, tex_list = None, reaction_zone='cell', zone_tags_a = None, zone_tags_i = None, in_mem_tag = False):

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

        * zone tags only apply to reactions occurring in 'cell' or 'mem' regions. For 'mit', the activator and inhibitor
        concentrations are always assumed to be inside the mitochondria, so are always also 'mit.

        Returns
        --------
        activator_alpha, inhibitor_alpha     string expressions representing the activator and inhibitor coeffs
        alpha_tex                            string of combined activator_alpha, inhibitor_alpha in LaTeX format

        """

        # in case it's not supplied, initialize a list to hold LaTeX string expressions for parameters
        if tex_list is None:
            tex_list = []

        # if zone tag lists aren't supplied or are invalid entries, create defaults:
        zone_tags_a, zone_tags_i = self.default_zones(zone_tags_a, zone_tags_i, a_list, i_list)

        # initialize string expressions for product of all activators and inhibitors
        activator_alpha = "("
        inhibitor_alpha = "("

        activator_alpha_tex = ""
        inhibitor_alpha_tex = ""

        # initialize lists for special "direct terms" that are the concentrations unmodified:

        direct_terms_list = []
        direct_terms_tex = []

        direct_term_alpha = "("
        direct_term_alpha_tex = ""

        # set the appropriate value for the absence of an activator or inhibitor expression:
        if reaction_zone == 'cell':
            dl = "np.ones(sim.cdl))"

        elif reaction_zone == 'mem':
            dl = "np.ones(sim.mdl))"

        elif reaction_zone == 'mit':
            dl = 'np.ones(sim.cdl))'

        elif reaction_zone == 'env':
            dl = "np.ones(sim.edl))"

        else:
            raise BetseSimConfigException("Reaction zone requested not an existing option!")


        if a_list is not None and a_list != 'None' and len(a_list) > 0:

            # initialize an empty list that will contain independently acting activator terms:
            independent_terms = []
            independent_terms_tex = []

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

                    Kmc = '({}*np.exp(-(sim.vm*{}*p.F)/(p.R*p.T)))'.format(Kmo, za)
                    Kme = '({}*np.exp((sim.vm*{}*p.F)/(p.R*p.T)))'.format(Kmo, za)

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

                        tex_name = name

                    elif zone_tag == 'env':

                        tex_name = name + '_{env}'

                        numo_string_a += "((self.env_concs['{}'][cells.map_cell2ecm]/{})**{})".format(name, Kme, n)
                        denomo_string_a += "(1 + (self.env_concs['{}'][cells.map_cell2ecm]/{})**{})".format(name, Kme, n)

                        direct_string_a += "self.env_concs['{}'][cells.map_cell2ecm]".format(name)

                elif reaction_zone == 'mem':

                    if zone_tag == 'cell':

                        numo_string_a += "((self.mem_concs['{}']/{})**{})".format(name, Kmc, n)
                        denomo_string_a += "(1 + (self.mem_concs['{}']/{})**{})".format(name, Kmc, n)

                        direct_string_a += "self.mem_concs['{}']".format(name)

                        tex_name = name

                    elif zone_tag == 'env':

                        tex_name = name + '_{env}'

                        numo_string_a += "((self.env_concs['{}'][cells.map_mem2ecm]/{})**{})".format(name, Kme, n)
                        denomo_string_a += "(1 + (self.env_concs['{}'][cells.map_mem2ecm]/{})**{})".format(name, Kme, n)

                        direct_string_a += "self.env_concs['{}'][cells.map_mem2ecm]".format(name)

                elif reaction_zone == 'mit':

                    tex_name = name + '_{mit}'

                    numo_string_a += "((self.mit_concs['{}']/{})**{})".format(name, Kmo, n)
                    denomo_string_a += "(1 + (self.mit_concs['{}']/{})**{})".format(name, Kmo, n)

                    direct_string_a += "self.mit_concs['{}']".format(name)

                elif reaction_zone == 'env':

                    if zone_tag == 'cell':

                        tex_name = name

                        numo_string_a += "1"

                        denomo_string_a += "1"

                        direct_string_a += "1"

                    elif zone_tag == 'env':

                        tex_name = name + '_{env}'

                        numo_string_a += "((self.env_concs['{}']/{})**{})".format(name, Kme, n)
                        denomo_string_a += "(1 + (self.env_concs['{}']/{})**{})".format(name, Kme, n)

                        direct_string_a += "self.env_concs['{}']".format(name)

                else:
                    raise BetseSimConfigException("You've asked for a reaction zone"
                                                   "that doesn't exist. Enable mitochondria or ensure all "
                                                   "reaction and transporter zones are 'cell' or 'env'.")

                numo_tex_a = r"\left(\frac{[%s]}{K_{%s}^{a}}\right)^{n_{%s}^{a}}" % (tex_name, tex_name, tex_name)
                denomo_tex_a = r"1+\left(\frac{[%s]}{K_{%s}^{a}}\right)^{n_{%s}^{a}}" % (tex_name, tex_name, tex_name)

                direct_tex_a =  r"\[%s]" % (tex_name)

                term = "(" + numo_string_a + "/" + denomo_string_a + ")"

                # write term as LaTeX expression:
                tex_term = r"\left(\frac{%s}{%s}\right)" % (numo_tex_a, denomo_tex_a)

                # write fixed parameter values to LaTeX----------
                kval = tex_val(Kmo)
                Ka_tex = "Ko_{%s}^{a} & =" % (tex_name)
                Ka_tex += kval

                nval = tex_val(n)
                n_a_tex = "n_{%s}^{a} & =" % (tex_name)
                n_a_tex += nval

                tex_list.append(Ka_tex)
                tex_list.append(n_a_tex)

                #---Update the final activator alpha term.
                if independence_tag is False and direct_tag is False:
                    # include the new activator term as a multiplier:
                    activator_alpha += term
                    activator_alpha_tex += tex_term

                elif independence_tag is True:
                    # include the term in the list; we will add them to the final coefficient at the end:
                    independent_terms.append(term)
                    independent_terms_tex.append(tex_term)

                elif direct_tag is True:
                    # include the direct term (no modifications) to the string as an activator
                    direct_terms_list.append(direct_string_a)
                    direct_terms_tex.append(direct_tex_a)


                if i < len(a_list) - 1 and direct_tag is False and independence_tag is False:

                    activator_alpha += "*"
                    activator_alpha_tex += r"\,"


                elif i >= len(a_list) -1: # if we've reached the end of the activators list:

                    # add any independently acting terms from the list:
                    for indt, indtex in zip(independent_terms, independent_terms_tex):

                        if activator_alpha == "(": # if nothing has been added to the activator alpha string
                            activator_alpha = '(1'
                            activator_alpha_tex = '1'

                        elif activator_alpha.endswith('*'): # otherwise, if it's a running list, remove the * operator
                            activator_alpha = activator_alpha[:-1]

                        activator_alpha += "+" + indt
                        activator_alpha_tex += "+" + indtex

                    if activator_alpha.endswith("*"):
                        activator_alpha = activator_alpha[:-1]

                    # cap things off with a final parens:
                    activator_alpha += ")"
                    activator_alpha_tex += r"\,"
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

                    Kmc = '({}*np.exp(-(sim.vm*{}*p.F)/(p.R*p.T)))'.format(Kmo, zi)
                    Kme = '({}*np.exp((sim.vm*{}*p.F)/(p.R*p.T)))'.format(Kmo, zi)

                else:
                    Kmc = '{}'.format(Kmo)
                    Kme = '{}'.format(Kmo)



                # initialize string expressions for activator and inhibitor, respectively
                numo_string_i = ""
                denomo_string_i = ""
                direct_string_i = ""

                if reaction_zone == 'cell':

                    if zone_tag == 'cell':

                        tex_name = name

                        numo_string_i += "1"
                        denomo_string_i += "(1 + (self.cell_concs['{}']/{})**{})".format(name, Kmc, n)

                        direct_string_i += "-self.cell_concs['{}']".format(name)

                    elif zone_tag == 'env':

                        tex_name = name + '_{env}'

                        numo_string_i += "1"
                        denomo_string_i += "(1 + (self.env_concs['{}'][cells.map_cell2ecm]/{})**{})".format(name, Kme, n)
                        direct_string_i += "-self.env_concs['{}'][cells.map_cell2ecm]".format(name)

                elif reaction_zone == 'mem':

                    if zone_tag == 'cell':

                        tex_name = name

                        numo_string_i += "1"
                        denomo_string_i += "(1 + (self.mem_concs['{}']/{})**{})".format(name, Kmc, n)

                        direct_string_i += "-self.mem_concs['{}']".format(name)

                    elif zone_tag == 'env':

                        tex_name = name + '_{env}'

                        numo_string_i += "1"
                        denomo_string_i += "(1 + (self.env_concs['{}'][cells.map_mem2ecm]/{})**{})".format(name, Kme, n)

                        direct_string_i += "-self.env_concs['{}'][cells.map_mem2ecm]".format(name)

                elif reaction_zone == 'mit':

                    tex_name = name + '_{mit}'

                    numo_string_i += "1"
                    denomo_string_i += "(1 + (self.mit_concs['{}']/{})**{})".format(name, Kmo, n)

                    direct_string_i += "-self.mit_concs['{}']".format(name)

                elif reaction_zone == 'env':

                    if zone_tag == 'cell':

                        tex_name = name

                        numo_string_i += "1"

                        denomo_string_i += "1"

                        direct_string_i += "1"

                    elif zone_tag == 'env':

                        tex_name = name + '_{env}'

                        numo_string_i += "1"
                        denomo_string_i += "(1 + (self.env_concs['{}']/{})**{})".format(name, Kme, n)
                        direct_string_i += "-self.env_concs['{}']".format(name)

                numo_tex_i = "1"
                denomo_tex_i = r"1+\left(\frac{[%s]}{K_{%s}^{i}}\right)^{n_{%s}^{i}}" % (tex_name, tex_name, tex_name)

                direct_tex_i =  r"[%s]" % (tex_name)

                term = "(" + numo_string_i + "/" + denomo_string_i + ")"

                tex_term = r"\left(\frac{%s}{%s}\right)" % (numo_tex_i, denomo_tex_i)

                # write fixed parameter values to LaTeX----------
                # FIXME differentiate between Kme and Kmc here
                kval = tex_val(Kmo)
                Ki_tex = "Ko_{%s}^{i} & =" % (tex_name)
                Ki_tex += kval

                nval = tex_val(n)
                n_i_tex = "n_{%s}^{i} & =" % (tex_name)
                n_i_tex += nval

                tex_list.append(Ki_tex)
                tex_list.append(n_i_tex)

                #-------------------------------------------------

                if direct_tag is False:
                    inhibitor_alpha += term
                    inhibitor_alpha_tex += tex_term

                elif direct_tag is True:

                    direct_terms_list.append(direct_string_i)
                    direct_terms_tex.append(direct_tex_i)

                if i < len(i_list) - 1:

                    inhibitor_alpha += "*"
                    inhibitor_alpha_tex += r"\,"

                elif i >= len(i_list) -1:   # else if we've reached the end:

                    if inhibitor_alpha.endswith("*"):

                        inhibitor_alpha = inhibitor_alpha[:-1]

                    inhibitor_alpha += ")"
                    inhibitor_alpha_tex += r"\,"

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

        # alpha = activator_alpha + '*' + inhibitor_alpha
        # finalize the LaTeX math string:
        alpha_tex = activator_alpha_tex + r"\," + inhibitor_alpha_tex


        if len(direct_terms_list) > 0: # if there's anything in the direct terms list

            for q, (dit, ditex) in enumerate(zip(direct_terms_list, direct_terms_tex)):

                direct_term_alpha += dit
                direct_term_alpha_tex += ditex

                if q < len(direct_terms_list) - 1:
                    direct_term_alpha += '+'
                    direct_term_alpha_tex += '+'

                else:
                    direct_term_alpha += ")"  # finish it off with a parens or space:
                    direct_term_alpha_tex += r"\,"

            alpha = alpha + "*" + direct_term_alpha
            alpha_tex = alpha_tex + r"\," + direct_term_alpha_tex


        return alpha, alpha_tex, tex_list

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

        self.growth_targets_cell = cells.cell_i

        self.gating_mod = 1.0

    def transport(self, sim, cells, p):
        """
        Transports the molecule across the membrane,
        through gap junctions, and if p.sim_ECM is true,
        through extracellular spaces and the environment.

        """


        self.c_env, self.c_cells, self.cc_at_mem, self.f_mem, self.f_gj, fenvx, fenvy = stb.molecule_mover(sim,
                                                                self.c_env,
                                                                self.c_cells,
                                                                cells, p,
                                                                z=self.z,
                                                                Dm = self.Dm,
                                                                Do = self.Do,
                                                                Dgj = self.Dgj,
                                                                c_bound = self.c_bound,
                                                                ignoreECM = self.ignore_ECM_pump,
                                                                smoothECM = p.smooth_concs,
                                                                ignoreTJ = self.ignoreTJ,
                                                                ignoreGJ = self.ignoreGJ,
                                                                rho = sim.rho_channel,
                                                                cmems=self.cc_at_mem)

    def updateC(self, flux, sim, cells, p):
        """

        General updater for a flux defined on membranes and updating concentrations in
        cells and environment.

        """
        self.c_cells, self.cc_at_mem, self.c_env = stb.update_Co(sim, self.c_cells, self.cc_at_mem,
                                                                self.c_env, flux, cells, p, ignoreECM=True)

    def update_intra(self, sim, cells, p):

        cav = self.c_cells[cells.mem_to_cells]  # concentration at cell centre
        cmi = self.cc_at_mem  # concentration at membrane
        z = self.z  # charge of ion
        Do = self.Do  # diffusion constant of ion

        cp = (cav + cmi) / 2  # concentration at midpoint between cell centre and membrane
        cg = (cmi - cav) / cells.R_rads  # concentration gradients

        # calculate normal component of microtubules at membrane:

        # umx, umy = sim.mtubes.mtubes_to_cell(cells, p)
        # umtn = umx[cells.mem_to_cells]*cells.mem_vects_flat[:, 2] + umy[cells.mem_to_cells]*cells.mem_vects_flat[:, 3]

        umtn = sim.mtubes.mtubes_x*cells.mem_vects_flat[:, 2] + sim.mtubes.mtubes_y*cells.mem_vects_flat[:, 3]

        cfluxo = (-Do*cg + ((Do*p.q*cp*z)/(p.kb*sim.T))*sim.Ec + umtn*self.u_mt*cp)

        # as no net mass must leave this intracellular movement, make the flux divergence-free:
        cflux = stb.single_cell_div_free(cfluxo, cells)

        # calculate the actual concentration at membranes by unpacking to concentration vectors:
        self.cc_at_mem = cmi + cflux*(cells.mem_sa/cells.mem_vol)*p.dt

        # smooth the concentration:
        self.cc_at_mem = sim.smooth_weight_mem*self.cc_at_mem + sim.smooth_weight_o*cav

        # # calculate the divergence so we can update the centre concentration:
        # divJ = np.dot(cells.M_sum_mems, -cflux*cells.mem_sa)/cells.cell_vol
        # self.c_cells = self.c_cells + divJ*p.dt

        # deal with the fact that our coarse diffusion model may leave some sub-zero concentrations:
        indsZ = (self.cc_at_mem < 0.0).nonzero()

        if len(indsZ[0]):

            raise BetseSimInstabilityException("Network concentration " + self.name + " on membrane below zero! Your simulation has"
                                               " become unstable.")

        # update the main matrices:
        self.flux_intra = cflux*1

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

                self.c_cells, self.c_env, flux = stb.molecule_pump(sim, self.c_cells, self.c_env,
                                                                     cells, p, Df=self.Do, z=self.z,
                                                                     pump_into_cell=self.pump_to_cell,
                                                                     alpha_max=self.pump_max_val, Km_X=self.pump_Km,
                                                                     Km_ATP=1.0, met = met_vect,
                                                                     ignoreECM = self.ignore_ECM_pump,
                                                                     rho = sim.rho_pump)
                if p.metabolism_enabled:
                    # update ATP concentrations after pump action:
                    sim.metabo.update_ATP(flux, sim, cells, p)

            else:

                self.c_cells, self.c_env, flux = stb.molecule_transporter(sim, self.c_cells, self.c_env,
                    cells, p, Df=self.Do, z=self.z, pump_into_cell=self.pump_to_cell, alpha_max=self.pump_max_val,
                    Km_X=self.pump_Km, Keq= 1.0, ignoreECM = self.ignore_ECM_pump, rho = sim.rho_pump)

    def gating(self, sim, cells, p):
        """
        Uses the molecule concentration to open an ion channel in the cell membranes.

        """

        # update membrane permeability if dye targets an ion channel:
        if self.use_gating_ligand:

            for ion_tag in self.gating_ion:

                # calculate any activators and/or inhibitor effects:
                if self.gating_extracell is False:

                    Dm_mod_mol = sim.rho_channel*self.gating_max_val*tb.hill(self.cc_at_mem,
                                                                            self.gating_Hill_K,self.gating_Hill_n)

                else:

                    if p.sim_ECM is False:

                        Dm_mod_mol = self.gating_max_val * tb.hill(self.c_env, self.gating_Hill_K, self.gating_Hill_n)

                    else:

                        # print(self.c_env.mean())

                        Dm_mod_mol = self.gating_max_val * tb.hill(self.c_env[cells.map_mem2ecm],
                                                                   self.gating_Hill_K, self.gating_Hill_n)

                # obtain concentration of ion inside and out of the cell, as well as its charge z:
                c_mem = sim.cc_cells[ion_tag][cells.mem_to_cells]

                if p.sim_ECM is True:
                    c_env = sim.cc_env[ion_tag][cells.map_mem2ecm]

                else:
                    c_env = sim.cc_env[ion_tag]

                IdM = np.ones(sim.mdl)

                z_ion = sim.zs[ion_tag]*IdM

                # membrane diffusion constant of the channel:
                Dchan = sim.rho_channel*Dm_mod_mol*self.gating_mod

                # calculate specific ion flux contribution for this channel:
                chan_flx = stb.electroflux(c_env, c_mem, Dchan, p.tm*IdM, z_ion, sim.vm, sim.T, p, rho=sim.rho_channel)

                # chan_flx = stb.electroflux(c_env, c_mem,
                #                        sim.Dm_cells[ion_tag] + Dchan, IdM*p.tm, z_ion, sim.vm, sim.T, p,
                #                        rho=sim.rho_channel)

                # update the sim flux keeper to ensure this contributes to net current:
                sim.fluxes_mem[ion_tag] = sim.fluxes_mem[ion_tag] + chan_flx

                # update ion concentrations in cell and ecm:
                sim.cc_cells[ion_tag], sim.cc_at_mem[ion_tag], sim.cc_env[ion_tag] = stb.update_Co(sim,
                                                                     sim.cc_cells[ion_tag], sim.cc_at_mem[ion_tag],
                                                                      sim.cc_env[ion_tag], chan_flx, cells, p,
                                                                      ignoreECM=False)

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

        # remove items from the mem concentration lists:
        ccmem2 = np.delete(self.cc_at_mem, target_inds_mem)
        # reassign the new data vector to the object:
        self.cc_at_mem = ccmem2[:]

        # remove cells from the cell gradient concentration list:
        cgx2 = np.delete(self.flux_intra, target_inds_mem)
        # reassign the new data vector to the object:
        self.cc_grad_x = cgx2[:]

        # # remove cells from the mems concentration list:
        # cmems2 = np.delete(self.c_mems, target_inds_mem)
        # # reassign the new data vector to the object:
        # self.c_mems = cmems2[:]

        if self.simple_growth is True:

            gmfc = np.delete(self.growth_mod_function_cells, target_inds_cell)
            self.growth_mod_function_cells = gmfc[:]

            # if len(self.growth_mod_function_cells) == 0:
            #     self.growth_mod_function_cells = 1

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
                self.c_env[:] = self.change_bounds_target*effector_MorphEnv + self.c_envo*(1-effector_MorphEnv)

            elif p.sim_ECM is True:
                self.c_bound = self.change_bounds_target*effector_MorphEnv + self.c_envo*(1-effector_MorphEnv)
                # self.c_bound = self.conc_MorphEnv*effector_MorphEnv + self.c_envo*(1-effector_MorphEnv)

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
        Create 1D plot of concentration in cell and environment for a single
        cell (params plot cell) as a function of simulation time.
        """

        c_cells = [arr[p.plot_cell] for arr in self.c_cells_time]
        plt.figure()
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
            plt.figure()
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

        # fig, ax, cb = viz.plotPrettyPolyData(self.c_mems,
        #     sim, cells, p,
        #     number_cells=p.enumerate_cells,
        #     clrAutoscale=self.plot_autoscale,
        #     clrMin=self.plot_min,
        #     clrMax=self.plot_max,
        #     clrmap=p.default_cm)

        # fig, ax, cb = viz.plotPolyData(sim, cells, p,
        #                                zdata=self.c_cells, number_cells=p.enumerate_cells, clrmap=p.default_cm,
        #                                clrMin=self.plot_min, clrMax=self.plot_max, clrAutoscale=self.plot_autoscale)

        fig, ax, cb = viz.plotPrettyPolyData(self.cc_at_mem,
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

        if p.smooth_level == 0.0:
            dyeEnv = gaussian_filter(self.c_env.reshape(cells.X.shape), 1.0)
        else:
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

    #FIXME: Ideally, this method should be refactored to comply with the
    #new pipeline API.
    @type_check
    def anim_cells(self, phase: SimPhase) -> None:
        """
        Create 2D animation of cell concentration.
        """

        #FIXME: To support GUI modification, refactor this class to access the
        #underlying YAML-based subconfiguration.
        conf = SimConfVisualCellsNonYAML(
            is_color_autoscaled=self.plot_autoscale,
            color_min=self.plot_min,
            color_max=self.plot_max)

        AnimFlatCellsTimeSeries(
            phase=phase,
            conf=conf,
            time_series=self.c_cells_time,
            label=self.name + '_cells',
            figure_title='Cytosolic ' + self.name,
            colorbar_title='Concentration [mmol/L]')

    #FIXME: Ideally, this method should be refactored to comply with the
    #new pipeline API.
    @type_check
    def anim_env(self, phase: SimPhase) -> None:
        """
        Create 2D animation of env concentration.
        """

        #FIXME: To support GUI modification, refactor this class to access the
        #underlying YAML-based subconfiguration.
        conf = SimConfVisualCellsNonYAML(
            is_color_autoscaled=self.plot_autoscale,
            color_min=self.plot_min,
            color_max=self.plot_max)

        env_time_series = [
            env.reshape(phase.cells.X.shape) for env in self.c_env_time]
        AnimEnvTimeSeries(
            phase=phase,
            conf=conf,
            time_series=env_time_series,
            label=self.name + '_env',
            figure_title='Environmental ' + self.name,
            colorbar_title='Concentration [mmol/L]')

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
            plt.figure()
            ax = plt.subplot(111)
            ax.plot(sim.time, r_rate)
            ax.set_xlabel('Time [s]')
            ax.set_ylabel('Rate [mM/s]')
            ax.set_title('Rate of ' + self.name + ' in cell ' + str(p.plot_cell))

        elif self.reaction_zone == 'mit' and self.mit_enabled is True:

            r_rate = [arr[p.plot_cell] for arr in self.rate_time]
            plt.figure()
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
            self.transporter_targets_env = cells.map_mem2ecm

    def plot_1D(self, sim, cells, p, saveImagePath):

        if len(self.flux_time) > 0:

            fig, ax, cb = viz.plotPrettyPolyData(self.flux_time[-1],
                sim, cells, p,
                number_cells=p.enumerate_cells,
                clrAutoscale=True,
                clrMin=0,
                clrMax=1,
                clrmap=p.default_cm)

            tit = 'Final Rate ' + self.name

            fig.suptitle(tit, fontsize=14, fontweight='bold')
            ax.set_xlabel('Spatial distance [um]')
            ax.set_ylabel('Spatial distance [um]')
            cb.set_label('Flux [mmol/s]')

            if p.autosave is True:
                savename = saveImagePath + 'TransporterRate2D_' + self.name + '.png'
                plt.savefig(savename, format='png', transparent=True)



            if len(self.flux_time[0]) == sim.cdl:
                r_rate = [arr[p.plot_cell] for arr in self.flux_time]
            else:
                mem_i = cells.cell_to_mems[p.plot_cell][0]
                r_rate = [arr[mem_i] for arr in self.flux_time]

            plt.figure()
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

    # FIXME something is wrong with channels -- biased along y axis.

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

        elif ion_string == 'K':

            self.dummy_dyna.maxDmK = max_val
            self.dummy_dyna.targets_vgK = self.channel_targets_mem
            class_string = vgk

        elif ion_string == 'Cl':

            self.dummy_dyna.maxDmCl = max_val
            self.dummy_dyna.targets_vgCl = self.channel_targets_mem
            class_string = vgcl

        elif ion_string == 'Ca':

            self.dummy_dyna.maxDmCa = max_val
            self.dummy_dyna.targets_vgCa = self.channel_targets_mem
            class_string = vgca

        elif ion_string == 'Fun':

            self.dummy_dyna.maxDmFun = max_val
            self.dummy_dyna.targets_vgFun = self.channel_targets_mem
            class_string = vgfun

        elif ion_string == 'Cat':

            self.dummy_dyna.maxDmCat = max_val
            self.dummy_dyna.targets_vgCat = self.channel_targets_mem
            class_string = vgcat

        else:

            raise BetseSimConfigException("Substance-modulated ion type not available. "
                                           "Valid choices: Na, K, Ca, NaP, Kir, and Fun")

            # create the desired voltage gated sodium channel instance:

        self.channel_core = getattr(class_string,type_string)()

        # if p.run_sim is True:
        #     # initialize the channel object
        self.channel_core.init(self.dummy_dyna, sim, cells, p)

    def update_channel(self, sim, cells, p):
        self.dummy_dyna.tissueProfiles(sim, cells, p)  # initialize all tissue profiles
        self.init_channel(self.channel_class, self.channel_type, self.channelMax, sim, cells, p)

    def plot_2D(self, sim, cells, p, saveImagePath):

        flux_chan = self.channel_core.chan_flux

        fig, ax, cb = viz.plotPrettyPolyData(flux_chan,
            sim, cells, p,
            number_cells=p.enumerate_cells,
            clrAutoscale=True,
            clrMin=0,
            clrMax=1,
            clrmap=p.default_cm)

        tit = 'Final Flux of Channel ' + self.name

        fig.suptitle(tit, fontsize=14, fontweight='bold')
        ax.set_xlabel('Spatial distance [um]')
        ax.set_ylabel('Spatial distance [um]')
        cb.set_label('Flux [mmol/s]')

        if p.autosave is True:
            savename = saveImagePath + 'ChannelFlux2D_' + self.name + '.png'
            plt.savefig(savename, format='png', transparent=True)

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

            raise BetseSimConfigException("You have requested a "
                                           "sim modulator that is not "
                                           "available. Available choices "
                                           "are: 'gj', 'Na/K-ATPase', 'H/K-ATPase', "
                                           "and 'V-ATPase' ")

def rgba2hex(rgb_col, alpha_val):
    """
    Convert an rgb tuple into a hex color
    code, with a desired alpha value added.

    The format of the returned hex color code
    works with GraphViz.

    Parameters
    -----------
    rgb_col     tuple of rgb alpha scaled 0 to 1
    alpha_val   float of alpha value, scaled 0 to 1

    Returns
    ----------
    hex_string   Hex color code with alpha value added.
                 Code string created for GraphViz color
                 specifications.

    """

    # convert the tuple into a list:
    rgb_col = list(rgb_col)

    # add the alpha value:
    rgb_col[-1] = alpha_val

    # convert the fractions into 255 color values:
    rgb_conv = [rg * 255 for rg in rgb_col]

    # write the string:
    hex_list = [format(int(rg), '02X') for rg in rgb_conv]

    # hex_string = "#" + hex_list[-1] + hex_list[0] + hex_list[1] + hex_list[2]
    hex_string = "#" + hex_list[0] + hex_list[1] + hex_list[2] + hex_list[3]

    return hex_string

def tex_val(v):
    """
    Checks a float value to see if
    it needs to be written in decimal or
    scientific notation, and returns
    a string representing the float in
    LaTeX math formalism.

    Parameters
    -------------
    v:         Float

    Returns
    --------
    v_str      String expressing float in normal or sci notation

    """

    if v != 0.0:

        v_check = int(math.log10(abs(v)))

    else:
        v_check = 0

    if v_check >= 3 or v_check <= -2:

        v_str = "%.2e" % (v)

    else:

        v_str = "%.2f" % (v)

    return v_str
