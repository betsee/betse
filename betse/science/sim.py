#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import copy
import os
import os.path
import time
import matplotlib.pyplot as plt
import numpy as np
from random import shuffle
from scipy import interpolate as interp
from scipy.ndimage.filters import gaussian_filter
from betse.util.io.log import logs
from betse.science import filehandling as fh
from betse.science import finitediff as fd
from betse.science import toolbox as tb
from betse.science.plot.anim.anim import AnimCellsWhileSolving
from betse.science import sim_toolbox as stb
from betse.science.tissue.channels_o import Gap_Junction
from betse.science.tissue.handler import TissueHandler
from betse.science.physics.ion_current import get_current
from betse.science.physics.flow import getFlow
from betse.science.physics.deform import (
    getDeformation, timeDeform, implement_deform_timestep)
from betse.science.physics.move_channels import eosmosis
from betse.science.physics.pressures import electro_F, getHydroF, osmotic_P

class Simulator(object):
    '''
    Contains the main routines used in the simulation of networked cell
    bioelectrical activity. For efficiency, all methods are implemented in
    terms of Numpy-based linear algebra.

    Methods
    -------
    baseInit_all(cells,p)           Prepares core data structures necessary for initialization and sim runs

    update_V(cells,p,t)             Gets charge densities in cells and environment
                                      and calculates respective voltages.

    update_C(ion_i,flux, cells, p)     Updates concentration of ion with index
                                            ion_i in cell and environment for a flux leaving the cell.

    acid_handler(cells,p,t)        Updates H+ concentrations in cell and
                                            environment, which are further influenced by the bicarbonate buffer action.
                                            Also, if included, runs the HKATPase pump and if included, runs the
                                            VATPase pump.

    update_gj(cells,p,t,i)                  Calculates the voltage gradient between two cells, the gating character of
                                            gap junctions, and updates concentration for ion 'i' and voltage of each
                                            cell after electrodiffusion of ion 'i' between gap junction connected cells.

    update_ecm(cells,p,t,i)                 Updates the environmental spaces by calculating electrodiffusive transport
                                            of ion 'i'.


    get_Efield(cells,p)                     Calculates electric fields in cells and environment.

    get_Bfield(cells,p)                     Calculates magnetic fields in cells and environment.

    get_current(cells,p)                    Calculates currents in cells and environment.

    getFlow(cells,p)                        Calculates electroosmotic flows in cells and environment.

    eosmosis(cells,p)                       Calculates lateral movement of membrane pumps and channels via tangential
                                            forces exerted by endogenous electric fields and electroosmotic flows.

    get_ion(label)                          Supply the ion name as a string input ('Na', 'K', 'Ca', etc) and it outputs
                                            the sim index of that ion type.

    initDenv(cells,p)                       Initializes the environmental diffusion matrix and corresponding weight
                                            matrices, including tight and adherin junctions.

    Attributes
    ----------
    _anim_cells_while_solving : AnimCellsWhileSolving
        In-place animation of cell voltage as a function of time plotted during
        (rather than after) simulation modelling if both requested and a
        simulation is currently being modelled _or_ `None` otherwise.
    vcell_time : np.ndarray
        Voltage at the inner membrane surface of each cell as a function of
        time.
    venv_time : np.ndarray
        Voltage at the outer membrane surface of each cell as a function of
        time.
    vm : np.ndarray
        Transmembrane voltage of each cell for the current time step.
    vm_time : np.ndarray
        Transmembrane voltage of each cell as a function of time.
    vm_Matrix : np.ndarray
        Transmembrane voltage of each cell as a function of time, resampled for
        use in smooth visualization (e.g., streamplots).
    '''

    def __init__(self, p: 'Parameters'):

        #FIXME: Define all other instance attributes as well.
        # For safety, defined subsequently accessed instance attributes.
        self._anim_cells_while_solving = None

        #FIXME: Defer until later. To quote the "simrunner" module, which
        #explicitly calls this public method:
        #   "Reinitialize save and load directories in case params defines new
        #    ones for this sim."
        #Hence, this method should instead be called as the first statement in
        #both the run_loop_no_ecm() and run_loop_with_ecm() methods.
        self.fileInit(p)

    def fileInit(self, p):
        '''
        Initializes the pathnames of top-level files and directories comprising
        the BETSE cache for subsequent initialization and simulation runs.

        This method currently implicitly assigns file names, but will
        (hopefully) permit caller-specified pathnames at some point.
        '''

        # Make the BETSE-specific cache directory if not found.
        betse_cache_dir = os.path.expanduser(p.init_path)
        os.makedirs(betse_cache_dir, exist_ok=True)

        sim_cache_dir = os.path.expanduser(p.sim_path)
        os.makedirs(sim_cache_dir, exist_ok=True)

        # Define data paths for saving an initialization and simulation run:
        self.savedInit = os.path.join(betse_cache_dir, p.init_filename)
        self.savedSim = os.path.join(sim_cache_dir, p.sim_filename)

    def baseInit_all(self, cells, p):
        """
        Creates a host of initialized data matrices for the main simulation,
        including intracellular and environmental concentrations, voltages, and specific
        diffusion constants.

        This method is only done once per world-creation, and therefore contains crucial
        parameters, such as the type of ions included in the simulation, which can't be
        changed after running an initialization.

        """

        i = -1  # dynamic index

        self.mdl = len(cells.mem_i)  # mems-data-length
        self.cdl = len(cells.cell_i)  # cells-data-length

        if p.sim_ECM is True:  # set environnment data length
            self.edl = len(cells.xypts)

        else:
            self.edl = len(cells.mem_i)

        # load in the object that mathematically handles individual gap junction functionality:
        self.gj_funk = Gap_Junction(p)

        # Identity matrix to easily make matrices out of scalars
        self.id_mems = np.ones(self.mdl)

        self.cc_cells = []  # cell concentrations at cell centres
        self.cc_mems = []  # cell concentrations at membranes initialized
        self.cc_er = []  # endoplasmic reticulum ion concentrations in each cell
        self.cc_env = []  # environmental concentrations initialized

        self.zs = []  # ion valence state initialized
        self.z_er = []  # ion valence states of er ions
        self.z_array = []  # ion valence array matched to cell points

        self.z_array_er = []
        self.Dm_cells = []  # membrane diffusion constants initialized
        self.Dm_er = []  # a list of endoplasmic reticulum membrane state
        self.D_free = []  # a list of single-valued free diffusion constants for each ion
        self.D_gj = []  # an array of diffusion constants for gap junctions
        self.movingIons = []  # moving ions indices
        self.ionlabel = {}  # dictionary to hold ion label names
        self.molar_mass = []

        self.T = p.T  # set the base temperature for the simulation

        self.flx_gj_i = np.zeros(len(cells.nn_i))  # flux matrix across gj for individual ions
        self.fluxes_gj = []

        self.J_gj_x = np.zeros(len(cells.nn_i))  # total current in the gj network
        self.J_gj_y = np.zeros(len(cells.nn_i))  # total current in the gj network

        self.J_cell_x = np.zeros(len(cells.cell_i))
        self.J_cell_y = np.zeros(len(cells.cell_i))

        # Membrane current data structure initialization
        self.flx_mem_i = np.zeros(len(cells.mem_i))
        self.fluxes_mem = []
        self.I_mem = np.zeros(len(cells.mem_i))  # total current across membranes

        self.P_cells = np.zeros(self.cdl)  # initialize pressure in cells

        self.v_cell = np.zeros(self.mdl)  # initialize intracellular voltage
        self.vm = np.zeros(self.mdl)     # initialize vmem

        if p.sim_ECM is True:  # special items specific to simulation of extracellular spaces only:

            # vectors storing separate cell and env voltages
            self.v_env = np.zeros(len(cells.xypts))


            self.z_array_env = []  # ion valence array matched to env points
            self.D_env = []  # an array of diffusion constants for each ion defined on env grid
            self.c_env_bound = []  # moving ion concentration at global boundary
            self.Dtj_rel = []  # relative diffusion constants for ions across tight junctions

            # # Initialize membrane thickness:
            # self.tm = np.zeros(len(cells.mem_i))
            # self.tm[:] = p.tm

            # initialize environmental fluxes and current data stuctures:
            self.flx_env_i = np.zeros(self.edl)
            self.fluxes_env_x = []
            self.fluxes_env_y = []
            self.I_env = np.zeros(len(cells.xypts))  # total current in environment

        else:  # items specific to simulation *without* extracellular spaces:
            # Initialize environmental volume:
            self.envV = np.zeros(self.mdl)
            self.envV[:] = p.vol_env

            # # Initialize membrane thickness:
            # self.tm = np.zeros(self.mdl)
            # self.tm[:] = p.tm

        if p.fluid_flow is True:
            # Electroosmosis Initialization:

            # initialize vectors for electroosmosis in the cell collection wrt each gap junction (note data type!):
            self.u_cells_x = np.zeros(self.cdl)
            self.u_cells_y = np.zeros(self.cdl)

            if p.sim_ECM is True:
                # initialize vectors for env flow (note enhanced data type!):
                self.u_env_x = np.zeros(cells.X.shape)
                self.u_env_y = np.zeros(cells.X.shape)


        if p.deformation is True:
            # initialize vectors for potential deformation:
            self.d_cells_x = np.zeros(self.cdl)
            self.d_cells_y = np.zeros(self.cdl)


        if p.gj_flux_sensitive is True:

            self.gj_rho = np.zeros(len(cells.nn_i))

        else:

            self.gj_rho = 0


        ion_names = list(p.ions_dict.keys())

        for name in ion_names:  # go through ion list/dictionary and initialize all sim structures

            if p.ions_dict[name] == 1:

                if name != 'H':

                    i = i+1 # update the dynamic index

                    str1 = 'i' + name  # create the ion index

                    setattr(self, str1, i)  # dynamically add this field to the object

                    self.ionlabel[vars(self)[str1]] = p.ion_long_name[name]

                    if name != 'P':
                        self.movingIons.append(vars(self)[str1])

                    # cell concentration for the ion
                    str_mems = 'c' + name + '_mems'
                    setattr(self, str_mems, np.zeros(self.mdl))
                    vars(self)[str_mems][:] = p.cell_concs[name]


                    # environmental concentration for the ion
                    str_env = 'c' + name + '_env'

                    setattr(self, str_env, np.zeros(self.edl))

                    vars(self)[str_env][:] = p.env_concs[name]

                    # base transmembrane diffusion for each ion
                    str_Dm = 'Dm' + name

                    setattr(self, str_Dm, np.zeros(self.mdl))

                    vars(self)[str_Dm][:] = p.mem_perms[name]

                    # base gap junction (intercellular) diffusion for each ion
                    str_Dgj = 'Dgj' + name

                    setattr(self, str_Dgj, np.zeros(len(cells.nn_i)))
                    vars(self)[str_Dgj][:] = p.free_diff[name]

                    # environmental diffusion for each ion:
                    if p.sim_ECM:
                        str_Denv = 'D' + name

                        setattr(self, str_Denv, np.zeros(len(cells.xypts)))
                        vars(self)[str_Denv][:] = p.free_diff[name]

                    # ion charge characteristic for intracellular:
                    str_z = 'z' + name

                    setattr(self, str_z, np.zeros(self.mdl))
                    vars(self)[str_z][:] = p.ion_charge[name]

                    if p.sim_ECM:  # ion charge characteristic for extracellular:
                        str_z2 = 'z2' + name

                        setattr(self, str_z2, np.zeros(len(cells.xypts)))
                        vars(self)[str_z2][:] = p.ion_charge[name]

                    self.cc_mems.append(vars(self)[str_mems])
                    self.cc_env.append(vars(self)[str_env])

                    self.zs.append(p.ion_charge[name])
                    self.molar_mass.append(p.molar_mass[name])
                    self.z_array.append(vars(self)[str_z])
                    self.Dm_cells.append(vars(self)[str_Dm])
                    self.D_gj.append(vars(self)[str_Dgj])
                    self.D_free.append(p.free_diff[name])

                    self.fluxes_gj.append(self.flx_gj_i)
                    self.fluxes_mem.append(self.flx_mem_i)

                    if p.sim_ECM:
                        self.c_env_bound.append(p.env_concs[name])
                        self.z_array_env.append(vars(self)[str_z2])
                        self.D_env.append(vars(self)[str_Denv])
                        self.Dtj_rel.append(p.Dtj_rel[name])
                        self.fluxes_env_x.append(self.flx_env_i)
                        self.fluxes_env_y.append(self.flx_env_i)

                    if name == 'Ca':
                        self.cCa_er = np.zeros(len(cells.cell_i))
                        self.cCa_er[:] = p.cCa_er

                        self.zCa_er = np.zeros(len(cells.cell_i))
                        self.zCa_er[:] = p.z_Ca

                        self.cc_er.append(self.cCa_er)
                        self.z_er.append(p.z_Ca)
                        self.z_array_er.append(self.zCa_er)

                    if name == 'M' and p.ions_dict['Ca'] == 1:
                        self.cM_er = np.zeros(len(cells.cell_i))
                        self.cM_er[:] = p.cCa_er

                        self.zM_er = np.zeros(len(cells.cell_i))
                        self.zM_er[:] = p.z_M

                        self.cc_er.append(self.cM_er)
                        self.z_er.append(p.z_M)
                        self.z_array_er.append(self.zM_er)

        # Do H+ separately as it's complicated by the bicarbonate buffer

        if p.ions_dict['H'] == 1:

            i = i + 1

            self.iH = i

            self.ionlabel[self.iH] = 'protons'

            self.movingIons.append(self.iH)

            # create concentrations of dissolved carbon dioxide (carbonic acid, non-dissociated):

            self.cHM_mems = np.zeros(self.mdl)
            self.cHM_mems[:] = 0.03 * p.CO2

            self.cHM_env = np.zeros(self.edl)

            self.cHM_env[:] = 0.03 * p.CO2

            self.cH_mems = np.zeros(self.mdl)
            self.cH_mems[:] = p.cH_cell

            # use Henderson-Hasselbach equation to obtain pH and cH concentrations:

            self.pH_cell = 6.1 + np.log10(self.cM_mems / self.cHM_mems)
            self.cH_mems = (10 ** (-self.pH_cell)) * 1e3  # units mmol/L

            self.pH_env = 6.1 + np.log10(self.cM_env / self.cHM_env)
            self.cH_env = (10 ** (-self.pH_env)) * 1e3  # units mmol/L

            # initialize diffusion constants

            DmH = np.zeros(self.mdl)

            DmH[:] = p.Dm_H

            self.zH = np.zeros(self.mdl)
            self.zH[:] = p.z_H

            # gap junction diffusion constant for H+
            DgjH = np.zeros(len(cells.nn_i))
            DgjH[:] = p.free_diff['H']

            if p.sim_ECM is True:
                self.zH2 = np.zeros(len(cells.xypts))
                self.zH2[:] = p.z_H

                # environmental diffusion for H+
                DenvH = np.zeros(len(cells.xypts))
                DenvH[:] = p.free_diff['H']

                # add fixed boundary concentration of H+
                self.c_env_bound.append(p.env_concs['H'])

            # append items to main data vectors:
            self.cc_mems.append(self.cH_mems)
            self.cc_env.append(self.cH_env)

            self.zs.append(p.z_H)
            self.molar_mass.append(p.M_H)
            self.z_array.append(self.zH)

            self.Dm_cells.append(DmH)
            self.D_gj.append(DgjH)
            self.D_free.append(p.Do_H)

            self.fluxes_gj.append(self.flx_gj)
            self.fluxes_mem.append(self.flx_mem_i)

            if p.sim_ECM is True:
                self.z_array_env.append(self.zH2)
                self.D_env.append(DenvH)
                self.Dtj_rel.append(p.Dtj_rel['H'])
                self.fluxes_env_x.append(self.flx_env_i)
                self.fluxes_env_y.append(self.flx_env_i)

        # -------------------------------------------------------------------------------------------------------

        if p.ions_dict['Ca'] == 1:  # initialize the endoplasmic reticulum
            # Define the diffusion matrix for the endoplasmic reticulum:
            self.Dm_er = np.zeros((2, len(cells.cell_i)))
            self.Dm_er[0, :] = p.Dm_Ca
            self.Dm_er[1, :] = p.Dm_M

            self.v_er = np.zeros(len(cells.cell_i))


            # initialize a time-zero vmem vector:

        self.vm_to = np.zeros(self.mdl)  # FIXME is this used anywhere?

        # convert all data structures to Numpy arrays:
        self.cc_mems = np.asarray(self.cc_mems)
        self.cc_env = np.asarray(self.cc_env)

        self.zs = np.asarray(self.zs)
        self.z_array = np.asarray(self.z_array)

        self.Dm_cells = np.asarray(self.Dm_cells)

        self.D_free = np.asarray(self.D_free)
        self.D_gj = np.asarray(self.D_gj)
        self.molar_mass = np.asarray(self.molar_mass)

        self.fluxes_gj = np.asarray(self.fluxes_gj)
        self.fluxes_mem = np.asarray(self.fluxes_mem)

        if p.ions_dict['Ca'] == 1:  # items specific for Calcium dynamics
            self.cc_er = np.asarray(self.cc_er)
            self.z_array_er = np.asarray(self.z_array_er)
            self.Dm_er = np.asarray(self.Dm_er)

        if p.sim_ECM is True:  # items specific for extracellular spaces simulation:

            self.z_array_env = np.asarray(self.z_array_env)
            self.D_env = np.asarray(self.D_env)
            self.fluxes_env_x = np.asarray(self.fluxes_env_x)
            self.fluxes_env_y = np.asarray(self.fluxes_env_y)

            # boundary conditions for voltages:

            # voltage (scheduled dynamics might vary these values)
            self.bound_V = {}
            self.bound_V['T'] = 0
            self.bound_V['B'] = 0
            self.bound_V['L'] = 0
            self.bound_V['R'] = 0

            # initialize the environmental diffusion matrix:
            self.initDenv(cells, p)

        # gap junction specific arrays:
        self.id_gj = np.ones(len(cells.mem_i))  # identity array for gap junction indices...
        self.gjopen = np.ones(len(cells.mem_i))  # holds gap junction open fraction for each gj
        self.gjl = np.zeros(len(cells.mem_i))  # gj length for each gj
        self.gjl[:] = cells.gj_len

        # produce the cell concentration array from the mems array:
        for arr in self.cc_mems:

            c_cells = np.dot(cells.M_sum_mems,arr)/cells.num_mems
            self.cc_cells.append(c_cells)

        self.cc_cells = np.asarray(self.cc_cells)

    def init_tissue(self, cells: 'Cells', p: 'Parameters') -> None:
        '''
        Prepares data structures pertaining to tissue profiles, dynamic
        activity, and optional methods such as electroosmotic fluid,
        which can be changed in between an initialization and simulation
        run.

        This method is called at the start of all simulations.
        '''

        # load in the gap junction dynamics object:
        self.gj_funk = Gap_Junction(p)

        self.J_gj_x = np.zeros(len(cells.mem_i))
        self.J_gj_y = np.zeros(len(cells.mem_i))

        if p.sim_ECM is True:
            #  Initialize diffusion constants for the extracellular transport:
            self.initDenv(cells,p)

        self.dyna = TissueHandler(self, cells, p)   # create the tissue dynamics object
        self.dyna.tissueProfiles(self, cells, p)  # initialize all tissue profiles

        if p.sim_ECM is True:
            # create a copy-base of the environmental junctions diffusion constants:
            self.D_env_base = copy.copy(self.D_env)

            # initialize current vectors
            self.J_env_x = np.zeros(len(cells.xypts))
            self.J_env_y = np.zeros(len(cells.xypts))

        # Initialize an array structure that will hold user-scheduled changes to membrane permeabilities:
        Dm_cellsA = np.asarray(self.Dm_cells)
        Dm_cellsER = np.copy(self.Dm_er)

        # if tb.emptyDict(p.scheduled_options) is False or tb.emptyDict(p.vg_options) is False or p.Ca_dyn is True:
        self.Dm_base = np.copy(Dm_cellsA) # make a copy that will serve as the unaffected values base

        # if tb.emptyDict(p.scheduled_options) is False:
        self.Dm_scheduled = np.copy(Dm_cellsA)
        self.Dm_scheduled[:] = 0

        # if tb.emptyDict(p.vg_options) is False:
            # Initialize an array structure that will hold dynamic voltage-gated channel changes to mem permeability:
        self.Dm_vg = np.copy(Dm_cellsA)
        self.Dm_vg[:] = 0

        # if p.Ca_dyn is True:
        # Initialize an array structure that will hold dynamic calcium-gated channel changes to mem perms:
        self.Dm_cag = np.copy(Dm_cellsA)
        self.Dm_cag[:] = 0

        self.Dm_stretch = np.copy(Dm_cellsA)   # array for stretch activated ion channels...
        self.Dm_stretch[:] = 0

        self.Dm_er_base = np.copy(Dm_cellsER)

        self.Dm_er_CICR = np.copy(Dm_cellsER)
        self.Dm_er_CICR[:] = 0

        self.dcc_ER = []

        self.Dm_morpho = np.copy(Dm_cellsA)
        self.Dm_morpho[:] = 0

        self.P_mod = np.copy(self.P_cells[:])
        self.P_base = np.copy(self.P_cells[:])

        # Noise Initialization ---------------------
        # add channel noise to the model:
        self.channel_noise_factor = np.random.random(self.mdl)

        self.Dm_cells[self.iK] = (p.channel_noise_level * self.channel_noise_factor + 1) * self.Dm_cells[self.iK]

        if p.dynamic_noise is True:
            # add a random walk on protein concentration to generate dynamic noise:
            self.protein_noise_factor = p.dynamic_noise_level * (np.random.random(self.mdl) - 0.5)

            if p.ions_dict['P'] == 1:
                self.cc_mems[self.iP] = self.cc_mems[self.iP] * (1 + self.protein_noise_factor)

        #--Blocks initialization--------------------

        if p.global_options['gj_block'] != 0:

            self.gj_block = np.ones(len(cells.mem_i))   # initialize the gap junction blocking vector to ones

        else:

            self.gj_block = 1

        # initialize dynamic pump blocking vectors:

        if p.global_options['NaKATP_block'] != 0:

            self.NaKATP_block = np.ones(self.mdl)  # initialize NaKATP blocking vector

        else:
            self.NaKATP_block = 1

        if p.HKATPase_dyn is True and p.global_options['HKATP_block'] != 0:
            self.HKATP_block = np.ones(self.mdl)  # initialize HKATP blocking vector
        else:
            self.HKATP_block = 1

        if p.VATPase_dyn is True and p.global_options['VATP_block'] != 0:
            self.VATP_block = np.ones(self.mdl)  # initialize HKATP blocking vector
        else:
            self.VATP_block = 1

        # -----auxiliary molecules initialization -------------------------

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:

            self.cIP3_mems = np.zeros(self.mdl)  # initialize a vector to hold IP3 concentrations
            self.cIP3_mems[:] = p.cIP3_to                 # set the initial concentration of IP3 from params file

            self.cIP3_cells = np.zeros(self.cdl)  # initialize a vector to hold IP3 concentrations
            self.cIP3_cells[:] = p.cIP3_to                 # set the initial concentration of IP3 from params file

            self.IP3_flux_x_gj = np.zeros(len(cells.nn_i))
            self.IP3_flux_y_gj = np.zeros(len(cells.nn_i))
            self.IP3_flux_mem = np.zeros(len(cells.mem_i))

            if p.sim_ECM is True:

                self.cIP3_env = np.zeros(self.edl)
                self.cIP3_env[:] = p.cIP3_to_env

                self.cIP3_flux_env_x = np.zeros(self.edl)
                self.cIP3_flux_env_y = np.zeros(self.edl)

            elif p.sim_ECM is False:
                self.cIP3_env = np.zeros(self.mdl)     # initialize IP3 concentration of the environment
                self.cIP3_env[:] = p.cIP3_to_env

        if p.voltage_dye is True:

            self.cDye_mems = np.zeros(self.mdl)   # initialize voltage sensitive dye array for mem regions
            self.cDye_mems[:] = p.cDye_to_cell

            self.cDye_cells = np.zeros(self.cdl)   # initialize voltage sensitive dye array for cell centroids
            self.cDye_cells[:] = p.cDye_to_cell

            self.Dye_flux_x_gj = np.zeros(len(cells.nn_i))
            self.Dye_flux_y_gj = np.zeros(len(cells.nn_i))
            self.Dye_flux_mem = np.zeros(len(cells.mem_i))

            self.c_dye_bound = p.cDye_to    # concentration of dye at the global boundaries

            self.cDye_env = np.zeros(self.mdl)
            self.cDye_env[:] = p.cDye_to

            if p.Dye_target_channel != 'None':

                # get the ion index of the target ion:
                self.dye_target = self.get_ion(p.Dye_target_channel)

                # make a copy of the appropriate ion list
                self.Dm_mod_dye = np.copy(self.Dm_cells[self.dye_target][:])

            else:

                self.dye_target = None

            if p.sim_ECM is True:

                self.cDye_env = np.zeros(self.edl)
                self.cDye_env[:] = p.cDye_to

                self.Dye_flux_env_x = np.zeros(self.edl)
                self.Dye_flux_env_y = np.zeros(self.edl)

            else:
                self.Dye_env = np.zeros(self.mdl)     # initialize Dye concentration in the environment
                self.Dye_env[:] = p.cDye_to

        #-----dynamic creation/anhilation of large Laplacian matrix computators!------------------

        if p.fluid_flow is True: # If at any time fluid flow is true, initialize the flow vectors to zeros
            # initialize data structures for flow:
            self.u_cells_x = np.zeros(self.cdl)
            self.u_cells_y = np.zeros(self.cdl)

            self.u_gj_x = np.zeros(self.mdl)
            self.u_gj_y = np.zeros(self.mdl)

        if p.calc_J is True:

            if cells.lapGJ_P_inv is None:

                # make a laplacian and solver for discrete transfers on closed, irregular cell network
                logs.log_info('Creating cell network Poisson solver...')
                cells.graphLaplacian(p)

            if p.sim_ECM is True and cells.lapENV_P_inv is None:

                logs.log_info('Creating environmental Poisson solver...')
                bdic = {'N': 'flux', 'S': 'flux', 'E': 'flux', 'W': 'flux'}
                cells.lapENV_P, cells.lapENV_P_inv = cells.grid_obj.makeLaplacian(bound=bdic)

                cells.lapENV_P = None  # get rid of the non-inverse matrix as it only hogs memory...

        elif p.fluid_flow is False and p.deformation is False and p.calc_J is False:
            # if there's no physics calcs required, we're running an init
            # get rid of all matrices rather than carry them around as large dead-weights

            if cells.lapGJinv is not None:
                cells.lapGJinv = None
                cells.lapGJ_P_inv = None
                cells.lapGJ = None
                cells.lapGJ_P = None

            if p.sim_ECM is True and cells.lapENV_P_inv is not None:
                cells.lapENVinv = None
                cells.lapENV_P_inv = None


        if p.run_sim is True:   # if we're running a simulation:

            if p.fluid_flow is True or p.deformation is True or p.calc_J is True:  # if user desires fluid flow:

                if p.fluid_flow is True or p.deformation is True:

                    if p.deformation is True:
                        cells.deform_tools(p)

                    # initialize data structures for flow:
                    self.u_cells_x = np.zeros(self.cdl)
                    self.u_cells_y = np.zeros(self.cdl)

                    self.u_gj_x = np.zeros(self.mdl)
                    self.u_gj_y = np.zeros(self.mdl)

                    # force electrostatic pressure:
                    p.deform_electro = True

                if cells.lapGJinv is None:

                    # make a laplacian and solver for discrete transfers on closed, irregular cell network
                    logs.log_info('Creating cell network Poisson solver...')
                    cells.graphLaplacian(p)

                    if p.td_deform is False: # if time-dependent deformation is not required

                        cells.lapGJ = None
                        cells.lapGJ_P = None       # null out the non-inverse matrices -- we don't need them

                    elif p.td_deform is True and cells.lapGJ is None:  # unforunately, we need to create them all again

                        logs.log_info('Creating cell network Poisson solver...')
                        cells.graphLaplacian(p)


                if p.sim_ECM is True:

                    if cells.lapENV_P_inv is None:

                        logs.log_info('Creating environmental Poisson solver...')
                        bdic = {'N': 'flux', 'S': 'flux', 'E': 'flux', 'W': 'flux'}
                        cells.lapENV_P, cells.lapENV_P_inv = cells.grid_obj.makeLaplacian(bound=bdic)

                        cells.lapENV_P = None  # get rid of the non-inverse matrix as it only hogs memory...

                    self.u_env_x = np.zeros(cells.X.shape)
                    self.u_env_y = np.zeros(cells.X.shape)


            elif p.fluid_flow is False and p.deformation is False and p.calc_J is False:
                # if there's no physics calcs required
                # get rid of all matrices rather than carry them around as large dead-weights

                if cells.lapGJ_P_inv is not None:
                    cells.lapGJinv = None
                    cells.lapGJ_P_inv = None
                    cells.lapGJ = None
                    cells.lapGJ_P = None

                if p.sim_ECM is True and cells.lapENV_P_inv is not None:
                    cells.lapENVinv = None
                    cells.lapENV_P_inv = None

        # if simulating electrodiffusive movement of membrane pumps and channels:-------------

        if p.sim_eosmosis is True:

            if cells.gradMem is None:
                cells.eosmo_tools(p)

            self.rho_pump = np.ones(len(cells.mem_i))
            self.rho_channel = np.ones(len(cells.mem_i))

        else:
            self.rho_pump = 1  # else just define it as identity.
            self.rho_channel = 1

        # Initialize all user-specified interventions and dynamic channels.
        self.dyna.runAllInit(self,cells,p)

    def run_sim_core(self, cells: 'Cells', p: 'Parameters') -> None:
        '''
        Runs and saves the current simulation phase (e.g. init or sim phases)

        This method drives the time-loop for the main simulation, including
        gap-junction connections and calls to all dynamic channels and entities (sim phase
        only).

        '''

        self.init_tissue(cells, p)  # Initialize all structures used for gap junctions, ion channels, and other dynamics

        # Reinitialize all time-data structures
        self.clear_storage(cells, p)

        # get the net, unbalanced charge and corresponding voltage in each cell to initialize values of voltages:
        self.update_V(cells, p)

        self.vm_to = np.copy(self.vm)  # create a copy of the original voltage

        if p.Ca_dyn is True:
            self.cc_er_to = np.copy(self.cc_er)

        # Display and/or save an animation during solving and calculate:
        #
        # * "tt", a time-steps vector appropriate for the current run.
        # * "tsamples", that vector resampled to save data at fewer times.
        tt, tsamples = self._plot_loop(cells, p)

        do_once = True  # a variable to time the loop only once

        for t in tt:  # run through the loop

            # start the timer to approximate time for the simulation
            if do_once is True:
                loop_measure = time.time()

            # Reinitialize flux storage device.
            self.fluxes_mem.fill(0)

            # Calculate the change in the voltage derivative.
            self.dvm = (self.vm - self.vm_to) / p.dt

            self.vm_to = np.copy(self.vm)  # reassign the history-saving vm

            if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:
                self.dcc_ER = (self.cc_er - self.cc_er_to) / p.dt

                self.cc_er_to = np.copy(self.cc_er)

            # Calculate the values of scheduled and dynamic quantities (e.g.
            # ion channel multipliers).
            if p.run_sim is True:
                self.dyna.runAllDynamics(self, cells, p, t)

            # -----------------PUMPS-------------------------------------------------------------------------------------

            if p.sim_ECM is True:
                # run the Na-K-ATPase pump:
                fNa_NaK, fK_NaK, self.rate_NaKATP = stb.pumpNaKATP(
                    self.cc_mems[self.iNa],
                    self.cc_env[self.iNa][cells.map_mem2ecm],
                    self.cc_mems[self.iK],
                    self.cc_env[self.iK][cells.map_mem2ecm],
                    self.vm,
                    self.T,
                    p,
                    self.NaKATP_block,
                )

            else:

                fNa_NaK, fK_NaK, self.rate_NaKATP = stb.pumpNaKATP(
                            self.cc_mems[self.iNa],
                            self.cc_env[self.iNa],
                            self.cc_mems[self.iK],
                            self.cc_env[self.iK],
                            self.vm,
                            self.T,
                            p,
                            self.NaKATP_block,
                        )

            # modify the fluxes by electrodiffusive membrane redistribution factor and add fluxes to storage:
            self.fluxes_mem[self.iNa] = self.rho_pump * fNa_NaK
            self.fluxes_mem[self.iK] = self.rho_pump * fK_NaK

            # update the concentrations of Na and K in cells and environment:
            self.cc_mems[self.iNa][:],  self.cc_env[self.iNa][:] =  stb.update_Co(self, self.cc_mems[self.iNa][:],
                                                                        self.cc_env[self.iNa][:],fNa_NaK, cells, p)

            self.cc_mems[self.iK][:], self.cc_env[self.iK][:] = stb.update_Co(self, self.cc_mems[self.iK][:],
                self.cc_env[self.iK][:], fK_NaK, cells, p)

            # update the concentrations intra-cellularly:
            self.cc_mems[self.iNa][:],  self.cc_cells[self.iNa][:], _ = \
                stb.update_intra(self, cells, self.cc_mems[self.iNa][:],
                                 self.cc_cells[self.iNa][:],
                                 self.D_free[self.iNa],
                                 self.zs[self.iNa], p)

            self.cc_mems[self.iK][:],  self.cc_cells[self.iK][:], _ = \
                stb.update_intra(self, cells, self.cc_mems[self.iK][:],
                                 self.cc_cells[self.iK][:],
                                 self.D_free[self.iK],
                                 self.zs[self.iK], p)


            # recalculate the net, unbalanced charge and voltage in each cell:
            self.update_V(cells, p)


            # ----------------ELECTRODIFFUSION---------------------------------------------------------------------------

            # electro-diffuse all ions (except for proteins, which don't move) across the cell membrane:

            shuffle(self.movingIons)

            for i in self.movingIons:

                if p.sim_ECM is True:

                    f_ED = stb.electroflux(self.cc_env[i][cells.map_mem2ecm], self.cc_mems[i],
                        self.Dm_cells[i], p.tm, self.zs[i], self.vm, self.T, p,
                        rho=self.rho_channel)

                else:

                    f_ED = stb.electroflux(self.cc_env[i],self.cc_mems[i],self.Dm_cells[i],p.tm,self.zs[i],
                                    self.vm,self.T,p,rho=self.rho_channel)


                self.fluxes_mem[i] = self.fluxes_mem[i] + f_ED

                # update ion concentrations in cell and ecm:

                self.cc_mems[i][:], self.cc_env[i][:] = stb.update_Co(self, self.cc_mems[i][:],
                    self.cc_env[i][:], f_ED, cells, p)

                # update the ion concentration intra-cellularly:
                self.cc_mems[i][:], self.cc_cells[i][:], _ = \
                    stb.update_intra(self, cells, self.cc_mems[i][:],
                        self.cc_cells[i][:],
                        self.D_free[i],
                        self.zs[i], p)

                # update flux between cells due to gap junctions
                self.update_gj(cells, p, t, i)

                if p.sim_ECM:
                    #update concentrations in the extracellular spaces:
                    self.update_ecm(cells, p, t, i)

            # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V(cells, p)

            # ----transport and handling of special ions---------------------------------------------------------------

            if p.ions_dict['Ca'] == 1:

                self.ca_handler(cells, p)

            if p.ions_dict['H'] == 1:
                self.acid_handler(cells, p)

            if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:

                self.cIP3_mems, self.cIP3_env, _, _, _, _, _ = stb.molecule_mover(self, self.cIP3_mems,
                    self.cIP3_env, cells, p, p.z_IP3, p.Dm_IP3,
                    p.Do_IP3, 0)

                # update concentrations intracellularly:
                self.cIP3_mems, self.cIP3_cells, _ = \
                    stb.update_intra(self, cells, self.cIP3_mems[:],
                        self.cIP3_cells[:],
                        p.Do_IP3,
                        p.z_IP3, p)

            # if p.voltage_dye=1 electrodiffuse voltage sensitive dye between cell and environment
            if p.voltage_dye == 1:

                if p.pump_Dye is True:

                    self.cDye_mems, self.cDye_env, fx = stb.molecule_pump(self, self.cDye_mems, self.cDye_env,
                                                                    cells, p, Df=p.Do_Dye, z=p.z_Dye,
                                                                    pump_into_cell=p.pump_Dye_in,
                                                                    alpha_max=p.pump_Dye_alpha, Km_X=1.0e-3,
                                                                    Km_ATP=1.0)

                self.cDye_mems, self.cDye_env, _, _, _, _, _ = stb.molecule_mover(self,self.cDye_mems,
                                                                    self.cDye_env,cells, p, z=p.z_Dye, Dm = p.Dm_Dye,
                                                                    Do = p.Do_Dye, c_bound = self.c_dye_bound,
                                                                    ignoreECM = True)

                # update concentrations intracellularly:
                self.cDye_mems, self.cDye_cells, _ = \
                    stb.update_intra(self, cells, self.cDye_mems, self.cDye_cells, p.Do_Dye, p.z_Dye, p)

            # dynamic noise handling-----------------------------------------------------------------------------------

            if p.dynamic_noise == 1 and p.ions_dict['P'] == 1:
                # add a random walk on protein concentration to generate dynamic noise:
                self.protein_noise_factor = p.dynamic_noise_level * (np.random.random(self.mdl) - 0.5)
                self.cc_mems[self.iP] = self.cc_mems[self.iP] * (1 + self.protein_noise_factor)

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V(cells, p)

            #-----forces, fields, and flow-----------------------------------------------------------------------------

            # main update of voltage in cells:
            # self.update_V(cells, p)  # Main update

            self.get_Efield(cells, p)  # FIXME update to also get a v_cell gradient

            # store charge in time vector:
            rho_cells_ave = np.dot(cells.M_sum_mems, self.rho_cells)/cells.num_mems
            self.charge_cells_time.append(rho_cells_ave)

            if p.sim_ECM:

                if p.smooth_level > 0.0:
                    # smooth the charge out as a derivative will be taken on it:
                    rho_env_sm = gaussian_filter(self.rho_env[:].reshape(cells.X.shape),p.smooth_level).ravel()

                else:
                    rho_env_sm = self.rho_env[:]

                # as flux is done in terms of env-grid squares, correct the volume density of charge:
                rho_env_sm = (cells.true_ecm_vol/cells.ecm_vol)*rho_env_sm

                self.charge_env_time.append(rho_env_sm)

            # get the current
            if p.calc_J:
                get_current(self, cells, p)


            # get forces from any hydrostatic (self.P_Cells) pressure:
            getHydroF(self,cells, p)

            # calculate specific forces and pressures:

            if p.deform_osmo is True:
                osmotic_P(self,cells, p)

            if p.deform_electro is True:
                electro_F(self,cells, p)

            if p.fluid_flow is True and p.run_sim is True:

                self.run_sim = True

                getFlow(self,cells, p)

            # if desired, electroosmosis of membrane channels
            if p.sim_eosmosis is True and p.run_sim is True:

                self.run_sim = True

                eosmosis(self,cells, p)  # modify membrane pump and channel density according to Nernst-Planck

            if p.deformation is True and p.run_sim is True:

                self.run_sim = True

                if p.td_deform is False:

                    getDeformation(self,cells, t, p)

                elif p.td_deform is True:

                    timeDeform(self,cells, t, p)

            stb.check_v(self.vm)

            # ---------time sampling and data storage---------------------------------------------------

            if t in tsamples:

                self.write2storage(t, cells, p)  # write data to time storage vectors

                # Update the currently displayed and/or saved animation.
                self._replot_loop(p)

            # Get time for loop and estimate total time for simulation.
            if do_once is True:
                loop_time = time.time() - loop_measure

                if p.run_sim is True:
                    time_estimate = round(loop_time * p.sim_tsteps, 2)
                else:
                    time_estimate = round(loop_time * p.init_tsteps, 2)
                logs.log_info("This run should take approximately " + str(time_estimate) + ' s to compute...')
                do_once = False

        # Find embedded functions that can't be pickled...
        fh.safe_pickle(self, p)

        cells.points_tree = None

        # Explicitly close the prior animation to conserve memory.
        self._deplot_loop()

        # save the init or sim and report results of potential interest to the user.
        self.save_and_report(cells, p)

        plt.close()
        logs.log_info('Simulation completed successfully.')

    #.................{  INITIALIZERS & FINALIZERS  }............................................
    def clear_storage(self, cells, p):
        """
        Re-initializes time storage vectors at the begining of a sim or init.

        """


        # clear mass flux storage vectors:
        self.fluxes_gj  = np.zeros(self.fluxes_gj.shape)
        self.fluxes_mem = np.zeros(self.fluxes_mem.shape)

        self.cc_time = []  # data array holding the concentrations at time points
        self.cc_env_time = [] # data array holding environmental concentrations at time points

        self.dd_time = []  # data array holding membrane permeabilites at time points
        self.vm_time = []  # data array holding voltage at time points
        self.vm_ave_time = []   # data array holding average vm (averaged to cell centres) FIXME temp hack
        self.vm_GHK_time = [] # data array holding GHK vm estimates
        self.dvm_time = []  # data array holding derivative of voltage at time points
        self.vcell_time = []
        self.time = []     # time values of the simulation
        self.gjopen_time = []   # stores the fractional gap junction open state at each time
        self.osmo_P_delta_time = []  # osmotic pressure difference between cell interior and exterior as func of time

        self.I_mem_time = []    # initialize membrane current time vector

        self.vm_Matrix = []    # initialize matrices for resampled data sets (used in smooth plotting and streamlines)
        self.I_gj_x_time = []
        self.I_gj_y_time = []
        self.I_tot_x_time = []
        self.I_tot_y_time = []

        self.efield_gj_x_time = []   # matrices storing smooth electric field in gj connected cells
        self.efield_gj_y_time = []

        self.P_cells_time = []
        self.u_cells_x_time = []
        self.u_cells_y_time = []

        self.rho_cells_time = []
        self.charge_cells_time =[]

        self.F_electro_time = []
        self.F_electro_x_time = []
        self.F_electro_y_time = []
        self.F_hydro_x_time = []
        self.F_hydro_y_time = []
        self.F_hydro_time = []
        self.P_electro_time = []

        self.rate_NaKATP_time =[]

        if p.deformation is True and p.run_sim is True:
            self.ecm_verts_unique_to = cells.ecm_verts_unique[:] # make a copy of original ecm verts as disp ref point

            self.cell_centres_time = []
            self.mem_mids_time = []
            self.maskM_time = []
            self.mem_edges_time = []
            self.cell_verts_time = []

            self.dx_cell_time = []
            self.dy_cell_time = []

            self.dx_time = []
            self.dy_time = []

            self.phi = np.zeros(self.mdl)
            self.phi_time = []

        if p.voltage_dye is True:
            self.cDye_time = []    # retains voltage-sensitive dye concentration as a function of time

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:
            self.cIP3_time = []    # retains IP3 concentration as a function of time

        if p.sim_eosmosis is True:
            self.rho_channel_time = []
            self.rho_pump_time = []

        if p.Ca_dyn is True:
            self.cc_er_time = []
            self.cc_er_to = np.copy(self.cc_er[:])

        if p.sim_ECM is True:

            # clear flux storage vectors for environment
            self.fluxes_env_x = np.zeros(self.fluxes_env_x.shape)
            self.fluxes_env_y = np.zeros(self.fluxes_env_y.shape)


            self.venv_time = []

            self.charge_env_time = []

            self.cc_er_time = []   # retains er concentrations as a function of time
            self.cIP3_time = []    # retains cellular ip3 concentrations as a function of time
            self.cIP3_env_time = []

            self.efield_ecm_x_time = []   # matrices storing smooth electric field in ecm
            self.efield_ecm_y_time = []

            # initialize time-storage vectors for electroosmotic data:
            self.u_env_x_time = []
            self.u_env_y_time = []

            self.rho_pump_time = []    # store pump and channel states as function of time...
            self.rho_channel_time = []

            vm_dato = np.zeros(len(cells.mem_i))
            dat_grid_vm = stb.vertData(vm_dato,cells,p)
            self.vm_Matrix.append(dat_grid_vm[:])

            self.cDye_env_time = []

    def write2storage(self,t,cells,p):

        if p.GHK_calc is True:
                stb.ghk_calculator(self,cells,p)
                self.vm_GHK_time.append(self.vm_GHK) # data array holding GHK vm estimates

        # add the new concentration and voltage data to the time-storage matrices:
        self.efield_gj_x_time.append(self.E_gj_x[:])

        self.efield_gj_y_time.append(self.E_gj_y[:])

        concs = np.copy(self.cc_mems[:])
        self.cc_time.append(concs)
        concs = None

        envsc = np.copy(self.cc_env[:])
        self.cc_env_time.append(envsc)
        envsc = None

        ddc = np.copy(self.Dm_cells[:])
        ddc.tolist()
        self.dd_time.append(ddc)
        ddc = None

        self.I_gj_x_time.append(self.J_gj_x[:])
        self.I_gj_y_time.append(self.J_gj_y[:])

        self.I_mem_time.append(self.I_mem[:])

        self.vm_time.append(self.vm[:])

        v_cells_ave = np.dot(cells.M_sum_mems, self.v_cell[:]) / cells.num_mems
        self.vcell_time.append(v_cells_ave)

        # FIXME this is a temporary hack until I have a chance to redo plotting for intracellular space transport...
        vm_ave = np.dot(cells.M_sum_mems,self.vm)/cells.num_mems
        self.vm_ave_time.append(vm_ave)

        self.dvm_time.append(self.dvm[:])

        # FIXME this is a temporary hack until I have a chance to redo plotting for intracellular space transport...
        rho_cells_ave = np.dot(cells.M_sum_mems,self.rho_cells)/cells.num_mems
        self.rho_cells_time.append(rho_cells_ave)

        self.rate_NaKATP_time.append(self.rate_NaKATP[:])
        self.P_cells_time.append(self.P_cells[:])
        self.F_hydro_x_time.append(self.F_hydro_x[:])
        self.F_hydro_y_time.append(self.F_hydro_y[:])

        if p.deform_osmo is True:
            self.osmo_P_delta_time.append(self.osmo_P_delta[:])

        if p.deform_electro is True:
            self.F_electro_time.append(self.F_electro[:])
            self.F_electro_x_time.append(self.F_electro_x[:])
            self.F_electro_y_time.append(self.F_electro_y[:])

            self.P_electro_time.append(self.P_electro[:])

        if p.deformation is True and p.run_sim is True:
            implement_deform_timestep(self,cells, t, p)


        if p.fluid_flow is True and p.run_sim is True:
            self.u_cells_x_time.append(self.u_cells_x[:])
            self.u_cells_y_time.append(self.u_cells_y[:])

        if p.sim_eosmosis is True and p.run_sim is True:
            self.rho_channel_time.append(self.rho_channel[:])
            self.rho_pump_time.append(self.rho_pump[:])

        if p.v_sensitive_gj is True:
            self.gjopen_time.append(self.gjopen[:])

        else:
            self.gjopen_time.append(self.gjopen)

        self.time.append(t)

        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:
            self.cIP3_time.append(self.cIP3_mems[:])

        if p.voltage_dye ==1:
            # FIXME this is a temporary hack until I have a chance to redo plotting for intracellular space transport...
            cDye_ave = np.dot(cells.M_sum_mems, self.cDye_mems) / cells.num_mems
            self.cDye_time.append(cDye_ave)

        if p.Ca_dyn == 1 and p.ions_dict['Ca']==1:
            self.cc_er_time.append(np.copy(self.cc_er[:]))

        if p.sim_ECM is True:

            self.efield_ecm_x_time.append(self.E_env_x[:])

            self.efield_ecm_y_time.append(self.E_env_y[:])

            self.I_tot_x_time.append(self.J_env_x[:])
            self.I_tot_y_time.append(self.J_env_y[:])

            # ecmsc = np.copy(self.cc_env[:])
            # ecmsc.tolist()
            # self.cc_env_time.append(ecmsc)
            # ecmsc = None

            # FIXME this is a temporary hack until I have a chance to redo plotting for intracellular space transport...


            self.venv_time.append(self.v_env[:])

            if p.fluid_flow is True and p.run_sim is True:
                self.u_env_x_time.append(self.u_env_x[:])
                self.u_env_y_time.append(self.u_env_y[:])

            # calculate interpolated verts and midpoint data for Vmem:
            dat_grid_vm = stb.vertData(self.vm[:],cells,p)

            self.vm_Matrix.append(dat_grid_vm[:])

            if p.voltage_dye ==1:
                self.cDye_env_time.append(self.cDye_env[:])

    def save_and_report(self,cells,p):

        # save the init or sim:
        if p.run_sim is False:
            datadump = [self, cells, p]
            fh.saveSim(self.savedInit, datadump)
            message_1 = 'Initialization run saved to' + ' ' + p.init_path
            logs.log_info(message_1)
        else:
            datadump = [self, cells, p]
            fh.saveSim(self.savedSim, datadump)
            message_2 = 'Simulation run saved to' + ' ' + p.sim_path
            logs.log_info(message_2)

        # report final output to user:
        for i in range(0,len(self.ionlabel)):
            endconc = np.round(np.mean(self.cc_time[-1][i]),6)
            label = self.ionlabel[i]
            concmess = 'Final average cytoplasmic concentration of'+ ' '+ label + ': '
            logs.log_info(concmess + str(endconc) + ' mmol/L')

        for i in range(0,len(self.ionlabel)):
            endconc = np.round(np.mean(self.cc_env_time[-1][i]),6)
            label = self.ionlabel[i]
            concmess = 'Final environmental concentration of'+ ' '+ label + ': '
            logs.log_info(concmess + str(endconc) + ' mmol/L')

        final_vmean = 1000*np.round(np.mean(self.vm_time[-1]),6)
        vmess = 'Final average cell Vmem of ' + ': '
        logs.log_info(vmess + str(final_vmean) + ' mV')

        if p.GHK_calc is True:
            final_vmean_GHK = 1000*np.round(np.mean(self.vm_GHK_time[-1]),6)
            vmess = 'Final average cell Vmem calculated using GHK: ' + ': '
            logs.log_info(vmess + str(final_vmean_GHK) + ' mV')

        if p.ions_dict['H'] == 1:
            final_pH = -np.log10(1.0e-3*np.mean((self.cc_time[-1][self.iH])))
            logs.log_info('Final average cell pH ' + str(np.round(final_pH, 2)))

            final_pH_env = -np.log10(np.mean(1.0e-3*(self.cc_env_time[-1][self.iH])))
            logs.log_info('Final environmental pH ' + str(np.round(final_pH_env, 2)))


        if p.scheduled_options['IP3'] != 0 or p.Ca_dyn is True:
            IP3_env_final = np.mean(self.cIP3_env)
            IP3_cell_final = np.mean(self.cIP3)
            logs.log_info('Final IP3 concentration in the environment: ' + str(np.round(IP3_env_final, 6)) + ' mmol/L')
            logs.log_info('Final average IP3 concentration in cells: ' + str(np.round(IP3_cell_final, 6)) + ' mmol/L')

        if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:

            endconc_er = np.round(np.mean(self.cc_er[0]),6)
            label = self.ionlabel[self.iCa]
            concmess = 'Final average ER concentration of'+ ' '+ label + ': '
            logs.log_info(concmess + str(endconc_er) + ' mmol/L')

        if p.voltage_dye ==1:
            dye_env_final = np.mean(self.cDye_env)
            dye_cell_final = np.mean(self.cDye_mems)
            logs.log_info('Final average morphogen concentration in the environment: ' + str(np.round(dye_env_final, 6))
                          + ' mmol/L')
            logs.log_info('Final average morphogen concentration in cells: ' + str(np.round(dye_cell_final, 6)) +
                             ' mmol/L')

    def sim_info_report(self,cells,p):

        logs.log_info('This world contains ' + str(cells.cell_number) + ' cells.')
        logs.log_info('Each cell has an average of ' + str(round(cells.average_nn, 2)) + ' nearest-neighbours.')

        logs.log_info('You are running the ion profile: ' + p.ion_profile)

        logs.log_info('Ions in this simulation: ' + str(self.ionlabel))
        logs.log_info(
            'If you have selected features using other ions, '
            'they will be ignored.')

        logs.log_info('Considering extracellular spaces: ' + str(p.sim_ECM))

        if p.sim_ECM:
            logs.log_info('Cells per env grid square ' + str(round(cells.ratio_cell2ecm, 2)))

        logs.log_info('Electroosmotic fluid flow: ' + str(p.fluid_flow))
        logs.log_info('Ion pump and channel electodiffusion in membrane: ' + str(p.sim_eosmosis))
        logs.log_info('Force-induced cell deformation: ' + str(p.deformation))
        logs.log_info('Osmotic pressure: ' + str(p.deform_osmo))
        logs.log_info('Electrostatic pressure: ' + str(p.deform_electro))

    # ................{ CORE DOERS & GETTERS }...............................................
    def get_Vall(self, cells, p) -> (np.ndarray, np.ndarray, np.ndarray):
        """
        Calculates transmembrane voltage (Vmem) and voltages inside the cell
        and extracellular space, assuming the cell is a capacitor consisting of
        two concentric spheres.

        Parameters
        --------
        cells:
        p:

        Returns
        --------
        vm : np.ndarray
            Transmembrane voltage.
        v_cell : np.ndarray
            Voltage at the inner membrane surface.
        v_env : np.ndarray
            Voltage at the outer membrane surface.
        """

        if p.sim_ECM is False:
            # if we're not modelling extracellular spaces, but assuming a perfectly mixing environment,
            # assume the cell is a concentric spherical capacitor with opposite but equal magnitude
            # charge inside the cell and out
            # calculate the voltage for the system using the measured capacitance of the cell membrane:

            # the vm at each membrane segment of a cell:
            v_cello = (1/(p.cm*cells.mem_sa))*self.rho_cells*cells.mem_vol

            v_cell_aveo = np.dot(cells.M_sum_mems, v_cello) / cells.num_mems

            v_env = 0

            # use finite volume method to integrate each region:
            # values at centroid mids:
            vcell_at_mids = (v_cello + v_cell_aveo[cells.mem_to_cells]) / 2

            # finite volume integral of membrane pie-box values:
            v_cell = np.dot(cells.M_int_mems, v_cello) + (1 / 2) * vcell_at_mids
            v_cell_ave = (1 / 2) * v_cell_aveo + np.dot(cells.M_sum_mems, vcell_at_mids) / (2 * cells.num_mems)

            vm = v_cell

        else:

            Qcells = self.rho_cells * cells.mem_vol

            Qecm = self.rho_env[cells.envInds_inClust] * cells.true_ecm_vol[cells.envInds_inClust]

            # concatenate the cell and ecm charge vectors to the maxwell capacitance vector:
            Q_max_vect = np.hstack((Qcells,Qecm))

           # original solver in terms of pseudo-inverse matrix:
            v_max_vect = np.dot(cells.M_max_cap_inv, Q_max_vect)

            # separate voltages for cells and ecm spaces
            v_cell = v_max_vect[cells.mem_range_a:cells.mem_range_b]
            v_ecm = v_max_vect[cells.env_range_a:cells.env_range_b]

            # use finite volume method to integrate each region of intracellular voltage:
            # values at centroid mids:
            v_cell_ave = np.dot(cells.M_sum_mems, v_cell) / cells.num_mems

            vcell_at_mids = (v_cell + v_cell_ave[cells.mem_to_cells]) / 2

            # finite volume integral of membrane pie-box values:
            v_cell = np.dot(cells.M_int_mems, v_cell) + (1 / 2) * vcell_at_mids
            v_cell_ave = (1 / 2) * v_cell_ave + np.dot(cells.M_sum_mems, vcell_at_mids) / (2 * cells.num_mems)

            # define the full environmental voltage:
            v_env = np.zeros(len(cells.xypts))
            v_env[cells.envInds_inClust] = v_ecm

            # finite volume integration of voltage:
            v_env = np.dot(cells.gridInt, v_env)

            # optional smoothing of voltage using gaussian:
            if p.smooth_level > 0.0:

                v_env = gaussian_filter(v_env.reshape(cells.X.shape),p.smooth_level).ravel()

            # v_cell = v_cell_ave[cells.mem_to_cells]

            # # set the conditions for the global boundaries:
            v_env[cells.bBot_k] = self.bound_V['B']
            v_env[cells.bTop_k] = self.bound_V['T']
            v_env[cells.bL_k] = self.bound_V['L']
            v_env[cells.bR_k] = self.bound_V['R']

            # calculate the vm
            vm = v_cell - v_env[cells.map_mem2ecm]


        return vm, v_cell, v_cell_ave, v_env

    def update_V(self,cells,p):

        if p.sim_ECM is True:
            # get the charge in cells and the environment:
            self.rho_cells = stb.get_charge_density(self.cc_mems, self.z_array, p)
            self.rho_env = stb.get_charge_density(self.cc_env, self.z_array_env, p)
            # self.rho_env = gaussian_filter(self.rho_env.reshape(cells.X.shape),p.smooth_level).ravel()
            # self.rho_env[cells.inds_env] = 0 # assumes charge screening in the bulk env

            self.vm, self.v_cell, self.v_cell_ave, self.v_env = self.get_Vall(cells,p)
            # self.v_env[cells.inds_env] = 0  # assumes charge screening in the bulk env

        else:
             self.rho_cells = stb.get_charge_density(self.cc_mems, self.z_array, p)
             self.vm, self.v_cell, self.v_cell_ave, _ = self.get_Vall(cells,p)

    def acid_handler(self,cells,p):

        # electrofuse the H+ ion between the cytoplasm and the environment
        if p.sim_ECM is True:

            # Electrofuse the H+ ion between the cytoplasm and the ecms.
            f_H1 = stb.electroflux(
                self.cc_env[self.iH][cells.map_mem2ecm],
                self.cc_mems[self.iH][cells.mem_to_cells],
                self.Dm_cells[self.iH],
                p.tm,
                self.zs[self.iH],
                self.vm,
                self.T,
                p,
            )

            self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + f_H1


        else:

            f_H1 = stb.electroflux(self.cc_env[self.iH],self.cc_mems[self.iH],self.Dm_cells[self.iH],p.tm,
                self.zs[self.iH],self.vm,self.T,p)

            self.fluxes_mem[self.iH] = f_H1[cells.mem_to_cells]

        # update H+ in cells and environment, first in absence of bicarbonate buffering:
        self.cc_mems[self.iH][:], self.cc_env[self.iH][:] = stb.update_Co(self, self.cc_mems[self.iH][:],
            self.cc_env[self.iH][:], f_H1, cells, p)


        # run the bicarbonate buffer to ensure realistic concentrations and pH in cell and environment:
        self.cc_mems[self.iH], _, self.cc_mems[self.iM], self.pH_cell = stb.bicarbonate_buffer(
                                                                                        self.cc_mems[self.iH],
                                                                                        self.cHM_mems,
                                                                                        self.cc_mems[self.iM],
                                                                                        p)

        self.cc_env[self.iH], _, self.cc_env[self.iM], self.pH_env = stb.bicarbonate_buffer(
                                                                                            self.cc_env[self.iH],
                                                                                            self.cHM_env,
                                                                                            self.cc_env[self.iM],
                                                                                            p)

        # recalculate the net, unbalanced charge and voltage in each cell:
        self.update_V(cells,p)

        if p.HKATPase_dyn == 1 and p.run_sim is True: # if there's an H,K ATPase pump

            if p.sim_ECM is True:

                f_H2, f_K2 = stb.pumpHKATP(self.cc_mems[self.iH][cells.mem_to_cells],
                    self.cc_env[self.iH][cells.map_mem2ecm],
                    self.cc_mems[self.iK][cells.mem_to_cells], self.cc_env[self.iK][cells.map_mem2ecm],
                    self.vm, self.T, p, self.HKATP_block)

                # modify fluxes by any uneven redistribution of pump location:
                f_H2 = self.rho_pump * f_H2
                f_K2 = self.rho_pump * f_K2

                self.HKATPase_rate = f_H2[:]

                self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + f_H2
                self.fluxes_mem[self.iK] = self.fluxes_mem[self.iK] + f_K2


            else:

                # if HKATPase pump is desired, run the H-K-ATPase pump:
                f_H2, f_K2 = stb.pumpHKATP(self.cc_mems[self.iH],self.cc_env[self.iH],self.cc_mems[self.iK],
                    self.cc_env[self.iK],self.vm,self.T,p,self.HKATP_block)

                # store fluxes for this pump:
                self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + self.rho_pump * f_H2[cells.mem_to_cells]
                self.fluxes_mem[self.iK] = self.fluxes_mem[self.iK] + self.rho_pump * f_K2[cells.mem_to_cells]

            # update the concentration in cells (assume environment steady and constant supply of ions)
            # update bicarbonate instead of H+, assuming buffer action holds:

            # calculate the update to K+ in the cell and environment:

            self.cc_mems[self.iK][:], self.cc_env[self.iK][:] = stb.update_Co(self, self.cc_mems[self.iK][:],
                self.cc_env[self.iK][:], f_K2, cells, p)

            # Update the anion (bicarbonate) concentration instead of H+, assuming bicarb buffer holds:
            self.cc_mems[self.iM][:], self.cc_env[self.iM][:] = stb.update_Co(self, self.cc_mems[self.iM][:],
                self.cc_env[self.iM][:], -f_H2, cells, p)


            # Calculate the new pH and H+ concentrations:
            # run the bicarbonate buffer to ensure realistic concentrations and pH in cell and environment:
            self.cc_mems[self.iH], _, self.cc_mems[
                self.iM], self.pH_cell = stb.bicarbonate_buffer(
                self.cc_mems[self.iH],
                self.cHM_mems,
                self.cc_mems[self.iM],
                p)

            self.cc_env[self.iH], _, self.cc_env[self.iM], self.pH_env = stb.bicarbonate_buffer(
                self.cc_env[self.iH],
                self.cHM_env,
                self.cc_env[self.iM],
                p)

            # update concentrations intracellularly:
            self.cc_mems[self.iH][:], self.cc_cells[self.iH][:], _ = \
                stb.update_intra(self, cells, self.cc_mems[self.iH][:],
                    self.cc_cells[self.iH][:],
                    self.D_free[self.iH],
                    self.zs[self.iH], p)

            # update concentrations intracellularly:
            self.cc_mems[self.iM][:], self.cc_cells[self.iM][:], _ = \
                stb.update_intra(self, cells, self.cc_mems[self.iM][:],
                    self.cc_cells[self.iM][:],
                    self.D_free[self.iM],
                    self.zs[self.iM], p)

            # recalculate the net, unbalanced charge and voltage in each cell:
            self.update_V(cells,p)

        if p.VATPase_dyn == 1 and p.run_sim is True:  # if there's a V-ATPase pump

            if p.sim_ECM is True:

                # if HKATPase pump is desired, run the H-K-ATPase pump:
                f_H3 = stb.pumpVATP(self.cc_mems[self.iH][cells.mem_to_cells], self.cc_env[self.iH][cells.map_mem2ecm],
                    self.vm, self.T, p, self.VATP_block)

                # modify flux by any uneven redistribution of pump location:
                f_H3 = self.rho_pump * f_H3

                self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + f_H3

            else:

                # if HKATPase pump is desired, run the H-K-ATPase pump:
                f_H3 = stb.pumpVATP(self.cc_mems[self.iH],self.cc_env[self.iH],self.vm,self.T,p,self.VATP_block)
                self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + self.rho_pump * f_H3[cells.mem_to_cells]

            # Update the anion (bicarbonate) concentration instead of H+, assuming bicarb buffer holds:
            self.cc_mems[self.iM][:], self.cc_env[self.iM][:] = stb.update_Co(self, self.cc_mems[self.iM][:],
                self.cc_env[self.iM][:], -f_H3, cells, p)

            # Calculate the new pH and H+ concentration:
             # run the bicarbonate buffer to ensure realistic concentrations and pH in cell and environment:

            self.cc_mems[self.iH], _, self.cc_mems[self.iM], self.pH_cell = stb.bicarbonate_buffer(
                 self.cc_mems[self.iH],
                 self.cHM_mems,
                 self.cc_mems[self.iM],
                 p)

            self.cc_env[self.iH], _, self.cc_env[self.iM], self.pH_env = stb.bicarbonate_buffer(
                 self.cc_env[self.iH],
                 self.cHM_env,
                 self.cc_env[self.iM],
                 p)


            # recalculate the net, unbalanced charge and voltage in each cell:
            self.update_V(cells,p)

    def ca_handler(self,cells,p):

        if p.sim_ECM is True:

            # run Ca-ATPase

            f_CaATP = stb.pumpCaATP(self.cc_mems[self.iCa][cells.mem_to_cells],
                self.cc_env[self.iCa][cells.map_mem2ecm],
                self.vm, self.T, p)

            f_CaATP = self.rho_pump * f_CaATP

            # add Ca++ flux to storage:
            self.fluxes_mem[self.iCa] = f_CaATP

        else:

            # run Ca-ATPase

            f_CaATP = stb.pumpCaATP(self.cc_mems[self.iCa], self.cc_env[self.iCa], self.vm, self.T, p)

            # store the transmembrane flux for this ion
            self.fluxes_mem[self.iCa] = self.rho_pump*f_CaATP[cells.mem_to_cells]

        # update calcium concentrations in cell and ecm:
        self.cc_mems[self.iCa][:], self.cc_env[self.iCa][:] = stb.update_Co(self, self.cc_mems[self.iCa][:],
            self.cc_env[self.iCa][:], f_CaATP, cells, p)

        # update concentrations intracellularly:
        self.cc_mems[self.iCa][:], self.cc_cells[self.iCa][:], _ = \
            stb.update_intra(self, cells, self.cc_mems[self.iCa][:],
                self.cc_cells[self.iCa][:],
                self.D_free[self.iCa],
                self.zs[self.iCa], p)

        # recalculate the net, unbalanced charge and voltage in each cell:
        self.update_V(cells, p)

        if p.Ca_dyn == 1:  # do endoplasmic reticulum handling  # FIXME not right!


            # run Ca-ATPase pump to the endoplasmic reticulum:

            f_Ca_ER_pump = stb.pumpCaER(
                self.cc_er[0],
                self.cc_mems[self.iCa],
                self.v_er,
                self.T,
                p,
            )


            # electrodiffusion of ions between cell and endoplasmic reticulum
            f_Ca_ER = \
                stb.electroflux(self.cc_mems[self.iCa], self.cc_er[0], self.Dm_er[0], p.tm, self.z_er[0],
                    self.v_er, self.T, p)

            # Electrodiffusion of charge compensation anion
            f_M_ER = \
                stb.electroflux(self.cc_mems[self.iM], self.cc_er[1], self.Dm_er[1], p.tm, self.z_er[1],
                    self.v_er, self.T, p)

            f_calcium = f_Ca_ER + f_Ca_ER_pump

            # update calcium concentrations in the ER and cell:
            self.cc_er[0] = self.cc_er[0] + f_calcium * ((cells.cell_sa) / (p.ER_vol * cells.cell_vol)) * p.dt

            self.cc_mems[self.iCa] = self.cc_mems[self.iCa] - f_calcium * (
                cells.cell_sa / cells.cell_vol) * p.dt

            # update anion concentration in the ER and cell:
            self.cc_er[1] = self.cc_er[1] + f_M_ER * ((cells.cell_sa) / (p.ER_vol * cells.cell_vol)) * p.dt
            self.cc_mems[self.iM] = self.cc_mems[self.iM] - f_M_ER * (cells.cell_sa / cells.cell_vol) * p.dt

            # recalculate the net, unbalanced charge and voltage in each cell:
            self.update_V(cells, p)

            # get the charge in the endoplasmic reticulum
            q_er = stb.get_charge(self.cc_er, self.z_array_er, p.ER_vol * cells.cell_vol, p)

            # get the voltage across the endoplasmic reticulum
            self.v_er = stb.get_volt(q_er, p.ER_sa * cells.cell_sa, p)

    def update_gj(self,cells,p,t,i):

        # calculate voltage difference (gradient*len_gj) between gj-connected cells:

        self.vgj = self.v_cell[cells.nn_i]- self.v_cell[cells.mem_i]

        if p.v_sensitive_gj is True:
            # determine the open state of gap junctions:
            gmin = self.gj_funk.gmin

            alpha_gj = self.gj_funk.alpha_gj
            beta_gj = self.gj_funk.beta_gj_p

            vgj = np.abs(self.vgj)

            dgjopen_dt = (1 - self.gjopen)*alpha_gj(vgj) - (self.gjopen - gmin)*beta_gj(vgj)

            self.gjopen = self.gjopen + 1e3*dgjopen_dt*p.dt

            # threshold to ensure 0 to 1 status
            inds_gj_over = (self.gjopen > 1.0).nonzero()
            self.gjopen[inds_gj_over] = 1.0

            inds_gj_under = (self.gjopen < 0.0).nonzero()
            self.gjopen[inds_gj_under] = 0.0

            self.gjopen = self.gj_block*self.gjopen

        else:
            self.gjopen = self.gj_block*np.ones(len(cells.mem_i))


        # voltage gradient:
        grad_vgj = self.vgj/(cells.gj_len)


        # concentration gradient for ion i:
        conc_mem = self.cc_mems[i]
        grad_cgj = (conc_mem[cells.nn_i] - conc_mem[cells.mem_i])/(cells.gj_len)

        # midpoint concentration:
        c = (conc_mem[cells.nn_i] + conc_mem[cells.mem_i])/2

        # electroosmotic fluid velocity at gap junctions:
        if p.fluid_flow is True:
            ux = self.u_gj_x
            uy = self.u_gj_y

            # get component of fluid tangent to gap junctions
            ugj = ux*cells.mem_vects_flat[:,2] + uy*cells.mem_vects_flat[:,3]

        else:
            ugj = 0

        # calculate nernst-planck flux tangent to gap junctions:
        fgj = stb.nernst_planck_vector(c,grad_cgj,grad_vgj, ugj,
            p.gj_surface*self.gjopen*self.D_gj[i],self.zs[i],self.T,p)

        delta_cc = (-fgj*cells.mem_sa)/cells.mem_vol

        self.cc_mems[i] = self.cc_mems[i] + p.dt * delta_cc

        # # update concentrations intracellularly:
        # self.cc_mems[i][:], self.cc_cells[i][:], _ = \
        #     stb.update_intra(self, cells, self.cc_mems[i][:],
        #         self.cc_cells[i][:],
        #         self.D_free[i],
        #         self.zs[i], p)

        # use finite volume method to integrate each region:
        # values at centroid mids:
        # c_at_mids = (cc_memso + self.cc_cells[i][cells.mem_to_cells]) / 2
        #
        # # finite volume integral of membrane pie-box values:
        # self.cc_mems[i] = np.dot(cells.M_int_mems, cc_memso) + (1 / 2) * c_at_mids
        # self.cc_cells[i] = (1 / 2) * self.cc_cells[i] + np.dot(cells.M_sum_mems, c_at_mids) / (2 * cells.num_mems)

        self.fluxes_gj[i] = fgj  # store gap junction flux for this ion

    def update_ecm(self,cells,p,t,i):

        if p.closed_bound is True:
            btag = 'closed'

        else:
            btag = 'open'

        # make v_env and cc_env into 2d matrices
        cenv = self.cc_env[i][:]

        v_env = self.v_env[:].reshape(cells.X.shape)

        # enforce voltage at boundary:
        v_env[:,0] = self.bound_V['L']
        v_env[:,-1] = self.bound_V['R']
        v_env[0,:] = self.bound_V['B']
        v_env[-1,:] = self.bound_V['T']

        cenv = cenv.reshape(cells.X.shape)

        # prepare concentrations and diffusion constants for MACs grid format
        # by resampling the values at the u v coordinates of the flux:
        cenv_x = np.zeros(cells.grid_obj.u_shape)
        cenv_y = np.zeros(cells.grid_obj.v_shape)

        # create the proper shape for the concentrations and state appropriate boundary conditions::
        cenv_x[:,1:] = cenv

        cenv_y[1:,:] = cenv

        if p.closed_bound is True: # insulation boundary conditions

            cenv_x[:,0] = cenv_x[:,1]
            cenv_x[:,-1] = cenv_x[:,-2]
            cenv_x[0,:] = cenv_x[1,:]
            cenv_x[-1,:] = cenv_x[-2,:]

            cenv_y[0,:] = cenv_y[1,:]
            cenv_y[-1,:] = cenv_y[-2,:]
            cenv_y[:,0] = cenv_y[:,1]
            cenv_y[:,-1] = cenv_y[:,-2]

        else:   # open and electrically grounded boundary conditions
            cenv_x[:,0] =  self.c_env_bound[i]
            cenv_x[:,-1] =  self.c_env_bound[i]
            cenv_x[0,:] =  self.c_env_bound[i]
            cenv_x[-1,:] =  self.c_env_bound[i]

            cenv_y[0,:] =  self.c_env_bound[i]
            cenv_y[-1,:] =  self.c_env_bound[i]
            cenv_y[:,0] =  self.c_env_bound[i]
            cenv_y[:,-1] =  self.c_env_bound[i]

            # try adding it to the second layer to get proper diffusion in...

            cenv_x[:,1] =  self.c_env_bound[i]
            cenv_x[:,-2] =  self.c_env_bound[i]
            cenv_x[1,:] =  self.c_env_bound[i]
            cenv_x[-2,:] =  self.c_env_bound[i]

            cenv_y[1,:] =  self.c_env_bound[i]
            cenv_y[-2,:] =  self.c_env_bound[i]
            cenv_y[:,1] =  self.c_env_bound[i]
            cenv_y[:,-2] =  self.c_env_bound[i]

        # calculate gradients in the environment
        grad_V_env_x, grad_V_env_y = cells.grid_obj.grid_gradient(v_env,bounds='closed')

        grad_cc_env_x, grad_cc_env_y = cells.grid_obj.grid_gradient(cenv,bounds=btag)

        # calculate fluxes for electrodiffusive transport:
        if p.fluid_flow is True:

            uenvx = np.zeros(cells.grid_obj.u_shape)
            uenvy = np.zeros(cells.grid_obj.v_shape)

            uenvx[:,1:] = self.u_env_x
            uenvy[1:,:] = self.u_env_y

            if p.closed_bound is False:

                uenvx[:,0] = uenvx[:,1]
                uenvx[:,-1]= uenvx[:,-2]
                uenvx[0,:] = uenvx[1,:]
                uenvx[-1,:] = uenvx[-2,:]

                uenvy[:,0] = uenvy[:,1]
                uenvy[:,-1]= uenvy[:,-2]
                uenvy[0,:] = uenvy[1,:]
                uenvy[-1,:] = uenvy[-2,:]

            else:

                uenvx[:,0] = 0
                uenvx[:,-1]= 0
                uenvx[0,:] = 0
                uenvx[-1,:] = 0

                uenvy[:,0] = 0
                uenvy[:,-1]= 0
                uenvy[0,:] = 0
                uenvy[-1,:] = 0


        else:
            uenvx = 0
            uenvy = 0

        # distance over which electric field is active is affected by electrolyte screening. Therefore
        # we define a field modulation in the bulk, which is assumed to be a fraction of the Debye length,
        # which is assumed to be about 1 nm:

        field_mod = (1e-9/cells.delta)*cells.mems_per_envSquare.mean()

        f_env_x, f_env_y = stb.np_flux_special(cenv_x,cenv_y,grad_cc_env_x,grad_cc_env_y,
            field_mod*grad_V_env_x, field_mod*grad_V_env_y, uenvx,uenvy,self.D_env_u[i],self.D_env_v[i],
            self.zs[i],self.T,p)

        # slow fluxes, if desired by user:
        f_env_x = p.env_delay_const*f_env_x
        f_env_y = p.env_delay_const*f_env_y

        if p.smooth_level > 0.0:

            f_env_x = gaussian_filter(f_env_x,p.smooth_level)  # smooth out the flux terms  #FIXME might not need these later
            f_env_y = gaussian_filter(f_env_y, p.smooth_level)  # smooth out the flux terms

        if p.closed_bound is False:

            f_env_x[:,0] = f_env_x[:,1]
            f_env_x[:,-1]= f_env_x[:,-2]
            f_env_x[0,:] = f_env_x[1,:]
            f_env_x[-1,:] = f_env_x[-2,:]

            f_env_y[:,0] = f_env_y[:,1]
            f_env_y[:,-1]= f_env_y[:,-2]
            f_env_y[0,:] = f_env_y[1,:]
            f_env_y[-1,:] = f_env_y[-2,:]

        else:

            f_env_x[:,0] = 0
            f_env_x[:,-1]= 0
            f_env_x[0,:] = 0
            f_env_x[-1,:] = 0

            f_env_y[:,0] = 0
            f_env_y[:,-1]= 0
            f_env_y[0,:] = 0
            f_env_y[-1,:] = 0

        # calculate the divergence of the total flux, which is equivalent to the total change per unit time

        d_fenvx = -(f_env_x[:,1:] - f_env_x[:,0:-1])/cells.delta
        d_fenvy = -(f_env_y[1:,:] - f_env_y[0:-1,:])/cells.delta

        delta_c = d_fenvx + d_fenvy

        #-----------------------
        cenv = cenv + delta_c*p.dt

        # finite volume integrator:
        # cenv = np.dot(cells.gridInt, cenv.ravel())
        # cenv = cenv.reshape(cells.X.shape)

        if p.closed_bound is True:
            # Neumann boundary condition (flux at boundary)
            # zero flux boundaries for concentration:

            cenv[:,-1] = cenv[:,-2]
            cenv[:,0] = cenv[:,1]
            cenv[0,:] = cenv[1,:]
            cenv[-1,:] = cenv[-2,:]

        elif p.closed_bound is False:
            # if the boundary is open, set the concentration at the boundary
            # open boundary
            cenv[:,-1] = self.c_env_bound[i]
            cenv[:,0] = self.c_env_bound[i]
            cenv[0,:] = self.c_env_bound[i]
            cenv[-1,:] = self.c_env_bound[i]

            cenv[:,-2] = self.c_env_bound[i]
            cenv[:,1] = self.c_env_bound[i]
            cenv[1,:] = self.c_env_bound[i]
            cenv[-2,:] = self.c_env_bound[i]

        self.cc_env[i] = cenv.ravel()

        fenvx = (f_env_x[:,1:] + f_env_x[:,0:-1])/2
        fenvy = (f_env_y[1:,:] + f_env_y[0:-1,:])/2

        self.fluxes_env_x[i] = fenvx.ravel()  # store ecm junction flux for this ion
        self.fluxes_env_y[i] = fenvy.ravel()  # store ecm junction flux for this ion

    def get_Efield(self,cells,p):

         # calculate voltage difference (gradient*len_gj) between gj-connected cells:
        if p.sim_ECM is True:

            vmem = self.v_cell

            self.Egj = - (vmem[cells.nn_i]- vmem[cells.mem_i])/cells.gj_len

            # in the environment:
            venv = self.v_env.reshape(cells.X.shape)
            genv_x, genv_y = fd.gradient(venv, cells.delta)

            self.E_env_x = -genv_x.ravel()*cells.ave2ecmV
            self.E_env_y = -genv_y.ravel()*cells.ave2ecmV

            if p.smooth_level > 0.0:

                self.E_env_x = gaussian_filter(self.E_env_x.reshape(cells.X.shape),p.smooth_level)
                self.E_env_y = gaussian_filter(self.E_env_y.reshape(cells.X.shape),p.smooth_level)

        else:

            self.Egj = - (self.vm[cells.nn_i] - self.vm[cells.mem_i])/cells.gj_len

        # self.Egj = self.Egj*cells.ave2cellV # scale from cell space volume to whole cell volume

        # get x and y components of the electric field:
        self.E_gj_x = cells.mem_vects_flat[:,2]*self.Egj
        self.E_gj_y = cells.mem_vects_flat[:,3]*self.Egj

    def get_ion(self,label):
        """
        Given a string input, returns the simulation index of the appropriate ion.

        """

        if label == 'Na':

            ion = self.iNa

        elif label == 'K':

            ion = self.iK

        elif label == 'Ca':

            ion = self.iCa

        elif label == 'Cl':

            ion = self.iCl

        else:
            logs.warning('Oops! Morphogen gated ion channel target not found!')
            ion = []

        return ion

    def initDenv(self,cells,p):

        self.D_env_u = np.zeros((self.D_env.shape[0],cells.grid_obj.u_shape[0],cells.grid_obj.u_shape[1]))
        self.D_env_v = np.zeros((self.D_env.shape[0],cells.grid_obj.v_shape[0],cells.grid_obj.v_shape[1]))

        for i, dmat in enumerate(self.D_env):

            if p.env_type is False: # if air surrounds, first set everything to zero and add in cluster data...
                self.D_env[i][:] = 0
            # for all cells and mems in the cluster, set the internal diffusion constant for adherens junctions:
            dummyMems = np.ones(len(cells.mem_i))*self.D_free[i]*p.D_adh

            # get a list of all membranes for boundary cells:
            neigh_to_bcells,_,_ = tb.flatten(cells.cell_nn[cells.bflags_cells])

            all_bound_mem_inds = cells.cell_to_mems[cells.bflags_cells]
            interior_bound_mem_inds = cells.cell_to_mems[neigh_to_bcells]
            interior_bound_mem_inds,_,_ = tb.flatten(interior_bound_mem_inds)
            all_bound_mem_inds, _ ,_ = tb.flatten(all_bound_mem_inds)

            # set external membrane of boundary cells to the diffusion constant of tight junctions:
            dummyMems[all_bound_mem_inds] = self.D_free[i]*p.D_tj*self.Dtj_rel[i]
            dummyMems[interior_bound_mem_inds] = self.D_free[i]*p.D_tj*self.Dtj_rel[i]
            dummyMems[cells.bflags_mems] = self.D_free[i]

            # interp the membrane data to an ecm grid, fill values correspond to environmental diffusion consts:
            if p.env_type is True:
                Denv_o = interp.griddata((cells.mem_vects_flat[:,0],cells.mem_vects_flat[:,1]),dummyMems,
                    (cells.X,cells.Y),method='nearest',fill_value=self.D_free[i])

            else:
                Denv_o = interp.griddata((cells.mem_vects_flat[:,0],cells.mem_vects_flat[:,1]),dummyMems,
                    (cells.X,cells.Y),method='nearest',fill_value=0)

            Denv_o = Denv_o.ravel()
            Denv_o[cells.inds_env] = self.D_free[i]

            # create an ecm diffusion grid filled with the environmental values
            self.D_env[i] = Denv_o

        # create a matrix that weights the relative transport efficiency in the world space:
        D_env_weight = self.D_env[self.iP]/self.D_env[self.iP].max()
        self.D_env_weight = D_env_weight.reshape(cells.X.shape)
        self.D_env_weight_base = np.copy(self.D_env_weight)

        for i, dmat in enumerate(self.D_env):

            if p.env_type is True:

                self.D_env_u[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                    (cells.grid_obj.u_X,cells.grid_obj.u_Y),method='nearest',fill_value = self.D_free[i])

                self.D_env_v[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                    (cells.grid_obj.v_X,cells.grid_obj.v_Y),method='nearest',fill_value=self.D_free[i])

            else:
                self.D_env_u[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                    (cells.grid_obj.u_X,cells.grid_obj.u_Y),method='nearest',fill_value = 0)

                self.D_env_v[i] = interp.griddata((cells.xypts[:,0],cells.xypts[:,1]),dmat.ravel(),
                    (cells.grid_obj.v_X,cells.grid_obj.v_Y),method='nearest',fill_value = 0)

        self.D_env_weight_u = self.D_env_u[self.iP]/self.D_env_u[self.iP].max()

        self.D_env_weight_v = self.D_env_v[self.iP]/self.D_env_v[self.iP].max()

        if p.closed_bound is True:  # set full no slip boundary condition at exterior bounds

            self.D_env_weight_u[:,0] = 0
            self.D_env_weight_u[:,-1] = 0
            self.D_env_weight_u[0,:] = 0
            self.D_env_weight_u[-1,:] = 0

            self.D_env_weight_v[:,0] = 0
            self.D_env_weight_v[:,-1] = 0
            self.D_env_weight_v[0,:] = 0
            self.D_env_weight_v[-1,:] = 0


    # ..................{ PLOTTERS                           }..................

    def _plot_loop(self, cells: 'Cells', p: 'Parameters') -> (np.ndarray, set):
        '''
        Display and/or save an animation during solving if requested _and_
        calculate data common to solving both with and without extracellular
        spaces.

        Returns
        --------
        tt : np.ndarray
            Time-steps vector appropriate for the current run.
        tsamples : set
            Time-steps vector resampled to save data at substantially fewer
            times. The length of this vector governs the number of frames in
            plotted animations, for example.
        '''


        # Human-readable type and maximum number of steps of the current run.
        if p.run_sim is False:
            figure_type_label = 'Initializing'
            loop_type_label = 'initialization'
            loop_time_step_max = p.init_tsteps
        else:
            figure_type_label = 'Simulating'
            loop_type_label = 'simulation'
            loop_time_step_max = p.sim_tsteps

        # Maximum number of seconds simulated by the current run.
        loop_seconds_max = loop_time_step_max * p.dt

        # Time-steps vector appropriate for the current run.
        tt = np.linspace(0, loop_seconds_max, loop_time_step_max)

        #FIXME: Refactor into a for loop calling the range() builtin. Sunsets!
        # Resample this vector to save data at substantially fewer times.
        tsamples = set()
        i = 0
        while i < len(tt) - p.t_resample:
            i += p.t_resample
            tsamples.add(tt[i])

        # Log this run.
        logs.log_info(
            'Your {} is running from 0 to {:.1f} seconds of in-world time '
            'in {} time-steps ({} sampled).'.format(
            loop_type_label,
            loop_seconds_max,
            len(tt),
            len(tsamples),
        ))

        # If displaying and/or saving an animation during solving, do so.
        if p.plot_while_solving:
            self._anim_cells_while_solving = AnimCellsWhileSolving(
                sim=self, cells=cells, p=p,
                type='Vmem',
                figure_title='Vmem while {}'.format(
                    figure_type_label),
                colorbar_title='Voltage [mV]',
                is_color_autoscaled=p.autoscale_Vmem,
                color_min=p.Vmem_min_clr,
                color_max=p.Vmem_max_clr,
            )
        # Else, nullify the object encapsulating this animation for safety.
        else:
            self._anim_cells_while_solving = None

        return tt, tsamples

    def _replot_loop(self, p: 'Parameters') -> None:
        '''
        Update the currently displayed and/or saved animation during solving
        with the results of the most recently solved time step, if requested.
        '''

        # Update this animation for the "last frame" if desired, corresponding
        # to the results of the most recently solved time step.
        if p.plot_while_solving:
            self._anim_cells_while_solving.plot_frame(frame_number=-1)

    def _deplot_loop(self) -> None:
        '''
        Explicitly close the previously displayed and/or saved animation if
        any _or_ noop otherwise.

        To conserve memory, this method nullifies and hence garbage collects
        both this animation and this animation's associated figure.
        '''

        if self._anim_cells_while_solving is not None:
            self._anim_cells_while_solving.close()
            self._anim_cells_while_solving = None



#-----------------------------------------------------------------------------------------------------------------------
# WASTELANDS
#-----------------------------------------------------------------------------------------------------------------------
#FIXME: I don't quite grok our usage of "sim.run_sim". This undocumented
#attribute appears to be internally set by the Simulator.run_phase_sans_ecm()
#method. That makes sense; however, what's the parallel "p.run_sim" attribute
#for, then?  Interestingly, the "SimRunner" class sets "p.run_sim" as follows:
#
#* To "False" if an initialization is being performed.
#* To "True" if a simulation is being performed.
#
#This doesn't seem quite ideal, however. Ideally, there would exist one and only
#one attribute whose value is an instance of a multi-state "PhaseEnum" class
#rather than two binary boolean attributes. Possible enum values might include:
#
#* "PhaseEnum.seed" when seeding a new cluster.
#* "PhaseEnum.init" when initializing a seeded cluster.
#* "PhaseEnum.sim" when simulation an initialized cluster.
#
#This attribute would probably exist in the "Simulator" class -- say, as
#"sim.phase". In light of that, consider the following refactoring:
#
#* Define a new "PhaseEnum" class in the "sim" module with the above attributes.
#* Define a new "Simulator.phase" attribute initialized to None.
#* Replace all existing uses of the "p.run_sim" and "sim.run_sim" booleans with
#  "sim.phase" instead. Note that only the:
#  * "SimRunner" class sets "p.run_sim".
#  * "Simulator" class sets "sim.run_sim".
#
#Note also the "plot_type" parameter passed to the pipeline.plot_all() function
#*AND* seemingly duplicate "p.plot_type" attribute, which should probably
#receive similar treatment. Wonder temptress at the speed of light and the
#sound of love!


# def update_C(self,ion_i,flux,cells,p):
#
#     c_cells = self.cc_cells[ion_i][:]
#     c_env = self.cc_env[ion_i][:]
#
#     if p.sim_ECM is True:
#
#         delta_cells = flux*(cells.mem_sa/cells.mem_vol)
#
#         # interpolate the flux from mem points to the ENV GRID:
#         flux_env = -flux[cells.map_mem2ecm]
#
#         # save values at the cluster boundary:
#         bound_vals = flux_env[cells.ecm_bound_k]
#
#         # set the values of the global environment to zero:
#         flux_env[cells.inds_env] = 0
#
#         # finally, ensure that the boundary values are restored:
#         flux_env[cells.ecm_bound_k] = bound_vals
#
#         # Now that we have a nice, neat interpolation of flux from cell membranes, multiply by the cell2ecm ratio,
#         # which yields number of cells per unit env grid square, and then by cell surface area, and finally
#         # by the volume of the env grid square, to get the mol/s change in concentration (divergence):
#
#         delta_env = (flux_env * cells.memSa_per_envSquare) / cells.true_ecm_vol
#
#         # update the concentrations
#         self.cc_cells[ion_i] = c_cells + delta_cells*p.dt
#         self.cc_env[ion_i] = c_env + delta_env*p.dt
#
#
#     else:
#
#         delta_cells = flux* (cells.mem_sa/cells.mem_vol)
#         delta_env = -flux*(cells.mem_sa/p.vol_env)
#
#         self.cc_cells[ion_i] = c_cells + delta_cells * p.dt
#
#         c_env = c_env + delta_env * p.dt
#         # assume auto-mixing of environmental concentrations:
#         self.cc_env[ion_i][:] = c_env.mean()




