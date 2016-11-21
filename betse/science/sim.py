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
from betse.exceptions import BetseSimInstabilityException
from betse.util.io.log import logs
from betse.science import filehandling as fh
from betse.science import finitediff as fd
from betse.science import toolbox as tb
from betse.science import sim_toolbox as stb
from betse.science.chemistry.networks import MasterOfNetworks
from betse.science.chemistry.metabolism import  MasterOfMetabolism
from betse.science.chemistry.gene import MasterOfGenes
from betse.science.organelles.endo_retic import EndoRetic
from betse.science.physics.ion_current import get_current
from betse.science.physics.flow import getFlow
from betse.science.physics.deform import (
    getDeformation, timeDeform, implement_deform_timestep)
from betse.science.physics.move_channels import eosmosis
from betse.science.physics.pressures import electro_F, getHydroF, osmotic_P
from betse.science.tissue.channels.gap_junction import Gap_Junction
from betse.science.tissue.handler import TissueHandler
from betse.science.visual.anim.anim import AnimCellsWhileSolving


# FIXME update plotting

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

    Attributes (Counts)
    ----------
    cdl : int
        Number of cells in this simulated cluster.
    mdl : int
        Number of cell membranes in this simulated cluster.

    Attributes (Voltage)
    ----------
    vcell_time : ndarray
        Voltage at the inner membrane surface of each cell as a function of
        time.
    venv_time : ndarray
        Voltage at the outer membrane surface of each cell as a function of
        time.
    vm : ndarray
        One-dimensional Numpy array of length the number of cell membranes such
        that each element is the transmembrane voltage spatially situated
        across the cell membrane indexed by that element for the current time
        step.
    vm_time : ndarray
        Two-dimensional Numpy array of the transmembrane voltage across all
        cell membranes, whose:
        . First dimension indexes each simulation time step.
        . Second dimension indexes each cell membrane in this cluster such that
          each element is the transmembrane voltage spatially situated across
          the cell membrane indexed by that element for the current time step.
        Equivalently, this array is the concatenation of all :attr:`vm` arrays
        over all time steps.

    Attributes (Current Density)
    ----------
    I_cell_x_time : list
        Two-dimensional list whose:
        * First dimension indexes each simulation time step.
        * Second dimension indexes cells, whose length is the number of cells
          and each element is the X component of the intracellular current
          density vector spatially situated at the center of each cell for this
          time step.
    I_cell_y_time : list
        Two-dimensional list whose:
        * First dimension indexes each simulation time step.
        * Second dimension indexes cells, whose length is the number of cells
          and each element is the Y component of the intracellular current
          density vector spatially situated at the center of each cell for this
          time step.
    I_tot_x_time : list
        Two-dimensional list whose:
        * First dimension indexes each simulation time step.
        * Second dimension indexes square grid spaces, whose length is the
          number of grid spaces in either dimension and each element is the X
          component of the **total current density vector** (i.e., vector of
          both intra- _and_ extracellular current densities) spatially situated
          at the center of each grid space for this time step.
    I_tot_y_time : list
        Two-dimensional list whose:
        * First dimension indexes each simulation time step.
        * Second dimension indexes square grid spaces, whose length is the
          number of grid spaces in either dimension and each element is the Y
          component of the **total current density vector** (i.e., vector of
          both intra- _and_ extracellular current densities) spatially situated
          at the center of each grid space for this time step.

    Attributes (Electric Field)
    ----------
    E_gj_x : ndarray
        One-dimensional Numpy array indexing each simulated cell membrane, such
        that each element is the X component of the intracellular electric
        field vector spatially situated across the gap junction to which the
        current membrane connects: specifically, the difference of the X
        component of the voltage situated at this membrane with that of the
        voltage situated at the gap junction-connected membrane adjacent to
        this membrane, divided by the length in meters of this gap junction.
    E_gj_y : ndarray
        One-dimensional Numpy array indexing each simulated cell membrane, such
        that each element is the Y component of the intracellular electric
        field vector defined as for the corresponding :attr:`E_gj_X` array.
    efield_gj_x_time : list
        Two-dimensional list whose:
        * First dimension indexes each simulation time step.
        * Second dimension indexes each simulated cell membrane, such that each
          element is the X component of the intracellular electric field vector
          defined as for the corresponding :attr:`E_gj_x` array.
        Equivalently, this array is the concatenation of all :attr:`E_gj_x`
        arrays for all time steps.
    efield_gj_y_time : list
        Two-dimensional list whose:
        * First dimension indexes each simulation time step.
        * Second dimension indexes each simulated cell membrane, such that each
          element is the Y component of the intracellular electric field vector
          defined as for the corresponding :attr:`E_gj_y` array.
        Equivalently, this array is the concatenation of all :attr:`E_gj_y`
        arrays for all time steps.
    '''

    def __init__(self, p):

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

        self.ignore_ecm = False

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

        # initialize all extra substances related objects to None, to be filled in if desired later
        self.molecules = None
        self.metabo = None
        self.met_concs = None
        self.grn = None

        #FIXME: Consider shifting these attributes from this class into the
        #more appropriate "Cells" class.
        self.mdl = len(cells.mem_i)  # mems-data-length
        self.cdl = len(cells.cell_i)  # cells-data-length

        if p.sim_ECM is True:  # set environnment data length
            self.edl = len(cells.xypts)
        else:
            self.edl = self.mdl

        # initialize two extra arrays that allow additional substances (defined outside of sim) to affect Vmem:
        self.extra_rho_cells = np.zeros(self.cdl)
        self.extra_rho_env = np.zeros(self.edl)

        self.extra_rho_cells_o = np.zeros(self.cdl)

        self.vgj = np.zeros(self.mdl)

        self.gj_block = 1 # will update this according to user preferences in self.init_tissue()

        # Identity matrix to easily make matrices out of scalars
        self.id_mems = np.ones(self.mdl)

        self.cc_cells = []  # cell concentrations at cell centres
        self.cc_er = []  # endoplasmic reticulum ion concentrations in each cell
        self.cc_env = []  # environmental concentrations initialized

        self.zs = []  # ion valence state initialized
        self.z_er = []  # ion valence states of er ions
        self.z_array_er = []
        self.z_array = []  # ion valence array matched to cell points

        self.Dm_cells = []  # membrane diffusion constants initialized
        self.Dm_er = []  # a list of endoplasmic reticulum membrane state
        self.D_free = []  # a list of single-valued free diffusion constants for each ion
        self.D_gj = []  # an array of diffusion constants for gap junctions
        self.movingIons = []  # moving ions indices
        self.ionlabel = {}  # dictionary to hold ion label names
        self.molar_mass = []

        self.T = p.T  # set the base temperature for the simulation

        self.Jn = np.zeros(self.mdl)  # net normal current density at membranes
        self.J_weight = np.ones(self.mdl)  # moments of current passing through membranes
        self.rho_mems = np.zeros(self.mdl)  # membrane charge
        self.dvm = np.zeros(self.mdl)   # rate of change of Vmem
        self.drho = np.zeros(self.cdl)  # rate of change of charge in cell
        self.Phi_env = np.zeros(cells.X.shape)
        self.Pol_mem = np.zeros(self.mdl)  # polarization density at the membrane

        # Current averaged to cell centres:
        self.J_cell_x = np.zeros(self.cdl)
        self.J_cell_y = np.zeros(self.cdl)

        # Membrane current data structure initialization
        self.flx_mem_i = np.zeros(self.mdl)
        self.fluxes_mem = []
        self.I_mem = np.zeros(self.mdl)  # total current across membranes

        self.P_cells = np.zeros(self.cdl)  # initialize pressure in cells

        self.v_cell = np.zeros(self.cdl)  # initialize intracellular voltage
        # self.vm = -50.0e-3*np.ones(self.mdl)     # initialize vmem
        self.vm = np.zeros(self.mdl)     # initialize vmem
        self.rho_cells = np.zeros(self.cdl)

        self.E_gj_x = np.zeros(self.mdl)   # electric field components across gap junctions
        self.E_gj_y = np.zeros(self.mdl)


        if p.sim_ECM is True:  # special items specific to simulation of extracellular spaces only:

            # vectors storing separate cell and env voltages
            self.v_env = np.zeros(len(cells.xypts))
            self.rho_env = np.zeros(len(cells.xypts))

            self.z_array_env = []  # ion valence array matched to env points
            self.D_env = []  # an array of diffusion constants for each ion defined on env grid
            self.c_env_bound = []  # moving ion concentration at global boundary
            self.Dtj_rel = []  # relative diffusion constants for ions across tight junctions

            self.E_env_x = np.zeros(self.edl).reshape(cells.X.shape)
            self.E_env_y = np.zeros(self.edl).reshape(cells.X.shape)


        else:  # items specific to simulation *without* extracellular spaces:
            # Initialize environmental volume:
            self.envV = np.zeros(self.mdl)
            self.envV[:] = p.vol_env


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


        ion_names = list(p.ions_dict.keys())

        i = -1  # dynamic index

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
                    str_cells = 'c' + name + 'cells'
                    setattr(self, str_cells, np.zeros(self.cdl))
                    vars(self)[str_cells][:] = p.cell_concs[name]

                    # environmental concentration for the ion:
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

                    setattr(self, str_z, np.zeros(self.cdl))
                    vars(self)[str_z][:] = p.ion_charge[name]

                    if p.sim_ECM:  # ion charge characteristic for extracellular:
                        str_z2 = 'z2' + name

                        setattr(self, str_z2, np.zeros(len(cells.xypts)))
                        vars(self)[str_z2][:] = p.ion_charge[name]

                    self.cc_cells.append(vars(self)[str_cells])
                    self.cc_env.append(vars(self)[str_env])

                    self.zs.append(p.ion_charge[name])
                    self.molar_mass.append(p.molar_mass[name])
                    self.z_array.append(vars(self)[str_z])
                    self.Dm_cells.append(vars(self)[str_Dm])
                    self.D_gj.append(vars(self)[str_Dgj])
                    self.D_free.append(p.free_diff[name])

                    self.fluxes_mem.append(self.flx_mem_i)

                    if p.sim_ECM:
                        self.c_env_bound.append(p.env_concs[name])
                        self.z_array_env.append(vars(self)[str_z2])
                        self.D_env.append(vars(self)[str_Denv])
                        self.Dtj_rel.append(p.Dtj_rel[name])

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

            # self.movingIons.append(self.iH)

            # create concentration arrays of dissolved carbon dioxide (carbonic acid, non-dissociated):
            self.cHM_cells = np.zeros(self.cdl)
            self.cHM_cells[:] = 0.03 * p.CO2

            self.cHM_env = np.zeros(self.edl)
            self.cHM_env[:] = 0.03 * p.CO2

            self.cH_cells, self.pH_cell = stb.bicarbonate_buffer(self.cHM_cells, self.cc_cells[self.iM])
            self.cH_env, self.pH_env = stb.bicarbonate_buffer(self.cHM_env, self.cc_env[self.iM])

            # initialize diffusion constants
            DmH = np.zeros(self.mdl)
            DmH[:] = p.Dm_H

            self.zH = np.zeros(self.cdl)
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
                p.env_concs['H'] = self.cH_env.mean()
                self.c_env_bound.append(p.env_concs['H'])

            # append items to main data vectors:
            self.cc_cells.append(self.cH_cells)
            self.cc_env.append(self.cH_env)

            self.zs.append(p.z_H)
            self.molar_mass.append(p.M_H)
            self.z_array.append(self.zH)

            self.Dm_cells.append(DmH)
            self.D_gj.append(DgjH)
            self.D_free.append(p.Do_H)

            self.fluxes_mem.append(self.flx_mem_i)

            if p.sim_ECM is True:
                self.z_array_env.append(self.zH2)
                self.D_env.append(DenvH)
                self.Dtj_rel.append(p.Dtj_rel['H'])

        # -------------------------------------------------------------------------------------------------------

        # convert all data structures to Numpy arrays:
        self.cc_cells = np.asarray(self.cc_cells)
        self.cc_env = np.asarray(self.cc_env)

        self.zs = np.asarray(self.zs)
        self.z_array = np.asarray(self.z_array)

        self.Dm_cells = np.asarray(self.Dm_cells)

        self.D_free = np.asarray(self.D_free)
        self.D_gj = np.asarray(self.D_gj)
        self.molar_mass = np.asarray(self.molar_mass)

        self.fluxes_mem = np.asarray(self.fluxes_mem)

        if p.sim_ECM is True:  # items specific for extracellular spaces simulation:

            self.z_array_env = np.asarray(self.z_array_env)
            self.D_env = np.asarray(self.D_env)

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

        # GJ fluxes storage vector:
        self.fluxes_gj = np.copy(self.fluxes_mem)

    def init_tissue(self, cells, p):
        '''
        Prepares data structures pertaining to tissue profiles, dynamic
        activity, and optional methods such as electroosmotic fluid,
        which can be changed in between an initialization and simulation
        run.

        This method is called at the start of all simulations.
        '''

        # # load in the gap junction dynamics object:
        if p.v_sensitive_gj:
            self.gj_funk = Gap_Junction(self, cells, p)

        if p.sim_ECM is True:
            #  Initialize diffusion constants for the extracellular transport:
            self.initDenv(cells,p)

            # re-init global boundary fixed concentrations:

            for key, val in p.ions_dict.items():

                if val == 1 and key != 'H':
                    ion_i = self.get_ion(key)
                    # print("resetting c_env from ", self.c_env_bound[ion_i], 'to ', p.cbnd[key], "for ", key)
                    self.c_env_bound[ion_i] = p.cbnd[key]

                elif val ==1 and key == 'H':
                    self.c_env_bound[self.iH] = p.env_concs['H']

        self.dyna = TissueHandler(self, cells, p)   # create the tissue dynamics object
        self.dyna.tissueProfiles(self, cells, p)  # initialize all tissue profiles

        if p.sim_ECM is True:
            # create a copy-base of the environmental junctions diffusion constants:
            self.D_env_base = copy.copy(self.D_env)

            # initialize current vectors
            self.J_env_x = np.zeros(len(cells.xypts))
            self.J_env_y = np.zeros(len(cells.xypts))

            flx_env_ix = np.zeros(self.edl)
            flx_env_iy = np.zeros(self.edl)

            self.fluxes_env_x = [flx_env_ix for n in self.zs]
            self.fluxes_env_y = [flx_env_iy for n in self.zs]

            self.fluxes_env_x = np.asarray(self.fluxes_env_x)
            self.fluxes_env_y = np.asarray(self.fluxes_env_x)

            # self.J_TJ = np.zeros(self.mdl)  # tight junction current density

        # Initialize an array structure that will hold user-scheduled changes to membrane permeabilities:
        Dm_cellsA = np.asarray(self.Dm_cells)

        self.Dm_base = np.copy(Dm_cellsA) # make a copy that will serve as the unaffected values base

        # if tb.emptyDict(p.scheduled_options) is False:
        self.Dm_scheduled = np.copy(Dm_cellsA)
        self.Dm_scheduled[:] = 0

        # Initialize an array structure that will hold dynamic calcium-gated channel changes to mem perms:
        self.Dm_cag = np.copy(Dm_cellsA)
        self.Dm_cag[:] = 0

        self.Dm_stretch = np.copy(Dm_cellsA)   # array for stretch activated ion channels...
        self.Dm_stretch[:] = 0

        self.Dm_morpho = np.copy(Dm_cellsA)
        self.Dm_morpho[:] = 0

        self.Dm_custom = np.copy(Dm_cellsA)  # array for customized ion channels
        self.Dm_custom[:] = 0

        self.P_mod = np.copy(self.P_cells[:])
        self.P_base = np.copy(self.P_cells[:])

        self.div_u_osmo = np.zeros(self.cdl)

        # Noise Initialization ---------------------
        # add channel noise to the model:
        self.channel_noise_factor = np.random.random(self.mdl)

        self.Dm_cells[self.iK] = (p.channel_noise_level * self.channel_noise_factor + 1) * self.Dm_cells[self.iK]


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

        # initialize additional pump blocks:
        self.CaATP_block = np.ones(self.mdl)  # initialize CaATP blocking vector
        self.NaCaExch_block = np.ones(self.mdl)  # initialize CaATP blocking vector

        # initialize calcium dynamics if desired:
        if p.ions_dict['Ca'] == 1 and p.Ca_dyn is True:
            self.endo_retic = EndoRetic(self, cells, p)

        else:
            self.endo_retic = None


        # -----auxiliary molecules initialization -------------------------

        # create and initialize the auxiliary-molecules handler for this simulation:
        #(only do these initializations if they haven't been done yet)
        if p.molecules_enabled and self.molecules is None:

            # self.molecules = MasterOfMolecules(self, cells, p.molecules_config,p)
            self.molecules = MasterOfNetworks(self, cells, p.molecules_config, p,
                mit_enabled=p.mol_mit_enabled)

            if p.mol_mit_enabled is True:
                self.molecules.mit_enabled = True

            if p.reactions_enabled:

                self.molecules.read_reactions(p.reactions_config, self, cells, p)
                self.molecules.write_reactions()

                if self.molecules.mit_enabled is True:
                    self.molecules.write_reactions_mit()
                    self.molecules.create_reaction_matrix_mit()

            self.molecules.create_reaction_matrix()

            if p.transporters_enabled:

                self.molecules.read_transporters(p.transporters_config, self, cells, p)
                self.molecules.write_transporters(self, cells,p)

            if p.channels_enabled:

                self.molecules.read_channels(p.channels_config, self, cells, p)

            if p.modulators_enabled:

                self.molecules.read_modulators(p.modulators_config, self, cells, p)


        elif p.molecules_enabled and self.molecules is not None:
        # don't declare a whole new object, but re-read in parts that user may have changed:
            self.molecules.tissue_init(self, cells, p.molecules_config,p)

            if p.reactions_enabled:
                self.molecules.read_reactions(p.reactions_config, self, cells, p)
                self.molecules.write_reactions()

            self.molecules.create_reaction_matrix()

            if p.transporters_enabled:
                self.molecules.read_transporters(p.transporters_config, self, cells, p)
                self.molecules.write_transporters(self, cells, p)

            if p.channels_enabled:
                self.molecules.read_channels(p.channels_config, self, cells, p)

            if p.modulators_enabled:
                self.molecules.read_modulators(p.modulators_config, self, cells, p)

        #-----metabolism initialization -----------------------------------
        if p.metabolism_enabled and self.metabo is None:

            logs.log_info("Initializing metabolism...")

            # create an instance of the metabolism simulator
            self.metabo = MasterOfMetabolism(p)
            # read in the configuration settings for the metabolism simulator:
            self.metabo.read_metabo_config(self, cells, p)

            # create a dictionary pointing to key metabolic molecules used in sim: ATP, ADP and Pi:
            self.met_concs = {'cATP': self.metabo.core.cell_concs['ATP'][cells.mem_to_cells],
                              'cADP': self.metabo.core.cell_concs['ADP'][cells.mem_to_cells],
                              'cPi': self.metabo.core.cell_concs['Pi'][cells.mem_to_cells]}

        elif p.metabolism_enabled and self.metabo is not None:

            logs.log_info("Reinitializing the metabolism reaction network for simulation...")

            # re-read the config file again and reassign everything except for concentrations,
            #  to capture any user updates:
            self.metabo.reinitialize(self, cells, p)


        #-----gene regulatory network initialization-------------------------
        if p.grn_enabled and self.grn is None:

            logs.log_info("Initializing gene regulatory network...")

            # create an instance of the gene network simulator
            self.grn = MasterOfGenes(p)
            # read in the configuration settings for the metabolism simulator:
            self.grn.read_gene_config(self, cells, p)

        elif p.grn_enabled and self.grn is not None:

            logs.log_info("Reinitializing the gene regulatory network for simulation...")

            # re-read the config file again and reassign everything except for concentrations,
            #  to capture any user updates:
            self.grn.reinitialize(self, cells, p)


        #-----dynamic creation/anhilation of large Laplacian matrix computators!------------------
        if p.deform_osmo is True:
            # if considering osmotic water fluxes, initialize the divergence matrix:
            self.div_u_osmo = np.zeros(self.cdl)
            self.u_net = np.zeros(self.mdl)

        if p.fluid_flow is True: # If at any time fluid flow is true, initialize the flow vectors to zeros
            # initialize data structures for flow:
            self.u_cells_x = np.zeros(self.cdl)
            self.u_cells_y = np.zeros(self.cdl)

            self.u_gj_x = np.zeros(self.mdl)
            self.u_gj_y = np.zeros(self.mdl)

            if p.sim_ECM is True and cells.lapENV_P_inv is None:

                # initialize flow vectors:
                self.u_env_x = np.zeros(cells.X.shape)
                self.u_env_y = np.zeros(cells.X.shape)

                logs.log_info('Creating environmental Poisson solver for fluids...')
                bdic = {'N': 'flux', 'S': 'flux', 'E': 'flux', 'W': 'flux'}
                cells.lapENV_P, cells.lapENV_P_inv = cells.grid_obj.makeLaplacian(bound=bdic)

                cells.lapENV_P = None  # get rid of the non-inverse matrix as it only hogs memory...

        if p.deformation is True:  # if user desires deformation:

                cells.deform_tools(p)
                # create a copy of cells world, to apply deformations to for visualization purposes only:
                self.cellso = copy.deepcopy(cells)

                if p.td_deform is True and cells.lapGJ is None or cells.lapGJ_P is None:

                    # make a laplacian and solver for discrete transfers on closed, irregular cell network
                    logs.log_info('Creating cell network Poisson solver...')
                    cells.graphLaplacian(p)

        else:

            self.cellso = cells

        # if simulating electrodiffusive movement of membrane pumps and channels:-------------
        if p.sim_eosmosis is True:

            if cells.gradMem is None:
                logs.log_info("Creating tools for self-electrodiffusion of membrane pumps and channels.")
                cells.eosmo_tools(p)

            self.rho_pump = np.ones(len(cells.mem_i))
            self.rho_channel = np.ones(len(cells.mem_i))

        else:
            self.rho_pump = 1  # else just define it as identity.
            self.rho_channel = 1


        # Initialize core user-specified interventions:
        self.dyna.runAllInit(self,cells,p)

    def run_sim_core(self, cells, p):
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

        # Display and/or save an animation during solving and calculate:
        #
        # * "tt", a time-steps vector appropriate for the current run.
        # * "tsamples", that vector resampled to save data at fewer times.
        tt, tsamples = self._plot_loop(self.cellso, p)

        # Exception raised if this simulation goes unstable, permitting us to safely
        # handle this instability (e.g., by saving simulation results).
        exception_instability = None

        do_once = True  # a variable to time the loop only once

        try:
            for t in tt:  # run through the loop
                # start the timer to approximate time for the simulation
                if do_once is True:
                    loop_measure = time.time()

                # Reinitialize flux storage devices:
                self.fluxes_mem.fill(0)
                self.fluxes_gj.fill(0)

                if p.sim_ECM is True:

                    self.fluxes_env_x.fill(0)
                    self.fluxes_env_y.fill(0)

                # Calculate the values of scheduled and dynamic quantities (e.g.
                # ion channel multipliers).
                if p.run_sim is True:
                    self.dyna.runAllDynamics(self, cells, p, t)

                # -----------------PUMPS-------------------------------------------------------------------------------------

                if p.sim_ECM is True:
                    # run the Na-K-ATPase pump:
                    fNa_NaK, fK_NaK, self.rate_NaKATP = stb.pumpNaKATP(
                        self.cc_cells[self.iNa][cells.mem_to_cells],
                        self.cc_env[self.iNa][cells.map_mem2ecm],
                        self.cc_cells[self.iK][cells.mem_to_cells],
                        self.cc_env[self.iK][cells.map_mem2ecm],
                        self.vm,
                        self.T,
                        p,
                        self.NaKATP_block,
                        met = self.met_concs
                    )

                else:

                    fNa_NaK, fK_NaK, self.rate_NaKATP = stb.pumpNaKATP(
                                self.cc_cells[self.iNa][cells.mem_to_cells],
                                self.cc_env[self.iNa],
                                self.cc_cells[self.iK][cells.mem_to_cells],
                                self.cc_env[self.iK],
                                self.vm,
                                self.T,
                                p,
                                self.NaKATP_block,
                                met = self.met_concs
                            )


                # modify pump flux with any lateral membrane diffusion effects:
                fNa_NaK = self.rho_pump*fNa_NaK
                fK_NaK = self.rho_pump*fK_NaK

                if p.cluster_open is False:
                    fNa_NaK[cells.bflags_mems] = 0
                    fK_NaK[cells.bflags_mems] = 0

                # fNa_NaK, _, _ = cells.zero_div_cell(fNa_NaK, rho=0.0, bc=0.0, open_bounds=p.cluster_open)
                # fK_NaK, _, _ = cells.zero_div_cell(fK_NaK, rho=0.0, bc=0.0, open_bounds=p.cluster_open)

                # modify the fluxes by electrodiffusive membrane redistribution factor and add fluxes to storage:
                self.fluxes_mem[self.iNa] = self.fluxes_mem[self.iNa]  + fNa_NaK
                self.fluxes_mem[self.iK] = self.fluxes_mem[self.iK] + fK_NaK

                if p.metabolism_enabled:
                    # update ATP concentrations after pump action:
                    self.metabo.update_ATP(fNa_NaK, self, cells, p)

                # update the concentrations of Na and K in cells and environment:
                self.cc_cells[self.iNa],  self.cc_env[self.iNa] =  stb.update_Co(self, self.cc_cells[self.iNa],
                                                                            self.cc_env[self.iNa],fNa_NaK, cells, p,
                                                                            ignoreECM = self.ignore_ecm)

                self.cc_cells[self.iK], self.cc_env[self.iK] = stb.update_Co(self, self.cc_cells[self.iK],
                    self.cc_env[self.iK], fK_NaK, cells, p, ignoreECM = self.ignore_ecm)


                # recalculate the net, unbalanced charge and voltage in each cell:
                # self.update_V(cells, p)

                # ----------------ELECTRODIFFUSION---------------------------------------------------------------------------

                # electro-diffuse all ions (except for proteins, which don't move) across the cell membrane:

                shuffle(self.movingIons)

                for i in self.movingIons:

                    IdM = np.ones(self.mdl)

                    if p.sim_ECM is True:

                        f_ED = stb.electroflux(self.cc_env[i][cells.map_mem2ecm], self.cc_cells[i][cells.mem_to_cells],
                            self.Dm_cells[i], IdM*p.tm, self.zs[i]*IdM, self.vm, self.T, p,
                            rho=self.rho_channel)

                    else:

                        f_ED = stb.electroflux(self.cc_env[i],self.cc_cells[i][cells.mem_to_cells],
                                               self.Dm_cells[i],IdM*p.tm,self.zs[i]*IdM,
                                        self.vm,self.T,p,rho=self.rho_channel)

                    if p.cluster_open is False:
                        f_ED[cells.bflags_mems] = 0

                    # f_ED, _, _ = cells.zero_div_cell(f_ED, rho=0.0, bc=0.0, open_bounds=p.cluster_open)

                    # add membrane flux to storage
                    self.fluxes_mem[i] = self.fluxes_mem[i] + f_ED

                    # update ion concentrations in cell and ecm:
                    self.cc_cells[i], self.cc_env[i] = stb.update_Co(self, self.cc_cells[i],
                        self.cc_env[i], f_ED, cells, p, ignoreECM = self.ignore_ecm)

                    # update flux between cells due to gap junctions
                    self.update_gj(cells, p, t, i)

                    if p.sim_ECM:
                        #update concentrations in the extracellular spaces:
                        self.update_ecm(cells, p, t, i)

                    # ensure no negative concentrations:
                    stb.no_negs(self.cc_cells[i])



                # ----transport and handling of special ions---------------------------------------------------------------

                if p.ions_dict['Cl'] == 1:

                    self.cl_handler(cells, p)


                if p.ions_dict['Ca'] == 1:

                    self.ca_handler(cells, p)


                if p.ions_dict['H'] == 1:

                    self.acid_handler(cells, p)


                # update the molecules handler-----------------------------------------------------------------
                if p.molecules_enabled:

                    if p.transporters_enabled:

                        self.molecules.run_loop_transporters(t, self, self.molecules, cells, p)

                    if p.channels_enabled:

                        self.molecules.run_loop_channels(self, self.molecules, cells, p)

                    if p.modulators_enabled:

                        self.molecules.run_loop_modulators(self, self.molecules, cells, p)

                    self.molecules.run_loop(t, self, cells, p)

                # update metabolic handler----------------------------------------------------------------------

                if p.metabolism_enabled:

                    if self.metabo.transporters:
                        self.metabo.core.run_loop_transporters(t, self, self.metabo.core, cells, p)

                    if self.metabo.channels:

                        self.metabo.core.run_loop_channels(self, self.metabo.core, cells, p)

                    if self.metabo.modulators:

                        self.metabo.core.run_loop_modulators(self, self.metabo.core, cells, p)

                    self.metabo.core.run_loop(t, self, cells, p)

                # update gene regulatory network handler--------------------------------------------------------

                if p.grn_enabled:

                    if self.grn.transporters:
                        self.grn.core.run_loop_transporters(t, self, self.grn.core, cells, p)

                    if self.grn.channels and p.run_sim is True:
                        self.grn.core.run_loop_channels(self, self.grn.core, cells, p)

                    if self.grn.modulators:
                        self.grn.core.run_loop_modulators(self, self.grn.core, cells, p)

                    # update the main gene regulatory network:
                    self.grn.core.run_loop(t, self, cells, p)


                # dynamic noise handling-----------------------------------------------------------------------------------

                if p.dynamic_noise == 1 and p.ions_dict['P'] == 1:

                    # add a random walk on protein concentration to generate dynamic noise:
                    self.protein_noise_flux = p.dynamic_noise_level * (np.random.random(self.mdl) - 0.5)

                    # update the concentration of P in cells and environment:
                    self.cc_cells[self.iP], self.cc_env[self.iP] = stb.update_Co(self, self.cc_cells[self.iP],
                        self.cc_env[self.iP], self.protein_noise_flux, cells, p, ignoreECM = self.ignore_ecm)

                    # recalculate the net, unbalanced charge and voltage in each cell:
                    # self.update_V(cells, p)

                #-----forces, fields, and flow-----------------------------------------------------------------------------

                # get forces from any hydrostatic (self.P_Cells) pressure:
                getHydroF(self,cells, p)

                # calculate specific forces and pressures:

                if p.deform_osmo is True:
                    osmotic_P(self,cells, p)

                if p.deform_electro is True:
                    electro_F(self,cells, p)

                if p.fluid_flow is True:

                    self.run_sim = True

                    getFlow(self,cells, p)

                # if desired, electroosmosis of membrane channels
                if p.sim_eosmosis is True:

                    self.run_sim = True

                    eosmosis(self,cells, p)  # modify membrane pump and channel density according to Nernst-Planck

                if p.deformation is True:

                    self.run_sim = True

                    if p.td_deform is False:

                        getDeformation(self,cells, t, p)

                    elif p.td_deform is True:

                        timeDeform(self,cells, t, p)

                # recalculate the net, unbalanced charge and voltage in each cell:
                self.update_V(cells, p)

                # check for NaNs in voltage and stop simulation if found:
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
        except BetseSimInstabilityException as exception:
            exception_instability = exception

        #--------------
        cells.points_tree = None

        # Explicitly close the prior animation to conserve memory.
        self._deplot_loop()

        # Save this initialization or simulation and report results of
        # potential interest to the user.
        self.save_and_report(cells, p)

        # If the simulation did not go unstable, inform the user of success.
        if exception_instability is None:
            logs.log_info('Simulation completed successfully.')
        # Else, inform the user of this instability and re-raise the previously
        # raised exception to preserve the exact cause of this instability.
        else:
            logs.log_error('Simulation prematurely halted due to instability.')
            raise exception_instability

    #.................{  INITIALIZERS & FINALIZERS  }............................................
    def clear_storage(self, cells, p):
        """
        Re-initializes time storage vectors at the begining of a sim or init.

        """

        # clear mass flux storage vectors:
        self.fluxes_mem = np.zeros(self.fluxes_mem.shape)
        self.fluxes_gj = np.zeros(self.fluxes_gj.shape)

        self.Ax_time = []
        self.Ay_time = []

        self.Bz_time = []

        self.cc_time = []  # data array holding the concentrations at time points
        self.cc_env_time = [] # data array holding environmental concentrations at time points

        self.dd_time = []  # data array holding membrane permeabilites at time points
        self.vm_time = []  # data array holding voltage at time points
        self.vm_ave_time = []   # data array holding average vm (averaged to cell centres)
        self.vm_GHK_time = [] # data array holding GHK vm estimates
        self.time = []     # time values of the simulation
        self.gjopen_time = []   # stores the fractional gap junction open state at each time
        self.osmo_P_delta_time = []  # osmotic pressure difference between cell interior and exterior as func of time

        self.I_mem_time = []    # initialize membrane current time vector

        self.I_cell_x_time = []
        self.I_cell_y_time = []
        self.I_tot_x_time = []
        self.I_tot_y_time = []

        self.efield_gj_x_time = []   # matrices storing smooth electric field in gj connected cells
        self.efield_gj_y_time = []

        self.P_cells_time = []
        self.u_cells_x_time = []
        self.u_cells_y_time = []

        self.rho_cells_time = []

        self.F_electro_x_time = []
        self.F_electro_y_time = []

        #FIXME: Consider removing the following two list attributes, which no
        #longer appear to be initialized anywhere. Alternately, consider
        #appending the "sim.F_hydro_x" array onto the "sim.F_hydro_x_time"
        #array in the "betse.science.physics.pressures" submodule, which
        #defines the the "sim.F_hydro_x" array; likewise for the Y components.
        self.F_hydro_x_time = []
        self.F_hydro_y_time = []

        self.P_electro_time = []

        self.rate_NaKATP_time =[]

        self.pol_x_time = []  # polarization vectors
        self.pol_y_time = []

        self.Pol_tot_x_time = []
        self.Pol_tot_y_time = []

        self.Pol_tot_time = []

        self.vm_ave_time = []

        if p.deformation is True:
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

        if p.molecules_enabled:

            self.molecules.clear_cache()

        if p.metabolism_enabled:

            self.metabo.core.clear_cache()

        if p.grn_enabled:

            self.grn.core.clear_cache()

        if p.sim_eosmosis is True:
            self.rho_channel_time = []
            self.rho_pump_time = []

        if p.Ca_dyn is True:
            self.endo_retic.clear_cache()

        if p.sim_ECM is True:

            self.venv_time = []

            self.charge_env_time = []

            self.efield_ecm_x_time = []   # matrices storing smooth electric field in ecm
            self.efield_ecm_y_time = []

            # initialize time-storage vectors for electroosmotic data:
            self.u_env_x_time = []
            self.u_env_y_time = []

            self.rho_pump_time = []    # store pump and channel states as function of time...
            self.rho_channel_time = []

    def write2storage(self,t,cells,p):

        if p.GHK_calc is True:
                stb.ghk_calculator(self,cells,p)
                self.vm_GHK_time.append(self.vm_GHK) # data array holding GHK vm estimates

        # add the new concentration and voltage data to the time-storage matrices:
        self.efield_gj_x_time.append(self.E_gj_x[:])
        self.efield_gj_y_time.append(self.E_gj_y[:])

        concs = np.copy(self.cc_cells[:])
        self.cc_time.append(concs)
        concs = None

        envsc = np.copy(self.cc_env[:])
        self.cc_env_time.append(envsc)
        envsc = None

        ddc = np.copy(self.Dm_cells[:])
        ddc.tolist()
        self.dd_time.append(ddc)
        ddc = None

        self.I_cell_x_time.append(self.J_cell_x[:])
        self.I_cell_y_time.append(self.J_cell_y[:])

        self.I_mem_time.append(self.I_mem[:])

        self.vm_time.append(self.vm[:])

        # vm_ave = np.dot(cells.M_sum_mems,self.vm)/cells.num_mems
        # self.vm_ave_time.append(vm_ave)

        self.rho_cells_time.append(self.rho_cells[:])

        self.rate_NaKATP_time.append(self.rate_NaKATP[:])
        self.P_cells_time.append(self.P_cells[cells.mem_to_cells])

        if p.deform_osmo is True:
            self.osmo_P_delta_time.append(self.osmo_P_delta[cells.mem_to_cells])

        if p.deformation is True:

            # make a copy of cells to apply deformation to:
            # self.cellso = copy.deepcopy(cells)
            implement_deform_timestep(self,self.cellso, t, p)
            self.dx_cell_time.append(self.d_cells_x[:])
            self.dy_cell_time.append(self.d_cells_y[:])

        if p.fluid_flow is True:
            self.u_cells_x_time.append(self.u_cells_x[:])
            self.u_cells_y_time.append(self.u_cells_y[:])

        if p.sim_eosmosis is True:
            self.rho_channel_time.append(self.rho_channel[:])
            self.rho_pump_time.append(self.rho_pump[:])

        self.gjopen_time.append(self.gjopen[:])

        self.time.append(t)

        if p.molecules_enabled:

            self.molecules.write_data(self, cells, p)
            self.molecules.report(self, p)

        if p.metabolism_enabled:

            self.metabo.core.write_data(self, cells, p)
            self.metabo.core.report(self, p)

        if p.grn_enabled:

            self.grn.core.write_data(self, cells, p)
            self.grn.core.report(self, p)

        if p.Ca_dyn == 1 and p.ions_dict['Ca']==1:

            self.endo_retic.write_cache(self)

        if p.sim_ECM is True:

            self.efield_ecm_x_time.append(self.E_env_x[:])

            self.efield_ecm_y_time.append(self.E_env_y[:])

            self.I_tot_x_time.append(self.J_env_x[:])
            self.I_tot_y_time.append(self.J_env_y[:])

            ecmsc = np.copy(self.cc_env[:])
            ecmsc.tolist()
            self.cc_env_time.append(ecmsc)
            ecmsc = None

            self.venv_time.append(self.v_env)

            if p.fluid_flow is True:
                self.u_env_x_time.append(self.u_env_x[:])
                self.u_env_y_time.append(self.u_env_y[:])

        # polarizations:

        self.pol_x_time.append(self.pol_cell_x*1)
        self.pol_y_time.append(self.pol_cell_y*1)

        self.Pol_tot_x_time.append(self.Pol_x)
        self.Pol_tot_y_time.append(self.Pol_y)

        self.Pol_tot_time.append(self.Pol_tot)

        self.vm_ave_time.append(self.vm_ave)

    def save_and_report(self,cells,p):

        # save the init or sim:

        # get rid of the extra copy of cells
        if p.deformation:
            cells = copy.deepcopy(self.cellso)

        self.cellso = None

        if p.run_sim is False:
            datadump = [self, cells, p]
            fh.saveSim(self.savedInit, datadump)
            logs.log_info('Initialization saved to "%s".', p.init_path)
        else:
            datadump = [self, cells, p]
            fh.saveSim(self.savedSim, datadump)
            logs.log_info('Simulation saved to "%s".', p.sim_path)

        # Report final output to the user.
        for i in range(0, len(self.ionlabel)):
            endconc = np.round(np.mean(self.cc_time[-1][i]),6)
            logs.log_info(
                'Final average cytoplasmic concentration of %s: %g mmol/L',
                self.ionlabel[i], endconc)

        for i in range(0,len(self.ionlabel)):
            endconc = np.round(np.mean(self.cc_env_time[-1][i]),6)
            logs.log_info(
                'Final environmental concentration of %s: %g mmol/L',
                self.ionlabel[i], endconc)

        final_vmean = 1000*np.round(np.mean(self.vm_time[-1]),6)
        logs.log_info(
            'Final average cell Vmem: %g mV', final_vmean)

        if p.GHK_calc is True:
            final_vmean_GHK = 1000*np.round(np.mean(self.vm_GHK_time[-1]),6)
            logs.log_info(
                'Final average cell Vmem calculated using GHK: %s mV',
                final_vmean_GHK)

        if p.ions_dict['H'] == 1:
            final_pH = -np.log10(1.0e-3*np.mean((self.cc_time[-1][self.iH])))
            logs.log_info(
                'Final average cell pH: %g', np.round(final_pH, 2))

            final_pH_env = -np.log10(np.mean(1.0e-3*(self.cc_env_time[-1][self.iH])))
            logs.log_info(
                'Final environmental pH: %g', np.round(final_pH_env, 2))

        if p.molecules_enabled:
            self.molecules.report(self, p)

        if p.metabolism_enabled:
            self.metabo.core.report(self, p)

        if p.grn_enabled:
            self.grn.core.report(self, p)

        if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:
            ver_mean = np.round(1.0e3*self.endo_retic.Ver.mean(),3)
            logs.log_info("Final Ver: %g mV", ver_mean)

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

        if p.molecules_enabled:
            logs.log_info(
                'Auxiliary molecules are enabled from '
                '"biochemistry" section of main config file.')

        if p.metabolism_enabled:
            logs.log_info(
                'Metabolism is being simulated using file: %s',
                p.metabo_config_filename)

        if p.grn_enabled:
            logs.log_info(
                'A gene regulatory network is being simulated using file: %s',
                p.grn_config_filename)

    #.................{  DOOERs & GETTERS  }............................................

    def update_V(self,cells,p):

        # save the voltage as a placeholder:
        vmo = self.vm*1
        rhoo = self.rho_cells*1

        # get the charge density in the cells:
        self.rho_cells = stb.get_charge_density(self.cc_cells, self.z_array, p) + self.extra_rho_cells

        # get the currents and in-cell and environmental voltages:
        get_current(self, cells, p)

        # update charge at membranes:
        # self.rho_mems = self.rho_mems - self.Jn*p.dt

        # CHOICE 1: update Vmem in terms of current across each membrane segment:
        # dv = - (1/p.cm)*self.Jn*p.dt
        # self.vm = self.vm + dv


        # CHOICE 2: calculate vmem via average rho in cell; GJ calculations are with respect to specific Vmem:
        Q_cell = self.rho_cells * (cells.cell_vol / cells.cell_sa)

        # self.vm = (1 / p.cm) * Q_cell[cells.mem_to_cells]
        # self.vm = (1 / p.cm) * Q_cell[cells.mem_to_cells] - (1/p.cm)*self.Jn*p.dt*p.cell_polarizability

        self.vm = (1 / p.cm)*Q_cell[cells.mem_to_cells] - (1/p.cm)*self.Jn*p.dt*p.cell_polarizability
        # self.vm = self.vm + (1/p.cm)*self.dPmem*p.dt*1.0e-3
        # self.vm = self.vm - (1/p.cm)*self.Jn*p.dt*p.cell_polarizability


        if p.sim_ECM is True:  # handle polarization induced by external electric field:

            # correct the membrane voltage with the potentially finite environment voltage:
            self.vm = self.vm - p.env_modulator*self.v_env[cells.map_mem2ecm]


        # calculate the derivative of Vmem:
        self.dvm = (self.vm - vmo)/p.dt
        self.drho = (self.rho_cells - rhoo)/p.dt

        # average vm:
        self.vm_ave = np.dot(cells.M_sum_mems, self.vm)/cells.num_mems

        # polarization vectors at the membranes:
        pol_mem_x = self.vm*p.cm*p.tm*cells.mem_sa*cells.mem_vects_flat[:,2]
        pol_mem_y = self.vm*p.cm*p.tm*cells.mem_sa*cells.mem_vects_flat[:,3]

        # calculate polarization vectors for individual cells:
        self.pol_cell_x = np.dot(cells.M_sum_mems, pol_mem_x)/cells.num_mems
        self.pol_cell_y = np.dot(cells.M_sum_mems, pol_mem_y)/cells.num_mems

        # polarization density for whole cluster:
        self.Pol_x = np.sum(pol_mem_x)*(self.mdl/np.sum(cells.cell_vol))

        self.Pol_y = np.sum(pol_mem_y)*(self.mdl/np.sum(cells.cell_vol))

        self.Pol_tot = np.sqrt(self.Pol_x**2 + self.Pol_y**2)


        # try calculating a vector potential based on current density:

        self.Ax = -1.23e-6*np.dot(cells.lapGJ_P_inv, self.J_cell_x/cells.geom_weight)
        self.Ay = -1.23e-6*np.dot(cells.lapGJ_P_inv, self.J_cell_y/cells.geom_weight)

        self.Ax_time.append(self.Ax)
        self.Ay_time.append(self.Ay)

        _, _, self.Bz = cells.curl(self.Ax, self.Ay, 0)

        self.Bz_time.append(self.Bz)

    def acid_handler(self,cells,p):

        IdM = np.ones(self.mdl)

        # run the bicarbonate buffer to ensure realistic concentrations and pH in cell and environment:
        self.cc_cells[self.iH], self.pH_cell = stb.bicarbonate_buffer(self.cHM_cells, self.cc_cells[self.iM])
        self.cc_env[self.iH], self.pH_env = stb.bicarbonate_buffer(self.cHM_env, self.cc_env[self.iM])

        if p.HKATPase_dyn == 1: # if there's an H,K ATPase pump

            if p.sim_ECM is True:

                f_H2, f_K2 = stb.pumpHKATP(self.cc_cells[self.iH][cells.mem_to_cells],
                    self.cc_env[self.iH][cells.map_mem2ecm],
                    self.cc_cells[self.iK][cells.mem_to_cells], self.cc_env[self.iK][cells.map_mem2ecm],
                    self.vm, self.T, p, self.HKATP_block, met = self.met_concs)

                # modify fluxes by any uneven redistribution of pump location:
                f_H2 = self.rho_pump * f_H2
                f_K2 = self.rho_pump * f_K2

                if p.cluster_open is False:
                    f_H2[cells.bflags_mems] = 0
                    f_K2[cells.bflags_mems] = 0

                self.HKATPase_rate = f_H2[:]

                self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + f_H2
                self.fluxes_mem[self.iK] = self.fluxes_mem[self.iK] + f_K2

            else:

                f_H2, f_K2 = stb.pumpHKATP(self.cc_cells[self.iH][cells.mem_to_cells],
                                           self.cc_env[self.iH],self.cc_cells[self.iK][cells.mem_to_cells],
                                           self.cc_env[self.iK],self.vm,self.T,p,
                                           self.HKATP_block, met = self.met_concs)

                # modulate by asymmetric pump factor if lateral membrane movements are being considered:
                f_H2 = self.rho_pump*f_H2
                f_K2 = self.rho_pump*f_K2

                if p.cluster_open is False:
                    f_H2[cells.bflags_mems] = 0
                    f_K2[cells.bflags_mems] = 0

                # capture rate:
                self.HKATPase_rate = f_H2[:]

                # store fluxes for this pump:
                self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + f_H2
                self.fluxes_mem[self.iK] = self.fluxes_mem[self.iK] + f_K2


            if p.metabolism_enabled:
                # update ATP concentrations after pump action:
                self.metabo.update_ATP(f_H2, self, cells, p)


            # update the concentration in cells (assume environment steady and constant supply of ions)
            # update bicarbonate instead of H+, assuming buffer action holds:

            # calculate the update to K+ in the cell and environment:

            self.cc_cells[self.iK][:], self.cc_env[self.iK][:] = stb.update_Co(self, self.cc_cells[self.iK][:],
                self.cc_env[self.iK][:], f_K2, cells, p, ignoreECM = self.ignore_ecm)

            # Update the anion (bicarbonate) concentration instead of H+, assuming bicarb buffer holds:
            self.cc_cells[self.iM][:], self.cc_env[self.iM][:] = stb.update_Co(self, self.cc_cells[self.iM][:],
                self.cc_env[self.iM][:], -f_H2, cells, p, ignoreECM = self.ignore_ecm)


            # Calculate the new pH and H+ concentrations:
            # run the bicarbonate buffer to ensure realistic concentrations and pH in cell and environment:
            self.cc_cells[self.iH], self.pH_cell = stb.bicarbonate_buffer(self.cHM_mems,self.cc_cells[self.iM])

            self.cc_env[self.iH], self.pH_env = stb.bicarbonate_buffer(self.cHM_env,self.cc_env[self.iM])


        if p.VATPase_dyn == 1:  # if there's a V-ATPase pump

            if p.sim_ECM is True:

                # if HKATPase pump is desired, run the H-K-ATPase pump:
                f_H3 = stb.pumpVATP(self.cc_cells[self.iH][cells.mem_to_cells],
                                    self.cc_env[self.iH][cells.map_mem2ecm],
                    self.vm, self.T, p, self.VATP_block, met = self.met_concs)

                # modify flux by any uneven redistribution of pump location:
                f_H3 = self.rho_pump * f_H3

                if p.cluster_open is False:

                    f_H3[cells.bflags_mems] = 0

                self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + f_H3

            else:

                # if HKATPase pump is desired, run the H-K-ATPase pump:
                f_H3 = stb.pumpVATP(self.cc_cells[self.iH][cells.mem_to_cells],
                                    self.cc_env[self.iH],self.vm,self.T,p,self.VATP_block,
                                    met = self.met_concs)

                if p.cluster_open is False:
                    f_H3[cells.bflags_mems] = 0

                self.fluxes_mem[self.iH] = self.fluxes_mem[self.iH] + self.rho_pump * f_H3

            # Update the anion (bicarbonate) concentration instead of H+, assuming bicarb buffer holds:
            self.cc_cells[self.iM], self.cc_env[self.iM][:] = stb.update_Co(self, self.cc_cells[self.iM],
                self.cc_env[self.iM], -f_H3, cells, p, ignoreECM= self.ignore_ecm)

            if p.metabolism_enabled:
                # update ATP concentrations after pump action:
                self.metabo.update_ATP(f_H3, self, cells, p)

            # Calculate the new pH and H+ concentration:
             # run the bicarbonate buffer to ensure realistic concentrations and pH in cell and environment:

            self.cc_cells[self.iH], self.pH_cell = stb.bicarbonate_buffer(self.cHM_cells, self.cc_cells[self.iM])

            self.cc_env[self.iH], self.pH_env = stb.bicarbonate_buffer(self.cHM_env, self.cc_env[self.iM])

    def cl_handler(self, cells, p):

        if p.NaKCl_exch_dyn is True:

            cNai = self.cc_cells[self.iNa][cells.mem_to_cells]
            cKi = self.cc_cells[self.iK][cells.mem_to_cells]
            cCli = self.cc_cells[self.iCl][cells.mem_to_cells]

            if p.sim_ECM is True:
                cNao = self.cc_env[self.iNa][cells.map_mem2ecm]
                cKo = self.cc_env[self.iK][cells.map_mem2ecm]
                cClo = self.cc_env[self.iCl][cells.map_mem2ecm]

            else:

                cNao = self.cc_env[self.iNa]
                cKo = self.cc_env[self.iK]
                cClo = self.cc_env[self.iCl]

            f_Na, f_K, f_Cl = stb.exch_NaKCl(cNai,cNao,cKi,cKo,cCli,cClo,self.vm,self.T,p)

            # FIXME where does this add to fluxes_mems?

            # update concentrations of Na, K and Cl in cells and environment:

            self.cc_cells[self.iNa][:], self.cc_env[self.iNa][:] = stb.update_Co(self, self.cc_cells[self.iNa][:],
                self.cc_env[self.iNa][:], f_Na, cells, p, ignoreECM=self.ignore_ecm)

            self.cc_cells[self.iK][:], self.cc_env[self.iK][:] = stb.update_Co(self, self.cc_cells[self.iK][:],
                self.cc_env[self.iK][:], f_K, cells, p, ignoreECM=self.ignore_ecm)

            self.cc_cells[self.iCl][:], self.cc_env[self.iCl][:] = stb.update_Co(self, self.cc_cells[self.iCl][:],
                self.cc_env[self.iCl][:], f_Cl, cells, p, ignoreECM=self.ignore_ecm)


        if p.ClK_symp_dyn is True:

            cKi = self.cc_cells[self.iK][cells.mem_to_cells]
            cCli = self.cc_cells[self.iCl][cells.mem_to_cells]

            if p.sim_ECM is True:
                cKo = self.cc_env[self.iK][cells.map_mem2ecm]
                cClo = self.cc_env[self.iCl][cells.map_mem2ecm]

            else:
                cKo = self.cc_env[self.iK]
                cClo = self.cc_env[self.iCl]

            f_K, f_Cl = stb.symp_ClK(cKi, cKo, cCli, cClo, self.vm, self.T, p)

            # FIXME add to flux_mems

            # update concentrations of K and Cl in cells and environment:
            self.cc_cells[self.iK][:], self.cc_env[self.iK][:] = stb.update_Co(self, self.cc_cells[self.iK][:],
                self.cc_env[self.iK][:], f_K, cells, p, ignoreECM=self.ignore_ecm)

            self.cc_cells[self.iCl][:], self.cc_env[self.iCl][:] = stb.update_Co(self, self.cc_cells[self.iCl][:],
                self.cc_env[self.iCl][:], f_Cl, cells, p, ignoreECM=self.ignore_ecm)

    def ca_handler(self,cells,p):

        if p.sim_ECM is True:

            # run Ca-ATPase

            f_CaATP = stb.pumpCaATP(self.cc_cells[self.iCa][cells.mem_to_cells],
                self.cc_env[self.iCa][cells.map_mem2ecm],
                self.vm, self.T, p, self.CaATP_block, met = self.met_concs)

            f_CaATP = self.rho_pump * f_CaATP

            if p.NaCa_exch_dyn is True:
                # run Na Ca exchanger
                f_NaEx, f_CaEx = stb.exch_NaCa(self.cc_cells[self.iNa][cells.mem_to_cells],
                        self.cc_env[self.iNa][cells.map_mem2ecm],
                        self.cc_cells[self.iCa][cells.mem_to_cells],
                        self.cc_env[self.iCa][cells.map_mem2ecm],
                        self.vm, self.T, p, self.NaCaExch_block)

            else:
                f_NaEx = np.zeros(self.mdl)
                f_CaEx = np.zeros(self.mdl)

        else:

            # run Ca-ATPase

            f_CaATP = stb.pumpCaATP(self.cc_cells[self.iCa][cells.mem_to_cells],
                                    self.cc_env[self.iCa], self.vm, self.T, p,
                                    self.CaATP_block, met = self.met_concs)

            if p.NaCa_exch_dyn is True:
                # run Na Ca exchanger
                f_NaEx, f_CaEx = stb.exch_NaCa(self.cc_cells[self.iNa][cells.mem_to_cells],
                        self.cc_env[self.iNa],
                        self.cc_cells[self.iCa][cells.mem_to_cells],
                        self.cc_env[self.iCa],
                        self.vm, self.T, p)

            else:
                f_NaEx = np.zeros(self.mdl)
                f_CaEx = np.zeros(self.mdl)

        self.rate_CaATP = f_CaATP

        if p.cluster_open is False:
            f_CaATP[cells.bflags_mems] = 0
            f_NaEx[cells.bflags_mems] = 0
            f_CaEx[cells.bflags_mems] = 0

        # store the transmembrane flux for this ion
        self.fluxes_mem[self.iCa] = self.fluxes_mem[self.iCa]  + self.rho_pump*(f_CaATP + f_CaEx)
        self.fluxes_mem[self.iNa] = self.fluxes_mem[self.iNa] + self.rho_pump*(f_NaEx)


        if p.metabolism_enabled:
            # update ATP concentrations after pump action:
            self.metabo.update_ATP(f_CaATP, self, cells, p)

        # update calcium concentrations in cell and ecm:

        self.cc_cells[self.iCa], self.cc_env[self.iCa] = stb.update_Co(self, self.cc_cells[self.iCa],
            self.cc_env[self.iCa], f_CaATP + f_CaEx, cells, p, ignoreECM = True)

        if p.NaCa_exch_dyn:


            self.cc_cells[self.iNa], self.cc_env[self.iNa] = stb.update_Co(self, self.cc_cells[self.iNa],
                self.cc_env[self.iNa], f_NaEx, cells, p, ignoreECM = True)


        if p.Ca_dyn == 1:  # do endoplasmic reticulum handling

            self.endo_retic.update(self, cells, p)

    def update_gj(self,cells,p,t,i):

        # calculate voltage difference (gradient*len_gj) between gj-connected cells:

        self.vgj = self.vm[cells.nn_i]- self.vm[cells.mem_i]

        # self.vgj = self.vm_ave[cells.mem_to_cells][cells.nn_i] - self.vm_ave[cells.mem_to_cells][cells.mem_i]

        self.Egj = -self.vgj/cells.gj_len

        self.E_gj_x = self.Egj*cells.mem_vects_flat[:,2]
        self.E_gj_y = self.Egj*cells.mem_vects_flat[:,3]

        if p.v_sensitive_gj is True:

            # run the gap junction dynamics object to update gj open state of sim:
            self.gj_funk.run(self, cells, p)

        else:
            self.gjopen = self.gj_block*np.ones(len(cells.mem_i))

        # concentration gradient for ion i:
        conc_mem = self.cc_cells[i][cells.mem_to_cells]

        grad_cgj = (conc_mem[cells.nn_i] - conc_mem[cells.mem_i])/(cells.gj_len)

        gcx = grad_cgj*cells.mem_vects_flat[:, 2]
        gcy = grad_cgj*cells.mem_vects_flat[:, 3]

        # midpoint concentration:
        c = (conc_mem[cells.nn_i] + conc_mem[cells.mem_i])/2

        # electroosmotic fluid velocity at gap junctions:
        if p.fluid_flow is True:
            ux = self.u_gj_x
            uy = self.u_gj_y

        else:

            ux = 0
            uy = 0


        fgj_x, fgj_y = stb.nernst_planck_flux(c, gcx, gcy, -self.E_gj_x,
                                          -self.E_gj_y, ux, uy,
                                              p.gj_surface*self.gjopen*self.D_gj[i],
                                              self.zs[i],
                                              self.T, p)

        fgj_X = fgj_x*cells.mem_vects_flat[:,2] + fgj_y*cells.mem_vects_flat[:,3]

        # enforce zero flux at outer boundary:
        fgj_X[cells.bflags_mems] = 0.0

        # fgj_X, _, _ = cells.zero_div_cell(fgj_X, rho=0.0, bc=0.0, open_bounds=False)

        # divergence calculation for individual cells (finite volume expression)
        delta_cco = np.dot(cells.M_sum_mems, -fgj_X*cells.mem_sa) / cells.cell_vol

        # Calculate the final concentration change:
        self.cc_cells[i] = self.cc_cells[i] + p.dt*delta_cco

        # self.fluxes_mem[i] = self.fluxes_mem[i] - fgj_X  # store gap junction flux for this ion
        self.fluxes_gj[i] = self.fluxes_gj[i] + fgj_X  # store gap junction flux for this ion

    def update_ecm(self,cells,p,t,i):

        cenv = self.cc_env[i]
        cenv = cenv.reshape(cells.X.shape)

        if p.smooth_level > 0.0:
            cenv = gaussian_filter(cenv, p.smooth_level)

        v_env = self.v_env.reshape(cells.X.shape)

        # enforce voltage at boundary:
        v_env[:,0] = self.bound_V['L']
        v_env[:,-1] = self.bound_V['R']
        v_env[0,:] = self.bound_V['B']
        v_env[-1,:] = self.bound_V['T']

        cenv[:,0] =  self.c_env_bound[i]
        cenv[:,-1] =  self.c_env_bound[i]
        cenv[0,:] =  self.c_env_bound[i]
        cenv[-1,:] =  self.c_env_bound[i]
        #
        # print(self.bound_V['L'])

        gcx, gcy = fd.gradient(cenv, cells.delta)

        gvx, gvy = fd.gradient(v_env, cells.delta)

        self.E_env_x = -gvx
        self.E_env_y = -gvy

        fx, fy = stb.nernst_planck_flux(cenv, gcx, gcy, -self.E_env_x, -self.E_env_y, 0, 0,
                                        self.D_env[i].reshape(cells.X.shape), self.zs[i],
                                        self.T, p)


        div_fa = fd.divergence(-fx, -fy, cells.delta, cells.delta)

        self.fluxes_env_x[i] = self.fluxes_env_x[i] + fx.ravel()  # store ecm junction flux for this ion
        self.fluxes_env_y[i] = self.fluxes_env_y[i] + fy.ravel()  # store ecm junction flux for this ion

        cenv = cenv + div_fa * p.dt

        self.cc_env[i] = cenv.ravel()


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

        elif label == 'M':

            ion = self.iM

        elif label == 'H':

            ion = self.iH

        elif label == 'P':

            ion = self.iP

        else:
            logs.warning('Oops! Molecule gated ion channel target not found!')
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
            # dummyMems[cells.bflags_mems] = self.D_free[i]*p.D_tj*self.Dtj_rel[i]
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
            Denv_o[cells.map_cell2ecm][cells.bflags_cells] = self.D_free[i] * p.D_tj * self.Dtj_rel[i]

            # create an ecm diffusion grid filled with the environmental values
            self.D_env[i] = Denv_o[:]*1

            # self.D_env[i][cells.ecm_bound_k] = self.D_free[i] * p.D_tj * self.Dtj_rel[i]


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
    def _plot_loop(self, cells, p):
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

        #FIXME: Refactor these conditionally set local variables into
        #attributes of the "Parameters" object set instead by the
        #Parameters.set_time_profile() method.

        # Type and maximum number of steps of the current run.
        if p.run_sim is False:
            phase_verb = 'Initializing'
            phase_noun = 'initialization'
            phase_time_step_count = p.init_tsteps
        else:
            phase_verb = 'Simulating'
            phase_noun = 'simulation'
            phase_time_step_count = p.sim_tsteps

        #FIXME: Replace "phase_time_len" with "p.total_time", which is the same
        #exact value. (Duplication is bad for coding health.)

        # Total number of seconds simulated by the current run.
        phase_time_len = phase_time_step_count * p.dt

        # Time-steps vector appropriate for the current run.
        tt = np.linspace(0, phase_time_len, phase_time_step_count)

        tsamples = set()
        i = 0
        while i < len(tt) - p.t_resample:
            i += p.t_resample
            i = int(i)
            tsamples.add(tt[i])

        # Log this run.
        logs.log_info(
            'Your %s is running from 0 to %.2f s of in-world time '
            'in %d time steps (%d sampled).',
            phase_noun,
            phase_time_len,
            len(tt),
            len(tsamples),
        )

        # If plotting an animation during simulation computation, do so.
        if p.anim.is_while_sim:
            self._anim_cells_while_solving = AnimCellsWhileSolving(
                label='Vmem',
                figure_title='Vmem while {}'.format(phase_verb),
                colorbar_title='Voltage [mV]',
                color_min=p.Vmem_min_clr,
                color_max=p.Vmem_max_clr,
                is_color_autoscaled=p.autoscale_Vmem,
                sim=self, cells=cells, p=p,

                # The number of animation frames is the number of sampled time
                # steps, as the _replot_loop() method plotting each such frame
                # is called *ONLY* for each such step.
                time_step_count=len(tsamples),
            )
        # Else, nullify the object encapsulating this animation for safety.
        else:
            self._anim_cells_while_solving = None

        return tt, tsamples

    def _replot_loop(self, p):
        '''
        Update the currently displayed and/or saved animation during solving
        with the results of the most recently solved time step, if requested.
        '''

        # Update this animation for the "last frame" if desired, corresponding
        # to the results of the most recently solved time step.
        if p.anim.is_while_sim:
            self._anim_cells_while_solving.plot_frame(time_step=-1)

    def _deplot_loop(self):
        '''
        Explicitly close the previously displayed and/or saved animation if
        any _or_ noop otherwise.

        To conserve memory, this method nullifies and hence garbage collects
        both this animation and this animation's associated figure.
        '''

        # Close the previously displayed and/or saved animation if any.
        if self._anim_cells_while_solving is not None:
            self._anim_cells_while_solving.close()
            self._anim_cells_while_solving = None

        # For safety, close all remaining plots and animations as well.
        plt.close()

#-----------------------------------------------------------------------------------------------------------------------
# WASTELANDS
#-----------------------------------------------------------------------------------------------------------------------
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




# Calculate flux across TJ for this ion:

# IdM = np.ones(self.mdl)
#
# f_TJ = stb.electroflux(self.c_env_bound[i]*IdM, self.cc_env[i][cells.map_mem2ecm],
#                     self.D_env[i][cells.map_mem2ecm], IdM*p.rc, self.zs[i]*IdM, self.v_env[cells.map_mem2ecm],
#                        self.T, p)
#
#
# self.J_TJ = self.J_TJ - p.F*self.zs[i]*f_TJ


# #-----old way------------------------------------------------------------
#
#
# if p.closed_bound is True:
#     btag = 'closed'
#
# else:
#     btag = 'open'
#
# # make v_env and cc_env into 2d matrices
# cenv = self.cc_env[i]
#
# v_env = self.v_env.reshape(cells.X.shape)
#
# # enforce voltage at boundary:
# v_env[:,0] = self.bound_V['L']
# v_env[:,-1] = self.bound_V['R']
# v_env[0,:] = self.bound_V['B']
# v_env[-1,:] = self.bound_V['T']
#
# cenv = cenv.reshape(cells.X.shape)
#
# # if p.smooth_level > 0.0:
# #     cenv = gaussian_filter(cenv, p.smooth_level)
#
# # prepare concentrations and diffusion constants for MACs grid format
# # by resampling the values at the u v coordinates of the flux:
# cenv_x = np.zeros(cells.grid_obj.u_shape)
# cenv_y = np.zeros(cells.grid_obj.v_shape)
#
# # create the proper shape for the concentrations and state appropriate boundary conditions:
# cenv_x[:,1:] = cenv
#
# cenv_y[1:,:] = cenv
#
# if p.closed_bound is True: # insulation boundary conditions
#
#     cenv_x[:,0] = cenv_x[:,1]
#     cenv_x[:,-1] = cenv_x[:,-2]
#     cenv_x[0,:] = cenv_x[1,:]
#     cenv_x[-1,:] = cenv_x[-2,:]
#
#     cenv_y[0,:] = cenv_y[1,:]
#     cenv_y[-1,:] = cenv_y[-2,:]
#     cenv_y[:,0] = cenv_y[:,1]
#     cenv_y[:,-1] = cenv_y[:,-2]
#
# else:   # open and electrically grounded boundary conditions
#     cenv_x[:,0] =  self.c_env_bound[i]
#     cenv_x[:,-1] =  self.c_env_bound[i]
#     cenv_x[0,:] =  self.c_env_bound[i]
#     cenv_x[-1,:] =  self.c_env_bound[i]
#
#     cenv_y[0,:] =  self.c_env_bound[i]
#     cenv_y[-1,:] =  self.c_env_bound[i]
#     cenv_y[:,0] =  self.c_env_bound[i]
#     cenv_y[:,-1] =  self.c_env_bound[i]
#
#
# # calculate gradients in the environment
# grad_V_env_x, grad_V_env_y = cells.grid_obj.grid_gradient(v_env,bounds='closed')
#
# grad_cc_env_x, grad_cc_env_y = cells.grid_obj.grid_gradient(cenv,bounds=btag)
#
#
# grad_V_env_x[:,1:] = grad_V_env_x[:,1:]  + self.E_env_x
#
# grad_V_env_y[1:,:] = grad_V_env_y[1:,:] + self.E_env_y
#
# # calculate fluxes for electrodiffusive transport:
# if p.fluid_flow is True:
#
#     uenvx = np.zeros(cells.grid_obj.u_shape)
#     uenvy = np.zeros(cells.grid_obj.v_shape)
#
#     uenvx[:,1:] = self.u_env_x
#     uenvy[1:,:] = self.u_env_y
#
#     if p.closed_bound is False:
#
#         uenvx[:,0] = uenvx[:,1]
#         uenvx[:,-1]= uenvx[:,-2]
#         uenvx[0,:] = uenvx[1,:]
#         uenvx[-1,:] = uenvx[-2,:]
#
#         uenvy[:,0] = uenvy[:,1]
#         uenvy[:,-1]= uenvy[:,-2]
#         uenvy[0,:] = uenvy[1,:]
#         uenvy[-1,:] = uenvy[-2,:]
#
#     else:
#
#         uenvx[:,0] = 0
#         uenvx[:,-1]= 0
#         uenvx[0,:] = 0
#         uenvx[-1,:] = 0
#
#         uenvy[:,0] = 0
#         uenvy[:,-1]= 0
#         uenvy[0,:] = 0
#         uenvy[-1,:] = 0
#
#
# else:
#     uenvx = 0
#     uenvy = 0
#
# # distance over which electric field is active is affected by electrolyte screening. Therefore
# # we define a field modulation in the bulk, which is assumed to be a fraction of the Debye length,
# # which is assumed to be about 1 nm:
#
# # self.field_mod = (1e-9/p.cell_space)
# self.field_mod = 0.0
#
# f_env_x, f_env_y = stb.np_flux_special(cenv_x,cenv_y,grad_cc_env_x,grad_cc_env_y,
#     self.field_mod*grad_V_env_x, self.field_mod*grad_V_env_y, uenvx,uenvy,self.D_env_u[i],self.D_env_v[i],
#     self.zs[i],self.T,p)
#
#
# IdM = np.ones(self.mdl)
#
#
# # Calculate flux across TJ for this ion:
#
# f_TJ = stb.electroflux(self.c_env_bound[i]*IdM, self.cc_env[i][cells.map_mem2ecm],
#                     self.D_env[i][cells.map_mem2ecm], IdM*p.rc, self.zs[i]*IdM, self.v_env[cells.map_mem2ecm],
#                        self.T, p)
#
#
# self.J_TJ = self.J_TJ - p.F*self.zs[i]*f_TJ
#
#
# if p.closed_bound is False:
#
#     f_env_x[:,0] = f_env_x[:,1]
#     f_env_x[:,-1]= f_env_x[:,-2]
#     f_env_x[0,:] = f_env_x[1,:]
#     f_env_x[-1,:] = f_env_x[-2,:]
#
#     f_env_y[:,0] = f_env_y[:,1]
#     f_env_y[:,-1]= f_env_y[:,-2]
#     f_env_y[0,:] = f_env_y[1,:]
#     f_env_y[-1,:] = f_env_y[-2,:]
#
# else:
#
#     f_env_x[:,0] = 0
#     f_env_x[:,-1]= 0
#     f_env_x[0,:] = 0
#     f_env_x[-1,:] = 0
#
#     f_env_y[:,0] = 0
#     f_env_y[:,-1]= 0
#     f_env_y[0,:] = 0
#     f_env_y[-1,:] = 0
#
# # calculate the divergence of the total flux, which is equivalent to the total change per unit time
#
# d_fenvx = -(f_env_x[:,1:] - f_env_x[:,0:-1])/cells.delta
# d_fenvy = -(f_env_y[1:,:] - f_env_y[0:-1,:])/cells.delta
#
# delta_c = d_fenvx + d_fenvy
#
# #-----------------------
# cenv = cenv + delta_c*p.dt
#
#
# if p.closed_bound is True:
#     # Neumann boundary condition (flux at boundary)
#     # zero flux boundaries for concentration:
#
#     cenv[:,-1] = cenv[:,-2]
#     cenv[:,0] = cenv[:,1]
#     cenv[0,:] = cenv[1,:]
#     cenv[-1,:] = cenv[-2,:]
#
# elif p.closed_bound is False:
#     # if the boundary is open, set the concentration at the boundary
#     # open boundary
#     cenv[:,-1] = self.c_env_bound[i]
#     cenv[:,0] = self.c_env_bound[i]
#     cenv[0,:] = self.c_env_bound[i]
#     cenv[-1,:] = self.c_env_bound[i]
#
# self.cc_env[i] = cenv.ravel()
#
# fenvx = (f_env_x[:, 1:] + f_env_x[:, 0:-1]) / 2
# fenvy = (f_env_y[1:, :] + f_env_y[0:-1, :]) / 2
#
# self.fluxes_env_x[i] = self.fluxes_env_x[i] + fenvx.ravel()  # store ecm junction flux for this ion
# self.fluxes_env_y[i] = self.fluxes_env_y[i] + fenvy.ravel()  # store ecm junction flux for this ion
