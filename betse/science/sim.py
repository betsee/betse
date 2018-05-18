#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS                           }....................
import copy, time
import numpy as np
from betse.exceptions import BetseSimException, BetseSimUnstableException
from betse.science import filehandling as fh
from betse.science import sim_toolbox as stb
from betse.science.channels.gap_junction import Gap_Junction
from betse.science.chemistry.gene import MasterOfGenes
from betse.science.chemistry.molecules import MasterOfMolecules
from betse.science.config.confenum import SolverType
from betse.science.math import finitediff as fd
from betse.science.organelles.endo_retic import EndoRetic
from betse.science.physics.deform import (
    getDeformation, timeDeform, implement_deform_timestep)
from betse.science.physics.flow import getFlow
from betse.science.physics.ion_current import get_current
from betse.science.physics.pressures import osmotic_P
from betse.science.phase.phasecls import SimPhase
from betse.science.phase.phaseenum import SimPhaseKind
from betse.science.organelles.microtubules import Mtubes
from betse.science.tissue.tishandler import TissueHandler
from betse.science.visual.anim.animwhile import AnimCellsWhileSolving
from betse.util.io.log import logs
from betse.util.path import pathnames
from betse.util.type.contexts import noop_context
from betse.util.type.types import type_check, NoneType
from numpy import ndarray
from scipy.ndimage.filters import gaussian_filter

# ....................{ CLASSES                            }....................
class Simulator(object):
    '''
    Phase-specific simulator, simulating networked cell bioelectrical activity
    for a specific phase (e.g., seed, initialization, simulation) and storing
    the results of this activity.

    For efficiency and scalability, *all* simulation methods are implemented in
    terms of Numpy-based linear algebra.

    Methods
    -------
    update_V(cells,p,t)             Gets charge densities in cells and environment
                                      and calculates respective voltages.

    update_C(ion_i,flux, cells, p)     Updates concentration of ion with index
                                            ion_i in cell and environment for a flux leaving the cell.

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

    Attributes (Counts)
    ----------
    cdl : int
        Number of cells in this simulated cluster.
    edl : int
        Number of extracellular grid spaces in this simulated cluster if this
        simulation enables extracellular spaces *or* :attr:`mdl` otherwise.
    mdl : int
        Number of cell membranes in this simulated cluster.

    Attributes (Time)
    ----------
    time : list
        One-dimensional ordered list of all **sampled time steps** (i.e., time
        steps at which data is sampled, substantially reducing data storage).
        The length of this list governs the number of:
        * Frames in each exported animation.
        * Rows in each exported spreadsheet.

    Attributes (Current Density: Extracellular)
    ----------
    The following attributes are *always* computed when the full solver is
    enabled, regardless of whether extracellular spaces are enabled.

    J_env_x : ndarray
        One-dimensional Numpy array of the X components of all extracellular
        current densities for the current time step, whose:
        * Only dimension indexes each square environmental grid space (in
          either dimension) such that each item is the X component of the
          extracellular current density vector spatially situated at the centre
          of that space for this time step.
    J_env_y : ndarray
        One-dimensional Numpy array of the Y components of all extracellular
        current densities for the current time step, whose:
        * Only dimension indexes each square environmental grid space (in
          either dimension) such that each item is the Y component of the
          extracellular current density vector spatially situated at the centre
          of that space for this time step.
    I_tot_x_time : list
        Two-dimensional list of the X components of all extracellular current
        densities over all time steps, whose:
        * First dimension indexes each sampled time step.
        * Second dimension yields a one-dimensional Numpy array of the X
          components of all extracellular current densities for this time step,
          defined as for the corresponding :attr:`J_env_x` array.
        Equivalently, this list is the concatenation of all :attr:`J_env_x`
        arrays for all sampled time steps.
    I_tot_y_time : list
        Two-dimensional list of the Y components of all extracellular current
        densities over all time steps, whose:
        * First dimension indexes each sampled time step.
        * Second dimension yields a one-dimensional Numpy array of the Y
          components of all extracellular current densities for this time step,
          defined as for the corresponding :attr:`J_env_y` array.
        Equivalently, this list is the concatenation of all :attr:`J_env_y`
        arrays for all sampled time steps.

    Attributes (Current Density: Intracellular)
    ----------
    I_cell_x_time : list
        Two-dimensional list of the X components of all intracellular current
        densities, whose:
        * First dimension indexes each sampled time step.
        * Second dimension indexes each cell such that each item is the X
          component of the intracellular current density vector spatially
          situated at the center of that cell for this time step.
    I_cell_y_time : list
        Two-dimensional list of the Y components of all intracellular current
        densities, whose:
        * First dimension indexes each sampled time step.
        * Second dimension indexes each cell such that each item is the Y
          component of the intracellular current density vector spatially
          situated at the center of that cell for this time step.

    Attributes (Deformation)
    ----------
    The following attributes are defined *only* if this simulation enables
    deformations. If this is *not* the case for this simulation, these
    attributes remain undefined.

    cell_verts_time : ndarray
        Four-dimensional list whose:
        * First dimension indexes each sampled time step.
        * Second dimension indexes each cell such that each item is the
          matplotlib-compatible polygon patch for that cell defined as for the
          :attr:`betse.science.cells.Cells.cell_verts` array.
        Equivalently, this list is the concatenation of all
        :attr:`betse.science.cells.Cells.cell_verts` arrays for all sampled time
        steps.
    d_cells_x : ndarray
        One-dimensional Numpy array indexing each cell such that each item is
        the X component of the **total cellular displacement** (i.e., summation
        of all cellular deformations due to galvanotropic and osmotic pressure
        body forces) spatially situated at the centre of the cell indexed by
        that item for the current time step.
    d_cells_y : ndarray
        One-dimensional Numpy array indexing each cell such that each item is
        the Y component of the total cellular displacement spatially situated at
        the centre of the cell indexed by that item for the current time
        step.
    dx_cell_time : list
        Two-dimensional list of the X components of all cellular deformations,
        for all sampled time steps, whose:
        * First dimension indexes each sampled time step.
        * Second dimension indexes each cell such that each item is the X
          component of the total deformation for that cell defined as for the
          :attr:`d_cells_x` array.
        Equivalently, this list is the concatenation of all :attr:`d_cells_x`
        arrays for all sampled time steps.
    dy_cell_time : list
        Two-dimensional list of the Y components of all cellular deformations
        for all sampled time steps, whose:
        * First dimension indexes each sampled time step.
        * Second dimension indexes each cell such that each item is the Y
          component of the total deformation for that cell defined as for the
          :attr:`d_cells_y` array.
        Equivalently, this list is the concatenation of all :attr:`d_cells_y`
        arrays for all sampled time steps.

    Attributes (Electric Field: Extracellular)
    ----------
    The following attributes are defined *only* if this simulation enables
    extracellular spaces. If this is *not* the case for this simulation, these
    attributes remain undefined.

    E_env_x : ndarray
        Two-dimensional Numpy array of the X components of the extracellular
        electric field for the current time step, whose:
        * First dimension indexes each row of the extracellular spaces grid.
        * Second dimension indexes each column of the extracellular spaces grid
          such that each item is the X component of the extracellular
          electric field spatially situated at the centre of the extracellular
          grid space corresponding to the current row and column.
    E_env_y : ndarray
        Two-dimensional Numpy array of the Y components of the extracellular
        electric field for the current time step, whose:
        * First dimension indexes each row of the extracellular spaces grid.
        * Second dimension indexes each column of the extracellular spaces grid
          such that each item is the Y component of the extracellular
          electric field spatially situated at the centre of the extracellular
          grid space corresponding to the current row and column.
    efield_ecm_x_time : list
        Three-dimensional list of the X components of the extracellular
        electric fields over all time steps, whose:
        * First dimension indexes each sampled time step.
        * Second dimension yields a two-dimensional Numpy array of the X
          components of the extracellular electric field for this time step,
          defined as for the corresponding :attr:`E_env_x` array.
        Equivalently, this list is the concatenation of all :attr:`E_env_x`
        arrays for all sampled time steps.
    efield_ecm_y_time : list
        Three-dimensional list of the Y components of the extracellular
        electric fields over all time steps, whose:
        * First dimension indexes each sampled time step.
        * Second dimension yields a two-dimensional Numpy array of the Y
          components of the extracellular electric field for this time step,
          defined as for the corresponding :attr:`E_env_y` array.
        Equivalently, this list is the concatenation of all :attr:`E_env_y`
        arrays for all sampled time steps.

    Attributes (Electric Field: Intracellular)
    ----------
    E_gj_x : ndarray
        One-dimensional Numpy array indexing each cell membrane such that each
        item is the X component of the intracellular electric field vector
        for the current time step spatially situated across the gap junction to
        which the membrane indexed by that item connects: specifically, the
        difference of the X component of the voltage situated at this membrane
        with that of the voltage situated at the gap junction-connected membrane
        adjacent to this membrane, divided by the length in meters of this gap
        junction.
    E_gj_y : ndarray
        One-dimensional Numpy array indexing each cell membrane such that each
        item is the Y component of the intracellular electric field vector
        for the current time step defined as for the corresponding
        :attr:`E_gj_X` array.
    efield_gj_x_time : list
        Two-dimensional list of the X components of the intracellular electric
        fields for all time steps, whose:
        * First dimension indexes each sampled time step.
        * Second dimension indexes each cell membrane such that each item is
          the X component of the intracellular electric field vector defined as
          for the :attr:`E_gj_x` array.
        Equivalently, this list is the concatenation of all :attr:`E_gj_x`
        arrays for all sampled time steps.
    efield_gj_y_time : list
        Two-dimensional list of the Y components of the intracellular electric
        fields for all time steps, whose:
        * First dimension indexes each sampled time step.
        * Second dimension indexes each cell membrane such that each item is
          the Y component of the intracellular electric field vector defined as
          for the :attr:`E_gj_y` array.
        Equivalently, this list is the concatenation of all :attr:`E_gj_y`
        arrays for all sampled time steps.

    Attributes (Ion)
    ----------
    cc_cells : ndarray
        Two-dimensional Numpy array of all cellular ion concentrations for the
        current time step, whose:
        #. First dimension indexes each ion enabled by the current ion profile.
        #. Second dimension indexes each cell such that each item is the
           concentration of that ion in that cell.
    cc_time : list
        Three-dimensional list of all cellular ion concentrations for all time
        steps, whose:
        #. First dimension indexes each sampled time step.
        #. Second dimension indexes each ion such that each item is the array
           of all cellular concentrations of that ion for this time step,
           defined as for the :attr:`cc_cells` array.
        Equivalently, this list is the concatenation of all :attr:`cc_cells`
        arrays for all sampled time steps.

    Attributes (Ion: Index)
    ----------
    The following ion indices are dynamically defined by the
    :meth:`init_core` method.

    iCa : int
        0-based index of the calcium ion if enabled by this simulation's ion
        profile *or* undefined otherwise.
    iCl : int
        0-based index of the chloride ion if enabled by this simulation's ion
        profile *or* undefined otherwise.
    iH : int
        0-based index of the hydrogen ion if enabled by this simulation's ion
        profile *or* undefined otherwise.
    iNa : int
        0-based index of the sodium ion if enabled by this simulation's ion
        profile *or* undefined otherwise.
    iK : int
        0-based index of the potassium ion if enabled by this simulation's ion
        profile *or* undefined otherwise.
    iM : int
        0-based index of the bioarbonate ion if enabled by this simulation's ion
        profile *or* undefined otherwise.
    iP : int
        0-based index of all anionic proteins if enabled by this simulation's
        ion profile *or* undefined otherwise.

    Attributes (Microtubules)
    ----------
    mtubes : Mtubes
        Object encapsulating all microtubules for the current time step.
    mtubes_x_time : list
        Two-dimensional list whose:
        * First dimension indexes each sampled time step.
        * Second dimension indexes each cell membrane such that each item is
          the X component of the microtubule unit vector spatially situated at
          the midpoint of that membrane for this time step.
    mtubes_y_time : list
        Two-dimensional list whose:
        * First dimension indexes each sampled time step.
        * Second dimension indexes each cell membrane such that each item is
          the Y component of the microtubule unit vector spatially situated at
          the midpoint of that membrane for this time step.

    Attributes (Pressure)
    ----------
    P_cells : ndarray
        One-dimensional Numpy array indexing each cell such that each item is
        the **total cellular pressure** (i.e., summation of the mechanical and
        osmotic cellular pressure) spatially situated at the centre of the cell
        indexed by that item for the current time step.
    P_cells_time : list
        Two-dimensional list of all **total cellular pressures** (i.e.,
        summation of all mechanical and osmotic cellular pressures) for all time
        steps, whose:
        * First dimension indexes each sampled time step.
        * Second dimension indexes each cell such that each item is the
          total pressure for that cell defined as for the :attr:`P_cells` array.
        Equivalently, this list is the concatenation of all :attr:`P_cells`
        arrays for all sampled time steps.

    Attributes (Voltage: Extracellular)
    ----------
    v_env : ndarray
        One-dimensional Numpy array of all **extracellular voltages** (i.e.,
        voltages at the outer membrane surfaces of all cells) at the current
        time step, indexing each environmental grid space such that each item
        is the extracellular voltage spatially situated at the centre of that
        grid space.
    venv_time : list
        Two-dimensional list of all extracellular voltages over all time steps,
        whose:
        #. First dimension indexes each sampled time step.
        #. Second dimension indexes each environmental grid space such that each
           item is the extracellular voltage for that grid space at this time
           step, defined as for the :attr:`v_env` array.
        Equivalently, this array is the concatenation of all :attr:`v_env`
        arrays over all time steps.

    Attributes (Voltage: Intracellular)
    ----------
    vcell_time : ndarray
        Voltage at the inner membrane surface of each cell as a function of
        time.

    Attributes (Voltage: Transmembrane)
    ----------
    vm : ndarray
        One-dimensional Numpy array of all transmembrane voltages across all
        cell membranes at the current time step, indexing each cell membrane
        such that each item is the transmembrane voltage spatially situated
        across that cell membrane.
    vm_time : list
        Two-dimensional list of all transmembrane voltages across all cell
        membranes over all sampled time steps, whose:
        #. First dimension indexes each sampled time step.
        #. Second dimension indexes each cell membrane such that each item is
           the transmembrane voltage for that cell at this time step, defined as
           for the :attr:`vm` array.
        Equivalently, this array is the concatenation of all :attr:`vm` arrays
        over all sampled time steps.
    vm_ave : ndarray
        One-dimensional Numpy array indexing each cell such that each item is
        the transmembrane voltage spatially situated at the centre of the cell
        indexed by that item for the current sampled time step.
    vm_ave_time : list
        Two-dimensional list of all transmembrane voltages averaged from all
        cell membranes onto cell centres over all sampled time steps, whose:
        . First dimension indexes each sampled time step.
        . Second dimension indexes each cell such that each item is the
          transmembrane voltage spatially situated at the centre of the cell
          indexed by that item for this time step.
        Equivalently, this array is the concatenation of all :attr:`vm_ave`
        arrays over all time steps.
    '''

    # ..................{ INITIALIZORS                      }..................
    # Avoid circular import dependencies.
    @type_check
    def __init__(self, p: 'betse.science.parameters.Parameters') -> None:
        '''
        Initialize this simulation.

        Parameters
        ----------
        p : betse.science.parameters.Parameters
            Current simulation configuration.
        '''

        #FIXME: Default all other instance attributes as well.
        #FIXME: Do we still need the "ignore_ecm" flag? It's never set
        #elsewhere, which is a bit... awkward. Clarion winds of change, unveil!

        # Default all remaining attributes.
        self.ignore_ecm = True

        #FIXME: Shift these variables into the "Parameters" class by:
        #* Renaming "self.savedInit" to "p.init_pickle_filename".
        #* Renaming "self.savedSim"  to "p.sim_pickle_filename".

        # Define data paths for saving an initialization and simulation run:
        self.savedInit = pathnames.join(
            p.init_pickle_dirname, p.init_pickle_basename)
        self.savedSim  = pathnames.join(
            p.sim_pickle_dirname, p.sim_pickle_basename)


    @type_check
    def init_core(self, phase: SimPhase) -> None:
        '''
        Prepare core data structures required by the passed phase in a
        general-purpose manner applicable to *all* possible phases.

        This method initializes core computational matrices -- including those
        concerning intracellular and environmental concentrations, voltages,
        specific diffusion constants, and types of ions included in the
        simulation. This method is performed only once per seed (i.e., cell
        cluster creation) and thus contains crucial parameters, which cannot be
        changed after running an initialization.

        Parameters
        --------
        phase : SimPhase
            Current simulation phase.
        '''

        # Log this attempt.
        logs.log_info('Initializing core simulation matrices...')

        # Localize frequently referenced phase variables for convenience.
        p = phase.p
        cells = phase.cells

        #FIXME: Eliminate this crude hack by refactoring all references to
        #"sim.dyna" throughout the codebase to "phase.dyna" instead. Naturally,
        #this will require refactoring all methods referencing "sim.dyna" to
        #accept a "phase: SimPhase" parameter. After doing so, remove this line.
        #Praise be to the multifoliate rose!
        self.dyna = phase.dyna

        # initialize all extra substances related objects to None, to be filled in if desired later
        self.molecules = None
        self.metabo = None
        self.met_concs = None
        self.grn = None

        self.mdl = len(cells.mem_i)  # mems-data-length
        self.cdl = len(cells.cell_i)  # cells-data-length

        if p.is_ecm:  # set environnment data length
            self.edl = len(cells.xypts)
        else:
            self.edl = self.mdl

        # initialize extra arrays that allow additional substances (defined outside of sim) to affect Vmem:
        self.extra_rho_cells = np.zeros(self.cdl)
        self.extra_rho_mems = np.zeros(self.mdl)
        self.extra_rho_env = np.zeros(self.edl)
        self.extra_J_mem = np.zeros(self.mdl)

        self.extra_Jenv_x = np.zeros(self.edl)
        self.extra_Jenv_y = np.zeros(self.edl)

        self.extra_rho_cells_o = np.zeros(self.cdl)

        self.vgj = np.zeros(self.mdl)


        self.gj_block = 1 # will update this according to user preferences in self.init_dynamics()

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

        self.Jn = np.zeros(self.mdl)  # net normal transmembrane current density at membranes
        self.Jmem = np.zeros(self.mdl)   # current across membrane from intra- to extra-cellular space
        self.Jgj = np.zeros(self.mdl)     # total normal current across GJ coupled membranes from intra- to intra- space
        self.Jc = np.zeros(self.mdl)   # current across membrane from intracellular conductivity
        self.Ec = np.zeros(self.mdl)   # intracellular electric field at membrane
        self.dvm = np.zeros(self.mdl)   # rate of change of Vmem
        self.drho = np.zeros(self.cdl)  # rate of change of charge in cell

        # Current averaged to cell centres:
        self.J_cell_x = np.zeros(self.cdl)
        self.J_cell_y = np.zeros(self.cdl)

        # Membrane current data structure initialization
        self.flx_mem_i = np.zeros(self.mdl)
        self.fluxes_mem = []
        self.I_mem = np.zeros(self.mdl)  # total current across membranes

        self.P_cells = np.zeros(self.cdl)  # initialize pressure in cells

        self.v_cell = np.zeros(self.cdl)  # initialize intracellular (long range) voltage
        self.vm = np.zeros(self.mdl)     # initialize Vmem at membranes
        self.vm_ave = np.zeros(self.cdl)  # initialize average Vmem at cell centres
        self.rho_cells = np.zeros(self.cdl)

        self.E_gj_x = np.zeros(self.mdl)   # electric field components across gap junctions (micro-field)
        self.E_gj_y = np.zeros(self.mdl)

        self.E_cell_x = np.zeros(self.mdl)   # intracellular electric field components (macro-field)
        self.E_cell_y = np.zeros(self.mdl)

        self.TJ_modulator = None # register the tight junction modulator field to a 'None" object

        self.c_env_bound = []  # moving ion concentration at global boundary

        if p.is_ecm:  # special items specific to simulation of extracellular spaces only:
            # vectors storing separate cell and env voltages
            self.v_env = np.zeros(self.edl)   # voltage in the full environment
            self.rho_env = np.zeros(self.edl)  # charge in the full environment
            self.z_array_env = []  # ion valence array matched to env points
            self.D_env = []  # an array of diffusion constants for each ion defined on env grid
            self.Dtj_rel = []  # relative diffusion constants for ions across tight junctions

            self.E_env_x = np.zeros(cells.X.shape)  # electric field in environment, x component
            self.E_env_y = np.zeros(cells.X.shape)  # electric field in environment, y component

        else:  # items specific to simulation *without* extracellular spaces:
            # Initialize environmental volume:
            self.envV = np.zeros(self.mdl)
            self.envV[:] = p.vol_env
            self.v_env = np.zeros(len(cells.xypts))
            self.rho_env = np.zeros(len(cells.xypts))

        ion_names = list(p.ions_dict.keys())

        i = -1  # dynamic index

        # Go through ion list/dictionary and initialize all sim structures.
        for name in ion_names:
            # If this ion is enabled...
            if p.ions_dict[name] == 1:
                i = i+1 # update the dynamic index

                str1 = 'i' + name  # create the ion index

                # Dynamically add this field to the object.
                setattr(self, str1, i)

                self.ionlabel[vars(self)[str1]] = p.ion_long_name[name]

                # if name != 'P':
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
                if p.is_ecm:
                    str_Denv = 'D' + name

                    setattr(self, str_Denv, np.zeros(len(cells.xypts)))
                    vars(self)[str_Denv][:] = p.free_diff[name]

                # ion charge characteristic for intracellular:
                str_z = 'z' + name

                setattr(self, str_z, np.zeros(self.cdl))
                vars(self)[str_z][:] = p.ion_charge[name]

                if p.is_ecm:  # ion charge characteristic for extracellular:
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

                self.c_env_bound.append(p.env_concs[name])

                if p.is_ecm:

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

        if p.is_ecm:  # items specific for extracellular spaces simulation:
            self.z_array_env = np.asarray(self.z_array_env)
            self.D_env = np.asarray(self.D_env)

            # boundary conditions for voltages:
            # voltage (scheduled dynamics might vary these values)
            self.bound_V = {}
            self.bound_V['T'] = 0
            self.bound_V['B'] = 0
            self.bound_V['L'] = 0
            self.bound_V['R'] = 0

            # redo environmental protein handling so that it's only present in the cell cluster:
            # self.c_env_bound[self.iM] += 1*self.c_env_bound[self.iP]
            # self.c_env_bound[self.iP] = 0.0
            #
            # new_P = np.zeros(self.edl)
            # new_P[cells.map_mem2ecm] = self.cc_env[self.iP][cells.map_mem2ecm]
            # new_P = fd.integrator(new_P.reshape(cells.X.shape), 0.5).ravel()
            #
            # change_M = np.ones(self.edl)*self.cc_env[self.iP]
            # change_M[cells.map_mem2ecm] = 0.0
            # change_M = fd.integrator(change_M.reshape(cells.X.shape), 0.5).ravel()
            #
            # self.cc_env[self.iP] = new_P
            # self.cc_env[self.iM] += change_M

            # initialize the environmental diffusion matrix:
            self.initDenv(cells, p)

        # gap junction specific arrays:
        self.id_gj = np.ones(len(cells.mem_i))  # identity array for gap junction indices...
        self.gjopen = np.ones(len(cells.mem_i))*cells.gj_default_weights  # holds gap junction open fraction for each gj
        self.gjl = np.zeros(len(cells.mem_i))  # gj length for each gj
        self.gjl[:] = cells.gj_len

        # GJ fluxes storage vector:
        self.fluxes_gj = np.copy(self.fluxes_mem)

        self.gj_funk = None  # initialize this to None; set in init_dynamics

        # Initialize matrices to store concentration gradient information for each ion:
        self.fluxes_intra = np.zeros(self.fluxes_mem.shape)

        self.cc_at_mem = np.asarray([
            cc[cells.mem_to_cells] for cc in self.cc_cells])

        # load in the microtubules object:
        self.mtubes = Mtubes(self, cells, p)

        self.u_cells_x = 0.0
        self.u_cells_y = 0.0

        self.u_env_x = 0.0
        self.u_env_y = 0.0


    @type_check
    def init_dynamics(self, phase: SimPhase) -> None:
        '''
        Prepare tissue-centric data structures required by the passed phase in a
        general-purpose manner applicable to *all* possible phases.

        This method initializes core computational matrices -- including those
        concerning tissue, cut, and boundary profiles, dynamic activities, and
        optional simulation features safely modifiable between the
        initialization and simulation phases (e.g., electroosmotic fluid flow).

        Parameters
        --------
        phase : SimPhase
            Current simulation phase.
        '''

        # Log this attempt.
        logs.log_info('Initializing tissue and boundary profiles...')

        # Localize frequently referenced phase variables for convenience.
        p = phase.p
        cells = phase.cells

        # smoothing weights for membrane and central values:
        nfrac = p.smooth_cells
        self.smooth_weight_mem = ((nfrac*cells.num_mems[cells.mem_to_cells] -1)/(nfrac*cells.num_mems[cells.mem_to_cells]))
        self.smooth_weight_o = 1/(nfrac*cells.num_mems[cells.mem_to_cells])

        # # load in the gap junction dynamics object:
        if p.v_sensitive_gj:
            self.gj_funk = Gap_Junction(self, cells, p)

        # Initialize diffusion constants for the extracellular transport.
        self.initDenv(cells, p)

        # Re-initialize global boundary fixed concentrations.
        # if p.is_ecm:
        # For each possible ion...
        for key, val in p.ions_dict.items():
            # If this ion is enabled...
            if val == 1:

                ion_i = self.get_ion(key)
                # print("resetting c_env from ", self.c_env_bound[ion_i], 'to ', p.cbnd[key], "for ", key)

                if p.cbnd is not None:
                    self.c_env_bound[ion_i] = p.cbnd[key]

        # Initialize all tissue profiles.
        phase.dyna.tissueProfiles(self, cells, p)

        if p.is_ecm:
            # create a copy-base of the environmental junctions diffusion constants:
            self.D_env_base = copy.copy(self.D_env)

            # initialize current vectors
            self.J_env_x = np.zeros(len(cells.xypts))
            self.J_env_y = np.zeros(len(cells.xypts))

            self.fluxes_env_x = np.zeros((len(self.zs), self.edl))
            self.fluxes_env_y = np.zeros((len(self.zs), self.edl))

        # # Initialize an array structure that will hold user-scheduled changes to membrane permeabilities:
        Dm_cellsA = np.asarray(self.Dm_cells)

        self.Dm_base = np.copy(Dm_cellsA) # make a copy that will serve as the unaffected values base

        self.Dm_scheduled = np.copy(Dm_cellsA)
        self.Dm_scheduled[:] = 0

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

        # initialize additional pump blocks:
        self.CaATP_block = np.ones(self.mdl)  # initialize CaATP blocking vector


        # initialize calcium dynamics if desired:
        if p.ions_dict['Ca'] == 1 and p.Ca_dyn is True:
            self.endo_retic = EndoRetic(self, cells, p)

        else:
            self.endo_retic = None


        # -----auxiliary molecules initialization -------------------------

        # create and initialize the auxiliary-molecules handler for this simulation:
        #(only do these initializations if they haven't been done yet)
        if p.molecules_enabled and self.molecules is None:

            logs.log_info("Initializing general network...")

            # create an instance of the metabolism simulator
            self.molecules = MasterOfMolecules(p)
            # read in the configuration settings for the metabolism simulator:
            self.molecules.read_mol_config(self, cells, p)


        elif p.molecules_enabled and self.molecules is not None:
        # don't declare a whole new object, but re-read in parts that user may have changed:
            logs.log_info("Reinitializing the general regulatory network for simulation...")

            # re-read the config file again and reassign everything except for concentrations,
            #  to capture any user updates:
            self.molecules.reinitialize(self, cells, p)


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

            if p.is_ecm is True:

                # initialize flow vectors:
                self.u_env_x = np.zeros(cells.X.shape)
                self.u_env_y = np.zeros(cells.X.shape)

                # logs.log_info('Creating environmental Poisson solver for fluids...')
                # bdic = {'N': 'flux', 'S': 'flux', 'E': 'flux', 'W': 'flux'}
                # cells.lapENV_P, cells.lapENV_P_inv = cells.grid_obj.makeLaplacian(bound=bdic)
                #
                # cells.lapENV_P = None  # get rid of the non-inverse matrix as it only hogs memory...

        if p.deformation:  # if user desires deformation:
            # initialize vectors for potential deformation:
            self.d_cells_x = np.zeros(self.cdl)
            self.d_cells_y = np.zeros(self.cdl)

            cells.deform_tools(p)

            #FIXME: "cellso" is pretty... wierd. Do we still need this and, if
            #so, can we better document why? If we do ultimately need this,
            #we'll probably want to privatize this attribute to the less
            #ambiguous name "_cells_deformed". External callers should never
            #access or care about this.

            # Deepy copy of the current cell cluster, to isolate deformations to
            # for visualization purposes only.
            self.cellso = copy.deepcopy(cells)

            #FIXME: "td_deform" is currently forced to "False", implying this
            #branch to currently reduce to a noop. Is this still desired?
            if p.td_deform and (cells.lapGJ is None or cells.lapGJ_P is None):
                # Make a laplacian and solver for discrete transfers on closed,
                # irregular cell network.
                logs.log_info('Creating cell network Poisson solver...')
                cells.graphLaplacian(p)
        else:
            self.cellso = cells

        # if simulating electrodiffusive movement of membrane pumps and channels:-------------
        # if p.sim_eosmosis is True:
        #
        #     # if cells.gradMem is None:
        #     #     logs.log_info("Creating tools for self-electrodiffusion of membrane pumps and channels.")
        #     #     cells.eosmo_tools(p)
        #
        #     self.rho_pump = np.ones(self.mdl)
        #     self.rho_channel = np.ones(self.mdl)
        #
        #     self.move_pumps_channels = MoveChannel(self, cells, p)

        # else:
        #     self.rho_pump = 1  # else just define it as identity.
        #     self.rho_channel = 1
        #
        #     self.move_pumps_channels = None
        # FIXME: eventually remove these entirely
        self.rho_pump = 1
        self.rho_channel = 1

        # Initialize core user-specified interventions.
        phase.dyna.runAllInit(self, cells, p)

        # update the microtubules dipole for the case user changed it between init and sim:
        self.mtubes.reinit(cells, p)

        # calculate an inverse electrical double layer based on internal and external concentrations:
        self.ko_env = (np.sqrt(np.dot((p.NAv * (p.q ** 2) * self.zs ** 2) / (p.er * p.eo * p.kb * p.T), self.cc_env))).mean()
        self.ko_cell = (np.sqrt(np.dot((p.NAv * (p.q ** 2) * self.zs ** 2) / (p.er * p.eo * p.kb * p.T), self.cc_cells))).mean()

        self.cedl_env = p.er*p.eo*self.ko_env
        self.cedl_cell = p.er*p.eo*self.ko_cell

        self.cedl = p.er*p.eo*self.ko_env  # electrical double layer capacitance

        # calculate a basic system conductivity:
        # self.sigma = np.asarray([((z**2)*p.q*p.F*cc*D)/(p.kb*p.T) for
        #                          (z, cc, D) in zip(self.zs, self.cc_cells, self.D_free)]).mean()

        self.sigma = (np.asarray([((z**2)*p.q*p.F*cc*D)/(p.kb*p.T) for
                                 (z, cc, D) in zip(self.zs, self.cc_cells, self.D_free)]).sum(axis = 0)).mean()

        # calculate specific maps of conductivity in cells and environment
        # conductivity of cells is assumed to be 0.1x lower than that of free environment:
        self.sigma_cell = np.asarray([((z ** 2) * p.q * p.F * cc * D * 0.1) / (p.kb * p.T) for (z, cc, D) in
                                 zip(self.zs, self.cc_cells, self.D_free)]).sum(axis=0)
        # conductivity map for environment:
        self.sigma_env = np.asarray(
            [((z ** 2) * p.q * p.F * cc * D) / (p.kb * p.T) for (z, cc, D) in zip(self.zs, self.cc_env, self.D_env)]).sum(
            axis=0)

        if p.is_ecm:

            # environmental conductivity matrix needs to be smoothed to assured simulation stability:
            self.sigma_env = gaussian_filter(self.sigma_env.reshape(cells.X.shape), 1).ravel()

        # calculate a geometric resistivity factor:
        self.rho_factor = cells.mem_sa/cells.R_rads

        # capacitance of gj
        # self.cgj = 1 / ((2 / self.cedl_cell) + (2 / p.cm))
        self.cgj = 1/(2/p.cm)

        # If this is the fast BETSE solver, initialize this solver.
        if p.solver_type is SolverType.FAST:
            self.fast_sim_init(cells, p)

    # ..................{ SOLVERS                           }..................
    @type_check
    def run_sim_core(self, phase: SimPhase) -> None:
        '''
        Perform the passed simulation phase (e.g., initialization, simulation),
        pickling the results to files defined by the configuration associated
        with this phase.

        Parameters
        --------
        phase : SimPhase
            Current simulation phase.
        '''

        # Initialize all structures used for gap junctions, ion channels, and
        # other dynamics.
        self.init_dynamics(phase)

        # Reinitialize all time-data structures
        self.clear_storage(phase.cells, phase.p)

        # Get the net, unbalanced charge and corresponding voltage in each cell
        # to initialize values of voltages.
        self.update_V(phase.cells, phase.p)

        # Calculate the following simulation phase-specific locals:
        #
        # * "time_steps", the array of all time steps for this phase.
        # * "time_steps_sampled", this array resampled to reduce data storage.
        # * "solver_context", the context manager intended to contextualize the
        #   core time loop for this phase.
        time_steps, time_steps_sampled, solver_context = self._plot_loop(phase)

        # Notify the caller of the range of work performed by this subcommand.
        # The phase.callbacks.progressed() callback is called exactly once for
        # each sampled time step, implying the maximum progress value to be
        # equal to the total number of sampled time steps.
        phase.callbacks.progress_ranged(progress_max=len(time_steps_sampled))

        # Exception raised if this simulation becomes unstable, enabling safe
        # handling of this instability (e.g., by saving simulation results).
        exception_instability = None

        # Attempt to...
        try:
            # If this is the full BETSE solver, set appropriate locals.
            if phase.p.solver_type is SolverType.FULL:
                solver_label = 'Full BETSE simulator'
                solver_method = self._run_sim_core_loop
            # Else if this is the fast BETSE solver, set appropriate locals.
            elif phase.p.solver_type is SolverType.FAST:
                solver_label = 'Fast (equivalent circuit) simulator'
                solver_method = self._run_fast_sim_core_loop
            # Else, this solver is unrecognized. Raise an exception.
            else:
                raise BetseSimException(
                    'Solver type "{}" unrecognized.'.format(
                        phase.p.solver_type))

            # Log this solver type.
            logs.log_info('Solver: %s in use.', solver_label)

            # Perform the time loop for this simulation phase.
            with solver_context:
                solver_method(
                    phase=phase,
                    time_steps=time_steps,
                    time_steps_sampled=time_steps_sampled,

                    #FIXME: Horrible hack. Ideally, we instead want to:
                    #
                    #* Define a "AnimCellsWhileSolvingNoop" subclass such that:
                    #  * The plot_frame() method simply reduces to "pass".
                    #  * Like the "AnimCellsWhileSolving" subclass, the
                    #    "AnimCellsWhileSolvingNoop" subclass should also
                    #    satisfy the context manager API (but by doing
                    #    nothing).
                    #* Restructure the "AnimCellsABC" class hierarchy to
                    #  support this subclass.
                    #* Refactor the _plot_loop() method to return an instance
                    #  of the "AnimCellsWhileSolvingNoop" subclass rather than
                    #  the noop() context manager.
                    anim_cells=(solver_context if isinstance(
                        solver_context, AnimCellsWhileSolving) else None),
                )
        # If this phase becomes computationally unstable...
        except BetseSimUnstableException as exception:
            # Log this instability *BEFORE* logging a report and reraising this
            # exception, improving readability.
            logs.log_error(
                'Simulation halted prematurely '
                'due to computational instability.')

            # Preserve this exception *BEFORE* writing results to disk and
            # reraising this exception. This guarantees access to results even
            # in the case of computational instability, preventing data loss.
            exception_instability = exception
        # If any other type of exception is raised, an unexpected fatal error
        # has occurred. In this case, these results are likely to be in an
        # inconsistent, nonsensical state and hence safely discarded.

        # Save this initialization or simulation and report results of
        # potential interest to the user.
        self.save_and_report(phase.cells, phase.p)

        # If the simulation went unstable, inform the user and reraise the
        # previously raised exception to preserve the underlying cause. To
        # avoid data loss, this exception is raised *AFTER* all pertinent
        # simulation data has been saved to disk.
        if exception_instability is not None:
            raise exception_instability

    # ..................{ SOLVERS ~ full                    }..................
    @type_check
    def _run_sim_core_loop(
        self,
        phase: SimPhase,
        time_steps: ndarray,
        time_steps_sampled: set,
        anim_cells: (AnimCellsWhileSolving, NoneType),
    ) -> None:
        '''
        Drive the time loop for the current simulation phase, including:

        * Gap-junction connections.
        * Calls to all dynamic channels and entities (when simulating).

        Parameters
        --------
        phase : SimPhase
            Current simulation phase.
        time_steps : ndarray
            One-dimensional Numpy array defining the time-steps vector for the
            current phase.
        time_steps_sampled : set
            Subset of the ``time_steps`` array whose elements are **sampled
            time steps** (i.e., time step at which to sample data,
            substantially reducing data storage). In particular, the length of
            this set governs the number of frames in each exported animation.
        anim_cells : (AnimCellsWhileSolving, NoneType)
            A mid-simulation animation of cell voltage as a function of time if
            enabled by this configuration *or* ``None`` otherwise.
        '''

        # Localize frequently accessed variables for efficiency when iterating.
        p = phase.p
        cells = phase.cells

        # True only on the first time step of this phase.
        is_time_step_first = True

        for t in time_steps:  # run through the loop
            # Start the timer to approximate time for the simulation.
            if is_time_step_first:
                loop_measure = time.time()

            # Reinitialize flux storage devices.
            self.fluxes_mem.fill(0)
            self.fluxes_gj.fill(0)

            if p.is_ecm:
                self.fluxes_env_x = np.zeros((len(self.zs), self.edl))
                self.fluxes_env_y = np.zeros((len(self.zs), self.edl))
                # self.Phi_vect = np.zeros((len(self.zs), self.edl))
                # self.conc_J_x = np.zeros(self.edl)
                # self.conc_J_y = np.zeros(self.edl)

            # Calculate the values of scheduled and dynamic quantities (e.g..
            # ion channel multipliers).
            if p.run_sim:
                phase.dyna.runAllDynamics(self, cells, p, t)

            # -----------------PUMPS-------------------------------------------------------------------------------------

            # have the pump run only if the rate constant is larger than 0.0 (so people can shut it off):

            if p.alpha_NaK == 0.0:
                self.rate_NaKATP = np.zeros(self.mdl)

            if p.alpha_NaK > 0.0:
                if p.is_ecm:
                    # run the Na-K-ATPase pump:
                    fNa_NaK, fK_NaK, self.rate_NaKATP = stb.pumpNaKATP(
                        self.cc_at_mem[self.iNa],
                        self.cc_env[self.iNa][cells.map_mem2ecm],
                        self.cc_at_mem[self.iK],
                        self.cc_env[self.iK][cells.map_mem2ecm],
                        self.vm,
                        self.T,
                        p,
                        self.NaKATP_block,
                        met = self.met_concs
                    )

                else:
                    fNa_NaK, fK_NaK, self.rate_NaKATP = stb.pumpNaKATP(
                                self.cc_at_mem[self.iNa],
                                self.cc_env[self.iNa],
                                self.cc_at_mem[self.iK],
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

                # modify the fluxes by electrodiffusive membrane redistribution factor and add fluxes to storage:
                self.fluxes_mem[self.iNa] +=  fNa_NaK
                self.fluxes_mem[self.iK] += fK_NaK

                # update the concentrations of Na and K in cells and environment:
                # self.cc_cells[self.iNa], self.cc_at_mem[self.iNa], self.cc_env[self.iNa] =  stb.update_Co(
                #                                                             self, self.cc_cells[self.iNa],
                #                                                             self.cc_at_mem[self.iNa],
                #                                                             self.cc_env[self.iNa],fNa_NaK, cells, p,
                #                                                             ignoreECM = self.ignore_ecm)
                #
                # self.cc_cells[self.iK], self.cc_at_mem[self.iK], self.cc_env[self.iK] = stb.update_Co(
                #                                                              self, self.cc_cells[self.iK],
                #                                                              self.cc_at_mem[self.iK],
                #                                                              self.cc_env[self.iK], fK_NaK,
                #                                                              cells, p, ignoreECM = self.ignore_ecm)


            # ----------------ELECTRODIFFUSION---------------------------------------------------------------------------

            # electro-diffuse all ions (except for proteins, which don't move) across the cell membrane:

            # shuffle(self.movingIons)

            for i in self.movingIons:

                IdM = np.ones(self.mdl)

                if p.is_ecm:

                    f_ED = stb.electroflux(self.cc_env[i][cells.map_mem2ecm], self.cc_at_mem[i],
                        self.Dm_cells[i], IdM*p.tm, self.zs[i]*IdM, self.vm, self.T, p,
                        rho=self.rho_channel)

                else:

                    f_ED = stb.electroflux(self.cc_env[i],self.cc_at_mem[i],
                                            self.Dm_cells[i],IdM*p.tm,self.zs[i]*IdM,
                                    self.vm,self.T,p,rho=self.rho_channel)

                if p.cluster_open is False:
                    f_ED[cells.bflags_mems] = 0

                # add membrane flux to storage
                self.fluxes_mem[i] += f_ED

                # # update ion concentrations in cell and ecm:
                # self.cc_cells[i], self.cc_at_mem[i], self.cc_env[i] = stb.update_Co(self, self.cc_cells[i],
                #                                                                     self.cc_at_mem[i],
                #                                                                     self.cc_env[i], f_ED,
                #                                                                     cells, p,
                #                                                                     ignoreECM = self.ignore_ecm)

                # update flux between cells due to gap junctions
                self.update_gj(cells, p, t, i)

                if p.is_ecm:
                    #update concentrations in the extracellular spaces:
                    self.update_ecm(cells, p, t, i)

                # update concentration gradient to estimate concentrations at membranes:
                self.update_intra(cells, p, i)



            # ----transport and handling of special ions------------------------------------------------------------

            if p.ions_dict['Ca'] == 1:
                self.ca_handler(cells, p)

            # if p.ions_dict['H'] == 1:
            #
            #     self.acid_handler(cells, p)

            # update the microtubules:------------------------------------------------------------------------------

            if p.use_microtubules:
                self.mtubes.update_mtubes(cells, self, p)

            # update the general molecules handler-----------------------------------------------------------------
            if p.molecules_enabled:

                self.molecules.core.clear_run_loop(self)

                if self.molecules.transporters:
                    self.molecules.core.run_loop_transporters(t, self, cells, p)

                if self.molecules.channels:
                    self.molecules.core.run_loop_channels(phase)

                if self.molecules.modulators:
                    self.molecules.core.run_loop_modulators(self, cells, p)

                self.molecules.core.run_loop(t, self, cells, p)

            # update gene regulatory network handler--------------------------------------------------------

            if p.grn_enabled:

                self.grn.core.clear_run_loop(self)

                if self.grn.transporters:
                    self.grn.core.run_loop_transporters(t, self, cells, p)

                if self.grn.channels:
                    self.grn.core.run_loop_channels(phase)

                if self.grn.modulators:
                    self.grn.core.run_loop_modulators(self, cells, p)

                # update the main gene regulatory network:
                self.grn.core.run_loop(t, self, cells, p)


            # dynamic noise handling-----------------------------------------------------------------------------------

            if p.dynamic_noise == 1 and p.ions_dict['P'] == 1:

                # add a random walk on protein concentration to generate dynamic noise:
                self.protein_noise_flux = p.dynamic_noise_level * (np.random.random(self.mdl) - 0.5)

                # update the concentration of P in cells and environment:
                self.cc_cells[self.iP], self.cc_at_mem[self.iP], self.cc_env[self.iP] = stb.update_Co(self,
                                                        self.cc_cells[self.iP], self.cc_at_mem[self.iP],
                                                        self.cc_env[self.iP], self.protein_noise_flux, cells,
                                                        p, ignoreECM = self.ignore_ecm)

            #-----forces, fields, and flow-----------------------------------------------------------------------------

            # calculate specific forces and pressures:

            if p.deform_osmo is True:

                osmotic_P(self,cells, p)

            if p.fluid_flow is True:

                # self.run_sim = True

                getFlow(self,cells, p)

            # if desired, electroosmosis of membrane channels
            # if p.sim_eosmosis is True:
            #
            #     # self.run_sim = True
            #     self.move_pumps_channels.run(self, cells, p)



            if p.deformation is True:

                # self.run_sim = True

                if p.td_deform is False:

                    getDeformation(self,cells, t, p)

                elif p.td_deform is True:

                    timeDeform(self,cells, t, p)

            # Use fluxes to update all concentrations in the cells
            self.update_all_concs(cells, p)

            # recalculate the net, unbalanced charge and voltage in each cell:
            self.update_V(cells, p)

            # check for NaNs in voltage and stop simulation if found:
            stb.check_v(self.vm)


            # ---------time sampling and data storage---------------------------------------------------
            # If this time step is sampled...
            if t in time_steps_sampled:
                # Notify the caller that an additional sampled time step has
                # been successfully simulated.
                phase.callbacks.progressed_next()

                # Write data to time storage vectors.
                self.write2storage(t, cells, p)

                # If animating this phase, display and/or save the next frame
                # of this animation. For simplicity, pass "-1" implying the
                # last frame and hence the results of the most recently solved
                # time step.
                if anim_cells is not None:
                    anim_cells.plot_frame(time_step=-1)

            # Get time for loop and estimate total time for simulation.
            if is_time_step_first:
                # Ignore this conditional on all subsequent time steps.
                is_time_step_first = False

                loop_time = time.time() - loop_measure

                # Estimated number of seconds to complete this phase.
                if p.run_sim:
                    time_estimate = round(loop_time * p.sim_tsteps, 2)
                else:
                    time_estimate = round(loop_time * p.init_tsteps, 2)

                # Log this estimate.
                logs.log_info(
                    'This run should take approximately %fs to compute...',
                    time_estimate)

    # ..................{ SOLVERS ~ fast                    }..................
    def fast_sim_init(self, cells, p):
        '''
        Special initialization required for fast (equivalent circuit) sims.
        '''

        self.rev_E_dic = {}  # dictionary of reversal potentials
        self.cbar_dic = {}
        sigma_gj = []  # temporary array of gap junction conductivities
        sigma_mem = []  # temporary array of "leak" membrane conductivities

        for ion_n, ion_i in p.ions_dict.items():  # for each ion
            if ion_i == 1:  # if it's used in the simulation
                ii = self.get_ion(ion_n)  # get the index
                ccell = self.cc_cells[ii].mean()

                if p.is_ecm:
                    cenv = self.cc_env[ii].mean()
                else:
                    cenv = self.cc_env[ii]

                cbar = (ccell + cenv) / 2
                self.cbar_dic[ion_n] = cbar

                # calculate the reversal potential using the Nernst Equation:
                revE = ((p.R * self.T) / (self.zs[ii] * p.F)) * np.log(cenv / ccell)

                # store the reversal potential in the dictionary
                self.rev_E_dic[ion_n] = revE

                # estimate gap junction conductivity
                # sigma_gj.append((self.D_free[ii] * p.gj_surface * p.F * ccell *self.zs[ii]**2) / (p.cell_space * p.R * p.T))
                # sigma_mem.append((self.Dm_cells[ii]*p.F*cbar*self.zs[ii]**2)/(p.tm*p.R*p.T))

                sigma_gj.append((self.D_free[ii]* p.q * p.gj_surface * p.F * ccell *self.zs[ii]**2) / (p.cell_space * p.kb * p.T))
                sigma_mem.append((self.Dm_cells[ii]*p.q*p.F*cbar*self.zs[ii]**2)/(p.tm*p.kb*p.T))

        # get the leak channel reversibility represented by the resting potential of the base membrane:
        stb.ghk_calculator(self, cells, p)

        # get the conversion for geometry of the cluster (required to convert to conductivity):
        # self.geo_conv = (cells.cell_sa/ np.dot(cells.M_sum_mems, cells.mem_sa))*cells.num_nn
        self.geo_conv = 1.0

        self.E_Leak = self.vm_GHK

        self.sigma_mem = sigma_mem

        self.cbar_all = np.mean([v for k, v in self.cbar_dic.items()])
        self.cbar_sum = np.sum([v.mean() for k, v in self.cbar_dic.items()])

        self.G_Leak = (np.dot(cells.M_sum_mems, sum(sigma_mem)*cells.mem_sa)/cells.cell_sa)*self.geo_conv

        # get the average gap junction conductivity:
        # self.G_gj = sum(sigma_gj)*self.geo_conv*(cells.mem_sa.mean()/cells.cell_sa.mean())
        self.G_gj = sum(sigma_gj) * self.geo_conv*(1/cells.num_mems)

        self.Emx = np.zeros(self.mdl)
        self.Emy = np.zeros(self.mdl)


    @type_check
    def _run_fast_sim_core_loop(
        self,
        phase: SimPhase,
        time_steps: ndarray,
        time_steps_sampled: set,
        anim_cells: (AnimCellsWhileSolving, NoneType),
    ) -> None:
        '''
        Drive the time loop for the simulation phase using equivalent circuit
        formalism.

        Parameters
        --------
        phase : SimPhase
            Current simulation phase.
        time_steps : ndarray
            One-dimensional Numpy array defining the time-steps vector for the
            current phase.
        time_steps_sampled : set
            Subset of the ``time_steps`` array whose elements are **sampled
            time steps** (i.e., time step at which to sample data,
            substantially reducing data storage). In particular, the length of
            this set governs the number of frames in each exported animation.
        anim_cells : (AnimCellsWhileSolving, NoneType)
            A mid-simulation animation of cell voltage as a function of time if
            enabled by this configuration *or* ``None`` otherwise.
        '''

        # Localize frequently-accessed variables for efficiency when iterating.
        p = phase.p
        cells = phase.cells

        # True only on the first time step of this phase.
        is_time_step_first = True

        for t in time_steps:  # run through the loop
            # Start the timer to approximate time for the simulation.
            if is_time_step_first:
                loop_measure = time.time()

            # Reinitialize flux storage devices.
            self.fluxes_mem.fill(0)
            self.fluxes_gj.fill(0)

            # Calculate the values of scheduled and dynamic quantities (e.g..
            # ion channel multipliers).
            if p.run_sim:
                phase.dyna.runAllDynamics(self, cells, p, t)

            # update the microtubules:------------------------------------------------------------------------------

            if p.use_microtubules:
                self.mtubes.update_mtubes(cells, self, p)

            # update the general molecules handler-----------------------------------------------------------------
            if p.molecules_enabled:

                self.molecules.core.clear_run_loop(self)

                if self.molecules.transporters:
                    self.molecules.core.run_loop_transporters(t, self, cells, p)

                if self.molecules.channels:
                    self.molecules.core.run_fast_loop_channels(phase)

                if self.molecules.modulators:
                    self.molecules.core.run_loop_modulators(self, cells, p)

                self.molecules.core.run_loop(t, self, cells, p)

            # update gene regulatory network handler--------------------------------------------------------

            if p.grn_enabled:

                self.grn.core.clear_run_loop(self)

                if self.grn.transporters:
                    self.grn.core.run_loop_transporters(t, self, cells, p)

                if self.grn.channels:
                    self.grn.core.run_fast_loop_channels(phase)

                if self.grn.modulators:
                    self.grn.core.run_loop_modulators(self, cells, p)

                # update the main gene regulatory network:
                self.grn.core.run_loop(t, self, cells, p)

            # Update gap junctions:
            self.vgj = self.vm_ave[cells.cell_nn_i[:, 1]] - self.vm_ave[cells.cell_nn_i[:, 0]]

            if p.v_sensitive_gj is True:

                # run the gap junction dynamics object to update gj open state of sim:
                self.gj_funk.run(self, cells, p)

            else:
                self.gjopen = self.gj_block*np.ones(len(cells.mem_i))*cells.gj_default_weights

            Jgj = self.G_gj*np.dot(cells.M_sum_mems, self.vgj)

            Jmem = np.dot(cells.M_sum_mems, self.extra_J_mem*cells.mem_sa)/cells.cell_sa

            self.vm_ave += p.dt*(1/p.cm)*(Jgj - Jmem - self.G_Leak*(self.vm_ave - self.E_Leak))

            self.vm = self.vm_ave[cells.mem_to_cells]

            # Currents:
            Jtot = -self.vgj*self.G_gj[cells.mem_to_cells] + self.extra_J_mem

            self.Jn = Jtot

            Jmx = self.Jn*cells.mem_vects_flat[:,2]
            Jmy = self.Jn*cells.mem_vects_flat[:,3]

            self.Emx, self.Emy = cells.single_cell_div_free(Jmx/(0.1*self.sigma_cell.mean()), Jmy/(0.1*self.sigma_cell.mean()))

            Jcx = self.Jn * cells.mem_vects_flat[:, 2]
            Jcy = self.Jn * cells.mem_vects_flat[:, 3]

            # average intracellular current to cell centres
            self.J_cell_x = np.dot(cells.M_sum_mems, Jcx * cells.mem_sa) / cells.cell_sa
            self.J_cell_y = np.dot(cells.M_sum_mems, Jcy * cells.mem_sa) / cells.cell_sa

            # intracellular electric field:
            self.E_cell_x = self.J_cell_x / (0.1 * self.sigma_cell)
            self.E_cell_y = self.J_cell_y / (0.1 * self.sigma_cell)

            # # calculate electric field in cells using net intracellular current and cytosol conductivity:
            # self.Emc = (self.E_cell_x[cells.mem_to_cells] * cells.mem_vects_flat[:, 2] +
            #            self.E_cell_y[cells.mem_to_cells] * cells.mem_vects_flat[:, 3])

            # check for NaNs in voltage and stop simulation if found:
            stb.check_v(self.vm_ave)

            # ---------time sampling and data storage---------------------------------------------------
            # If this time step is sampled...
            if t in time_steps_sampled:
                # Notify the caller that an additional sampled time step has
                # been successfully simulated.
                phase.callbacks.progressed_next()

                # Write data to time storage vectors.
                self.vm_time.append(self.vm * 1)

                # microtubules:
                self.mtubes_x_time.append(self.mtubes.mtubes_x * 1)
                self.mtubes_y_time.append(self.mtubes.mtubes_y * 1)

                self.I_cell_x_time.append(self.J_cell_x * 1)
                self.I_cell_y_time.append(self.J_cell_y * 1)

                self.efield_gj_x_time.append(self.E_cell_x[cells.mem_to_cells] * 1)
                self.efield_gj_y_time.append(self.E_cell_y[cells.mem_to_cells] * 1)

                self.gjopen_time.append(self.gjopen*1)

                self.time.append(t * 1)

                if p.molecules_enabled:
                    self.molecules.core.write_data(self, cells, p)
                    self.molecules.core.report(self, p)

                if p.grn_enabled:
                    self.grn.core.write_data(self, cells, p)
                    self.grn.core.report(self, p)

                self.vm_ave_time.append(self.vm_ave*1)

                # If animating this phase, display and/or save the next frame
                # of this animation. For simplicity, pass "-1" implying the
                # last frame and hence the results of the most recently solved
                # time step.
                if anim_cells is not None:
                    anim_cells.plot_frame(time_step=-1)

            # Get time for loop and estimate total time for simulation.
            if is_time_step_first:
                # Ignore this conditional on all subsequent time steps.
                is_time_step_first = False

                loop_time = time.time() - loop_measure

                # Estimated number of seconds to complete this phase.
                if p.run_sim:
                    time_estimate = round(loop_time * p.sim_tsteps, 2)
                else:
                    time_estimate = round(loop_time * p.init_tsteps, 2)

                # Log this estimate.
                logs.log_info(
                    'This run should take approximately %fs to compute...',
                    time_estimate)

    #.................{ FINALIZERS                  }..........................
    def clear_storage(self, cells, p):
        '''
        Re-initializes time storage vectors at the begining of a sim or init.
        '''

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

        self.vm_ave_time = []

        self.venv_time = []

        # microtubules:
        self.mtubes_x_time = []
        self.mtubes_y_time = []

        if p.deformation:
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

            self.molecules.core.clear_cache()

        if p.grn_enabled:

            self.grn.core.clear_cache()

        # if p.sim_eosmosis is True:
        #     self.rho_channel_time = []
        #     self.rho_pump_time = []

        if p.Ca_dyn is True:
            self.endo_retic.clear_cache()

        if p.is_ecm is True:

            self.charge_env_time = []

            self.efield_ecm_x_time = []   # matrices storing smooth electric field in ecm
            self.efield_ecm_y_time = []

            # initialize time-storage vectors for electroosmotic data:
            self.u_env_x_time = []
            self.u_env_y_time = []

            self.rho_pump_time = []    # store pump and channel states as function of time...
            self.rho_channel_time = []


    def write2storage(self,t,cells,p):
        '''
        Append each multidimensional Numpy array covering all time steps (e.g.,
        :attr:`cc_env_time`) with the corresponding Numpy array of fewer
        dimensions specific to the passed time step (e.g., :attr:`cc_env`).
        '''

        if p.GHK_calc:
            stb.ghk_calculator(self,cells,p)
            self.vm_GHK_time.append(self.vm_GHK) # data array holding GHK vm estimates

        # add the new concentration and voltage data to the time-storage matrices:
        self.efield_gj_x_time.append(self.E_gj_x*1)
        self.efield_gj_y_time.append(self.E_gj_y*1)

        concs = np.copy(self.cc_cells)
        self.cc_time.append(concs)
        concs = None

        envsc = np.copy(self.cc_env)
        self.cc_env_time.append(envsc)
        envsc = None

        ddc = np.copy(self.Dm_cells)
        ddc.tolist()
        self.dd_time.append(ddc)
        ddc = None

        self.I_cell_x_time.append(self.J_cell_x*1)
        self.I_cell_y_time.append(self.J_cell_y*1)

        self.I_mem_time.append(self.I_mem*1)

        self.vm_time.append(self.vm*1)

        self.rho_cells_time.append(self.rho_cells*1)
        self.rate_NaKATP_time.append(self.rate_NaKATP*1)
        self.P_cells_time.append(self.P_cells)

        self.venv_time.append(self.v_env * 1)

        if p.deform_osmo:
            self.osmo_P_delta_time.append(self.osmo_P_delta)

        # microtubules:
        self.mtubes_x_time.append(self.mtubes.mtubes_x*1)
        self.mtubes_y_time.append(self.mtubes.mtubes_y*1)

        if p.deformation:
            # make a copy of cells to apply deformation to:
            # self.cellso = copy.deepcopy(cells)
            implement_deform_timestep(self, self.cellso, t, p)
            self.dx_cell_time.append(self.d_cells_x*1)
            self.dy_cell_time.append(self.d_cells_y*1)

        if p.fluid_flow:
            self.u_cells_x_time.append(self.u_cells_x*1)
            self.u_cells_y_time.append(self.u_cells_y*1)

        # if p.sim_eosmosis:
        #     self.rho_channel_time.append(self.rho_channel*1)
        #     self.rho_pump_time.append(self.rho_pump*1)

        self.gjopen_time.append(self.gjopen*1)
        self.time.append(t)

        if p.molecules_enabled:
            self.molecules.core.write_data(self, cells, p)
            self.molecules.core.report(self, p)

        if p.grn_enabled:
            self.grn.core.write_data(self, cells, p)
            self.grn.core.report(self, p)

        if p.Ca_dyn == 1 and p.ions_dict['Ca'] == 1:
            self.endo_retic.write_cache(self)

        self.I_tot_x_time.append(self.J_env_x*1)
        self.I_tot_y_time.append(self.J_env_y*1)

        if p.is_ecm:
            self.efield_ecm_x_time.append(self.E_env_x*1)
            self.efield_ecm_y_time.append(self.E_env_y*1)

            if p.fluid_flow:
                self.u_env_x_time.append(self.u_env_x*1)
                self.u_env_y_time.append(self.u_env_y*1)

        self.vm_ave_time.append(self.vm_ave)

        # # magnetic vector potential:
        # self.Ax_time.append(self.Ax)
        # self.Ay_time.append(self.Ay)

        # magnetic field
        # self.Bz_time.append(self.Bz)


    def save_and_report(self, cells, p) -> None:
        '''
        Save the results of running the current phase (e.g., initialization,
        simulation).
        '''

        cells.points_tree = None

        #FIXME: What is this? Why do we need an extra copy of "cells", anyway?
        # get rid of the extra copy of cells
        if p.deformation:
            cells = copy.deepcopy(self.cellso)

        self.cellso = None

        if p.run_sim is False:
            datadump = [self, cells, p]
            fh.saveSim(self.savedInit, datadump)
            logs.log_info('Initialization saved to "%s".', p.init_pickle_dirname)
        else:
            datadump = [self, cells, p]
            fh.saveSim(self.savedSim, datadump)
            logs.log_info('Simulation saved to "%s".', p.sim_pickle_dirname)

        final_vmean = 1000 * np.round(np.mean(self.vm_time[-1]), 6)
        logs.log_info('Final average cell Vmem: %g mV', final_vmean)

        # If this is the full BETSE solver...
        if p.solver_type is SolverType.FULL:
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

            if p.GHK_calc:
                final_vmean_GHK = 1000*np.round(np.mean(self.vm_GHK_time[-1]),6)
                logs.log_info(
                    'Final average cell Vmem calculated using GHK: %s mV',
                    final_vmean_GHK)

        if p.molecules_enabled:
            self.molecules.core.report(self, p)

        if p.grn_enabled:
            self.grn.core.report(self, p)


    def sim_info_report(self,cells,p):

        logs.log_info('This world contains ' + str(cells.cell_number) + ' cells.')
        logs.log_info('Each cell has an average of ' + str(round(cells.average_nn, 2)) + ' nearest-neighbours.')

        logs.log_info(
            'You are running the ion profile: %s', str(p.ion_profile).lower())
        logs.log_info('Ions in this simulation: %s', str(self.ionlabel))
        logs.log_info(
            'If you have selected features using other ions, '
            'they will be ignored.')

        logs.log_info('Considering extracellular spaces: ' + str(p.is_ecm))

        if p.is_ecm:
            logs.log_info('Cells per env grid square ' + str(round(cells.ratio_cell2ecm, 2)))

        logs.log_info('Electroosmotic fluid flow: ' + str(p.fluid_flow))
        # logs.log_info('Ion pump and channel electodiffusion in membrane: ' + str(p.sim_eosmosis))
        logs.log_info('Force-induced cell deformation: ' + str(p.deformation))
        logs.log_info('Osmotic pressure: ' + str(p.deform_osmo))
        # logs.log_info('Electrostatic pressure: ' + str(p.deform_electro))

        if p.molecules_enabled:
            logs.log_info(
                'Auxiliary molecules and properties are enabled from '
                '"General Networks" section of main config file.')

        if p.grn_enabled:
            logs.log_info(
                'A gene regulatory network is being simulated using file: %s',
                p.grn_config_filename)

    #.................{  DOOERs & GETTERS  }............................................

    def update_V(self,cells,p):

        # save the voltage as a placeholder:
        vmo = self.vm*1

        # get the currents and in-cell and environmental voltages:
        get_current(self, cells, p)

        # conductivity of cells:
        self.sigma_cell = np.asarray([((z ** 2) * (p.F**2) * cc * D * 0.1) / (p.R * p.T) for (z, cc, D) in
                                     zip(self.zs, self.cc_cells, self.D_free)]).mean(axis=0)


        if p.cell_polarizability == 0.0:  # allow users to have "simple" case behaviour

            # change in charge density at the membrane:
            # Jm = np.dot(cells.M_sum_mems, self.Jn*cells.mem_sa)/cells.cell_sa
            # self.vm += -(1/p.cm)*Jm[cells.mem_to_cells]*p.dt

            # In terms of intra and extracellular charge:
            rho_surf = self.rho_cells * cells.diviterm

            self.vm = (1/p.cm)*rho_surf[cells.mem_to_cells]

            # if p.is_ecm:
            #
            #     self.vm += -self.v_env[cells.map_mem2ecm]

            # average vm:
            # self.vm_ave = np.dot(cells.M_sum_mems, self.vm*cells.mem_sa)/cells.cell_sa
            self.vm_ave = np.dot(cells.M_sum_mems, self.vm) / cells.num_mems

            self.E_cell_x = self.J_cell_x/(self.sigma_cell)
            self.E_cell_y = self.J_cell_y/(self.sigma_cell)

            # calculate electric field in cells using net intracellular current and cytosol conductivity:
            self.Emc = (self.E_cell_x[cells.mem_to_cells] * cells.mem_vects_flat[:, 2] +
                       self.E_cell_y[cells.mem_to_cells] * cells.mem_vects_flat[:, 3])


        else:

            # Calculate the central voltage value in the cell from total change in cell charge:
            rho_surf = self.rho_cells * cells.diviterm

            # if p.is_ecm:
            #     vm_o = (1/p.cm)*rho_surf[cells.mem_to_cells] - self.v_env[cells.map_mem2ecm]
            #
            # else:
            vm_o = (1/p.cm)*rho_surf[cells.mem_to_cells]

            self.vm = (self.vm - (p.dt/p.cm)*self.Jn +
                       ((p.dt*self.sigma_cell[cells.mem_to_cells]*vm_o)/(p.cm*cells.R_rads)))/(1 +
                       ((p.dt*self.sigma_cell[cells.mem_to_cells])/(p.cm*cells.R_rads)))

            # average vm:
            self.vm_ave = np.dot(cells.M_sum_mems, self.vm) / cells.num_mems

            # True cell radii:
            Rcells = cells.R_rads*(p.true_cell_size/p.cell_radius)
            # Rcells = cells.R_rads

            # # average intracellular electric field at cell centres:
            gE = (self.vm - self.vm_ave[cells.mem_to_cells]) /Rcells  # concentration gradients
            gEx = -gE * cells.mem_vects_flat[:, 2]
            gEy = -gE * cells.mem_vects_flat[:, 3]

            self.E_cell_x = np.dot(cells.M_sum_mems, gEx * cells.mem_sa) / cells.cell_sa
            self.E_cell_y = np.dot(cells.M_sum_mems, gEy * cells.mem_sa) / cells.cell_sa

            # calculate electric field in cells using net intracellular current and cytosol conductivity:
            self.Emc = (self.E_cell_x[cells.mem_to_cells] * cells.mem_vects_flat[:, 2] +
                       self.E_cell_y[cells.mem_to_cells] * cells.mem_vects_flat[:, 3])

        # calculate the derivative of Vmem:
        self.dvm = (self.vm - vmo)/p.dt

    def update_all_concs(self, cells, p):

        for i in self.movingIons:

            f_mem_i = self.fluxes_mem[i]

            f_gj_i = self.fluxes_gj[i]

            self.cc_cells[i], self.cc_at_mem[i], self.cc_env[i] = stb.update_Co(self, self.cc_cells[i],
                                                                                self.cc_at_mem[i],
                                                                                self.cc_env[i], f_mem_i,
                                                                                cells, p,
                                                                                ignoreECM = self.ignore_ecm)

            delta_cgj = np.dot(cells.M_sum_mems, -f_gj_i*cells.mem_sa) / cells.cell_vol

            self.cc_cells[i] +=  p.dt*delta_cgj

            # ensure no negative concentrations:
            stb.no_negs(self.cc_cells[i])

    def acid_handler(self, cells, p) -> None:
        '''
        Update H+ concentrations in both the cell cluster and environment,
        which are further influenced by the bicarbonate buffer action.

        This method additionally runs the HKATPase and VATPase pumps if enabled
        by the current simulation configuration.
        '''

        pass

    def ca_handler(self,cells,p):


        # operate default Ca2+ pumps only if the user-supplied rate is greater than zero, so user can shut them off.
        if p.alpha_Ca > 0.0:

            if p.is_ecm is True:

                # run Ca-ATPase

                f_CaATP = stb.pumpCaATP(self.cc_at_mem[self.iCa],
                    self.cc_env[self.iCa][cells.map_mem2ecm],
                    self.vm, self.T, p, self.CaATP_block, met = self.met_concs)

                f_CaATP = self.rho_pump * f_CaATP

            else:

                # run Ca-ATPase

                f_CaATP = stb.pumpCaATP(self.cc_at_mem[self.iCa],
                                        self.cc_env[self.iCa], self.vm, self.T, p,
                                        self.CaATP_block, met = self.met_concs)


            self.rate_CaATP = f_CaATP

            if p.cluster_open is False:
                f_CaATP[cells.bflags_mems] = 0

            # store the transmembrane flux for this ion
            self.fluxes_mem[self.iCa] += self.rho_pump*(f_CaATP)


        if p.Ca_dyn == 1:  # do endoplasmic reticulum handling

            self.endo_retic.update(self, cells, p)

    def update_gj(self,cells,p,t,i):

        # calculate voltage difference (gradient*len_gj) between gj-connected cells:

        self.vgj = self.vm[cells.nn_i]- self.vm[cells.mem_i]

        ## smooth the vgj:
        # vgj_ave = np.dot(cells.M_sum_mems, self.vgj * cells.mem_sa) / cells.cell_sa
        # self.vgj = self.smooth_weight_mem * self.vgj + vgj_ave[cells.mem_to_cells] * self.smooth_weight_o

        # store transjunctional electric field:
        self.Egj = -self.vgj/cells.gj_len

        # self.Egj = -self.vgj/cells.nn_len

        self.E_gj_x = self.Egj*cells.mem_vects_flat[:,2]
        self.E_gj_y = self.Egj*cells.mem_vects_flat[:,3]

        if p.v_sensitive_gj is True:

            # run the gap junction dynamics object to update gj open state of sim:
            self.gj_funk.run(self, cells, p)

        else:
            self.gjopen = self.gj_block*np.ones(len(cells.mem_i))*cells.gj_default_weights


        conc_mem = self.cc_at_mem[i]

        fgj_X = stb.electroflux(conc_mem[cells.mem_i],
                       conc_mem[cells.nn_i],
                       self.D_gj[i]*p.gj_surface*self.gjopen,
                       cells.gj_len*np.ones(self.mdl),
                       self.zs[i]*np.ones(self.mdl),
                       self.vgj,
                       p.T,
                       p,
                       rho=1
                       )

        # enforce zero flux at outer boundary:
        fgj_X[cells.bflags_mems] = 0.0


        self.fluxes_gj[i] = self.fluxes_gj[i] + fgj_X   # store gap junction flux for this ion

    def update_ecm(self,cells,p,t,i):

        cenv = self.cc_env[i]
        cenv = cenv.reshape(cells.X.shape)

        cenv[:,0] =  self.c_env_bound[i]
        cenv[:,-1] =  self.c_env_bound[i]
        cenv[0,:] =  self.c_env_bound[i]
        cenv[-1,:] =  self.c_env_bound[i]

        gcx, gcy = fd.gradient(cenv, cells.delta)

        if p.fluid_flow is True:

            ux = self.u_env_x
            uy = self.u_env_y

        else:

            ux = np.zeros(cells.X.shape)
            uy = np.zeros(cells.X.shape)


        denv = self.D_env[i].reshape(cells.X.shape)*self.TJ_modulator[i].reshape(cells.X.shape)

        # this equation assumes environmental transport is electrodiffusive--------------------------------------------:
        fx, fy = stb.nernst_planck_flux(cenv, gcx, gcy, -self.E_env_x, -self.E_env_y, ux, uy,
                                          denv, self.zs[i], self.T, p)


        self.fluxes_env_x[i] = fx.ravel()  # store ecm junction flux for this ion
        self.fluxes_env_y[i] = fy.ravel()  # store ecm junction flux for this ion

        # divergence of total flux:
        div_fa = fd.divergence(-fx, -fy, cells.delta, cells.delta)

        # update concentration in the environment:
        cenv = cenv + div_fa * p.dt

        if p.sharpness < 1.0:

            # smooth concentration in the environment:
            cenv = fd.integrator(cenv, sharp = p.sharpness)

        self.cc_env[i] = cenv.ravel()

    def update_intra(self, cells, p, i):

        cav = self.cc_cells[i][cells.mem_to_cells]  # concentration at cell centre
        # cmi = self.cc_at_mem[i]  # concentration at membrane
        # z = self.zs[i]    # charge of ion
        # Do = 0.1*self.D_free[i]  # diffusion constant of ion, assuming diffusion in cytoplasm is 10x slower than free
        #
        # cp = (cav + cmi)/2   # concentration at midpoint between cell centre and membrane
        # cg = (cmi - cav)/cells.R_rads  # concentration gradients
        #
        # # normal component of electric field at membranes:
        # En = (self.E_cell_x[cells.mem_to_cells] * cells.mem_vects_flat[:, 2] +
        #       self.E_cell_y[cells.mem_to_cells] * cells.mem_vects_flat[:, 3])

        # calculate normal component of microtubules at membrane:
        # umtn = self.mtubes.mtubes_x*cells.mem_vects_flat[:, 2] + self.mtubes.mtubes_y*cells.mem_vects_flat[:, 3]


        # if p.cell_polarizability != 0.0:
        #
        #     # cflux = np.zeros(self.mdl)
        #     # self.cc_at_mem[i] = cav*1
        #
        #     cfluxo = (-Do*cg + ((Do*p.q*cp*z)/(p.kb*self.T))*En)*p.cell_polarizability
        #
        #     # as no net mass must leave this intracellular movement, make the flux divergence-free:
        #     cflux = stb.single_cell_div_free(cfluxo, cells)
        #
        #     # update the concentration at membranes:
        #     # flux is positive as the field is internal to the cell, working in the opposite direction to transmem fluxes
        #     self.cc_at_mem[i] = cmi + cflux*(cells.mem_sa/cells.mem_vol)*p.dt
        #
        #     # smooth the concentration:
        #     self.cc_at_mem[i] = self.smooth_weight_mem*self.cc_at_mem[i] + cav*self.smooth_weight_o
        #
        # elif p.cell_polarizability == 0.0:

            # cflux = np.zeros(self.mdl)
            # self.cc_at_mem[i] = cav*1


        # deal with the fact that our coarse diffusion model may leave some sub-zero concentrations:
        # indsZ = (self.cc_at_mem[i] < 0.0).nonzero()
        #
        # if len(indsZ[0]):
        #
        #     raise BetseSimUnstableException("Ion concentration value on membrane below zero! Your simulation has"
        #                                        " become unstable.")
        #
        # # update the main matrices:
        # self.fluxes_intra[i] = cflux * 1

        # uncomment this to skip the above computational loop ---------------
        self.cc_at_mem[i] = cav*1

    def get_ion(self, ion_name: str) -> int:
        '''
        0-based index assigned to the ion with the passed name (e.g., ``Na``) if
        this ion is enabled by this simulation *or* the empty list otherwise.

        This index is guaranteed to uniquely (but arbitrarily) identify this ion
        with respect to this simulation.
        '''

        if ion_name == 'Na':
            ion = self.iNa
        elif ion_name == 'K':
            ion = self.iK
        elif ion_name == 'Ca':
            ion = self.iCa
        elif ion_name == 'Cl':
            ion = self.iCl
        elif ion_name == 'M':
            ion = self.iM
        elif ion_name == 'H':
            ion = self.iH
        elif ion_name == 'P':
            ion = self.iP
        else:
            #FIXME: Shouldn't a fatal exception be raised instead?
            logs.warning(
                'Oops! Molecule gated ion "%s" channel target not found!',
                ion_name)
            ion = []

        return ion

    def initDenv(self,cells,p):
        '''
        Initialize the environmental diffusion matrix and corresponding weight
        matrices, including tight and adherin junctions.
        '''

        if p.is_ecm is True:

            for i, dmat in enumerate(self.D_env):

                Denv_o = np.ones(self.edl) * self.D_free[i]

                # adherens junctions slow diffusion throughout the cell cluster:
                Denv_o[cells.envInds_inClust] = self.D_free[i]*p.D_adh

                # if p.env_type is True:
                Denv_o[cells.all_bound_mem_inds] = self.D_free[i]*p.D_tj*self.Dtj_rel[i]
                Denv_o[cells.interior_bound_mem_inds] = self.D_free[i] * p.D_tj * self.Dtj_rel[i]
                Denv_o[cells.ecm_inds_bound_cell] = self.D_free[i] * p.D_tj * self.Dtj_rel[i]

                # create an ecm diffusion grid filled with the environmental values
                self.D_env[i] = Denv_o*1.0

            # create a matrix that weights the relative transport efficiency in the world space:
            D_env_weight = self.D_env[self.iNa]/self.D_env[self.iNa].max()
            self.D_env_weight = D_env_weight.reshape(cells.X.shape)
            self.D_env_weight_base = np.copy(self.D_env_weight)

            self.TJ_modulator = np.ones(self.D_env.shape)


        else:

            Denv_o = np.ones(len(cells.xypts))

            # if p.env_type is True:
            Denv_o[cells.all_bound_mem_inds] = p.D_tj
            Denv_o[cells.interior_bound_mem_inds] = p.D_tj
            Denv_o[cells.ecm_inds_bound_cell] = p.D_tj

            Denv_o = gaussian_filter(Denv_o.reshape(cells.X.shape), 1.0)

            self.D_env = np.copy(self.D_free)


            # create a matrix that weights the relative transport efficiency in the world space:
            self.D_env_weight = Denv_o.reshape(cells.X.shape)
            self.D_env_weight_base = np.copy(self.D_env_weight)

        self.TJ_targets = np.hstack(
            (cells.all_bound_mem_inds, cells.interior_bound_mem_inds, cells.ecm_inds_bound_cell))

        self.Chi = np.ones(len(cells.xypts)) * p.er  # electrical susceptibility of pure water
        # self.Chi[cells.envInds_inClust] = 2.0e5  # electrical susceptibility of tissue

        self.Chi = gaussian_filter(self.Chi.reshape(cells.X.shape), 2)

    # ..................{ PLOTTERS                           }.................
    def _plot_loop(self, phase: SimPhase) -> tuple:
        '''
        Display and/or save an animation during solving if requested *and*
        calculate data common to solving both with and without extracellular
        spaces.

        Parameters
        --------
        phase : SimPhase
            Current simulation phase.

        Returns
        --------
        (ndarray, set, ContextManager)
            3-tuple ``(time_steps, time_steps_sampled, solver_context)`` where:

            * ``time_steps`` is a one-dimensional Numpy array defining the
              time-steps vector for the current phase.
            * ``time_steps_sampled`` is the subset of the ``time_steps`` array
              whose elements are **sampled time steps** (i.e., time steps at
              which to sample data, substantially reducing data storage). In
              particular, the length of this set is exactly equal to:

              * The maximum progress value to be passed to the
                :meth:`SimCallbacksAPI.progress_ranged` callback.
              * The number of frames exported from each animation.

            * ``solver_context`` is the context manager intended to wrap the
              core time loop for this phase. Specifically, this is either:
              * If the configuration for this phase enables non-blocking display
                and/or saving of one or more in-phase exports (e.g., plots), the
                context manager doing so.
              * Else, the empty context manager doing nothing.
        '''

        #FIXME: Refactor these conditionally set local variables into
        #attributes of the "Parameters" object set instead by the
        #Parameters.set_time_profile() method.

        if phase.kind is SimPhaseKind.INIT:
            phase_verb = 'Initializing'
            phase_noun = 'initialization'
            phase_time_step_count = phase.p.init_tsteps
        else:
            phase_verb = 'Simulating'
            phase_noun = 'simulation'
            phase_time_step_count = phase.p.sim_tsteps

        # Total number of seconds simulated by the current run.
        phase_time_len = phase_time_step_count * phase.p.dt

        # Time-steps vector appropriate for the current run.
        time_steps = np.linspace(0, phase_time_len, phase_time_step_count)

        time_steps_sampled = set()
        i = 0
        while i < len(time_steps) - phase.p.t_resample:
            i += phase.p.t_resample
            i = int(i)
            time_steps_sampled.add(time_steps[i])

        # Log this run.
        logs.log_info(
            'Your %s is running from 0 to %.2f s of in-world time '
            'in %d time steps (%d sampled).',
            phase_noun,
            phase_time_len,
            len(time_steps),
            len(time_steps_sampled),
        )

        # Mid-simulation animation of cell voltage as a function of time if
        # enabled by this configuration or None otherwise.
        solver_context = None

        # If this animation is enabled...
        if phase.p.anim.is_while_sim:
            #FIXME: This is terrible. Ideally, all deformations should already
            #be properly incorporated into the current cell cluster *WITHOUT*
            #requiring "cellso" shenanigans. (Wake me up when Utopia arrives.)

            # Deformed simulation phase, replacing the current cell cluster with
            # the deformed cell cluster *ONLY* for this animation. If
            # deformations are disabled, this is a noop; if deformations are
            # enabled, this is required to animate deformations while solving.
            phase_deformed = SimPhase(
                kind=phase.kind, sim=phase.sim, cells=self.cellso, p=phase.p)

            # Create this animation.
            solver_context = AnimCellsWhileSolving(
                phase=phase_deformed,
                conf=phase.p.anim.anim_while_sim,

                # Number of frames to animate, corresponding to the number of
                # sampled time steps. The plot_frame() method of this animation
                # will be called only for each such step.
                time_step_count=len(time_steps_sampled),

                # Animation metadata.
                label='Vmem',
                figure_title='Vmem while {}'.format(phase_verb),
                colorbar_title='Voltage [mV]',
            )
        else:
            solver_context = noop_context()

        # Return the 3-tuple of these objects to the caller.
        return time_steps, time_steps_sampled, solver_context
