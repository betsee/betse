#!/usr/bin/env python3
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME this module will load parameters from a yaml file!
# FIXME put *all* constants and options in here, including plotting colormaps, etc...

# Lodish H, Berk A, Zipursky SL, et al. Molecular Cell Biology. 4th edition. New York: W. H. Freeman;
# 2000. Section 15.4, Intracellular Ion Environment and Membrane Electric Potential.
# Available from: http://www.ncbi.nlm.nih.gov/books/NBK21627/

from betse.science import simconfig
from betse.util.path import paths
import numpy as np
import math
import matplotlib.cm as cm
import os

# define the basic class that holds variables
class Parameters(object):
    '''
    The object that stores all constants used in world-building, simulation, and
    plotting.
    '''
    def __init__(self):
        self.time_profile_init = 'initialize'        # choose time profile for initialization sim
        self.time_profile_sim = 'simulate_somatic'   # choice of 'simulate_excitable' or 'simulate_somatic'

        self.time4init = 10*60      # set the time for the initialization sim [s]
        self.time4sim = 5*60        # set total time for simulation [s]

        # File saving
        self.cache_path = os.path.expanduser("~/.betse/cache/100umInit/")  # world, inits, and sims are saved and read to/from this directory.
        self.sim_path = os.path.expanduser("~/.betse/cache/100umInit/sim_test") # folder to save unique simulation and data linked to init
        self.sim_results = os.path.expanduser("~/.betse/cache/100umInit/sim_test/results") # folder to auto-save results (graphs, images, animations)

        # Geometric constants and factors
        self.wsx = 100e-6  # the x-dimension of the world space [m] recommended range 50 to 1000 um
        self.wsy = 100e-6  # the y-dimension of the world space [m] recommended range 50 to 1000 um
        self.rc = 5e-6  # radius of single cell
        self.cell_height = 5.0e-6  # the height of a cell in the z-direction (for volume and surface area calculations)
        self.cell_space = 26.0e-9  # the true cell-cell spacing (width of extracellular space)
        self.nl = 0.8  # noise level for the lattice
        self.vol_env = 2*self.wsx*self.wsy*self.cell_height    # volume of the environmental space [m3]

        self.T = 310   # World temperature [K]

        # gap junction constants and network connectivity
        self.search_d =1.5     # distance to search for nearest neighbours (relative to cell diameter dc) min 1.0 max 5.0

        self.gj_vthresh = 60e-3              # cell-cell voltage threshhold at which gj close [V]
        self.gj_vgrad  = 30e-3               # the range over which gj goes from open to shut at threshold [V]

        # set ion profile to be used: 'basic' (4 ions), 'basic_Ca' (5 ions), 'animal' (7 ions), 'invertebrate' (7 ions)
        self.ion_profile = 'animal'

        # include full calcium dynamics in the situation (i.e. endoplasmic reticulum, etc)? Yes = 1, No =0
        self.Ca_dyn = 1

        # include HK-ATPase in the simulation? Yes =1, No = 0
        self.HKATPase_dyn = 0

        # include diffusion of a voltage sensitive dye? Yes = 1, No = 0
        self.voltage_dye = 0

        self.Dm_Dye = 1.0e-12  # voltage sensitive dye membrane diffusion coefficient [m2/s]
        self.Do_Dye = 1.0e-9   # gap junction diffusion constant of voltage-sensitive dye [m2/s]
        self.z_Dye = 1         # charge valence of dye
        self.cDye_to = 1.0e-3    # initial concentration of voltage sensitive dye in environment [mol/m3]

    #..................................................................................................................
        # default membrane diffusion constants: easy control of cell's base resting potential
        self.Dm_Na = 1.0e-18     # membrane diffusion constant sodium [m2/s]
        self.Dm_K = 15.0e-18      # membrane diffusion constant potassium [m2/s]
        self.Dm_Cl = 2.0e-18     # membrane diffusion constant chloride [m2/s]
        self.Dm_Ca = 1.0e-18     # membrane diffusion constant calcium [m2/s]
        self.Dm_H = 1.0e-16      # membrane diffusion constant hydrogen [m2/s]
        self.Dm_M = 1.0e-18     # membrane diffusion constant anchor ion [m2/s]
        self.Dm_P = 0.0        # membrane diffusion constant proteins [m2/s]

        #.............................Scheduled Interventions..........................................................

        # Schedule global changes to all cells in the collective:
        self.global_options = {'K_env':0,'Cl_env':0,'Na_env':0,'gj_block':0,'T_change':0,'NaKATP_block':0,
            'HKATP_block':0}
        # K_env, Cl_env, Na_env, T_change: [time on, time off, rate change, multiplier]
        # gj_block, NaKATP_block,HKATP_block, CaATP_block: [time on, time off, rate change]

        # cell to effect in scheduled intervention: (choices = 'none','all','random1','random50', [1,2,3])
        self.scheduled_targets = [0]    # targets do not affect: gj_block,T_change, NaKATP_block, HKATP_block,CaATP_block,CaER_block

        #self.ion_options specifications list is [time on, time off, rate of change, multiplier]
        self.scheduled_options = {'Na_mem':0,'K_mem':0,'Cl_mem':0,'Ca_mem':0,'K_env':0,'Cl_env':0,
            'IP3':[5,15,1,1e-4]}

        #...................................Voltage Gated Channels......................................................

        # cells to effect with voltage gated channels: (choices = 'none','all','random1','random50', [1,2,3])
        self.gated_targets = 'all'
        # self.vg_options specifications list for voltage gated ion channel options:
        vgNa = [1.0e-15,-50e-3,30e-3,-52e-3,5e-3,10e-3]  # [max Na mem diffusion m2/s, v on, v inactive, v deactivate,duration active (s), duration inactive]
        vgK = [0.5e-15, -20e-3,-75e-3,10.0e-3]           # [max K mem diffusion (m2/s), v on, v off, duration (s)]
        vgCa = [1.0e-15,-40e-3,40e-3,0.75e-3,200.0e-6]  # [maxCa mem diffusion m2/s, v on, v off, Ca2+ off mmol/L, Ca2+ reactivate]
        cagK = [2.0e-16,7.5e-4,3]                    # [maxK mem diffusion (m2/s), half-max Ca2+ for gating, hill coefficient]

        self.vg_options = {'Na_vg':0,'K_vg':0,'Ca_vg':0,'K_cag':cagK}

        # Calcium Dynamics: Calcium Induced Calcium Release (CICR) and Store Operated Calcium Entry (SOCE)..............


        ERstore_dyn = [5e-15,0.8,0.5]   # base dynamics of endoplasmic reticulum Ca2+ store:
                                        # [max diffusion m2/s, full Ca thresh mol/m3, empty Ca thresh mol/m3]
        ca_reg = []   # central concentration for Ca-act-Ca release [Ca mid, Ca width]
        #ca_reg =[]                  # leave this empty to have no Ca-influence on the dynamics

        ip3_reg = [1e-3,3.4]   # max Ca2+ diffusion constant through ER membrane, IP3 half-max, Hill coefficient
        #ip3_reg = []            # leave this empty to have no ip3-influence on the dynamics
        self.FMmod = 1              # frequency modulate ER response to IP3? 1 = yes, 0 = no
        self.ip3FM = 0.8            # degree to which ip3 affects Ca2+ frequency (higher more effect, max 0.9)


        cicr = [ERstore_dyn,ca_reg,ip3_reg]
        self.Ca_dyn_options = {'CICR':cicr}

        self.Dm_IP3 = 1.0e-18   # membrane diffusion constant of IP3
        self.Do_IP3 = 1.0e-5    # IP3 free diffusion constant [m2/s] (this is artificially high due to gj being artificially low)
        self.z_IP3 = -3        # charge valence of IP3
        self.cIP3_to = 1e-6     # initial value of IP3 in all cells
        self.cIP3_to_env = 1e-6  # initial value of IP3 in environment

        #..........................PLOTTING OPTIONS and OUTPUT..........................................................

         # Default colormap
        self.default_cm = cm.coolwarm   # options include cm.rainbow, cm.jet, cm.Blues, cm.Greens, see:
                                        # http://matplotlib.org/examples/color/colormaps_reference.html

        self.gj_cm = cm.bone           # colormap for plotting gj currents on top of default colormap

        self.plot_while_solving = True  # create a 2d plot of cell vmems while solution is taking place
        self.save_solving_plot = True   # save the 2d plot generated while solving (warning: will slow sim down!)

        self.enumerate_cells = False    # number cells on the static 2D maps with their simulation index (this can help
                                        # decide on the value of self.plot_cell

        self.plot_cell = 11             # State the cell index to use for single-cell time plots

        self.plot_single_cell_graphs = True # plot graphs of concentration and voltage in self.plot_cell with time

        self.showCells = True     # plots and ani are individual cell plots if True; as interpolated mesh data if False

        self.plot_vm2d = False                # 2d plot of final vmem ?
        self.plot_ca2d = False                # 2d plot of final cell calcium ?
        self.plot_ip32d = False               # 2d plot of final cIP3 ?
        self.plot_dye2d = False               # 2d plot of voltage sensitive dye in cell collective?

        self.createAnimations = True   # create all animations = True; turn off all animations = False

        # specify desired animations:
        self.ani_vm2d = True                # 2d animation of vmem with time?
        self.ani_ca2d = False                # 2d animation of cell calcium with time ?
        self.ani_ip32d = False               # 2d animation of cIP3 with time?
        self.ani_dye2d = False               # 2d animation of voltage sensitive dye in cell collective with time?
        self.ani_vmgj2d = False              # 2d animation of vmem with superimposed gj network showing current direction

        self.autosave = True           # autosave all still images to a results directory in the simulation folder
        self.saveAnimations = True    # save all animations as png sequences in animation-specific folders

        self.exportData = True        # export all stored data for the plot_cell to a csv text file

        self.clip = 20e-6



        # ........................Rarely changed constants and calculations.............................................

        # default free diffusion constants (cytoplasmic)
        self.Do_Na = 1.33e-9      # free diffusion constant sodium [m2/s]
        self.Do_K = 1.96e-9      # free diffusion constant potassium [m2/s]
        self.Do_Cl = 2.03e-9     # free diffusion constant chloride [m2/s]
        self.Do_Ca = 1.0e-10     # free diffusion constant calcium [m2/s]
        self.Do_H = 2.5e-9      # free diffusion constant hydrogen [m2/s]
        self.Do_M = 1.0e-9     # free diffusion constant mystery anchor ion [m2/s]
        self.Do_P = 5.0e-10      # free diffusion constant protein [m2/s]

        # pump parameters
        self.alpha_NaK = 1.0e-17 # maximum rate constant sodium-potassium ATPase [m3/mols] (range 1e-17 to 5e-16)
        self.halfmax_NaK = 12   # the free energy level at which pump activity is halved [kJ]
        self.slope_NaK = 24  # the energy window width of the NaK-ATPase pump [kJ]

        self.alpha_Ca = 5.0e-15 # pump rate for calcium ATPase in membrane [m3/mols] 2.0e-15
        self.alpha_CaER = 5.0e-14  # pump rate for calcium ATPase in endoplasmic reticulum
        self.halfmax_Ca = 12
        self.slope_Ca = 24

        self.alpha_HK = 5.0e-14  # pump rate for the H-K-ATPase
        self.halfmax_HK = 12
        self.slope_HK = 24

         # Endoplasmic reticulum
        self.ER_vol = 0.1                  # volume of endoplasmic reticulum as a fraction of cell volume
        self.ER_sa = 1.0                    # surface area of endoplasmic reticulum as a fraction of cell surface area

        # partial pressure dissolved CO2
        self.CO2 = 0.03*40
        self.bicarb = 25.0

        # charge states of ions
        self.z_Na = 1
        self.z_K = 1
        self.z_Cl = -1
        self.z_Ca = 2
        self.z_H = 1
        self.z_P = -1
        self.z_M = -1

        # fundamental constants
        self.F = 96485 # Faraday constant [J/V*mol]
        self.R = 8.314  # Gas constant [J/K*mol]

        self.deltaGATP = 20*self.R*self.T    # free energy released in ATP hydrolysis [J/mol]

        self.ac = 1e-6  # cell-cell separation for drawing
        self.scale_cell = 0.9          # the amount to scale cell membranes in from ecm edges (only affects drawing)
        self.cm = 0.022            # patch capacitance of cell membrane up to 0.022 [F/m2]
        self.tm = 7.5e-9           # thickness of cell membrane [m]
        self.cell_sides = 4      # minimum number of membrane domains per cell (must be >2)
        self.scale_alpha = 1.0   # the amount to scale (1/d_cell) when calculating the concave hull (boundary search)

        self.d_cell = self.rc * 2  # diameter of single cell
        self.nx = int(self.wsx / self.d_cell)  # number of lattice sites in world x index
        self.ny = int(self.wsy / self.d_cell)  # number of lattice sites in world y index
        self.wsx = self.wsx + 5 * self.nl * self.d_cell  # readjust the world size for noise
        self.wsy = self.wsy + 5 * self.nl * self.d_cell

        self.gjl = 2*self.tm + self.cell_space     # gap junction length


        self.um = 1e6    # multiplication factor to convert m to um

        # simplest ion ion_profile giving realistic results with minimal ions (Na+ & K+ focus):
        if self.ion_profile == 'basic':

            self.cNa_env = 145.0
            self.cK_env = 5.0
            self.cP_env = 9.0

            zs = [self.z_Na, self.z_K, self.z_P]

            conc_env = [self.cNa_env,self.cK_env, self.cP_env]
            self.cM_env, self.z_M_env = bal_charge(conc_env,zs)

            assert self.z_M_env == -1

            self.cNa_cell = 5.4
            self.cK_cell = 140.44
            self.cP_cell = 138.0

            conc_cell = [self.cNa_cell,self.cK_cell, self.cP_cell]

            self.cM_cell, self.z_M_cell = bal_charge(conc_cell,zs)

            assert self.z_M_cell == -1

            self.ions_dict = {'Na':1,'K':1,'Cl':0,'Ca':0,'H':0,'P':1,'M':1}


        if self.ion_profile == 'basic_Ca':

            self.cNa_env = 145.0
            self.cK_env = 5.0
            self.cCa_env = 1.0
            self.cP_env = 9.0

            zs = [self.z_Na, self.z_K, self.z_Ca, self.z_P]

            conc_env = [self.cNa_env,self.cK_env, self.cCa_env, self.cP_env]
            self.cM_env, self.z_M_env = bal_charge(conc_env,zs)

            assert self.z_M_env == -1

            self.cNa_cell = 5.4
            self.cK_cell = 140.44
            self.cCa_cell = 1.0e-3
            self.cP_cell = 138.0

            conc_cell = [self.cNa_cell,self.cK_cell, self.cCa_cell, self.cP_cell]

            self.cM_cell, self.z_M_cell = bal_charge(conc_cell,zs)

            assert self.z_M_cell == -1

            self.cCa_er = 0.9
            self.cM_er = self.cCa_er

            self.ions_dict = {'Na':1,'K':1,'Cl':0,'Ca':1,'H':0,'P':1,'M':1}

        # default environmental and cytoplasmic initial values mammalian cells
        if self.ion_profile == 'animal':

            self.cNa_env = 145.0
            self.cK_env = 5.0
            self.cCl_env = 105.0
            self.cCa_env = 1.0
            self.cH_env = 3.98e-5
            self.cP_env = 9.0

            zs = [self.z_Na, self.z_K, self.z_Cl, self.z_Ca, self.z_H, self.z_P]

            conc_env = [self.cNa_env,self.cK_env, self.cCl_env, self.cCa_env, self.cH_env, self.cP_env]
            self.cM_env, self.z_M_env = bal_charge(conc_env,zs)

            assert self.z_M_env == -1

            self.cNa_cell = 5.4
            self.cK_cell = 140.44
            self.cCl_cell = 6.0
            self.cCa_cell = 1.0e-3
            self.cH_cell = 6.31e-5
            self.cP_cell = 138.0

            conc_cell = [self.cNa_cell,self.cK_cell, self.cCl_cell, self.cCa_cell, self.cH_cell, self.cP_cell]
            self.cM_cell, self.z_M_cell = bal_charge(conc_cell,zs)

            assert self.z_M_cell == -1

            self.cCa_er = 0.9
            self.cM_er = - self.cCa_er

            self.ions_dict = {'Na':1,'K':1,'Cl':1,'Ca':1,'H':1,'P':1,'M':1}

         # default environmental and cytoplasm values invertebrate cells
        if self.ion_profile == 'invertebrate':
            self.cNa_env = 440.0
            self.cK_env = 20.0
            self.cCl_env = 460.0
            self.cCa_env = 10.0
            self.cH_env = 3.98e-5
            self.cP_env = 7.0

            zs = [self.z_Na, self.z_K, self.z_Cl, self.z_Ca, self.z_H, self.z_P]

            conc_env = [self.cNa_env,self.cK_env, self.cCl_env, self.cCa_env, self.cH_env, self.cP_env]
            self.cM_env, self.z_M_env = bal_charge(conc_env,zs)

            assert self.z_M_env == -1

            self.cNa_cell = 8.66
            self.cK_cell = 406.09
            self.cCl_cell = 45.56
            self.cCa_cell = 3.0e-4
            self.cH_cell = 6.31e-5
            self.cP_cell = 350.0

            conc_cell = [self.cNa_cell,self.cK_cell, self.cCl_cell, self.cCa_cell, self.cH_cell, self.cP_cell]
            self.cM_cell, self.z_M_cell = bal_charge(conc_cell,zs)

            assert self.z_M_cell == -1

            self.cCa_er = 3.0e-4
            self.cM_er = 3.0e-4

            self.ions_dict = {'Na':1,'K':1,'Cl':1,'Ca':1,'H':1,'P':1,'M':1}

    def load_yaml(self, config_filename: str):
        # Dictionary loaded from such YAML file.
        config = simconfig.load(config_filename)

        # Absolute path of the parent directory of such file.
        config_dirname = paths.get_dirname(config_filename)

        #FIXME: Right. These dictionary entries have probably changed.
        self.cache_path = paths.join(
            config_dirname, config['init']['cache file'])  # world, inits, and sims are saved and read to/from this directory.
        self.sim_path = paths.join(
            config_dirname, config['run']['cache file']) # folder to save unique simulation and data linked to init
        self.sim_results = paths.join(
            config_dirname, config['plot']['media dir']) # folder to auto-save results (graphs, images, animations)

        membrane = config['variables']['membrane diffusion']
        self.Dm_Na = membrane['Dm_Na']

    def set_time_profile(self,time_profile):

        if time_profile == 'simulate_somatic':

            self.dt = 5e-3    # Simulation step-size [s] recommended range 5e-3 to 1e-4 for regular sims; 5e-5 for neural
            self.sim_end = self.time4sim         # world time to end the simulation
            self.resamp = 0.1         # time to resample in world time

            self.sim_tsteps = self.sim_end/self.dt    # Number of timesteps for the simulation
            self.t_resample = self.resamp/self.dt         # resample the time vector every x steps
            self.method = 0            # Solution method. For 'Euler' = 0, for 'RK4' = 1.

            self.gj_radius = 1.0e-9              # effective radius of gap junctions connecting cells [m] (range 0 to 5.0 e-9 m)
            self.gjsa = math.pi*((self.gj_radius)**2)      # total gap junction surface area as fraction of cell surface area

        if time_profile == 'simulate_excitable':

            self.dt = 5e-5    # Simulation step-size [s] recommended range 5e-3 to 1e-4 for regular sims; 5e-5 for neural
            self.sim_end = self.time4sim         # world time to end the simulation
            self.resamp = 5e-4         # time to resample in world time

            self.sim_tsteps = self.sim_end/self.dt    # Number of timesteps for the simulation
            self.t_resample = self.resamp/self.dt         # resample the time vector every x steps
            self.method = 0            # Solution method. For 'Euler' = 0, for 'RK4' = 1.

            self.gj_radius = 5.0e-9              # effective radius of gap junctions connecting cells [m] (range 0 to 5.0 e-9 m)
            self.gjsa = math.pi*((self.gj_radius)**2)      # total gap junction surface area as fraction of cell surface area


        elif time_profile == 'initialize':

            self.dt = 1e-2    # Simulation step-size [s] recommended range 1e-2 to 1e-3 for regular sims; 5e-5 for neural
            self.init_end = self.time4init      # world time to end the initialization simulation time [s]
            self.resamp = 1.0         # time to resample in world time

            self.init_tsteps = self.init_end/self.dt # Number of timesteps for an initialization from scratch (range 50000 to 100000)
            self.t_resample = self.resamp/self.dt         # resample the time vector every x steps
            self.method = 0            # Solution method. For 'Euler' = 0, for 'RK4' = 1.

            self.gj_radius = 1.0e-9              # effective radius of gap junctions connecting cells [m] (range 0 to 5.0 e-9 m)
            self.gjsa = math.pi*((self.gj_radius)**2)      # total gap junction surface area as fraction of cell surface area


def bal_charge(concentrations,zs):

    q = 0

    for conc,z in zip(concentrations,zs):
        q = q+ conc*z

        to_zero = -q
        bal_conc = abs(to_zero)
        valance = np.sign(to_zero)

        assert bal_conc >= 0

    return bal_conc,valance

#FIXME: Is this actually used anywhere? If not, we should probably remove this;
#it's probably consuming a bit of memory.
params = Parameters()
