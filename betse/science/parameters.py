#!/usr/bin/env python3
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME this module will load parameters from a yaml file!
# FIXME put all constants in as fields of the structure then from science.parameters import Parameters as p wherever
#they're needed

# Lodish H, Berk A, Zipursky SL, et al. Molecular Cell Biology. 4th edition. New York: W. H. Freeman;
# 2000. Section 15.4, Intracellular Ion Environment and Membrane Electric Potential.
# Available from: http://www.ncbi.nlm.nih.gov/books/NBK21627/

import numpy as np

# define the basic class that holds variables
class Parameters(object):
    """
    For now, a very simple object that stores simulation constants.

    """
    def __init__(self,profile=None):

        self.dt = 1e-2    # Simulation step-size [s]

        # basic constants
        self.F = 96485 # Faraday constant [J/V*mol]
        self.R = 8.314  # Gas constant [J/K*mol]
        self.T = 310   # Temperature [K]

        # geometric constants and factors
        self.wsx = 100e-6  # the x-dimension of the world space [m] recommended range 50 to 1000 um
        self.wsy = 100e-6  # the y-dimension of the world space [m] recommended range 50 to 1000 um
        self.rc = 5e-6  # radius of single cell
        self.d_cell = self.rc * 2  # diameter of single cell
        self.nx = int(self.wsx / self.d_cell)  # number of lattice sites in world x index
        self.ny = int(self.wsy / self.d_cell)  # number of lattice sites in world y index
        self.ac = 1e-6  # cell-cell separation for drawing
        self.nl = 0.8  # noise level for the lattice
        self.wsx = self.wsx + 5 * self.nl * self.d_cell  # readjust the world size for noise
        self.wsy = self.wsy + 5 * self.nl * self.d_cell
        self.vol_env = 1        # volume of the environmental space [m3]
        self.search_d =1.5     # distance to search for nearest neighbours (relative to cell diameter dc) min 1.0 max 5.0
        self.scale_cell = 0.9          # the amount to scale cell membranes in from ecm edges (only affects drawing)
        self.cell_sides = 4      # minimum number of membrane domains per cell (must be >2)
        self.scale_alpha = 1.0   # the amount to scale (1/d_cell) when calculating the concave hull (boundary search)
        self.cell_height = 5.0e-6  # the height of a cell in the z-direction (for volume and surface area calculations)
        self.cell_space = 26.0e-9  # the true cell-cell spacing (width of extracellular space)
        self.cm = 0.010            # patch capacitance of cell membrane up to 0.022 [F/m2]
        self.tm = 7.5e-9           # thickness of cell membrane [m]
        self.um = 1e6    # multiplication factor to convert m to um

        # diffusion constants
        self.Dm_Na = 1.0e-18     # membrane diffusion constant sodium [m2/s]
        self.Dm_K = 1.0e-16      # membrane diffusion constant potassium [m2/s]
        self.Dm_Cl = 1.0e-17     # membrane diffusion constant chloride [m2/s]
        self.Dm_Ca = 1.0e-18     # membrane diffusion constant calcium [m2/s]
        self.Dm_H = 1.0e-18      # membrane diffusion constant hydrogen [m2/s]
        self.Dm_M = 1.0e-18     # membrane diffusion constant anchor ion [m2/s]
        self.Dm_P = 0.0        # membrane diffusion constant proteins [m2/s]

        self.Do_Na = 1.0e-9      # free diffusion constant sodium [m2/s]
        self.Do_K = 1.0e-9      # free diffusion constant potassium [m2/s]
        self.Do_Cl = 1.0e-9     # free diffusion constant chloride [m2/s]
        self.Do_Ca = 1.0e-9     # free diffusion constant calcium [m2/s]
        self.Do_H = 1.0e-9      # free diffusion constant hydrogen [m2/s]
        self.Do_M = 1.0e-9     # free diffusion constant mystery anchor ion [m2/s]
        self.Do_P = 5.0e-9      # free diffusion constant protein [m2/s]

        # charge states of ions
        self.z_Na = 1
        self.z_K = 1
        self.z_Cl = -1
        self.z_Ca = 2
        self.z_H = 1
        self.z_P = -1
        self.z_M = -1

        zs = [self.z_Na, self.z_K, self.z_Cl, self.z_Ca, self.z_H, self.z_P]

        # default environmental and initial values mammalian cells and plasma
        if profile == 'mammalian' or profile == None:
            self.cNa_env = 145.0
            self.cK_env = 5.0
            self.cCl_env = 105.0
            self.cCa_env = 1.0
            self.cH_env = 4.0e-8
            self.cP_env = 9.0

            conc_env = [self.cNa_env,self.cK_env, self.cCl_env, self.cCa_env, self.cH_env, self.cP_env]
            self.cM_env, self.z_M_env = bal_charge(conc_env,zs)

            self.cNa_cell = 17.0
            self.cK_cell = 131.0
            self.cCl_cell = 6.0
            self.cCa_cell = 1.0e-6
            self.cH_cell = 6.3e-8
            self.cP_cell = 138.0

            conc_cell = [self.cNa_cell,self.cK_cell, self.cCl_cell, self.cCa_cell, self.cH_cell, self.cP_cell]
            self.cM_cell, self.z_M_cell = bal_charge(conc_cell,zs)

         # default environmental and initial values invertebrate cells and plasma
        if profile == 'invertebrate':
            self.cNa_env = 440.0
            self.cK_env = 20.0
            self.cCl_env = 560.0
            self.cCa_env = 10.0
            self.cH_env = 4.0e-8
            self.cP_env = 7.0

            conc_env = [self.cNa_env,self.cK_env, self.cCl_env, self.cCa_env, self.cH_env, self.cP_env]
            self.cM_env, self.z_M_env = bal_charge(conc_env,zs)

            self.cNa_cell = 50.0
            self.cK_cell = 400.0
            self.cCl_cell = 75.0
            self.cCa_cell = 3.0e-4
            self.cH_cell = 6.3e-8
            self.cP_cell = 350.0

            conc_cell = [self.cNa_cell,self.cK_cell, self.cCl_cell, self.cCa_cell, self.cH_cell, self.cP_cell]
            self.cM_cell, self.z_M_cell = bal_charge(conc_cell,zs)

        # pump parameters
        self.deltaGATP = 50e3    # free energy released in ATP hydrolysis [J/mol]
        self.alpha_NaK = 5.0e-17 # rate constant sodium-potassium ATPase [m3/mols]  range 1.0e-9 to 1.0e-10 for dt =1e-2
        self.halfmax_NaK = 12   # the free energy level at which pump activity is halved [kJ]
        self.slope_NaK = 24  # the energy window width of the NaK-ATPase pump [kJ]


def bal_charge(concentrations,zs):

    q = 0

    for conc,z in zip(concentrations,zs):
        q = q+ conc*z

        to_zero = -q
        bal_conc = abs(to_zero)
        valance = np.sign(to_zero)

    return bal_conc,valance


params = Parameters(profile=None)
