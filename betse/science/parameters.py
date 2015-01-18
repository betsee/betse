#!/usr/bin/env python3
# Copyright 2014-2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

# FIXME this module will load parameters from a yaml file!
# FIXME put all constants in as fields of the structure then from science.parameters import Parameters as p wherever
#they're needed



# define the basic class that holds variables
class Parameters(object):
    """
    For now, a very simple object that stores simulation constants

    """
    def __init__(self):
        self.wsx = 100e-6  # the x-dimension of the world space [m] recommended range 50 to 1000 um
        self.wsy = 100e-6  # the y-dimension of the world space [m] recommended range 50 to 1000 um
        self.rc = 5e-6  # radius of single cell
        self.d_cell = self.rc * 2  # diameter of single cell
        self.nx = int(self.wsx / self.d_cell)  # number of lattice sites in world x index
        self.ny = int(self.wsy / self.d_cell)  # number of lattice sites in world y index
        self.ac = 1e-6  # cell-cell separation for drawing
        self.dc = self.rc * 2  # cell diameter
        self.nl = 0.8  # noise level for the lattice
        self.wsx = self.wsx + 5 * self.nl * self.d_cell  # readjust the world size for noise
        self.wsy = self.wsy + 5 * self.nl * self.d_cell
        self.search_d =1.5     # distance to search for nearest neighbours (relative to cell diameter dc) min 1.0 max 5.0
        self.scale_cell = 0.9          # the amount to scale cell membranes in from ecm edges (only affects drawing)
        self.cell_sides = 4      # minimum number of membrane domains per cell (must be >2)
        self.scale_alpha = 1.0   # the amount to scale (1/d_cell) when calculating the concave hull (boundary search)
        self.cell_height = 5.0e-6  # the height of a cell in the z-direction (for volume and surface area calculations)
        self.cell_space = 26.0e-9  # the true cell-cell spacing (width of extracellular space)
        self.um = 1e6    # multiplication factor to convert m to um
        self.F = 96485 # Faraday constant [J/V*mol]
        self.R = 8.314  # Gas constant [J/K*mol]


params = Parameters()