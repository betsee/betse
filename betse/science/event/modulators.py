#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

import numpy as np


def periodic(pc,cells,p):
    """
    Creates a time-dependent signal of a sinusoidal function.

    Parameters
    ----------
    pc                 Indices of cells (p.sim_ECM == False) or membranes (p.sim_ECM == True)
    cells              Instance of Cells module
    p                  Instance of parameters

    Returns
    ---------
    y              Data values from 0 to 1 defining a periodic signal.

    scalar         Null output corresponding to spatial variation.

    """

    scalar = 1

    y = lambda t: np.sin(t*np.pi*(p.periodic_properties['frequency']/1) + p.periodic_properties['phase'])**2

    return scalar, y

def f_sweep(pc,cells,p):
    """
    Creates a time-dependent sweep of a sinusoidal function through various frequencies.

    Parameters
    ----------
    pc                 Indices of cells (p.sim_ECM == False) or membranes (p.sim_ECM == True)
    cells              Instance of Cells module
    p                  Instance of parameters

    Returns
    ---------
    y              Data values from 0 to 1 defining a periodic signal.

    scalar         Null output corresponding to spatial variation.

    """

    scalar = 1

    if p.f_scan_properties['f slope'] is None and p.run_sim is True:

        p.f_scan_properties['f slope'] = (p.f_scan_properties['f stop'] -
                                         p.f_scan_properties['f start'])/p.sim_end

    f_slope = p.f_scan_properties['f slope']

    y = lambda t: np.sin(np.pi*(f_slope*t + p.f_scan_properties['f start'])*t)**2

    return scalar, y

def gradient_x(pc,cells,p):
    """
    Creates a spatial gradient along the x-axis from 0 to 1 over a patch
    of cells defined by indices.

    If p.sim_ECM is True, the data is defined and returned on membrane midpoints,
    otherwise cell centres are the defining structure.

    Parameters
    ----------
    pc                 Indices of cells (p.sim_ECM == False) or membranes (p.sim_ECM == True)
    cells              Instance of Cells module
    p                  Instance of parameters

    Returns
    ---------
    fx              Data values from 0 to 1 defining a gradient in the x-direciton.
                    If p.sim_ECM == True, length of fx is number of membranes, otherwise
                    length of fx is equal to cell number.

    dynamics        Null quantity corresponding to time dynamics.
    """

    if p.sim_ECM is False:

        fx = np.zeros(len(cells.cell_i))

        fx[pc] = np.abs(cells.cell_centres[pc,0] - p.gradient_x_properties['offset'])

        fx = (fx/fx.max())*p.gradient_x_properties['slope']

    else:

        fx = np.zeros(len(cells.mem_i))

        fx[pc] = np.abs(cells.mem_mids_flat[pc,0] - p.gradient_x_properties['offset'])

        fx = (fx/fx.max())*p.gradient_x_properties['slope']

    dynamics = lambda t: 1

    return fx, dynamics

def gradient_y(pc, cells,p):
    """
    Creates a spatial gradient along the y-axis from 0 to 1 over a patch
    of cells defined by indices.

    If p.sim_ECM is True, the data is defined and returned on membrane midpoints,
    otherwise cell centres are the defining structure.

    Parameters
    ----------
    pc                 Indices of cells (p.sim_ECM == False) or membranes (p.sim_ECM == True)
    cells              Instance of Cells module
    p                  Instance of parameters

    Returns
    ---------
    fy              Data values from 0 to 1 defining a gradient in the y-direciton.
                    If p.sim_ECM == True, length of fy is number of membranes, otherwise
                    length of fy is equal to cell number.
    """

    if p.sim_ECM is False:

        fy = np.zeros(len(cells.cell_i))

        fy[pc] = np.abs(cells.cell_centres[pc,1] - p.gradient_x_properties['offset'])

        fy = (fy/fy.max())*p.gradient_x_properties['slope']

    else:

        fy = np.zeros(len(cells.mem_i))

        fy[pc] = np.abs(cells.mem_mids_flat[pc,1] - p.gradient_x_properties['offset'])

        fy = (fy/fy.max())*p.gradient_x_properties['slope']

    dynamics = lambda t: 1

    return fy, dynamics

def gradient_r(pc, cells,p):

    """
        Creates a spatial gradient in a radial direction from 0 to 1 over a patch
        of cells defined by indices.

        If p.sim_ECM is True, the data is defined and returned on membrane midpoints,
        otherwise cell centres are the defining structure.

        Parameters
        ----------
        pc                 Indices of cells (p.sim_ECM == False) or membranes (p.sim_ECM == True)
        cells              Instance of Cells module
        p                  Instance of parameters

        Returns
        ---------
        r              Data values from 0 to 1 defining a gradient in the radial direction.
                        If p.sim_ECM == True, length of r is number of membranes, otherwise
                        length of r is equal to cell number.

    """
    if p.sim_ECM is False:

        fx = np.zeros(len(cells.cell_i))
        fy = np.zeros(len(cells.cell_i))

        fx[pc] = cells.cell_centres[pc,0] - cells.centre[0] - p.gradient_r_properties['offset']
        fy[pc] = cells.cell_centres[pc,1] - cells.centre[1] - p.gradient_r_properties['offset']

        r = np.sqrt(fx**2 + fy**2)

        r = r/r.max()

        r = r*p.gradient_r_properties['slope']

    else:

        fx = np.zeros(len(cells.mem_i))
        fy = np.zeros(len(cells.mem_i))

        fx[pc] = cells.mem_mids_flat[pc,0] - cells.centre[0] - p.gradient_r_properties['offset']
        fy[pc] = cells.mem_mids_flat[pc,1] - cells.centre[1] - p.gradient_r_properties['offset']

        r = np.sqrt(fx**2 + fy**2)

        r = r/r.max()

        r = r*p.gradient_r_properties['slope']

    dynamics = lambda t: 1


    return r, dynamics