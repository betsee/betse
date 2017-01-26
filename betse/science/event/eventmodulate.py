#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
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

def gradient_x(pc, cells,p):
    """
    Creates a spatial gradient along the x-axis from 0 to 1 over a patch
    of cells defined by indices.

    Parameters
    ----------
    pc                 Indices of membranes or cells
    cells              Instance of Cells module
    p                  Instance of parameters
    targs              Indices of targets within mem or cell inds to apply gradient

    Returns
    ---------
    fx              Data values from 0 to 1 defining a gradient in the x-direciton.
                    If p.sim_ECM == True, length of fx is number of membranes, otherwise
                    length of fx is equal to cell number.

    dynamics        Null quantity corresponding to time dynamics.
    """

    if len(pc) == len(cells.mem_i):
        fx = np.zeros(len(cells.mem_i))
        x_vals_o = cells.mem_mids_flat[:, 0]

    elif len(pc) == len(cells.cell_i):
        fx = np.zeros(len(cells.mem_i))
        x_vals_o = cells.cell_centres[:, 0]

    min_grad = x_vals_o.min()

    # shift values to have a minimum of zero:
    x_vals = x_vals_o - min_grad

    # scale values to have a maximum of 1:
    x_vals = x_vals/x_vals.max()

    # shift the x_mid value by the offset:
    x_mid = 0.5 + p.gradient_x_properties['offset']

    n = p.gradient_x_properties['exponent']

    grad_slope = ((x_vals/x_mid)**n)/(1+((x_vals/x_mid)**n))

    fx = (grad_slope)*p.gradient_x_properties['slope']

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
    pc                 Indices of membranes
    cells              Instance of Cells module
    p                  Instance of parameters

    Returns
    ---------
    fy              Data values from 0 to 1 defining a gradient in the y-direciton.
                    If p.sim_ECM == True, length of fy is number of membranes, otherwise
                    length of fy is equal to cell number.
    """

    if len(pc) == len(cells.mem_i):
        fy = np.zeros(len(cells.mem_i))
        y_vals_o = cells.mem_mids_flat[:, 1]

    elif len(pc) == len(cells.cell_i):
        fy = np.zeros(len(cells.mem_i))
        y_vals_o = cells.cell_centres[:, 1]

    min_grad = y_vals_o.min()

    # shift values to have a minimum of zero:
    y_vals = y_vals_o - min_grad

    # scale values to have a maximum of 1:
    y_vals = y_vals / y_vals.max()

    # shift the y_mid value by the offset:
    y_mid = 0.5 + p.gradient_y_properties['offset']

    n = p.gradient_y_properties['exponent']

    grad_slope = ((y_vals / y_mid) ** n) / (1 + ((y_vals / y_mid) ** n))

    fy = (grad_slope) * p.gradient_y_properties['slope']

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
        pc                 Indices of membranes
        cells              Instance of Cells module
        p                  Instance of parameters

        Returns
        ---------
        r              Data values from 0 to 1 defining a gradient in the radial direction.
                        If p.sim_ECM == True, length of r is number of membranes, otherwise
                        length of r is equal to cell number.

    """

    fx = np.zeros(len(cells.mem_i))
    fy = np.zeros(len(cells.mem_i))

    fx[pc] = cells.mem_mids_flat[pc,0] - cells.centre[0] - p.gradient_r_properties['offset']
    fy[pc] = cells.mem_mids_flat[pc,1] - cells.centre[1] - p.gradient_r_properties['offset']

    r = np.sqrt(fx**2 + fy**2)

    r = r/r.max()

    r = r*p.gradient_r_properties['slope']

    dynamics = lambda t: 1


    return r, dynamics