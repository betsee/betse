#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS                            }....................
import numpy as np
from scipy import interpolate
from betse.lib.pil import pilnumpy
from betse.lib.pil.pilnumpy import ImageModeType
from betse.util.path import pathnames
from scipy.ndimage.filters import gaussian_filter
from betse.science.math import finitediff as fd

# ....................{ IMPORTS                            }....................
def periodic(pc,cells,p):
    """
    Creates a time-dependent signal of a sinusoidal function.

    Parameters
    ----------
    pc                 Indices of cells (p.is_ecm == False) or membranes (p.is_ecm == True)
    cells              Instance of Cells module
    p                  Instance of parameters

    Returns
    ---------
    y              Data values from 0 to 1 defining a periodic signal.
    scalar         Null output corresponding to spatial variation.
    """

    scalar = 1

    y = lambda t: np.sin(
        t*np.pi*(p.periodic_properties['frequency']/1) +
        p.periodic_properties['phase']
    ) ** 2

    return scalar, y

def f_sweep(pc,cells,p):
    """
    Creates a time-dependent sweep of a sinusoidal function through various
    frequencies.

    Parameters
    ----------
    pc                 Indices of cells (p.is_ecm == False) or membranes (p.is_ecm == True)
    cells              Instance of Cells module
    p                  Instance of parameters

    Returns
    ---------
    y              Data values from 0 to 1 defining a periodic signal.
    scalar         Null output corresponding to spatial variation.
    """

    scalar = 1

    if p.f_scan_properties['f slope'] is None and p.run_sim is True:
        p.f_scan_properties['f slope'] = (
            p.f_scan_properties['f stop'] -
            p.f_scan_properties['f start']) / p.sim_time_total

    f_slope = p.f_scan_properties['f slope']

    y = lambda t: np.sin(
        np.pi*(f_slope*t + p.f_scan_properties['f start'])*t)**2

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
    fx
        Data values from 0 to 1 defining a gradient in the x-direciton. If
        p.is_ecm == True, length of fx is number of membranes, otherwise length
        of fx is equal to cell number.
    dynamics
        Null quantity corresponding to time dynamics.
    """

    if len(pc) == len(cells.mem_i):
        x_vals_o = cells.mem_mids_flat[:, 0]
    elif len(pc) == len(cells.cell_i):
        x_vals_o = cells.cell_centres[:, 0]

    min_grad = x_vals_o.min()

    # shift values to have a minimum of zero:
    x_vals = x_vals_o - min_grad

    # scale values to have a maximum of 1:
    x_vals = x_vals/x_vals.max()

    # shift the x_mid value by the offset:
    x_mid = 0.5 + p.gradient_x_properties['x-offset']

    n = p.gradient_x_properties['exponent']
    grad_slope = ((x_vals/x_mid)**n)/(1+((x_vals/x_mid)**n))
    fx = (grad_slope)*p.gradient_x_properties['slope'] + p.gradient_x_properties['z-offset']

    dynamics = lambda t: 1

    return fx, dynamics

def gradient_y(pc, cells,p):
    """
    Creates a spatial gradient along the y-axis from 0 to 1 over a patch
    of cells defined by indices.

    If p.is_ecm is True, the data is defined and returned on membrane
    midpoints, otherwise cell centres are the defining structure.

    Parameters
    ----------
    pc                 Indices of membranes
    cells              Instance of Cells module
    p                  Instance of parameters

    Returns
    ---------
    fy
        Data values from 0 to 1 defining a gradient in the y-direciton. If
        p.is_ecm == True, length of fy is number of membranes, otherwise length
        of fy is equal to cell number.
    dynamics
        Null quantity corresponding to time dynamics.
    """

    if len(pc) == len(cells.mem_i):
        y_vals_o = cells.mem_mids_flat[:, 1]
    elif len(pc) == len(cells.cell_i):
        y_vals_o = cells.cell_centres[:, 1]

    min_grad = y_vals_o.min()

    # shift values to have a minimum of zero:
    y_vals = y_vals_o - min_grad

    # scale values to have a maximum of 1:
    y_vals = y_vals / y_vals.max()

    # shift the y_mid value by the offset:
    y_mid = 0.5 + p.gradient_y_properties['x-offset']

    n = p.gradient_y_properties['exponent']
    grad_slope = ((y_vals / y_mid) ** n) / (1 + ((y_vals / y_mid) ** n))
    fy = (grad_slope) * p.gradient_y_properties['slope'] + p.gradient_y_properties['z-offset']

    dynamics = lambda t: 1

    return fy, dynamics

def gradient_r(pc, cells,p):
    """
    Creates a spatial gradient in a radial direction from 0 to 1 over a patch
    of cells defined by indices.

    If p.is_ecm is True, the data is defined and returned on membrane
    midpoints, otherwise cell centres are the defining structure.

    Parameters
    ----------
    pc                 Indices of membranes
    cells              Instance of Cells module
    p                  Instance of parameters

    Returns
    ---------
    r
        Data values from 0 to 1 defining a gradient in the radial direction.  If
        p.is_ecm == True, length of r is number of membranes, otherwise length
        of r is equal to cell number.
    """

    fx = np.zeros(len(cells.mem_i))
    fy = np.zeros(len(cells.mem_i))

    fx[pc] = (
        cells.mem_mids_flat[pc,0] -
        cells.centre[0] -
        p.gradient_r_properties['x-offset'])
    fy[pc] = (
        cells.mem_mids_flat[pc,1] -
        cells.centre[1] -
        p.gradient_r_properties['x-offset'])

    r = np.sqrt(fx**2 + fy**2)

    r = r/r.max()

    r = r*p.gradient_r_properties['slope'] + p.gradient_r_properties['z-offset']

    dynamics = lambda t: 1

    return r, dynamics

def gradient_bitmap(pc, cells, p, bitmap_filename = None):

    """
    This modulator reads in a bitmap supplied by the user from the
    directory in params.

    Red is treated as positive, blue is treated as negative, and
    black is treated as zero.

    Parameters
    ------------

    pc:   array of cell centre or membrane midpoints
    cells:  BETSE cells object
    p:      BETSE parameters object

    """

    if bitmap_filename is None:  # default this to loading the 'gradient bitmap function' from params
        bitmap_filename = p.grad_bm_fn
        grad_bm_offset = np.max((p.grad_bm_offset, 0.0))

    else:
        grad_bm_offset = 0.0


    if len(pc) == len(cells.mem_i):
        xmap = cells.map_mem2ecm

    elif len(pc) == len(cells.cell_i):
        xmap = cells.map_cell2ecm

    # xmi = cells.xmin - 4*p.cell_radius
    # ymi = cells.ymin - 4*p.cell_radius
    # xma = cells.xmax + 4*p.cell_radius
    # yma = cells.ymax + 4*p.cell_radius

    xmi = cells.xmin
    ymi = cells.ymin
    xma = cells.xmax
    yma = cells.ymax

    xx = np.linspace(cells.xmin, cells.xmax, cells.X.shape[1])
    yy = np.linspace(cells.ymin, cells.ymax, cells.X.shape[0])

    fn1 = pathnames.join(p.conf_dirname, bitmap_filename)

    # Three-dimensional Numpy array of the RGB-ordered integer components of all
    # pixels loaded from this image.
    a1o = pilnumpy.load_image(filename=fn1, mode=ImageModeType.COLOR_RGB)

    a1 = np.asarray(a1o, dtype=np.float64)

    a1_F = (a1[:, :, 0] -  a1[:, :, 2]) / 255

    a1_F = np.flipud(a1_F)
    # a1_F = fd.integrator(a1_F, sharp=0.5) # smooth a little to avoid bizarre visual effects

    xa = np.linspace(xmi, xma, a1_F.shape[1])
    ya = np.linspace(ymi, yma, a1_F.shape[0])

    spline_F = interpolate.interp2d(xa, ya, a1_F, kind='linear', fill_value=0.0)
    fe = spline_F(xx, yy)

    fe = fd.integrator(fe, sharp=0.5) # smooth a little to avoid bizarre visual effects

    f = fe.ravel()[xmap]

    # f = (f/f.max()) + grad_bm_offset
    f = f + grad_bm_offset

    # indz = (f < 0.0).nonzero()
    # f[indz] = 0.0
    dynamics = lambda t:1

    return f, dynamics
