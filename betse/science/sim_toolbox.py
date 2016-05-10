#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# FIXME the NaKATPase has been changed, but all other pumps need updating
# FIXME all pumps should now take and return cATP, cADP and cP as parameters (or a "metabo" object)

import numpy as np
import numpy.ma as ma
from scipy import interpolate as interp
from scipy.ndimage.filters import gaussian_filter

from betse.exceptions import BetseExceptionSimulation


# Toolbox of functions used in the Simulator class to calculate key bioelectric properties.

def electroflux(cA,cB,Dc,d,zc,vBA,T,p,rho=1):
    """
    Electro-diffusion between two connected volumes. Note for cell work, 'b' is
    'inside', 'a' is outside, with a positive flux moving from a to b. The
    voltage is defined as Vb - Va (Vba), which is equivalent to Vmem.

    This function defaults to regular diffusion if Vba == 0.0.

    This function takes numpy matrix values as input. All inputs must be
    matrices of the same shape.

    This is the Goldman Flux/Current Equation (not to be confused with the
    Goldman Equation).

    Parameters
    ----------
    cA          concentration in region A [mol/m3] (out)
    cB          concentration in region B [mol/m3] (in)
    Dc          Diffusion constant of c  [m2/s]
    d           Distance between region A and region B [m]
    zc          valence of ionic species c
    vBA         voltage difference between region B (in) and A (out) = Vmem
    p           an instance of the Parameters class

    Returns
    --------
    flux        Chemical flux magnitude between region A and B [mol/s]

    """

    # modify the diffusion constant by the membrane density
    Dc = rho*Dc

    alpha = (zc*vBA*p.F)/(p.R*T)

    exp_alpha = np.exp(-alpha)

    deno = 1 - exp_alpha   # calculate the denominator for the electrodiffusion equation,..

    izero = (deno==0).nonzero()     # get the indices of the zero and non-zero elements of the denominator
    inotzero = (deno!=0).nonzero()

    # initialize data matrices to the same shape as input data
    flux = np.zeros(deno.shape)

    if len(deno[izero]):   # if there's anything in the izero array:
         # calculate the flux for those elements as standard diffusion [mol/m2s]:
        flux[izero] = -(Dc[izero]/d[izero])*(cB[izero] - cA[izero])

    if len(deno[inotzero]):   # if there's any indices in the inotzero array:

        # calculate the flux for those elements:
        flux[inotzero] = -((Dc[inotzero]*alpha[inotzero])/d[inotzero])*((cB[inotzero] -
                        cA[inotzero]*exp_alpha[inotzero])/deno[inotzero])

    return flux

def pumpNaKATP(cNai,cNao,cKi,cKo,Vm,T,p,block):

    """
    Parameters
    ----------
    cNai            Concentration of Na+ inside the cell
    cNao            Concentration of Na+ outside the cell
    cKi             Concentration of K+ inside the cell
    cKo             Concentration of K+ outside the cell
    Vm              Voltage across cell membrane [V]
    p               An instance of Parameters object


    Returns
    -------
    f_Na            Na+ flux (into cell +)
    f_K             K+ flux (into cell +)
    """

    deltaGATP_o = p.deltaGATP  # standard free energy of ATP hydrolysis reaction in J/(mol K)

    # At the moment the concentrations are fixed as the metabolism module isn't ready yet.
    cATP = p.cATP  # concentration of ATP in mmol/L
    cADP = p.cADP  # concentration of ADP in mmol/L
    cPi = p.cPi  # concentration of Pi in mmol/L

    # calculate the reaction coefficient Q:
    Qnumo = cADP * cPi * (cNao ** 3) * (cKi ** 2)
    Qdenomo = cATP * (cNai ** 3) * (cKo ** 2)

    # ensure no chance of dividing by zero:
    inds_Z = (Qdenomo == 0.0).nonzero()
    Qdenomo[inds_Z] = 1.0e-10

    Q = Qnumo / Qdenomo

    # calculate the equilibrium constant for the pump reaction:
    Keq = np.exp(-deltaGATP_o / (p.R * T) + ((p.F * Vm) / (p.R * T)))

    # calculate the reaction rate coefficient
    alpha = block * p.alpha_NaK * (1 - (Q / Keq))

    # calculate the enzyme coefficient:
    numo_E = ((cNai/p.KmNK_Na)**3) * ((cKo/p.KmNK_K)**2) * (cATP/p.KmNK_ATP)
    denomo_E = (1 + (cNai/p.KmNK_Na)**3)*(1+(cKo/p.KmNK_K)**2)*(1+(cATP/p.KmNK_ATP))

    f_Na = -alpha * (numo_E / denomo_E)  # flux as [mol/m2s]   scaled to concentrations Na in and K out

    f_K = -(2/3)*f_Na          # flux as [mol/m2s]

    return f_Na, f_K, -f_Na  # FIXME get rid of this return of extra -f_Na!!

def pumpCaATP(cCai,cCao,Vm,T,p):

    """
    Parameters
    ----------
    cCai            Concentration of Ca2+ inside the cell
    cCao            Concentration of Ca2+ outside the cell
    voli            Volume of the cell [m3]
    volo            Volume outside the cell [m3]
    Vm              Voltage across cell membrane [V]
    p               An instance of Parameters object


    Returns
    -------
    cCai2           Updated Ca2+ inside cell
    cCao2           Updated Ca2+ outside cell
    f_Ca            Ca2+ flux (into cell +)
    """


    deltaGATP_o = p.deltaGATP

    # At the moment the concentrations are fixed as the metabolism module isn't ready yet.
    cATP = p.cATP  # concentration of ATP in mmol/L
    cADP = p.cADP  # concentration of ADP in mmol/L
    cPi = p.cPi  # concentration of Pi in mmol/L

    # calculate the reaction coefficient Q:
    Qnumo = cADP * cPi * cCao
    Qdenomo = cATP * cCai

    # ensure no chance of dividing by zero:
    inds_Z = (Qdenomo == 0.0).nonzero()
    Qdenomo[inds_Z] = 1.0e-12

    Q = Qnumo / Qdenomo

    # calculate the equilibrium constant for the pump reaction:
    Keq = np.exp(-deltaGATP_o / (p.R * T) + 2*((p.F * Vm) / (p.R * T)))

    # calculate the reaction rate coefficient
    alpha = p.alpha_Ca * (1 - (Q / Keq))

    # calculate the enzyme coefficient:
    numo_E = (cCai/p.KmCa_Ca) * (cATP/p.KmCa_ATP)
    denomo_E = (1 + (cCai/p.KmCa_Ca)) * (1+ (cATP/p.KmCa_ATP))

    f_Ca = -alpha * (numo_E / denomo_E)  # flux as [mol/m2s]

    return f_Ca

def pumpCaER(cCai,cCao,Vm,T,p):  # FIXME this should be replaced and use only pumpCaATP, defined above!
    """
    Pumps calcium out of the cell and into the endoplasmic reticulum.
    Vm is the voltage across the endoplasmic reticulum membrane.

    """

    deltaGATP = 20*p.R*T

    delG_Ca = p.R*T*np.log(cCai/cCao) + 2*p.F*Vm
    delG_CaATP = deltaGATP - (delG_Ca)
    delG_pump = (delG_CaATP/1000)

    alpha = p.alpha_CaER*(delG_pump - p.halfmax_Ca)

    f_Ca  = alpha*(cCao)      #flux as [mol/s]

    return f_Ca

def pumpHKATP(cHi,cHo,cKi,cKo,Vm,T,p,block):

    """
    Parameters
    ----------
    cHi            Concentration of H+ inside the cell
    cHo            Concentration of H+ outside the cell
    cKi             Concentration of K+ inside the cell
    cKo             Concentration of K+ outside the cell
    voli            Volume of the cell [m3]
    volo            Volume outside the cell [m3]
    Vm              Voltage across cell membrane [V]
    p               An instance of Parameters object


    Returns
    -------
    cHi2           Updated H+ inside cell
    cHo2           Updated H+ outside cell
    cKi2            Updated K+ inside cell
    cKo2            Updated K+ outside cell
    f_Na            Na+ flux (into cell +)
    f_K             K+ flux (into cell +)
    """

    deltaGATP_o = p.deltaGATP

    # At the moment the concentrations are fixed as the metabolism module isn't ready yet.
    cATP = p.cATP  # concentration of ATP in mmol/L
    cADP = p.cADP  # concentration of ADP in mmol/L
    cPi = p.cPi  # concentration of Pi in mmol/L

    # calculate the reaction coefficient Q:
    Qnumo = cADP * cPi * (cHo) * (cKi)
    Qdenomo = cATP * (cHi) * (cKo)

    # ensure no chance of dividing by zero:
    inds_Z = (Qdenomo == 0.0).nonzero()
    Qdenomo[inds_Z] = 1.0e-10

    Q = Qnumo / Qdenomo

    # calculate the equilibrium constant for the pump reaction:
    Keq = np.exp(-deltaGATP_o / (p.R * T))

    # calculate the reaction rate coefficient
    alpha = block * p.alpha_HK * (1 - (Q / Keq))

    # calculate the enzyme coefficient:
    numo_E = (cHi / p.KmHK_H) * (cKo / p.KmHK_K) * (cATP / p.KmHK_ATP)
    denomo_E = (1 + (cHi / p.KmHK_H)) *(1+ (cKo / p.KmHK_K)) *(1+ (cATP / p.KmHK_ATP))


    f_H  = -alpha*(numo_E/denomo_E)      #flux as [mol/s], scaled by concentrations in and out

    f_K = -f_H          # flux as [mol/s]

    return f_H, f_K

def pumpVATP(cHi,cHo,Vm,T,p,block):

    deltaGATP_o = p.deltaGATP

    # At the moment the concentrations are fixed as the metabolism module isn't ready yet.
    cATP = 1.5  # concentration of ATP in mmol/L
    cADP = 0.2  # concentration of ADP in mmol/L
    cPi = 0.5  # concentration of Pi in mmol/L

    # calculate the reaction coefficient Q:
    Qnumo = cADP * cPi * cHo
    Qdenomo = cATP * cHi

    # ensure no chance of dividing by zero:
    inds_Z = (Qdenomo == 0.0).nonzero()
    Qdenomo[inds_Z] = 1.0e-10

    Q = Qnumo / Qdenomo

    # calculate the equilibrium constant for the pump reaction:
    Keq = np.exp(-deltaGATP_o / (p.R * T) + ((p.F * Vm) / (p.R * T)))

    # calculate the reaction rate coefficient
    alpha = block * p.alpha_V * (1 - Q / Keq)

    # calculate the enzyme coefficient:
    numo_E = (cHi/p.KmV_H) * (cATP/p.KmV_ATP)
    denomo_E = (1 + (cHi/p.KmV_H)) * (1 + (cATP/p.KmV_ATP))

    f_H = -alpha * (numo_E / denomo_E)  # flux as [mol/m2s]

    return f_H

def pumpATPSynth(cHi,cHo,cATP,cADP,cPi,Vm,T,p):
    """
    Calculates H+ flux into the mitochondrial matrix (+
    flow is into the mitochondria), which is equivalent
    to ATP created and ADP consumed in the ATP Synthase
    enzyme reaction of the mitochondria.

    Parameters
    ------------
    cHi      Hydrogen concentration in the mitocondrial matrix
    cHo      Hydrogen concentration outside of the mitochondria
    cATP     ATP concentration in the mitochondrial matrix
    cADP     ADP concentration in the mitochondrial matrix
    cPi      Phosphate concentration in the mitochondrial matrix
    Vm        Mitochondrial membrane potential
    T         Temperature
    p         Parameters object p

    Returns
    -------
    f_h      hydrogen flux into the matrix and ATP synthesis rate

    """

    deltaGATP_o = p.deltaGATP

    # calculate the reaction coefficient Q:
    Qnumo = cATP * (cHi ** 3)
    Qdenomo = cADP * cPi * (cHo ** 3)

    # ensure no chance of dividing by zero:
    inds_Z = (Qdenomo == 0.0).nonzero()
    Qdenomo[inds_Z] = 1.0e-6

    Q = Qnumo / Qdenomo

    # calculate the equilibrium constant for the pump reaction:
    Keq = np.exp(deltaGATP_o / (p.R * T) - 3 * ((p.F * Vm) / (p.R * T)))

    # calculate the reaction rate coefficient
    alpha = p.alpha_AS * (1 - Q / Keq)

    # calculate the enzyme coefficient:
    numo_E = (cHo / p.KmAS_H) * (cADP / p.KmAS_ADP) * (cPi / p.KmAS_P)
    denomo_E = (1 + (cHo / p.KmAS_H)) *(1+ (cADP / p.KmAS_ADP)) *(1 + (cPi / p.KmAS_P))

    f_H = alpha * (numo_E / denomo_E)  # flux as [mol/m2s]

    return f_H

def get_volt(q,sa,p):

    """
    Calculates the voltage for a net charge on a capacitor.

    Parameters
    ----------
    q           Net electrical charge [C]
    sa          Surface area [m2]

    Returns
    -------
    V               Voltage on the capacitive space holding charge

    """

    cap = sa*p.cm
    V = (1/cap)*q

    return V

def get_charge(concentrations,zs,vol,p):
    """
    Calculates the total charge in a space given a set of concentrations
    and their ionic charges, along with the space volume.

    Parameters
    ---------------
    concentrations       An array of arrays listing concentrations of different ions in multiple spaces [mol/m2].
    zs                   An array of ion charge states
    vol                  Volume of the spaces [m3]
    p                    Instance of Parameters module object

    Returns
    ----------------
    netcharge            The unbalanced charge in the space or spaces [C]

    """

    netcharge = np.sum(p.F*vol*concentrations*zs,axis=0)

    return netcharge

def get_charge_density(concentrations,zs,p):

    """
    Calculates the charge density given ion concentrations in an array of spaces.

    Parameters
    --------------
    concentrations:  an array of array of concentration of ions in spaces [mol/m3]
    zs:              valence of each ion
    p:               Parameters object instance

    Returns
    -------------
    netcharge     the net charge density in spaces C/m3
    """

    netcharge = np.sum(p.F*zs*concentrations, axis=0)

    return netcharge

def get_molarity(concentrations,p):

    q = 0

    for conc in concentrations:
        q = q+ conc

    netmolarity = q

    return netmolarity

def cell_ave(cells,vm_at_mem):

    """
    Averages Vmem over membrane domains to return a mean value for each cell

    Parameters
    ----------
    cells               An instance of the Cells module cells object
    vm_at_mem           Vmem at individual membrane domains


    Returns
    --------
    v_cell              Cell Vm averaged over the whole cell

    """

    v_cell = []

    for i in cells.cell_i:
        cellinds = (cells.mem_to_cells == i).nonzero()
        v_cell_array = vm_at_mem[cellinds]
        v_cell_ave = np.mean(v_cell_array)
        v_cell.append(v_cell_ave)

    v_cell = np.asarray(v_cell)

    return v_cell

def check_v(vm):
    """
    Does a quick check on Vmem values
    and displays error warning or exception if the value
    indicates the simulation is unstable.

    """


    isnans = np.isnan(vm)

    if isnans.any():  # if there's anything in the isubzeros matrix...
        raise BetseExceptionSimulation("Your simulation has become unstable. Please try a smaller time step,"
                                       "reduce gap junction radius, and/or reduce pump rate coefficients.")

def vertData(data, cells, p):
    """
    Interpolate data from midpoints to verts
    and sample it over a uniform grid.
    Produces data suitable for 2D mesh and
    streamline plots.

    Parameters
    -----------
    data          A numpy vector of data points on cell membrane mids
    cells         An instance of the Cells object
    p             An instance of the Parameters object

    Returns
    ----------
    dat_grid      THe data sampled on a uniform grid
    """

    # interpolate vmem defined on mem mids to cell vertices:
    verts_data = np.dot(data,cells.matrixMap2Verts)

    # amalgamate both mem mids and verts data into one stack:
    plot_data = np.hstack((data,verts_data))

    # interpolate the stack to the plotting grid:
    dat_grid = interp.griddata((cells.plot_xy[:,0],cells.plot_xy[:,1]),plot_data,(cells.Xgrid,cells.Ygrid),
                               method=p.interp_type,
                               fill_value=0)

    # # smooth out the data a bit:
    dat_grid = gaussian_filter(dat_grid,1)

    # get rid of values that bleed into the environment:
    # dat_grid = np.multiply(dat_grid,cells.maskM)

    dat_grid = ma.masked_array(dat_grid, np.logical_not(cells.maskM))

    return dat_grid

def nernst_planck_flux(c, gcx, gcy, gvx, gvy,ux,uy,D,z,T,p):
    """
     Calculate the flux component of the Nernst-Planck equation

     Parameters
     ------------

    c:     concentration
    gcx:   concentration gradient, x component
    gcy:   concentration gradient, y component
    gvx:   voltage gradient, x component
    gvy:   voltage gradient, y component
    ux:    fluid velocity, x component
    uy:    fluid velocity, y component
    D:     diffusion constant, D
    z:     ion charge
    T:     temperature
    p:     parameters object

    Returns
    --------
    fx, fx        mass flux in x and y directions
    """

    alpha = (D*z*p.q)/(p.kb*T)
    fx =  -D*gcx - alpha*gvx*c + ux*c
    fy =  -D*gcy - alpha*gvy*c + uy*c

    return fx, fy

def np_flux_special(cx,cy,gcx,gcy,gvx,gvy,ux,uy,Dx,Dy,z,T,p):
    """
     Calculate the flux component of the Nernst-Planck equation on a MACs grid

     Parameters
     ------------

    c:     concentration
    gcx:   concentration gradient, x component
    gcy:   concentration gradient, y component
    gvx:   voltage gradient, x component
    gvy:   voltage gradient, y component
    ux:    fluid velocity, x component
    uy:    fluid velocity, y component
    D:     diffusion constant, D
    z:     ion charge
    T:     temperature
    p:     parameters object

    Returns
    --------
    fx, fx        mass flux in x and y directions
    """

    alphax = (Dx*z*p.q)/(p.kb*T)
    alphay = (Dy*z*p.q)/(p.kb*T)

    fx =  -Dx*gcx - alphax*gvx*cx + cx*ux

    fy =  -Dy*gcy - alphay*gvy*cy + cy*uy

    return fx, fy

def no_negs(data):

    # ensure no NaNs:
    inds_nan = (np.isnan(data)).nonzero()
    data[inds_nan] = 0

    # ensure that data has no less than zero values:
    inds_neg = (data < 0).nonzero()
    data[inds_neg] = 0

    # if len(inds_nan[0]) > 0 or len(inds_neg[0]) > 0:
    #     loggers.log_info("Warning: invalid value (0 or Nan) found in concentration data.")

    return data

def bicarbonate_buffer(cH, cCO2, cHCO3, p):
    """
    This most amazing buffer handles influx of H+,
    HCO3-, H2CO3 (from dissolved carbon dioxide) to
    handle pH in real time.

    Uses the bicarbonate dissacociation reaction:

    H2CO3 ----> HCO3 + H

    Where all dissolved carbon dioxide is assumed
    converted to carbonic acid via carbonic anhydrase enzyme.

    """

    Q = (cH*cHCO3)/(cCO2)

    v_ph = p.vm_ph*(cCO2/(1+cCO2))*(1 - (Q/p.Keqm_ph))

    cCO2 = cCO2 - v_ph*p.dt
    cHCO3 = cHCO3 + v_ph*p.dt
    cH = cH + v_ph*p.dt

    pH = -np.log10(cH*1e-3)

    return cH, cCO2, cHCO3, pH


def ghk_calculator(sim, cells, p):
    """
    Uses simulation parameters in the Goldman (GHK) equation
    to calculate an alternative Vmem for validation purposes.

    """

    # begin by initializing all summation arrays for the cell network:
    sum_PmAnion_out = np.zeros(len(cells.cell_i))
    sum_PmAnion_in = np.zeros(len(cells.cell_i))
    sum_PmCation_out = np.zeros(len(cells.cell_i))
    sum_PmCation_in = np.zeros(len(cells.cell_i))

    for i, z in enumerate(sim.zs):

        ion_type = np.sign(z)

        if ion_type == -1:

            if p.sim_ECM is True:

                Dm = np.dot(cells.M_sum_mems, sim.Dm_cells[i]) / cells.num_mems

                sum_PmAnion_in = sum_PmAnion_in + Dm * sim.cc_cells[i] * (1 / p.tm)
                sum_PmAnion_out = sum_PmAnion_out + Dm * sim.cc_env[i][cells.map_cell2ecm] * (1 / p.tm)

            else:
                sum_PmAnion_in = sum_PmAnion_in + sim.Dm_cells[i] * sim.cc_cells[i] * (1 / p.tm)
                sum_PmAnion_out = sum_PmAnion_out + sim.Dm_cells[i] * sim.cc_env[i] * (1 / p.tm)

        if ion_type == 1:

            if p.sim_ECM is True:

                Dm = np.dot(cells.M_sum_mems, sim.Dm_cells[i]) / cells.num_mems

                sum_PmCation_in = sum_PmCation_in + Dm * sim.cc_cells[i] * (1 / p.tm)
                sum_PmCation_out = sum_PmCation_out + Dm * sim.cc_env[i][cells.map_cell2ecm] * (1 / p.tm)

            else:
                sum_PmCation_in = sum_PmCation_in + sim.Dm_cells[i] * sim.cc_cells[i] * (1 / p.tm)
                sum_PmCation_out = sum_PmCation_out + sim.Dm_cells[i] * sim.cc_env[i] * (1 / p.tm)

    sim.vm_GHK = ((p.R * sim.T) / p.F) * np.log(
        (sum_PmCation_out + sum_PmAnion_in) / (sum_PmCation_in + sum_PmAnion_out))
