#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# FIXME all pumps should now take and return cATP, cADP and cP as parameters (or a "metabo" object)

import numpy as np
import numpy.ma as ma
from scipy import interpolate as interp
from scipy.ndimage.filters import gaussian_filter

from betse.exceptions import BetseSimInstabilityException
from betse.science.math import finitediff as fd


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

    # # modify the diffusion constant by the membrane density
    # Dc = rho*Dc
    #
    # alpha = (zc*vBA*p.F)/(p.R*T)
    #
    # exp_alpha = np.exp(-alpha)
    #
    # deno = 1 - exp_alpha + 1e-15  # calculate the denominator with small term in case it's equal to zero
    #
    # # calculate the flux for those elements:
    # flux = -((Dc*alpha)/d)*((cB - cA*exp_alpha)/deno)

    # modify the diffusion constant by the membrane density
    # print("cB", cB.mean(), cB.max(), cB.min())


    # Dc = rho*Dc

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

    flux = flux*rho


    return flux

def pumpNaKATP(cNai,cNao,cKi,cKo,Vm,T,p,block, met = None):

    """
    Parameters
    ----------
    cNai            Concentration of Na+ inside the cell
    cNao            Concentration of Na+ outside the cell
    cKi             Concentration of K+ inside the cell
    cKo             Concentration of K+ outside the cell
    Vm              Voltage across cell membrane [V]
    p               An instance of Parameters object

    met             A "metabolism" vector containing concentrations of ATP, ADP and Pi


    Returns
    -------
    f_Na            Na+ flux (into cell +)
    f_K             K+ flux (into cell +)
    """

    deltaGATP_o = p.deltaGATP  # standard free energy of ATP hydrolysis reaction in J/(mol K)

    if met is None:

        # if metabolism vector not supplied, use singular defaults for concentrations
        cATP = p.cATP
        cADP = p.cADP
        cPi  = p.cPi

    else:

        cATP = met['cATP']  # concentration of ATP in mmol/L
        cADP = met['cADP']  # concentration of ADP in mmol/L
        cPi = met['cPi']  # concentration of Pi in mmol/L


    # calculate the reaction coefficient Q:
    Qnumo = cADP * cPi * (cNao ** 3) * (cKi ** 2)
    Qdenomo = cATP * (cNai ** 3) * (cKo ** 2)

    # ensure no chance of dividing by zero:
    inds_Z = (Qdenomo == 0.0).nonzero()
    Qdenomo[inds_Z] = 1.0e-15

    Q = Qnumo / Qdenomo


    # calculate the equilibrium constant for the pump reaction:
    Keq = np.exp(-(deltaGATP_o / (p.R * T) - ((p.F * Vm) / (p.R * T))))

    # calculate the enzyme coefficient:
    numo_E = ((cNai/p.KmNK_Na)**3) * ((cKo/p.KmNK_K)**2) * (cATP/p.KmNK_ATP)
    denomo_E = (1 + (cNai/p.KmNK_Na)**3)*(1+(cKo/p.KmNK_K)**2)*(1+(cATP/p.KmNK_ATP))

    fwd_co = numo_E/denomo_E

    f_Na = -3*block*p.alpha_NaK*fwd_co*(1 - (Q/Keq))  # flux as [mol/m2s]   scaled to concentrations Na in and K out

    f_K = -(2/3)*f_Na          # flux as [mol/m2s]

    return f_Na, f_K, -f_Na  # FIXME get rid of this return of extra -f_Na!!

def pumpCaATP(cCai,cCao,Vm,T,p, block, met = None):

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

    no_negs(cCai)
    no_negs(cCao)

    if met is None:

        # if metabolism vector not supplied, use singular defaults for concentrations
        cATP = p.cATP
        cADP = p.cADP
        cPi  = p.cPi

    else:

        cATP = met['cATP']  # concentration of ATP in mmol/L
        cADP = met['cADP']  # concentration of ADP in mmol/L
        cPi = met['cPi']  # concentration of Pi in mmol/L

    # calculate the reaction coefficient Q:
    Qnumo = cADP * cPi * cCao
    Qdenomo = cATP * cCai

    # ensure no chance of dividing by zero:
    inds_Z = (Qdenomo == 0.0).nonzero()
    Qdenomo[inds_Z] = 1.0e-16

    Q = Qnumo / Qdenomo

    # calculate the equilibrium constant for the pump reaction:
    Keq = np.exp(-(deltaGATP_o / (p.R * T) - 2*((p.F * Vm) / (p.R * T))))

    # calculate the enzyme coefficient for forward reaction:
    numo_E = (cCai/p.KmCa_Ca) * (cATP/p.KmCa_ATP)
    denomo_E = (1 + (cCai/p.KmCa_Ca)) * (1+ (cATP/p.KmCa_ATP))

    frwd = numo_E/denomo_E

    # calculate the enzyme coefficient for backward reaction:
    numo_Eb = (cCao/p.KmCa_Ca)
    denomo_Eb = (1 + (cCao/p.KmCa_Ca))

    bkwrd = numo_Eb/denomo_Eb

    f_Ca = -p.alpha_Ca*frwd*(1 - (Q/Keq))  # flux as [mol/m2s]

    return f_Ca

def pumpCaER(cCai,cCao,Vm,T,p):
    """
    Pumps calcium out of the cell and into the endoplasmic reticulum.
    Vm is the voltage across the endoplasmic reticulum membrane.

    """

    deltaGATP_o = p.deltaGATP

    cATP = p.cATP
    cADP = p.cADP
    cPi = p.cPi

    # calculate the reaction coefficient Q:
    Qnumo = cADP * cPi * cCai
    Qdenomo = cATP * cCao

    # ensure no chance of dividing by zero:
    inds_Z = (Qdenomo == 0.0).nonzero()
    Qdenomo[inds_Z] = 1.0e-16

    Q = Qnumo / Qdenomo

    # calculate the equilibrium constant for the pump reaction:
    Keq = np.exp(-deltaGATP_o / (p.R * T) - 2 * ((p.F * Vm) / (p.R * T)))


    # calculate the enzyme coefficient for forward reaction:
    numo_E = (cCao / p.KmCa_Ca) * (cATP / p.KmCa_ATP)
    denomo_E = (1 + (cCao / p.KmCa_Ca)) * (1 + (cATP / p.KmCa_ATP))

    frwd = numo_E / denomo_E

    # calculate the enzyme coefficient for backward reaction:
    numo_Eb = (cCai / p.KmCa_Ca)
    denomo_Eb = (1 + (cCai / p.KmCa_Ca))

    bkwrd = numo_Eb / denomo_Eb

    f_Ca = p.serca_max * frwd * (1 - (Q / Keq))  # flux as [mol/m2s]

    return f_Ca

def pumpHKATP(cHi,cHo,cKi,cKo,Vm,T,p,block, met = None):

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


    if met is None:

        # if metabolism vector not supplied, use singular defaults for concentrations
        cATP = p.cATP
        cADP = p.cADP
        cPi  = p.cPi

    else:

        cATP = met['cATP']  # concentration of ATP in mmol/L
        cADP = met['cADP']  # concentration of ADP in mmol/L
        cPi = met['cPi']  # concentration of Pi in mmol/L

    # calculate the reaction coefficient Q:
    Qnumo = cADP * cPi * (cHo) * (cKi)
    Qdenomo = cATP * (cHi) * (cKo)

    # ensure no chance of dividing by zero:
    inds_Z = (Qdenomo == 0.0).nonzero()
    Qdenomo[inds_Z] = 1.0e-15

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

def pumpHKATP_m(cMi,cMo,cKi,cKo,Vm,T,p,block, met = None):

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


    if met is None:

        # if metabolism vector not supplied, use singular defaults for concentrations
        cATP = p.cATP
        cADP = p.cADP
        cPi  = p.cPi

    else:

        cATP = met['cATP']  # concentration of ATP in mmol/L
        cADP = met['cADP']  # concentration of ADP in mmol/L
        cPi = met['cPi']  # concentration of Pi in mmol/L

    # calculate the reaction coefficient Q:
    Qnumo = cADP * cPi * (cMi) * (cKi)
    Qdenomo = cATP * (cMo) * (cKo)

    # ensure no chance of dividing by zero:
    inds_Z = (Qdenomo == 0.0).nonzero()
    Qdenomo[inds_Z] = 1.0e-15

    Q = Qnumo / Qdenomo

    # calculate the equilibrium constant for the pump reaction:
    Keq = np.exp(-deltaGATP_o / (p.R * T))

    # calculate the reaction rate coefficient
    alpha = block * p.alpha_HK * (1 - (Q / Keq))

    # calculate the enzyme coefficient:
    numo_E = (cMo / p.KmHK_H) * (cKo / p.KmHK_K) * (cATP / p.KmHK_ATP)
    denomo_E = (1 + (cMo / p.KmHK_H)) *(1+ (cKo / p.KmHK_K)) *(1+ (cATP / p.KmHK_ATP))

    f_M  = alpha*(numo_E/denomo_E)      #flux as [mol/s], scaled by concentrations in and out

    f_K = f_M          # flux as [mol/s]

    # print(block, f_M.mean())

    return f_M, f_K

def pumpVATP(cHi,cHo,Vm,T,p, block, met = None):

    deltaGATP_o = p.deltaGATP

    if met is None:

        # if metabolism vector not supplied, use singular defaults for concentrations
        cATP = p.cATP
        cADP = p.cADP
        cPi  = p.cPi

    else:

        cATP = met['cATP']  # concentration of ATP in mmol/L
        cADP = met['cADP']  # concentration of ADP in mmol/L
        cPi = met['cPi']  # concentration of Pi in mmol/L

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

def exch_NaCa(cNai, cNao, cCai, cCao, Vm, T, p, block):
    """
    The sodium-calcium exchanger, which removes one
    calcium from the cell in exchange for 3 sodium
    into the cell. Important for voltage gated
    calcium signals.

    3 Na_out + 1 Ca_in --> 3 Na_in + 1 Ca_out

    cNai:  Sodium inside the cell
    cNao:  Sodium outside the cell
    cCai:  Calcium inside the cell
    cCao:  Calcium outside the cell
    Vm:    Transmembrane voltage
    T:     Temperature
    p:     Instance of parameters

    Returns
    -------
    f_Na      Flux of sodium
    f_Ca      Flux of calcium

    """

    # calculate the reaction coefficient Q:
    Qnumo = (cNai ** 3) * (cCao)
    Qdenomo = (cNao ** 3) * (cCai)

    # ensure no chance of dividing by zero:
    inds_Z = (Qdenomo == 0.0).nonzero()
    Qdenomo[inds_Z] = 1.0e-10

    Q = Qnumo / Qdenomo

    # # calculate the equilibrium constant for the pump reaction:
    Keq = np.exp(-((p.F * Vm) / (p.R * T)))

    # calculate the reaction rate coefficient
    alpha = block*p.alpha_NaCaExch * (1 - (Q / Keq))

    # calculate the enzyme coefficient:
    numo_E = ((cNao / p.KmNC_Na) ** 3) * (cCai / p.KmNC_Ca)
    denomo_E = (1 + (cNao / p.KmNC_Na) ** 3) * (1 + (cCai / p.KmNC_Ca))

    f_Na = alpha * (numo_E / denomo_E)  # flux as [mol/m2s]

    f_Ca = -(1/ 3) * f_Na  # flux as [mol/m2s]

    return f_Na, f_Ca

def exch_NaKCl(cNai, cNao, cKi, cKo, cCli, cClo, Vm, T, p, block=1.0):
    """
    The Na-K-Cl cotransporter. Important for maintainig the correct
    orientation of the trans-epithelial potential but maintaining
    the correct salt balance and providing stimulus for Na/K-ATPase
    pumpt activity.

    Na_out + K_out + 2 Cl_out --> Na_in + K_in + 2 Cl_in

    cNai:  Sodium inside the cell
    cNao:  Sodium outside the cell

    cKi:  Potassium inside the cell
    cKo:  Potassium outside the cell

    cCli:  Chloride inside the cell
    cClo:  Chloride outside the cell

    Vm:    Transmembrane voltage
    T:     Temperature
    p:     Instance of parameters

    Returns
    -------
    f_Na      Flux of sodium
    f_K      Flux of potassium
    f_K      Flux of chloride

    """

    # calculate the reaction coefficient Q:
    Qnumo = (cNai) * (cKi)* (cCli**2)
    Qdenomo = (cNao) * (cKo)*(cClo**2)

    # ensure no chance of dividing by zero:
    inds_Z = (Qdenomo == 0.0).nonzero()
    Qdenomo[inds_Z] = 1.0e-10

    Q = Qnumo / Qdenomo

    # # calculate the equilibrium constant for the pump reaction:
    Keq = 1

    # calculate the reaction rate coefficient
    alpha = block*p.alpha_NaKCl * (1 - (Q / Keq))

    # calculate the enzyme coefficient:
    numo_E = (cNao / p.KmNaKCl_Na) * (cKo / p.KmNaKCl_K)* ((cClo / p.KmNaKCl_Cl)**2)
    denomo_E = (1 + (cNao / p.KmNaKCl_Na)) * (1 + (cKo / p.KmNaKCl_K)) * (1 + (cClo / p.KmNaKCl_Cl)**2)

    f_Na = alpha * (numo_E / denomo_E)  # flux as [mol/m2s]

    f_K = f_Na

    f_Cl = 2 * f_Na  # flux as [mol/m2s]

    # print(f_Na.mean())

    return f_Na, f_K, f_Cl

def symp_ClK(cKi, cKo, cCli, cClo, Vm, T, p, block=1.0):
    """
    The Na-K-Cl cotransporter. Important for maintainig the correct
    orientation of the trans-epithelial potential but maintaining
    the correct salt balance and providing stimulus for Na/K-ATPase
    pumpt activity.

    K_in + Cl_in --> K_out + Cl_out

    cKi:  Potassium inside the cell
    cKo:  Potassium outside the cell

    cCli:  Chloride inside the cell
    cClo:  Chloride outside the cell

    Vm:    Transmembrane voltage
    T:     Temperature
    p:     Instance of parameters

    Returns
    -------
    f_Na      Flux of sodium
    f_K      Flux of potassium
    f_K      Flux of chloride

    """

    # calculate the reaction coefficient Q:
    Qdenomo = (cKi)* (cCli)
    Qnumo = (cKo)*(cClo)

    # ensure no chance of dividing by zero:
    inds_Z = (Qdenomo == 0.0).nonzero()
    Qdenomo[inds_Z] = 1.0e-10

    Q = Qnumo / Qdenomo

    # # calculate the equilibrium constant for the pump reaction:
    Keq = 1

    # calculate the reaction rate coefficient
    alpha = block*p.alpha_ClK * (1 - (Q / Keq))

    # calculate the enzyme coefficient:
    numo_E = (cKi / p.KmClK_K)* (cCli / p.KmClK_Cl)
    denomo_E = (1 + (cKi / p.KmClK_K)) * (1 + (cCli / p.KmClK_Cl))

    f_K = - alpha * (numo_E / denomo_E)  # flux as [mol/m2s]

    f_Cl = f_K  # flux as [mol/m2s]

    # print(f_K.mean(), cCli.mean())

    return f_K, f_Cl

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
        raise BetseSimInstabilityException(
            "Your simulation has become unstable. Please try a smaller time step,"
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
    dat_grid = gaussian_filter(dat_grid,p.smooth_level)

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

def nernst_planck_vector(c, gc, gv,u,D,z,T,p):
    """
     Calculate the flux component of the Nernst-Planck equation
     along a directional gradient (e.g. gap junction)

     Parameters
     ------------

    c:     concentration
    gc:   concentration gradient
    gv:   voltage gradient
    u:    fluid velocity
    D:     diffusion constant, D
    z:     ion charge
    T:     temperature
    p:     parameters object

    Returns
    --------
    fx, fx        mass flux in x and y directions
    """

    alpha = (D*z*p.q)/(p.kb*T)
    f =  -D*gc - alpha*gv*c + u*c

    return f

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
    """
    This function screens an (concentration) array to
    ensure there are no NaNs and no negative values,
    crashing with an instability message if it finds any.

    """

    # ensure no NaNs:
    inds_nan = (np.isnan(data)).nonzero()

    # ensure that data has no less than zero values:
    inds_neg = (data < 0.0).nonzero()

    if len(inds_neg[0]) > 0:

        data[inds_neg] = 0.0 # add in a small bit to protect from crashing

    if len(inds_nan[0]) > 0:

        raise BetseSimInstabilityException(
            "Your simulation has become unstable. Please try a smaller time step,"
            "reduce gap junction radius, and/or reduce rate coefficients.")

    return data

def bicarbonate_buffer(cCO2, cHCO3):
    """
    This most amazing buffer handles influx of H+,
    HCO3-, H2CO3 (from dissolved carbon dioxide) to
    handle pH in real time.

    Uses the bicarbonate dissacociation reaction:

    H2CO3 ----> HCO3 + H

    Where all dissolved carbon dioxide is assumed
    converted to carbonic acid via carbonic anhydrase enzyme.

    """

    pH = 6.1 + np.log10(cHCO3/cCO2)

    cH = 10**(-pH)*1e3

    return cH, pH

def bicarb_reaction_buffer(cH, cHCO3, cCO2, p):

    Q = (cH*cHCO3)/(cCO2)

    v_ph = p.vm_ph*(cCO2/(1+cCO2))*(1 - (Q/p.Keqm_ph))

    cCO2 = cCO2 - v_ph*p.dt
    cHCO3 = cHCO3 + v_ph*p.dt
    cH = cH + v_ph*p.dt

    pH = -np.log10(cH*1e-3)

    # print("RUNNING BUFFFER")

    return cH, cHCO3, cCO2, pH

def ghk_calculator(sim, cells, p):
    """
    Uses simulation parameters in the Goldman (GHK) equation
    to calculate an alternative Vmem for validation purposes.

    """


    # FIXME the Goldman calculator should be altered to account for network pumps and channels!
    # begin by initializing all summation arrays for the cell network:
    sum_PmAnion_out = np.zeros(len(cells.cell_i))
    sum_PmAnion_in = np.zeros(len(cells.cell_i))
    sum_PmCation_out = np.zeros(len(cells.cell_i))
    sum_PmCation_in = np.zeros(len(cells.cell_i))

    for i, z in enumerate(sim.zs):

        # tag as anion or cation
        ion_type = np.sign(z)

        # average values from membranes or environment to cell centres:
        Dm = np.dot(cells.M_sum_mems, sim.Dm_cells[i]) / cells.num_mems
        conc_cells = sim.cc_cells[i]

        if p.sim_ECM is True:
            # average entities from membranes to the cell centres:
            conc_env = np.dot(cells.M_sum_mems, sim.cc_env[i][cells.map_mem2ecm]) / cells.num_mems

        else:

            conc_env = np.dot(cells.M_sum_mems, sim.cc_env[i]) / cells.num_mems

        if ion_type == -1:

            sum_PmAnion_in = sum_PmAnion_in + Dm * conc_cells * (1 / p.tm)
            sum_PmAnion_out = sum_PmAnion_out + Dm * conc_env * (1 / p.tm)


        if ion_type == 1:

            sum_PmCation_in = sum_PmCation_in + Dm * conc_cells * (1 / p.tm)
            sum_PmCation_out = sum_PmCation_out + Dm * conc_env * (1 / p.tm)


    # NaKrate = (np.dot(cells.M_sum_mems, sim.rate_NaKATP)/cells.num_mems)

    # sum together contributions for Na and K flux across the membrane:
    # NaKflux = NaKrate - (2/3)*NaKrate

    sim.vm_GHK = ((p.R * sim.T) / p.F) * np.log(
        (sum_PmCation_out + sum_PmAnion_in) / (sum_PmCation_in + sum_PmAnion_out))

def molecule_pump(sim, cX_cell_o, cX_env_o, cells, p, Df=1e-9, z=0, pump_into_cell =False, alpha_max=1.0e-8, Km_X=1.0,
                 Km_ATP=1.0, met = None, n=1, ignoreECM = True, rho = 1.0):


    """
    Defines a generic active transport pump that can be used to move
    a general molecule (such as serotonin or glutamate)
    into or out of the cell by active transport.

    Works on the basic premise of enzymatic pumps defined elsewhere:

    pump_out is True:

    cX_cell + cATP  -------> cX_env + cADP + cPi

    pump_out is False:

    cX_env + cATP  <------- cX_cell + cADP + cPi

    Parameters
    -------------
    cX_cell_o           Concentration of X in the cell         [mol/m3]
    cX_env_o            Concentration of X in the environment  [mol/m3]
    cells               Instance of cells
    p                   Instance of parameters
    z                   Charge of X
    pump_out            Is pumping out of cell (pump_out = True) or into cell (pump_out = False)?
    alpha_max           Maximum rate constant of pump reaction [mol/s]
    Km_X                Michaelis-Mentin 1/2 saturation value for X [mol/m3]
    Km_ATP              Michaelis-Mentin 1/2 saturation value for ATP [mol/m3]

    Returns
    ------------
    cX_cell_1     Updated concentration of X in cells
    cX_env_1      Updated concentration of X in environment
    f_X           Flux of X (into the cell +)

    """

    deltaGATP_o = p.deltaGATP  # standard free energy of ATP hydrolysis reaction in J/(mol K)

    if met is None:

        # if metabolism vector not supplied, use singular defaults for concentrations
        cATP = p.cATP
        cADP = p.cADP
        cPi  = p.cPi

    else:

        cATP = met['cATP']  # concentration of ATP in mmol/L
        cADP = met['cADP']  # concentration of ADP in mmol/L
        cPi = met['cPi']  # concentration of Pi in mmol/L

    if p.sim_ECM is True:

        cX_env = cX_env_o[cells.map_mem2ecm]

        cX_cell = cX_cell_o[cells.mem_to_cells]

    else:
        cX_env = cX_env_o[:]
        cX_cell = cX_cell_o[cells.mem_to_cells]

    if pump_into_cell is False:

        # active pumping of molecule from cell and into environment:
        # calculate the reaction coefficient Q:
        Qnumo = cADP * cPi * (cX_env**n)
        Qdenomo = cATP * (cX_cell**n)

        # ensure no chance of dividing by zero:
        inds_Z = (Qdenomo == 0.0).nonzero()
        Qdenomo[inds_Z] = 1.0e-10

        Q = Qnumo / Qdenomo

        # calculate the equilibrium constant for the pump reaction:
        Keq = np.exp(-deltaGATP_o / (p.R * sim.T) + ((n*z * p.F * sim.vm) / (p.R * sim.T)))

        # calculate the reaction rate coefficient
        alpha = alpha_max * (1 - (Q / Keq))

        # calculate the enzyme coefficient:
        numo_E = ((cX_cell / Km_X)**n) * (cATP / Km_ATP)
        denomo_E = (1 + (cX_cell / Km_X)**n) * (1 + (cATP / Km_ATP))

        f_X = -rho*alpha * (numo_E / denomo_E)  # flux as [mol/m2s]   scaled to concentrations Na in and K out

    else:

        # active pumping of molecule from environment and into cell:
        # calculate the reaction coefficient Q:
        Qnumo = cADP * cPi * cX_cell
        Qdenomo = cATP * cX_env

        # ensure no chance of dividing by zero:
        inds_Z = (Qdenomo == 0.0).nonzero()
        Qdenomo[inds_Z] = 1.0e-10

        Q = Qnumo / Qdenomo

        # calculate the equilibrium constant for the pump reaction:
        Keq = np.exp(-deltaGATP_o / (p.R * sim.T) - ((z * p.F * sim.vm) / (p.R * sim.T)))

        # calculate the reaction rate coefficient
        alpha = alpha_max * (1 - (Q / Keq))

        # calculate the enzyme coefficient:
        numo_E = (cX_env / Km_X) * (cATP / Km_ATP)
        denomo_E = (1 + (cX_env / Km_X)) * (1 + (cATP / Km_ATP))

        f_X = rho* alpha * (numo_E / denomo_E)  # flux as [mol/m2s]   scaled to concentrations Na in and K out

    if p.cluster_open is False:
        f_X[cells.bflags_mems] = 0

    cmems = cX_cell_o[cells.mem_to_cells]

    # update cell and environmental concentrations
    cX_cell_1, _, cX_env_1 = update_Co(sim, cX_cell_o, cmems, cX_env_o, f_X, cells, p, ignoreECM = ignoreECM)


    if p.sim_ECM is False:
        cX_env_1_temp = cX_env_1.mean()
        cX_env_1[:] = cX_env_1_temp

    return cX_cell_1, cX_env_1, f_X

def molecule_transporter(sim, cX_cell_o, cX_env_o, cells, p, Df=1e-9, z=0, pump_into_cell=False, alpha_max=1.0e-8,
        Km_X=1.0, Keq=1.0, n = 1.0, ignoreECM = True, rho = 1.0):


    """
    Defines a generic facillitated transporter that can be used to move
    a general molecule (such as glucose).

    ATP is not used for the transporter

    Works on the basic premise of enzymatic pumps defined elsewhere:

    pump_out is True:

    cX_cell  -------> cX_env

    pump_out is False:

    cX_env   <------- cX_cell

    Parameters
    -------------
    cX_cell_o           Concentration of X in the cell         [mol/m3]
    cX_env_o            Concentration of X in the environment  [mol/m3]
    cells               Instance of cells
    p                   Instance of parameters
    z                   Charge of X
    pump_out            Is pumping out of cell (pump_out = True) or into cell (pump_out = False)?
    alpha_max           Maximum rate constant of pump reaction [mol/s]
    Km_X                Michaelis-Mentin 1/2 saturation value for X [mol/m3]

    Returns
    ------------
    cX_cell_1     Updated concentration of X in cells
    cX_env_1      Updated concentration of X in environment
    f_X           Flux of X (into the cell +)

    """

    if p.sim_ECM is True:

        cX_env = cX_env_o[cells.map_mem2ecm]

        cX_cell = cX_cell_o[cells.mem_to_cells]

    else:
        cX_env = cX_env_o
        cX_cell = cX_cell_o[cells.mem_to_cells]

    if pump_into_cell is False:

        # active pumping of molecule from cell and into environment:
        # calculate the reaction coefficient Q:
        Qnumo = cX_env
        Qdenomo = cX_cell

        # ensure no chance of dividing by zero:
        inds_Z = (Qdenomo == 0.0).nonzero()
        Qdenomo[inds_Z] = 1.0e-15

        Q = Qnumo / Qdenomo

        # modify equilibrium constant by membrane voltage if ion is charged:
        Keq = Keq*np.exp((z * p.F * sim.vm) / (p.R * sim.T))

        # calculate the reaction rate coefficient
        alpha = alpha_max * (1 - (Q / Keq))

        # calculate the enzyme coefficient:
        numo_E = (cX_cell / Km_X)
        denomo_E = (1 + (cX_cell / Km_X))

        f_X = -rho*alpha * (numo_E / denomo_E)  # flux as [mol/m2s]   scaled to concentrations Na in and K out


    else:

        # active pumping of molecule from environment and into cell:
        # calculate the reaction coefficient Q:
        Qnumo = cX_cell
        Qdenomo = cX_env

        # ensure no chance of dividing by zero:
        inds_Z = (Qdenomo == 0.0).nonzero()
        Qdenomo[inds_Z] = 1.0e-15

        Q = Qnumo / Qdenomo

        # modify equilibrium constant by membrane voltage if ion is charged:
        Keq = Keq*np.exp(-(z * p.F * sim.vm) / (p.R * sim.T))

        # calculate the reaction rate coefficient
        alpha = alpha_max * (1 - (Q / Keq))

        # calculate the enzyme coefficient:
        numo_E = (cX_env / Km_X)
        denomo_E = (1 + (cX_env / Km_X))

        f_X = rho*alpha * (numo_E / denomo_E)  # flux as [mol/m2s]   scaled to concentrations Na in and K out

    if p.cluster_open is False:
        f_X[cells.bflags_mems] = 0

    cmems = cX_cell_o[cells.mem_to_cells]

    # update cell and environmental concentrations
    cX_cell_1, _, cX_env_1 = update_Co(sim, cX_cell_o, cmems, cX_env_o, f_X, cells, p, ignoreECM= ignoreECM)

    # next electrodiffuse concentrations around the cell interior:
    # cX_cell_1 = update_intra(sim, cells, cX_cell_1, Df, z, p)

    # ensure that there are no negative values
    # cX_cell_1 = no_negs(cX_cell_1)
    # cX_env_1 = no_negs(cX_env_1)


    if p.sim_ECM is False:
        cX_env_1_temp = cX_env_1.mean()
        cX_env_1[:] = cX_env_1_temp

    return cX_cell_1, cX_env_1, f_X

def molecule_mover(sim, cX_env_o, cX_cells, cells, p, z=0, Dm=1.0e-18, Do=1.0e-9, Dgj=1.0e-12, c_bound=1.0e-6,
                   ignoreECM = False, smoothECM = False, ignoreTJ = False, ignoreGJ = False, rho = 1, cmems = None):

    """
    Transports a generic molecule across the membrane,
    through gap junctions, and if p.sim_ECM is true,
    through extracellular spaces and the environment.

    Parameters
    -----------
    cX_cell_o           Concentration of molecule in the cytosol [mol/m3]
    cX_env_o            Concentration of molecule in the environment [mol/m3]
    cells               Instance of Cells
    p                   Instance of Parameters
    z                   Charge state of molecule
    Dm                  Membrane diffusion constant [m2/s]
    Do                  Free diffusion constant [m2/s]
    c_bound             Concentration of molecule at global bounds (required for sim_ECM True only)

    Returns
    -----------
    cX_cell_1         Updated concentration of molecule in the cell
    cX_env_1          Updated concentration of molecule in the environment

    """


    if p.sim_ECM is True:

        cX_env = cX_env_o[cells.map_mem2ecm]

    else:
        cX_env = cX_env_o[:]

    if cmems is None:

        cX_mems = cX_cells[cells.mem_to_cells]

    elif len(cmems) == sim.mdl:
        cX_mems = cmems

    Dm_vect = np.ones(len(cX_mems))*Dm

    # Transmembrane: electrodiffuse molecule X between cell and extracellular space------------------------------

    IdM = np.ones(sim.mdl)

    f_X_ED = electroflux(cX_env, cX_mems, Dm_vect, p.tm*IdM, z*IdM, sim.vm, sim.T, p, rho = rho)

    if p.cluster_open is False:
        f_X_ED[cells.bflags_mems] = 0

    # update concentrations due to electrodiffusion:

    cX_cells, cX_mems, cX_env_o = update_Co(sim, cX_cells, cX_mems, cX_env_o, f_X_ED, cells, p, ignoreECM = ignoreECM)

    # ------------------------------------------------------------
    if ignoreGJ is False:
        # Update dye concentration in the gj connected cell network:

        grad_cgj = (cX_mems[cells.nn_i] - cX_mems[cells.mem_i]) / cells.gj_len

        gcx = grad_cgj*cells.mem_vects_flat[:, 2]
        gcy = grad_cgj*cells.mem_vects_flat[:, 3]

        # midpoint concentration:
        cX_mids = (cX_mems[cells.nn_i] + cX_mems[cells.mem_i]) / 2

        # # electroosmotic fluid velocity:
        if p.fluid_flow is True:
            ux = sim.u_cells_x[cells.mem_to_cells]
            uy = sim.u_cells_y[cells.mem_to_cells]

        else:
            ux = 0
            uy = 0


        fgj_x, fgj_y = nernst_planck_flux(cX_mids, gcx, gcy, -sim.E_gj_x,
                                          -sim.E_gj_y, ux, uy,
                                          sim.gjopen*Dgj*sim.gj_block, z, sim.T, p)

        fgj_X = fgj_x*cells.mem_vects_flat[:,2] + fgj_y*cells.mem_vects_flat[:,3]


        # enforce zero flux at outer boundary:
        fgj_X[cells.bflags_mems] = 0.0

        # divergence calculation for individual cells (finite volume expression)
        delta_cco = np.dot(cells.M_sum_mems, -fgj_X*cells.mem_sa) / cells.cell_vol


        # Calculate the final concentration change (the acceleration effectively speeds up time):
        cX_cells = cX_cells + p.dt*delta_cco
        cX_mems = cX_mems - fgj_X*p.dt*(cells.mem_sa/cells.mem_vol)


    else:
        fgj_X = np.zeros(sim.mdl)


    #------------------------------------------------------------------------------------------------------------

    # Transport through environment, if p.sim_ECM is True-----------------------------------------------------

    if p.sim_ECM is True: #-----------------------------------------------------------------------------------------

        cenv = cX_env_o
        cenv = cenv.reshape(cells.X.shape)

        if smoothECM is True and p.smooth_level > 0.0:
            cenv = gaussian_filter(cenv, p.smooth_level, mode='constant', cval = c_bound)

        cenv[:, 0] = c_bound
        cenv[:, -1] = c_bound
        cenv[0, :] = c_bound
        cenv[-1, :] = c_bound

        denv = Do * np.ones(len(cells.xypts))
        denv = denv.reshape(cells.X.shape)*sim.D_env_weight

        gcx, gcy = fd.gradient(cenv, cells.delta)

        if p.fluid_flow is True:

            ux = sim.u_env_x.reshape(cells.X.shape)
            uy = sim.u_env_y.reshape(cells.X.shape)

        else:

            ux = 0.0
            uy = 0.0


        fx, fy = nernst_planck_flux(cenv, gcx, gcy, -sim.E_env_x, -sim.E_env_y, ux, uy,
                                        denv, z, sim.T, p)

        div_fa = fd.divergence(-fx, -fy, cells.delta, cells.delta)

        fenvx = fx
        fenvy = fy

        cenv = cenv + div_fa * p.dt

        cX_env_o = cenv.ravel()

    else:
        cX_env_temp = cX_env_o.mean()
        cX_env_o = np.zeros(sim.mdl)
        cX_env_o[:] = cX_env_temp
        fenvx = 0
        fenvy = 0

    return cX_env_o, cX_cells, cX_mems, f_X_ED, fgj_X, fenvx, fenvy

def update_Co(sim, cX_cell, cX_mem, cX_env, flux, cells, p, ignoreECM = True):
    """

    General updater for a concentration defined on
    cells and in environment.

    Parameters
    --------------
    cX_cell     Concentration inside cell defined on cell centres [mol/m3]
    cX_env      Concentration in extracellular spaces [mol/m3]
    flux        Flux of X, where into cell is + and component is normal to membranes  [mol/m2 s]
    cells       Instance of Cells
    p           Instance of Parameters

    Returns
    --------------
    cX_cell     Updated concentration in cell [mol/m3]
    cX_env      Updated concentration in environment [mol/m3]

    """

    # take the divergence of the flux for each enclosed cell:
    delta_cells = np.dot(cells.M_sum_mems, flux * cells.mem_sa) / cells.cell_vol

    # update cell concentration of substance:
    cX_cell = cX_cell + delta_cells * p.dt

    cX_mem = cX_mem + flux*(cells.mem_sa/cells.mem_vol)*p.dt

    if p.sim_ECM is True:

        flux_env = np.zeros(sim.edl)
        flux_env[cells.map_mem2ecm] = -flux

        # save values at the cluster boundary:
        bound_vals = flux_env[cells.ecm_bound_k]

        # set the values of the global environment to zero:
        flux_env[cells.inds_env] = 0

        # finally, ensure that the boundary values are restored:
        flux_env[cells.ecm_bound_k] = bound_vals

        delta_env = (flux_env * cells.memSa_per_envSquare) / cells.ecm_vol

        # update the environmental concentrations:
        cX_env = cX_env + delta_env * p.dt


    else:

        delta_env = -flux * (cells.mem_sa / p.vol_env)

        cX_env_o = cX_env + delta_env * p.dt

        # assume auto-mixing of environmental concentrations:
        cX_env = cX_env_o.mean()

    return cX_cell, cX_mem, cX_env

def div_env(flux, sim, cells, p):

    flux_env = np.zeros(sim.edl)
    flux_env[cells.map_mem2ecm] = -flux

    # save values at the cluster boundary:
    bound_vals = flux_env[cells.ecm_bound_k]

    # set the values of the global environment to zero:
    flux_env[cells.inds_env] = 0

    # finally, ensure that the boundary values are restored:
    flux_env[cells.ecm_bound_k] = bound_vals

    # Now that we have a nice, neat interpolation of flux from cell membranes, multiply by the
    # true membrane surface area in the square, and divide by the true ecm volume of the env grid square,
    # to get the mol/s change in concentration (divergence):

    # if sim.ignore_ecm is False:

        # delta_env = (flux_env * cells.memSa_per_envSquare) / cells.true_ecm_vol

    # else:

    delta_env = (flux_env * cells.memSa_per_envSquare) / cells.ecm_vol

    # if p.smooth_level > 0.0:
    #     delta_env = gaussian_filter(delta_env.reshape(cells.X.shape), p.smooth_level).ravel()

    return delta_env

def HH_Decomp(JJx, JJy, cells, bounds = None):

    if bounds is not None:

        Lb = bounds['L']
        Rb = bounds['R']
        Tb = bounds['T']
        Bb = bounds['B']

    else:
        Lb = 0.0
        Rb = 0.0
        Tb = 0.0
        Bb = 0.0

    # ----divergence-free component--------------------------------------
    # FIXME specify boundary condishs better!

    Jxr = -JJy.reshape(cells.X.shape)
    Jyr = JJx.reshape(cells.X.shape)

    divJr = fd.divergence(Jxr, Jyr, cells.delta, cells.delta)

    divJr[:, 0] = 0.0
    divJr[:, -1] = 0.0
    divJr[0, :] = 0.0
    divJr[-1, :] = 0.0

    AA = np.dot(cells.lapENVinv, -divJr.ravel())

    gAx, gAy = fd.gradient(AA.reshape(cells.X.shape), cells.delta)

    Fx = -gAy
    Fy = gAx
    #
    # F = np.sqrt(Fx ** 2 + Fy ** 2)

    # ----curl free component------------------------------------------

    divJd = fd.divergence(JJx.reshape(cells.X.shape), JJy.reshape(cells.X.shape), cells.delta, cells.delta)

    # set boundary conditions for normal component
    # enforce applied voltage condition at the boundary:
    divJd[:, 0] = Lb * (1 / cells.delta ** 2)
    divJd[:, -1] = Rb * (1 / cells.delta ** 2)
    divJd[0, :] = Bb * (1 / cells.delta ** 2)
    divJd[-1, :] = Tb * (1 / cells.delta ** 2)


    BB = np.dot(cells.lapENVinv, divJd.ravel())

    Gx, Gy = fd.gradient(BB.reshape(cells.X.shape), cells.delta)

    # G = np.sqrt(Gx ** 2 + Gy ** 2)

    return AA, Fx, Fy, BB, Gx, Gy

def div_free(Fxo, Fyo, cells):

    """

    Uses a "hidden potential" method to
    calculate a divergence-free field

    :param Fxo:  x-component of input field
    :param Fyo:  y-component of input field
    :param cells:  cells instance
    :return: Fx, Fy divergence-free field components (finite divergence at bounds)

    """
    # calculate divergence of vector field:
    divF = fd.divergence(Fxo.reshape(cells.X.shape), Fyo.reshape(cells.X.shape), cells.delta, cells.delta)

    # value of the correcting potenial:
    Phi = np.dot(cells.lapENVinv, divF.ravel())

    gPhix, gPhiy = fd.gradient(Phi.reshape(cells.X.shape), cells.delta)

    Fx = Fxo.reshape(cells.X.shape) - gPhix
    Fy = Fyo.reshape(cells.X.shape) - gPhiy

    return Fx, Fy, Phi

def single_cell_div_free(cfluxo, cells):
    """
    Averages a membrane-specific vector field to the cell centre
    in order to provide a divergence-free field wrt to a single cell.

    Parameters
    ------------
    cfluxo:            Component of non-divergence free field wrt membranes
    cells:             Cells object

    Returns
    -------
    cflux               Divergence-free field normal component to membrane
    """

    # as no flux is leaving the cell, cflux must be a divergence-free field:
    cfxo = np.dot(cells.M_sum_mems, cfluxo * cells.mem_vects_flat[:, 2] * cells.mem_sa) / cells.cell_sa
    cfyo = np.dot(cells.M_sum_mems, cfluxo * cells.mem_vects_flat[:, 3] * cells.mem_sa) / cells.cell_sa

    cflux = cfxo[cells.mem_to_cells]*cells.mem_vects_flat[:, 2] + cfyo[cells.mem_to_cells]*cells.mem_vects_flat[:,3]

    return cflux


#----------------------------------------------------------------------------------------------------------------
# WASTELANDS
#---------------------------------------------------------------------------------------------------------------
# def update_intra(sim, cells, cX_mems, cX_cells, D_x, zx, p):
#     """
#     Perform electrodiffusion on intracellular vertices
#     to update concentration around the cell interior
#     in response to interior voltage and concentration
#     gradients, in addition to any electroosmotic flows.
#
#     """
#
#     # x and y components of membrane tangent unit vectors
#     nx = cells.mem_vects_flat[:, 2]
#     ny = cells.mem_vects_flat[:, 3]
#
#     if p.fluid_flow is True and p.run_sim is True:
#         # get intracellular fluid flow vector
#         ux_mem = sim.u_cells_x[cells.mem_to_cells]
#         uy_mem = sim.u_cells_y[cells.mem_to_cells]
#
#         # component of fluid flow velocity normal to the membrane:
#         u = ux_mem*nx + uy_mem*ny
#
#     else:
#         u = 0
#
#     # get the gradient of rho concentration for each cell centre wrt to each membrane midpoint:
#     grad_c = (cX_mems - cX_cells[cells.mem_to_cells])/cells.chords
#     # grad_v = (sim.v_cell - sim.v_cell_ave[cells.mem_to_cells])/cells.chords
#
#     # field-modulate the grad_v to account for screening (assumes motion primarily near the double layer):
#     # grad_v = (1.0e-9/p.d_cell)*grad_v
#     # grad_v = p.field_modulation*grad_v
#
#     # obtain an average concentration at the pie-slice midpoints:
#     c_at_mids = (cX_mems + cX_cells[cells.mem_to_cells])/2
#
#     # flux_intra = - D_x * p.cell_delay_const * grad_c + u * c_at_mids
#     flux_intra = nernst_planck_vector(c_at_mids, grad_c, grad_v, u, D_x, zx, sim.T, p)
#
#     # net flux across inner centroid surfaces of each cell:
#     net_flux = flux_intra * (cells.mem_sa / 2)
#
#     # sum of the flux entering each centroid region:
#     flux_sum = np.dot(cells.M_sum_mems, net_flux)
#
#     # divergence of the flux wrt the centroid region:
#     divF_centroid = flux_sum/cells.centroid_vol
#
#     # divergence of the flux wrt each membrane region
#     divF_mems = net_flux/cells.mem_vol
#
#     # update concentrations with the appropriate divergence:
#     cX_mems = cX_mems + divF_mems * p.dt
#     cX_cells = cX_cells - divF_centroid * p.dt
#
#     # # use finite volume method to integrate each region:
#     # # values at centroid mids:
#     # c_at_mids = (cX_memso + cX_cellso[cells.mem_to_cells]) / 2
#     #
#     # # finite volume integral of membrane pie-box values:
#     # cX_mems = np.dot(cells.M_int_mems, cX_memso) + (1 / 2) * c_at_mids
#     # cX_cells = (1 / 2) * cX_cellso + np.dot(cells.M_sum_mems, c_at_mids) / (2 * cells.num_mems)
#
#     return cX_mems, cX_cells, net_flux

# MACS METHOD for environmental handling:

# if p.closed_bound is True:
#     btag = 'closed'
#
# else:
#     btag = 'open'
#
# # make v_env and cc_env into 2d matrices
# cenv = cX_env_o
# denv = Do * np.ones(len(cells.xypts))
#
# v_env = sim.v_env.reshape(cells.X.shape)
#
# v_env[:, 0] = sim.bound_V['L']
# v_env[:, -1] = sim.bound_V['R']
# v_env[0, :] = sim.bound_V['B']
# v_env[-1, :] = sim.bound_V['T']
#
# cenv = cenv.reshape(cells.X.shape)
#
# # prepare concentrations and diffusion constants for MACs grid format
# # by resampling the values at the u v coordinates of the flux:
# cenv_x = np.zeros(cells.grid_obj.u_shape)
# cenv_y = np.zeros(cells.grid_obj.v_shape)
#
# # create the proper shape for the concentrations and state appropriate boundary conditions::
# cenv_x[:, 1:] = cenv[:]
# cenv_x[:, 0] = cenv_x[:, 1]
# cenv_y[1:, :] = cenv[:]
# cenv_y[0, :] = cenv_y[1, :]
#
# if p.closed_bound is True:  # insulation boundary conditions
#
#     cenv_x[:, 0] = cenv_x[:, 1]
#     cenv_x[:, -1] = cenv_x[:, -2]
#     cenv_x[0, :] = cenv_x[1, :]
#     cenv_x[-1, :] = cenv_x[-2, :]
#
#     cenv_y[0, :] = cenv_y[1, :]
#     cenv_y[-1, :] = cenv_y[-2, :]
#     cenv_y[:, 0] = cenv_y[:, 1]
#     cenv_y[:, -1] = cenv_y[:, -2]
#
# else:  # open and electrically grounded boundary conditions
#     cenv_x[:, 0] = c_bound
#     cenv_x[:, -1] = c_bound
#     cenv_x[0, :] = c_bound
#     cenv_x[-1, :] = c_bound
#
#     cenv_y[0, :] = c_bound
#     cenv_y[-1, :] = c_bound
#     cenv_y[:, 0] = c_bound
#     cenv_y[:, -1] = c_bound
#
# denv = denv.reshape(cells.X.shape)
#
# denv_x = interp.griddata((cells.xypts[:, 0], cells.xypts[:, 1]), denv.ravel(),
#     (cells.grid_obj.u_X, cells.grid_obj.u_Y), method='nearest', fill_value=Do)
#
# denv_y = interp.griddata((cells.xypts[:, 0], cells.xypts[:, 1]), denv.ravel(),
#     (cells.grid_obj.v_X, cells.grid_obj.v_Y), method='nearest', fill_value=Do)
#
# if ignoreTJ is False:
#
#     denv_x = denv_x * sim.D_env_weight_u
#     denv_y = denv_y * sim.D_env_weight_v
#
# # calculate gradients in the environment
# grad_V_env_x, grad_V_env_y = cells.grid_obj.grid_gradient(v_env, bounds='closed')
#
# grad_cc_env_x, grad_cc_env_y = cells.grid_obj.grid_gradient(cenv, bounds=btag)
#
# # calculate fluxes for electrodiffusive transport in environment:
#
# if p.fluid_flow is True:
#
#     uenvx = np.zeros(cells.grid_obj.u_shape)
#     uenvy = np.zeros(cells.grid_obj.v_shape)
#
#     uenvx[:, 1:] = sim.u_env_x
#     uenvy[1:, :] = sim.u_env_y
#
#     if p.closed_bound is False:
#
#         uenvx[:, 0] = uenvx[:, 1]
#         uenvx[:, -1] = uenvx[:, -2]
#         uenvx[0, :] = uenvx[1, :]
#         uenvx[-1, :] = uenvx[-2, :]
#
#         uenvy[:, 0] = uenvy[:, 1]
#         uenvy[:, -1] = uenvy[:, -2]
#         uenvy[0, :] = uenvy[1, :]
#         uenvy[-1, :] = uenvy[-2, :]
#
#     else:
#
#         uenvx[:, 0] = 0
#         uenvx[:, -1] = 0
#         uenvx[0, :] = 0
#         uenvx[-1, :] = 0
#
#         uenvy[:, 0] = 0
#         uenvy[:, -1] = 0
#         uenvy[0, :] = 0
#         uenvy[-1, :] = 0
#
# else:
#     uenvx = 0
#     uenvy = 0
#
# field_mod = 1.0
#
# f_env_x_X, f_env_y_X = np_flux_special(cenv_x, cenv_y, grad_cc_env_x, grad_cc_env_y,
#     field_mod*grad_V_env_x, field_mod*grad_V_env_y, uenvx, uenvy, denv_x, denv_y, z, sim.T, p)
#
#
# # calculate the divergence of the total (negative) flux to obtain the total change per unit time:
# d_fenvx = -(f_env_x_X[:, 1:] - f_env_x_X[:, 0:-1]) / cells.delta
# d_fenvy = -(f_env_y_X[1:, :] - f_env_y_X[0:-1, :]) / cells.delta
#
# delta_c = d_fenvx + d_fenvy
#
# cenv = cenv + delta_c * p.dt
#
# if p.closed_bound is True:
#     # Neumann boundary condition (flux at boundary)
#     # zero flux boundaries for concentration:
#     cenv[:, -1] = cenv[:, -2]
#     cenv[:, 0] = cenv[:, 1]
#     cenv[0, :] = cenv[1, :]
#     cenv[-1, :] = cenv[-2, :]
#
# elif p.closed_bound is False:
#     # if the boundary is open, set the concentration at the boundary
#     # open boundary
#     cenv[:, -1] = c_bound
#     cenv[:, 0] = c_bound
#     cenv[0, :] = c_bound
#     cenv[-1, :] = c_bound
#
# if smoothECM is True:
#
#     cenv = gaussian_filter(cenv, p.smooth_level)
#
# # reshape the matrices into vectors:
# # self.v_env = self.v_env.ravel()
# cX_env_o = cenv.ravel()
#
# # average flux at the midpoint of the MACs grid:
# fenvx = (f_env_x_X[:, 1:] + f_env_x_X[:, 0:-1]) / 2
# fenvy = (f_env_y_X[1:, :] + f_env_y_X[0:-1, :]) / 2
