#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# FIXME all pumps should now take and return cATP, cADP and cP as parameters (or a "metabo" object)

import numpy as np
import numpy.ma as ma
from scipy import interpolate as interp
from scipy.ndimage.filters import gaussian_filter
from betse.science import finitediff as fd

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

    deno = 1 - exp_alpha + 1e-15  # calculate the denominator with small term in case it's equal to zero

    # calculate the flux for those elements:
    flux = -((Dc*alpha)/d)*((cB - cA*exp_alpha)/deno)

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
    numo_E = ((cNai/p.KmNK_Na)) * ((cKo/p.KmNK_K)) * (cATP/p.KmNK_ATP)
    denomo_E = (1 + (cNai/p.KmNK_Na))*(1+(cKo/p.KmNK_K))*(1+(cATP/p.KmNK_ATP))

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

    alpha = p.alpha_Ca*(delG_pump - p.halfmax_Ca)

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

    # print(block, f_H.mean())

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

        # tag as anion or cation
        ion_type = np.sign(z)

        # average values from membranes or environment to cell centres:
        Dm = np.dot(cells.M_sum_mems, sim.Dm_cells[i]) / cells.num_mems
        conc_cells = np.dot(cells.M_sum_mems, sim.cc_mems[i]) / cells.num_mems

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


    sim.vm_GHK = ((p.R * sim.T) / p.F) * np.log(
        (sum_PmCation_out + sum_PmAnion_in) / (sum_PmCation_in + sum_PmAnion_out))

def molecule_pump(sim, cX_cell_o, cX_env_o, cells, p, Df=1e-9, z=0, pump_into_cell =False, alpha_max=1.0e-8, Km_X=1.0,
                 Km_ATP=1.0):


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

    # At the moment the concentrations are fixed as the metabolism module isn't ready yet.
    cATP = p.cATP  # concentration of ATP in mmol/L
    cADP = p.cADP  # concentration of ADP in mmol/L
    cPi = p.cPi  # concentration of Pi in mmol/L

    if p.sim_ECM is True:

        cX_env = cX_env_o[cells.map_mem2ecm]

        cX_cell = cX_cell_o[:]

    else:
        cX_env = cX_env_o[:]
        cX_cell = cX_cell_o[:]

    if pump_into_cell is False:

        # active pumping of molecule from cell and into environment:
        # calculate the reaction coefficient Q:
        Qnumo = cADP * cPi * (cX_env)
        Qdenomo = cATP * (cX_cell)

        # ensure no chance of dividing by zero:
        inds_Z = (Qdenomo == 0.0).nonzero()
        Qdenomo[inds_Z] = 1.0e-10

        Q = Qnumo / Qdenomo

        # calculate the equilibrium constant for the pump reaction:
        Keq = np.exp(-deltaGATP_o / (p.R * sim.T) + ((z * p.F * sim.vm) / (p.R * sim.T)))

        # calculate the reaction rate coefficient
        alpha = alpha_max * (1 - (Q / Keq))

        # calculate the enzyme coefficient:
        numo_E = (cX_cell / Km_X) * (cATP / Km_ATP)
        denomo_E = (1 + (cX_cell / Km_X)) * (1 + (cATP / Km_ATP))

        f_X = -alpha * (numo_E / denomo_E)  # flux as [mol/m2s]   scaled to concentrations Na in and K out


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

        f_X = alpha * (numo_E / denomo_E)  # flux as [mol/m2s]   scaled to concentrations Na in and K out

    # update cell and environmental concentrations
    cX_cell_1, cX_env_1 = update_Co(sim, cX_cell_o, cX_env_o, f_X, cells, p)

    # next electrodiffuse concentrations around the cell interior:
    # cX_cell_1 = update_intra(sim, cells, cX_cell_1, Df, z, p)

    # ensure that there are no negative values
    cX_cell_1 = no_negs(cX_cell_1)
    cX_env_1 = no_negs(cX_env_1)

    if p.sim_ECM is True:
        # cX_env_1_temp = gaussian_filter(cX_env_1.reshape(cells.X.shape), p.smooth_level)
        # cX_env_1 = cX_env_1_temp.ravel()
        pass

    else:
        cX_env_1_temp = cX_env_1.mean()
        cX_env_1[:] = cX_env_1_temp

    return cX_cell_1, cX_env_1, f_X

def molecule_mover(sim, cX_mems_o, cX_env_o, cells, p, z=0, Dm=1.0e-18, Do=1.0e-9, c_bound=1.0e-6, ignoreECM = False):
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

        cX_mems = cX_mems_o[:]


    else:
        cX_env = cX_env_o[:]
        cX_mems = cX_mems_o[:]

    Dm_vect = np.ones(len(cX_mems))*Dm

    # Transmembrane: electrodiffuse molecule X between cell and extracellular space------------------------------

    f_X_ED = electroflux(cX_env, cX_mems, Dm_vect, p.tm, z, sim.vm, sim.T, p)

    # update concentrations due to electrodiffusion:

    cX_mems, cX_env_o = update_Co(sim, cX_mems, cX_env_o, f_X_ED, cells, p, ignoreECM = ignoreECM)

    # ------------------------------------------------------------

    # Update dye concentration in the gj connected cell network:

    # Intracellular voltage gradient:
    grad_vgj = sim.vgj / cells.gj_len

    grad_cgj = (cX_mems[cells.nn_i] - cX_mems[cells.mem_i]) / cells.gj_len

    # midpoint concentration:
    cX_mids = (cX_mems[cells.nn_i] + cX_mems[cells.mem_i]) / 2

    # electroosmotic fluid velocity:
    if p.fluid_flow is True:
        ux = sim.u_gj_x
        uy = sim.u_gj_y

        # get component of fluid tangent to gap junctions
        ugj = ux * cells.mem_vects_flat[:, 2] + uy * cells.mem_vects_flat[:, 3]

    else:
        ugj = 0

    fgj_X = nernst_planck_vector(cX_mids, grad_cgj, grad_vgj, ugj,
        p.gj_surface*sim.gjopen*Do, z, sim.T, p)

    # divergence calculation for individual cells (finite volume expression)
    delta_cc = (-fgj_X * cells.mem_sa) / cells.mem_vol

    cX_mems = cX_mems + p.dt * delta_cc

    #------------------------------------------------------------------------------------------------------------

    # electrodiffuse intracellular concentrations
    # cX_cell_1 = update_intra(sim, cells, cX_cell_1, Do, z, p)

    # Transport dye through environment, if p.sim_ECM is True-----------------------------------------------------

    if p.sim_ECM is True:

        if p.closed_bound is True:
            btag = 'closed'

        else:
            btag = 'open'

        # make v_env and cc_env into 2d matrices
        cenv = cX_env_o
        denv = Do * np.ones(len(cells.xypts))

        v_env = sim.v_env.reshape(cells.X.shape)

        v_env[:, 0] = sim.bound_V['L']
        v_env[:, -1] = sim.bound_V['R']
        v_env[0, :] = sim.bound_V['B']
        v_env[-1, :] = sim.bound_V['T']

        cenv = cenv.reshape(cells.X.shape)

        # prepare concentrations and diffusion constants for MACs grid format
        # by resampling the values at the u v coordinates of the flux:
        cenv_x = np.zeros(cells.grid_obj.u_shape)
        cenv_y = np.zeros(cells.grid_obj.v_shape)

        # create the proper shape for the concentrations and state appropriate boundary conditions::
        cenv_x[:, 1:] = cenv[:]
        cenv_x[:, 0] = cenv_x[:, 1]
        cenv_y[1:, :] = cenv[:]
        cenv_y[0, :] = cenv_y[1, :]

        if p.closed_bound is True:  # insulation boundary conditions

            cenv_x[:, 0] = cenv_x[:, 1]
            cenv_x[:, -1] = cenv_x[:, -2]
            cenv_x[0, :] = cenv_x[1, :]
            cenv_x[-1, :] = cenv_x[-2, :]

            cenv_y[0, :] = cenv_y[1, :]
            cenv_y[-1, :] = cenv_y[-2, :]
            cenv_y[:, 0] = cenv_y[:, 1]
            cenv_y[:, -1] = cenv_y[:, -2]

        else:  # open and electrically grounded boundary conditions
            cenv_x[:, 0] = c_bound
            cenv_x[:, -1] = c_bound
            cenv_x[0, :] = c_bound
            cenv_x[-1, :] = c_bound

            cenv_y[0, :] = c_bound
            cenv_y[-1, :] = c_bound
            cenv_y[:, 0] = c_bound
            cenv_y[:, -1] = c_bound

        denv = denv.reshape(cells.X.shape)

        denv_x = interp.griddata((cells.xypts[:, 0], cells.xypts[:, 1]), denv.ravel(),
            (cells.grid_obj.u_X, cells.grid_obj.u_Y), method='nearest', fill_value=p.Do_Dye)

        denv_y = interp.griddata((cells.xypts[:, 0], cells.xypts[:, 1]), denv.ravel(),
            (cells.grid_obj.v_X, cells.grid_obj.v_Y), method='nearest', fill_value=p.Do_Dye)

        denv_x = denv_x * sim.D_env_weight_u
        denv_y = denv_y * sim.D_env_weight_v

        # calculate gradients in the environment
        grad_V_env_x, grad_V_env_y = cells.grid_obj.grid_gradient(v_env, bounds='closed')

        grad_cc_env_x, grad_cc_env_y = cells.grid_obj.grid_gradient(cenv, bounds=btag)

        # calculate fluxes for electrodiffusive transport in environment:

        if p.fluid_flow is True:

            uenvx = np.zeros(cells.grid_obj.u_shape)
            uenvy = np.zeros(cells.grid_obj.v_shape)

            uenvx[:, 1:] = sim.u_env_x
            uenvy[1:, :] = sim.u_env_y

            if p.closed_bound is False:

                uenvx[:, 0] = uenvx[:, 1]
                uenvx[:, -1] = uenvx[:, -2]
                uenvx[0, :] = uenvx[1, :]
                uenvx[-1, :] = uenvx[-2, :]

                uenvy[:, 0] = uenvy[:, 1]
                uenvy[:, -1] = uenvy[:, -2]
                uenvy[0, :] = uenvy[1, :]
                uenvy[-1, :] = uenvy[-2, :]

            else:

                uenvx[:, 0] = 0
                uenvx[:, -1] = 0
                uenvx[0, :] = 0
                uenvx[-1, :] = 0

                uenvy[:, 0] = 0
                uenvy[:, -1] = 0
                uenvy[0, :] = 0
                uenvy[-1, :] = 0

        else:
            uenvx = 0
            uenvy = 0

        field_mod = (1e-9 / p.cell_space)

        f_env_x_X, f_env_y_X = np_flux_special(cenv_x, cenv_y, grad_cc_env_x, grad_cc_env_y,
            field_mod*grad_V_env_x, field_mod*grad_V_env_y, uenvx, uenvy, denv_x, denv_y, z, sim.T, p)

        # if p.smooth_level > 0.0:
        #     f_env_x_X = gaussian_filter(f_env_x_X,p.smooth_level)  # smooth out the flux terms
        #     f_env_y_X = gaussian_filter(f_env_y_X, p.smooth_level)  # smooth out the flux terms

        # calculate the divergence of the total (negative) flux to obtain the total change per unit time:
        d_fenvx = -(f_env_x_X[:, 1:] - f_env_x_X[:, 0:-1]) / cells.delta
        d_fenvy = -(f_env_y_X[1:, :] - f_env_y_X[0:-1, :]) / cells.delta

        delta_c = d_fenvx + d_fenvy

        cenv = cenv + delta_c * p.dt

        if p.closed_bound is True:
            # Neumann boundary condition (flux at boundary)
            # zero flux boundaries for concentration:
            cenv[:, -1] = cenv[:, -2]
            cenv[:, 0] = cenv[:, 1]
            cenv[0, :] = cenv[1, :]
            cenv[-1, :] = cenv[-2, :]

        elif p.closed_bound is False:
            # if the boundary is open, set the concentration at the boundary
            # open boundary
            cenv[:, -1] = c_bound
            cenv[:, 0] = c_bound
            cenv[0, :] = c_bound
            cenv[-1, :] = c_bound

            cenv[:, -2] = c_bound
            cenv[:, 1] = c_bound
            cenv[1, :] = c_bound
            cenv[-2, :] = c_bound

        # reshape the matrices into vectors:
        # self.v_env = self.v_env.ravel()
        cX_env_o = cenv.ravel()

        # average flux at the midpoint of the MACs grid:
        fenvx = (f_env_x_X[:, 1:] + f_env_x_X[:, 0:-1]) / 2
        fenvy = (f_env_y_X[1:, :] + f_env_y_X[0:-1, :]) / 2


    else:
        cX_env_temp = cX_env_o.mean()
        cX_env_o = np.zeros(len(cX_env_o))
        cX_env_o[:] = cX_env_temp
        fenvx = 0
        fenvy = 0

    return cX_mems, cX_env_o, f_X_ED, fgj_X, fenvx, fenvy

def update_Co(sim, cX_cell, cX_env, flux, cells, p, ignoreECM = False):
    """

    General updater for a concentration defined on
    cells and in environment.

    Parameters
    --------------
    cX_cell     Concentration inside cell  [mol/m3]
    cX_env      Concentration in extracellular spaces [mol/m3]
    flux        Flux of X, where into cell is +   [mol/m2 s]
    cells       Instance of Cells
    p           Instance of Parameters

    Returns
    --------------
    cX_cell     Updated concentration in cell [mol/m3]
    cX_env      Updated concentration in environment [mol/m3]

    """

    if p.sim_ECM is True:

        delta_cells = flux * (cells.mem_sa / cells.mem_vol)

        flux_env = np.zeros(sim.edl)
        flux_env[cells.map_mem2ecm] = -flux

        # save values at the cluster boundary:
        bound_vals = flux_env[cells.ecm_bound_k]

        # set the values of the global environment to zero:
        flux_env[cells.inds_env] = 0

        # finally, ensure that the boundary values are restored:
        flux_env[cells.ecm_bound_k] = bound_vals

        # Now that we have a nice, neat interpolation of flux from cell membranes, multiply by the,
        # true membrane surface area in the square, and divide by the true ecm volume of the env grid square,
        # to get the mol/s change in concentration (divergence):

        if ignoreECM is False:

            delta_env = (flux_env * cells.memSa_per_envSquare) / cells.true_ecm_vol

        else:

            delta_env = (flux_env * cells.memSa_per_envSquare) / cells.ecm_vol

        # update the concentrations:
        cX_cell = cX_cell + delta_cells * p.dt
        cX_env = cX_env + delta_env * p.dt

    else:

        delta_cells = flux * (cells.mem_sa / cells.mem_vol)
        delta_env = -flux * (cells.mem_sa / p.vol_env)

        cX_cell = cX_cell + delta_cells * p.dt

        cX_env_o = cX_env + delta_env * p.dt

        # assume auto-mixing of environmental concentrations:
        cX_env[:] = cX_env_o.mean()

    return cX_cell, cX_env

def update_intra(sim, cells, cX_mems, cX_cells, D_x, zx, p):
    """
    Perform electrodiffusion on intracellular vertices
    to update concentration around the cell interior
    in response to interior voltage and concentration
    gradients, in addition to any electroosmotic flows.

    """

    # x and y components of membrane tangent unit vectors
    nx = cells.mem_vects_flat[:, 2]
    ny = cells.mem_vects_flat[:, 3]


    if p.fluid_flow is True and p.run_sim is True:
        # get intracellular fluid flow vector
        ux_mem = sim.u_cells_x[cells.mem_to_cells]
        uy_mem = sim.u_cells_y[cells.mem_to_cells]

        # component of fluid flow velocity normal to the membrane:
        u = ux_mem*nx + uy_mem*ny

    else:
        u = 0

    # get the gradient of rho concentration for each cell centre wrt to each membrane midpoint:
    grad_c = (cX_mems - cX_cells[cells.mem_to_cells])/cells.chords
    grad_v = (sim.v_cell - sim.v_cell_ave[cells.mem_to_cells])/cells.chords

    # field-modulate the grad_v to account for screening (assumes motion primarily near the double layer):
    grad_v = (1.0e-9/p.d_cell)*grad_v

    # obtain an average concentration at the pie-slice midpoints:
    c_at_mids = (cX_mems + cX_cells[cells.mem_to_cells])/2

    # flux_intra = - D_x * p.cell_delay_const * grad_c + u * c_at_mids
    flux_intra = nernst_planck_vector(c_at_mids, grad_c, grad_v, u, D_x, zx, sim.T, p)

    # net flux across inner centroid surfaces of each cell:
    net_flux = flux_intra * (cells.mem_sa / 2)

    # sum of the flux entering each centroid region:
    flux_sum = np.dot(cells.M_sum_mems, net_flux)

    # divergence of the flux wrt the centroid region:
    divF_centroid = flux_sum/cells.centroid_vol

    # divergence of the flux wrt each membrane region
    divF_mems = net_flux/cells.mem_vol

    # update concentrations with the appropriate divergence:
    cX_mems = cX_mems + divF_mems * p.dt
    cX_cells = cX_cells - divF_centroid * p.dt

    # # use finite volume method to integrate each region:
    # # values at centroid mids:
    # c_at_mids = (cX_memso + cX_cellso[cells.mem_to_cells]) / 2
    #
    # # finite volume integral of membrane pie-box values:
    # cX_mems = np.dot(cells.M_int_mems, cX_memso) + (1 / 2) * c_at_mids
    # cX_cells = (1 / 2) * cX_cellso + np.dot(cells.M_sum_mems, c_at_mids) / (2 * cells.num_mems)

    return cX_mems, cX_cells, net_flux
