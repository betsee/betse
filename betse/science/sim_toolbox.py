#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

# ....................{ IMPORTS                            }....................
import numpy as np
import numpy.ma as ma
from scipy import interpolate as interp
from scipy.ndimage.filters import gaussian_filter
# from betse.science.math import toolbox as tb
from betse.exceptions import BetseSimUnstableException
from betse.science.math import finitediff as fd

# ....................{ UTILITIES                          }....................
# Toolbox of functions used in the Simulator class to calculate key bioelectric
# properties.

def electroflux(cA,cB,Dc,d,zc,vBA,T,p,rho=1):
    """
    Electro-diffusion between two connected volumes. Note for cell work, 'b' is
    'inside', 'a' is outside, with a positive flux moving from a to b. The
    voltage is defined as Vb - Va (Vba), which is equivalent to Vmem.

    This function defaults to regular diffusion if Vba == 0.0.

    This function takes numpy matrix values as input. All inputs must be
    matrices of the same shape.

    This is the Goldman Flux/Current Equation (not to be confused with the
    Goldman Equation). Note: the Nernst-Planck equation has been trialed in
    place of this, and it does not reproduce proper reversal potentials.

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

    # Reasonably small real number, preventing divide by zero errors.
    # Note that "betse.util.type.numeric.floats.FLOAT_MIN", the
    # smallest possible real number, is too small for this use case.
    FLOAT_NONCE = 1.0e-25

    vBA += FLOAT_NONCE

    zc += FLOAT_NONCE

    alpha = (zc*vBA*p.F)/(p.R*T)

    exp_alpha = np.exp(-alpha)

    deno = -np.expm1(-alpha)   # calculate the denominator for the electrodiffusion equation,..

    # calculate the flux for those elements:
    flux = -((Dc*alpha)/d)*((cB -cA*exp_alpha)/deno)*rho

    # flux = flux*rho

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

    cATP = p.cATP
    cADP = p.cADP
    cPi  = p.cPi

    # calculate the reaction coefficient Q:
    Qnumo = (cADP*1e-3)*(cPi*1e-3)*((cNao*1e-3)**3)*((cKi*1e-3)**2)
    Qdenomo = (cATP*1e-3)*((cNai*1e-3)**3)*((cKo*1e-3)** 2)

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

    cATP = p.cATP
    cADP = p.cADP
    cPi  = p.cPi
    #

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

#FIXME: For efficiency, all calls to this function should be replaced by calls
#to the dramatically faster betse.lib.numpy.nptest.die_if_nan() function. Hola!
def check_v(vm):
    """
    Does a quick check on Vmem values
    and displays error warning or exception if the value
    indicates the simulation is unstable.
    """

    isnans = np.isnan(vm)

    if isnans.any():  # if there's anything in the isubzeros matrix...
        raise BetseSimUnstableException(
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
    dat_grid = gaussian_filter(dat_grid,1)

    # get rid of values that bleed into the environment:
    # dat_grid = np.multiply(dat_grid,cells.maskM)

    dat_grid = ma.masked_array(dat_grid, np.logical_not(cells.maskM))

    return dat_grid

def nernst_planck_flux(c, gcx, gcy, gvx, gvy,ux,uy,D,z,T,p, mu = 0.0):
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
    fx =  -D*gcx - alpha*gvx*c + ux*c - mu*c*gvx
    fy =  -D*gcy - alpha*gvy*c + uy*c - mu*c*gvy

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

    #FIXME: For efficiency, this NaN check should be replaced by a call to the
    #dramatically faster betse.lib.numpy.nptest.die_if_nan() function. O arise!

    # Ensure no NaNs.
    inds_nan = (np.isnan(data)).nonzero()

    # Ensure that data has no negative values.
    inds_neg = (data < 0.0).nonzero()

    #FIXME: This seems to contradict the documentation. It also seems a bit
    #unsafe. Shouldn't this raise an exception rather than silently replace all
    #negative values with 0.0? Or maybe this is O.K.? Cloudy marshmallows!
    if len(inds_neg[0]) > 0:
        data[inds_neg] = 0.0 # add in a small bit to protect from crashing

    if len(inds_nan[0]) > 0:

        raise BetseSimUnstableException(
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

def ghk_calculator(sim, cells, p):
    """
    Uses simulation parameters in the Goldman (GHK) equation
    to calculate an alternative Vmem for validation purposes.

    """


    # FIXME the Goldman calculator must be altered to account for network pumps and channels!!
    # begin by initializing all summation arrays for the cell network:
    sum_PmAnion_out = []
    sum_PmAnion_in = []
    sum_PmCation_out = []
    sum_PmCation_in = []

    for i, z in enumerate(sim.zs):

        # tag as anion or cation
        ion_type = np.sign(z)

        # average values from membranes or environment to cell centres:
        Dm = np.dot(cells.M_sum_mems, sim.Dm_cells[i]) / cells.num_mems
        conc_cells = sim.cc_cells[i]

        if p.is_ecm is True:
            # average entities from membranes to the cell centres:
            conc_env = np.dot(cells.M_sum_mems, sim.cc_env[i][cells.map_mem2ecm]) / cells.num_mems

        else:

            conc_env = np.dot(cells.M_sum_mems, sim.cc_env[i]) / cells.num_mems

        if ion_type == -1:

            sum_PmAnion_in.append(Dm * conc_cells * (1 / p.tm))
            sum_PmAnion_out.append(Dm * conc_env * (1 / p.tm))


        if ion_type == 1:

            sum_PmCation_in.append(Dm * conc_cells * (1 / p.tm))
            sum_PmCation_out.append( Dm * conc_env * (1 / p.tm))

    if p.molecules_enabled:

        for name in sim.molecules.core.channels:
            obj = sim.molecules.core.channels[name]

            for ii, relP in zip(obj.channel_core.ions, obj.channel_core.rel_perm):

                ion_i = sim.get_ion(ii)
                zi = sim.zs[ion_i]
                conc_cells = sim.cc_cells[ion_i]

                # tag as anion or cation
                ion_type = np.sign(zi)

                if p.is_ecm is True:
                    # average entities from membranes to the cell centres:
                    conc_env = np.dot(cells.M_sum_mems, sim.cc_env[ion_i][cells.map_mem2ecm]) / cells.num_mems

                else:

                    conc_env = np.dot(cells.M_sum_mems, sim.cc_env[ion_i]) / cells.num_mems

                if obj.channel_core.DChan is not None:
                    Dmo = obj.channel_core.DChan*relP
                    Dm = np.dot(cells.M_sum_mems, Dmo) / cells.num_mems

                else:
                    Dm = 0.0

                if ion_type == -1:
                    sum_PmAnion_in.append(Dm * conc_cells * (1 / p.tm))
                    sum_PmAnion_out.append(Dm * conc_env * (1 / p.tm))

                if ion_type == 1:
                    sum_PmCation_in.append( Dm * conc_cells * (1 / p.tm))
                    sum_PmCation_out.append(Dm * conc_env * (1 / p.tm))

    if p.grn_enabled:

        for name in sim.grn.core.channels:
            obj = sim.grn.core.channels[name]

            for ii, relP in zip(obj.channel_core.ions, obj.channel_core.rel_perm):

                ion_i = sim.get_ion(ii)
                zi = sim.zs[ion_i]
                conc_cells = sim.cc_cells[ion_i]

                # tag as anion or cation
                ion_type = np.sign(zi)

                if p.is_ecm is True:
                    # average entities from membranes to the cell centres:
                    conc_env = np.dot(cells.M_sum_mems,
                                      sim.cc_env[ion_i][cells.map_mem2ecm]) / cells.num_mems

                else:

                    conc_env = np.dot(cells.M_sum_mems, sim.cc_env[ion_i]) / cells.num_mems

                if obj.channel_core.DChan is not None:
                    Dmo = obj.channel_core.DChan * relP
                    Dm = np.dot(cells.M_sum_mems, Dmo) / cells.num_mems

                else:
                    Dm = 0.0

                if ion_type == -1:
                    sum_PmAnion_in.append(Dm * conc_cells * (1 / p.tm))
                    sum_PmAnion_out.append(Dm * conc_env * (1 / p.tm))

                if ion_type == 1:
                    sum_PmCation_in.append(Dm * conc_cells * (1 / p.tm))
                    sum_PmCation_out.append(Dm * conc_env * (1 / p.tm))


    # NaKrate = (np.dot(cells.M_sum_mems, sim.rate_NaKATP)/cells.num_mems)

    # sum together contributions for Na and K flux across the membrane:
    # NaKflux = NaKrate - (2/3)*NaKrate

    sum_PmAnion_in_i = np.sum(sum_PmAnion_in, axis = 0)
    sum_PmAnion_out_i = np.sum(sum_PmAnion_out, axis=0)
    sum_PmCation_in_i = np.sum(sum_PmCation_in, axis=0)
    sum_PmCation_out_i = np.sum(sum_PmCation_out, axis=0)

    sim.vm_GHK = ((p.R * sim.T) / p.F) * np.log(
        (sum_PmCation_out_i + sum_PmAnion_in_i) / (sum_PmCation_in_i + sum_PmAnion_out_i))

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

    if p.is_ecm is True:

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


    if p.is_ecm is False:
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

    if p.is_ecm is True:

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


    if p.is_ecm is False:
        cX_env_1_temp = cX_env_1.mean()
        cX_env_1[:] = cX_env_1_temp

    return cX_cell_1, cX_env_1, f_X

def molecule_mover(sim, cX_env_o, cX_cells, cells, p, z=0, Dm=1.0e-18, Do=1.0e-9, Dgj=1.0e-12, Ftj = 1.0, c_bound=0.0,
                   ignoreECM = True, smoothECM = False, ignoreTJ = False, ignoreGJ = False, rho = 1, cmems = None,
                   time_dilation_factor = 1.0, update_intra = False, name = "Unknown", transmem = False,
                   mu_mem = 0.0, umt = 0.0):

    """
    Transports a generic molecule across the membrane,
    through gap junctions, and if p.is_ecm is true,
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
    Ftj                 Factor influencing relative diffusion of substance across tight junction barrier
    c_bound             Concentration of molecule at global bounds (required for is_ecm True only)

    Returns
    -----------
    cX_cell_1         Updated concentration of molecule in the cell
    cX_env_1          Updated concentration of molecule in the environment

    """


    if p.is_ecm is True:

        cX_env = cX_env_o[cells.map_mem2ecm]

    else:
        cX_env = cX_env_o[:]

    if cmems is None:

        cX_mems = cX_cells[cells.mem_to_cells]

    elif len(cmems) == sim.mdl:
        cX_mems = cmems

    Dm_vect = np.ones(len(cX_mems))*Dm

    if Dm != 0.0:  # if there is some finite membrane diffusivity, exchange between cells and env space:

        IdM = np.ones(sim.mdl)

        f_X_ED = electroflux(cX_env, cX_mems, Dm_vect, p.tm*IdM, z*IdM, sim.vm, sim.T, p, rho = rho)

        if p.cluster_open is False:
            f_X_ED[cells.bflags_mems] = 0

    else:

        f_X_ED = np.zeros(sim.mdl)

    # update concentrations due to electrodiffusion:
    cX_cells, cX_mems, cX_env_o = update_Co(sim, cX_cells, cX_mems, cX_env_o, f_X_ED, cells, p,
                                            ignoreECM=ignoreECM, update_at_mems = update_intra)

    # ------------------------------------------------------------
    if ignoreGJ is False:

        # Update gap junctions using the GHK-flux equation
        fgj_X = electroflux(cX_mems[cells.mem_i],
                       cX_mems[cells.nn_i],
                       Dgj*sim.gj_block*sim.gjopen,
                       cells.gj_len*np.ones(sim.mdl),
                       z*np.ones(sim.mdl),
                       sim.vgj,
                       p.T,
                       p,
                       rho=1
                       )

        # enforce zero flux at outer boundary:
        fgj_X[cells.bflags_mems] = 0.0

        delta_cco = np.dot(cells.M_sum_mems, -fgj_X * cells.mem_sa) / cells.cell_vol

        # Calculate the final concentration change (the acceleration effectively speeds up time):
        if update_intra is False: # do the GJ transfer assuming instant mixing in the cell:
            # divergence calculation for individual cells (finite volume expression)
            cX_cells = cX_cells + p.dt*delta_cco*time_dilation_factor
            cX_mems = cX_cells[cells.mem_to_cells]

        else: # only update the membranes; calculate cell centre as the average:
            cX_cells = cX_cells + p.dt * delta_cco * time_dilation_factor
            # the 3/4 is the volume of the pie-shaped wedge half-triangle alloted to specific membrane volume:
            cX_mems += -fgj_X*(cells.mem_sa/((3/4)*cells.mem_vol))*time_dilation_factor*p.dt

    else:
        fgj_X = np.zeros(sim.mdl)

    #----Motor protein transport-----------------------------------------------------------------------------
    if update_intra is False and ignoreGJ and transmem is True:

        nx = cells.mem_vects_flat[:, 2]
        ny = cells.mem_vects_flat[:, 3]

        # concentration gradient at cell centres:
        gcc, gccmx, gccmy = cells.gradient(cX_cells)

        # mean value of concentration between two cells:
        mcc = (cX_mems[cells.nn_i] + cX_mems[cells.mem_i])/2
        # mcc = cX_cells[cells.mem_to_cells]

        # mean value of u-field between cells, treating motor protein transport like convection field:
        umtx = (sim.mtubes.mtubes_x[cells.nn_i] + sim.mtubes.mtubes_x[cells.mem_i])/2
        umty = (sim.mtubes.mtubes_y[cells.nn_i] + sim.mtubes.mtubes_y[cells.mem_i])/2

        # umtx = sim.mtubes.uxmt[cells.mem_to_cells]
        # umty = sim.mtubes.uymt[cells.mem_to_cells]

        emtx = (sim.E_cell_x[cells.mem_to_cells][cells.nn_i] + sim.E_cell_x[cells.mem_to_cells][cells.nn_i])/2
        emty = (sim.E_cell_y[cells.mem_to_cells][cells.nn_i] + sim.E_cell_y[cells.mem_to_cells][cells.nn_i])/2

        # emtx = sim.E_cell_x[cells.mem_to_cells]
        # emty = sim.E_cell_y[cells.mem_to_cells]

        alpha = (Do*p.q*z)/(p.kb*sim.T)

        # bgrad_x = ((umtx[cells.bflags_mems]*umt + (mu_mem + alpha)*emtx[cells.bflags_mems])*mcc[cells.bflags_mems])/Do
        # bgrad_y = ((umty[cells.bflags_mems]*umt + (mu_mem + alpha)*emty[cells.bflags_mems])*mcc[cells.bflags_mems])/Do
        #
        # gccmx[cells.bflags_mems] = bgrad_x
        # gccmy[cells.bflags_mems] = bgrad_y

        flux_mtx = -Do*gccmx + umtx*mcc*umt + (mu_mem + alpha)*emtx*mcc
        flux_mty = -Do*gccmy + umty*mcc*umt + (mu_mem + alpha)*emty*mcc

        # flux_mtx[cells.bflags_mems] = -umtx[cells.bflags_mems]*mcc[cells.bflags_mems]*umt
        # flux_mty[cells.bflags_mems] = -umty[cells.bflags_mems]*mcc[cells.bflags_mems]*umt

        flux_mtn = flux_mtx*nx  + flux_mty*ny
        # flux_mtn = -Do*gcc + sim.mtubes.umtn*mcc*umt

        flux_mtn[cells.bflags_mems] = 0.0

        div_ccmt = -np.dot(cells.M_sum_mems, flux_mtn*cells.mem_sa)/cells.cell_vol

        # update cell concentration:
        cX_cells += div_ccmt*p.dt*time_dilation_factor

        cX_mems = cX_cells[cells.mem_to_cells]*1



    # Transport through environment, if p.is_ecm is True-----------------------------------------------------
    if p.is_ecm is True:

        env_check = len((cX_env_o != 0.0).nonzero()[0])

        if  env_check != 0.0 or c_bound > 1.0e-15:

            cenv = cX_env_o
            cenv = cenv.reshape(cells.X.shape)

            cenv[:, 0] = c_bound
            cenv[:, -1] = c_bound
            cenv[0, :] = c_bound
            cenv[-1, :] = c_bound

            denv_multiplier = np.ones(len(cells.xypts))

            if ignoreTJ is False:
                # if tight junction barrier applies, create a mask that defines relative strength of barrier:
                denv_multiplier = denv_multiplier.reshape(cells.X.shape)*sim.D_env_weight

                # at the cluster boundary, further modify the env diffusion map by a relative TJ diffusion factor:
                denv_multiplier.ravel()[sim.TJ_targets] = denv_multiplier.ravel()[sim.TJ_targets]*Ftj

            else:
                denv_multiplier = denv_multiplier.reshape(cells.X.shape)

            gcx, gcy = fd.gradient(cenv, cells.delta)

            if p.fluid_flow is True:

                ux = sim.u_env_x.reshape(cells.X.shape)
                uy = sim.u_env_y.reshape(cells.X.shape)

            else:

                ux = 0.0
                uy = 0.0

            fx, fy = nernst_planck_flux(cenv, gcx, gcy, -sim.E_env_x, -sim.E_env_y, ux, uy,
                                            denv_multiplier*Do, z, sim.T, p, mu = mu_mem)

            div_fa = fd.divergence(-fx, -fy, cells.delta, cells.delta)

            fenvx = fx
            fenvy = fy

            cenv = cenv + div_fa * p.dt*time_dilation_factor

            if p.sharpness < 1.0:

                cenv = fd.integrator(cenv, sharp = p.sharpness)

            cX_env_o = cenv.ravel()

        else:

            fenvx = np.zeros(sim.edl)
            fenvy = np.zeros(sim.edl)

    else:
        cX_env_temp = cX_env_o.mean()
        cX_env_o = np.zeros(sim.mdl)
        cX_env_o[:] = cX_env_temp
        fenvx = 0
        fenvy = 0

    # check for sub-zero concentrations:
    indsZm = (cX_mems < 0.0).nonzero()[0]


    indsZc = (cX_cells < 0.0).nonzero()[0]

    if len(indsZc) > 0:
        raise BetseSimUnstableException(
            "Network concentration of " + name + " in cells below zero! Your simulation has"
                                                   " become unstable.")

    if len(indsZm) > 0:
        raise BetseSimUnstableException(
            "Network concentration of " + name + " on membrane below zero! Your simulation has"
                                                   " become unstable.")

    indsZe = (cX_env_o < 0.0).nonzero()[0]

    if len(indsZe) > 0:
        print(name, c_bound)
        raise BetseSimUnstableException(
            "Network concentration of " + name + " in environment below zero! Your simulation has"
                                                   " become unstable.")

    return cX_env_o, cX_cells, cX_mems, f_X_ED, fgj_X, fenvx, fenvy

def update_Co(sim, cX_cell, cX_mem, cX_env, flux, cells, p, ignoreECM = True, update_at_mems = False):
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
    if update_at_mems is False: # treat cell mem and centre values as equal
        cX_cell = cX_cell + delta_cells * p.dt
        cX_mem = cX_cell[cells.mem_to_cells]

    else: # update only the membranes, treat cell centre as the average:
        cX_mem += flux*(cells.mem_sa/((3/4)*cells.mem_vol))*p.dt # the 3/4 is the volume of the pi wedge half-triangle
        cX_cell = cX_cell + delta_cells * p.dt


    if p.is_ecm is True:

        # delta_env = np.dot(cells.M_divmap_mem2ecm, -flux)
        delta_env = div_env(-flux, cells,p)

        # update the environmental concentrations:
        cX_env = cX_env + delta_env * p.dt


    else:

        delta_env = -flux * (cells.mem_sa / p.vol_env)

        cX_env_o = cX_env + delta_env * p.dt

        # assume auto-mixing of environmental concentrations:
        cX_env = cX_env_o.mean()

    return cX_cell, cX_mem, cX_env

def div_env(flux, cells, p):

    """
    Calculates the divergence of the normal net flux component defined on a set of cell membranes
    with respect to the environmental grid spaces.

    :param flux:
    :param cells:
    :return:
    """

    if p.fast_update_ecm:

        flux_env = np.zeros(len(cells.xypts))
        flux_env[cells.map_mem2ecm] = flux

        delta_env = (flux_env * cells.memSa_per_envSquare) / (cells.ecm_vol)

    else:
        # Method # 2:
        delta_env = np.dot(cells.M_divmap_mem2ecm, flux)/(p.cell_height*cells.delta**2)

        # use the "integrator" function to conservatively distribute this exchange to nearest neighbours of the env grid:
        # delta_env = fd.integrator(delta_env.reshape(cells.X.shape), sharp = 0.5).ravel()

    return delta_env

def HH_Decomp(JJx, JJy, cells, bounds = None, sigma = 1.0):

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

    divJd = fd.divergence(JJx.reshape(cells.X.shape)/sigma, JJy.reshape(cells.X.shape)/sigma, cells.delta, cells.delta)

    # set boundary conditions for normal component
    # enforce applied voltage condition at the boundary:
    divJd[:, 0] = -Lb * (1 / cells.delta ** 2)
    divJd[:, -1] = -Rb * (1 / cells.delta ** 2)
    divJd[0, :] = -Bb * (1 / cells.delta ** 2)
    divJd[-1, :] = -Tb * (1 / cells.delta ** 2)


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

def smooth_flux(Fxo, Fyo, cells):

    """
    Uses the Helmholtz-Hodge decomposition to smooth a vector field in a manner that does not change the physics.
    Fxo, Fyo must be arrays defined on the environmental grid

    """

    _, Frx, Fry, _, Fdx, Fdy = HH_Decomp(Fxo, Fyo, cells)

    # free current density at each membrane, "smoothed" using Helmholtz-Hodge decomposition and reconstruction:
    Fx = Frx + Fdx
    Fy = Fry + Fdy


    return Fx, Fy

def map_to_cells(Fx, Fy, cells, p, smoothing=1.0):
    """
    Takes a vector field defined on the environmental spaces and maps it to the cell grid.

    """
    if smoothing > 0.0:
        # smooth environmental electric field:
        F_env_x = gaussian_filter(Fx.reshape(cells.X.shape), smoothing).ravel()
        F_env_y = gaussian_filter(Fy.reshape(cells.X.shape), smoothing).ravel()

    else:

        F_env_x = Fx
        F_env_y = Fy

    # env electric field mapped to cells:
    Fx_atmem = F_env_x.ravel()[cells.map_mem2ecm]
    Fy_atmem = F_env_y.ravel()[cells.map_mem2ecm]

    # map env x to cell membrane component:
    Fmem = (Fx_atmem * cells.mem_vects_flat[:, 2] +
            Fy_atmem * cells.mem_vects_flat[:, 3])

    Fx_atcell = (np.dot(cells.M_sum_mems, Fx_atmem * cells.mem_sa) / cells.cell_sa)
    Fy_atcell = (np.dot(cells.M_sum_mems, Fy_atmem * cells.mem_sa) / cells.cell_sa)

    return Fx_atcell, Fy_atcell, Fmem

def get_conductivity(D, z, c, d, p):
    """
    Simple function returns conductivity (S/m)
    for diffusion constant D (m2/s), mean
    concentration c, and membrane thickness 'd'.

    """
    # (D * (z ** 2) * p.F * c) / (d * p.R * p.T)

    return (D * p.q *(z**2) * p.F * c) / (d * p.kb * p.T)

def rotate_field(gFxo, gFyo, angle):
    """
    Rotates a vector field gFxo, gFyo by a specified angle (in degrees)

    """
    # rotate the axis of the model:
    rotangle = (angle * np.pi) / 180.0

    gFx = gFxo * np.cos(rotangle) - gFyo * np.sin(rotangle)
    gFy = gFxo * np.sin(rotangle) + gFyo * np.cos(rotangle)

    return gFx, gFy



