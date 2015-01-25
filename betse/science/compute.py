#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

"""
A toolbox of workhorse functions and routines used in the main simulation. This is the matrix version, which
implements the simulation in terms of Numpy arrays.

"""
# FIXME implement stability safety threshhold parameter checks and loss-of-stability detection + error message
# FIXME think about testing python loop (with numba jit) versus numpy matrix versions regarding speed...

import numpy as np
from betse.science.parameters import params as p


class Simulator(object):


    def __init__(self,cells,trueVol=None,method=None):   # FIXME Sess says faster to use constants and not strings

        if method == 'Euler':
            self.mthd = 0
        elif method == 'RK4':
            self.mthd = 1
        elif method == None:
            self.mthd =0

        self.trueVol = trueVol
        self.initData(cells)


    def initData(self,cells):
        """
        Creates a host of initialized data matrices for the main simulation,
        including intracellular and environmental concentrations, voltages, and specific
        diffusion constants.

        """
        # initialization of data-arrays holding time-related information
        self.cc_time = []  # data array holding the concentrations at time points
        self.vm_time = []  # data array holding voltage at time points
        self.time = []     # time values of the simulation

        # whether to use the unique volume and surface area of cells or a single average value:
        if self.trueVol == None or self.trueVol==0:

            self.volcelli = np.mean(cells.cell_vol,axis=0)
            self.sacelli = np.mean(cells.cell_sa,axis=0)

            self.volcell = np.zeros(len(cells.cell_i))
            self.volcell[:]=self.volcelli

            self.sacell = np.zeros(len(cells.cell_i))
            self.sacell[:]=self.sacelli

        elif self.trueVol == 1:

            self.volcell = cells.cell_vol
            self.sacell = cells.cell_sa

        # Initialize cellular concentrations of ions:
        cNa_cells = np.zeros(len(cells.cell_i))
        cNa_cells[:]=p.cNa_cell

        cK_cells = np.zeros(len(cells.cell_i))
        cK_cells[:]=p.cK_cell

        cCl_cells = np.zeros(len(cells.cell_i))
        cCl_cells[:]=p.cCl_cell

        cCa_cells = np.zeros(len(cells.cell_i))
        cCa_cells[:]=p.cCa_cell

        cH_cells = np.zeros(len(cells.cell_i))
        cH_cells[:]=p.cH_cell

        cP_cells = np.zeros(len(cells.cell_i))
        cP_cells[:]=p.cP_cell

        cM_cells = np.zeros(len(cells.cell_i))
        cM_cells[:]=p.cM_cell

        # Initialize environmental ion concentrations:
        cNa_env = np.zeros(len(cells.cell_i))
        cNa_env[:]=p.cNa_env

        cK_env = np.zeros(len(cells.cell_i))
        cK_env[:]=p.cK_env

        cCl_env = np.zeros(len(cells.cell_i))
        cCl_env[:]=p.cCl_env

        cCa_env = np.zeros(len(cells.cell_i))
        cCa_env[:]=p.cCa_env

        cH_env = np.zeros(len(cells.cell_i))
        cH_env[:]=p.cH_env

        cP_env = np.zeros(len(cells.cell_i))
        cP_env[:]=p.cP_env

        cM_env = np.zeros(len(cells.cell_i))
        cM_env[:]=p.cM_env

        # Initialize membrane diffusion co-efficients:
        DmNa = np.zeros(len(cells.cell_i))
        DmNa[:] = p.Dm_Na

        DmK = np.zeros(len(cells.cell_i))
        DmK[:] = p.Dm_K

        DmCl = np.zeros(len(cells.cell_i))
        DmCl[:] = p.Dm_Cl

        DmCa = np.zeros(len(cells.cell_i))
        DmCa[:] = p.Dm_Ca

        DmH = np.zeros(len(cells.cell_i))
        DmH[:] = p.Dm_H

        DmP = np.zeros(len(cells.cell_i))
        DmP[:] = p.Dm_P

        DmM = np.zeros(len(cells.cell_i))
        DmM[:] = p.Dm_M

        # Initialize membrane thickness:
        self.tm = np.zeros(len(cells.cell_i))
        self.tm[:] = p.tm

        # Initialize environmental volume:
        self.envV = np.zeros(len(cells.cell_i))
        self.envV[:] = p.vol_env


        # Create vectors holding a range of ion-matched data
        self.cc_cells = [cNa_cells,cK_cells,cCl_cells,cCa_cells,cH_cells,cP_cells,cM_cells]  # cell concentrations
        self.cc_env = [cNa_env,cK_env,cCl_env,cCa_env,cH_env,cP_env,cM_env]   # environmental concentrations
        self.zs = [p.z_Na, p.z_K, p.z_Cl, p.z_Ca, p.z_H, p.z_P, p.z_M]   # matched ion valence state
        self.Dm_cells = [DmNa, DmK, DmCl,DmCa,DmH,DmP,DmM]              # matched membrane diffusion constants

        self.iNa=0     # indices to each ion for use in above arrays
        self.iK = 1
        self.iCl=2
        self.iCa = 3
        self.iH = 4
        self.iP = 5
        self.iM = 6

        self.movingIons = [self.iNa,self.iK,self.iCl,self.iCa,self.iH,self.iM]

        # Initialize gap-junction data
    def runSim(self,timesteps):
        """
        Drives the actual time-loop iterations for the simulation.
        """
        # create a time-steps vector:
        tt = np.linspace(0,timesteps*p.dt,timesteps)

        # report
        print('Your simulation is running from',0,'to',timesteps*p.dt,'seconds.')
        # FIXME would be nice to have a time estimate for the simulation

        for t in tt:   # run through the loop

            # get the net, unbalanced charge in each cell:
            q_cells = get_charge(self.cc_cells,self.zs,self.volcell)

            # calculate the voltage in the cell (which is also Vmem as environment is zero):
            vm = get_volt(q_cells,self.sacell)

            # run the Na-K-ATPase pump:  # FIXME would be nice to track ATP use
            self.cc_cells[self.iNa],self.cc_env[self.iNa],self.cc_cells[self.iK],self.cc_env[self.iK], fNa_NaK, fK_NaK =\
                pumpNaKATP(self.cc_cells[self.iNa],self.cc_env[self.iNa],self.cc_cells[self.iK],self.cc_env[self.iK],
                    self.volcell,self.envV,vm,method=self.mthd)

            # electro-diffuse all ions (except for proteins, which don't move!) across the cell membrane:

            for i in self.movingIons:

                self.cc_env[i],self.cc_cells[i],fNa = \
                    electrofuse(self.cc_env[i],self.cc_cells[i],self.Dm_cells[i],self.tm,self.sacell,
                        self.envV,self.volcell,self.zs[i],vm,method=self.mthd)

            # self.cc_env[self.iK],self.cc_cells[self.iK],fK =\
            #     electrofuse(self.cc_env[self.iK],self.cc_cells[self.iK],self.Dm_cells[self.iK],self.tm,self.sacell,
            #         self.envV,self.volcell,self.zs[self.iK],vm,method=self.mthd)
            #
            # self.cc_env[self.iCl],self.cc_cells[self.iCl],fCl = \
            #     electrofuse(self.cc_env[self.iCl],self.cc_cells[self.iCl],self.Dm_cells[self.iCl],self.tm,self.sacell,
            #         self.envV,self.volcell,self.zs[self.iCl],vm,method=self.mthd)
            #
            # self.cc_env[self.iCa],self.cc_cells[self.iCa],fCa =\
            #     electrofuse(self.cc_env[self.iCa],self.cc_cells[self.iCa],self.Dm_cells[self.iCa],self.tm,self.sacell,
            #         self.envV,self.volcell,self.zs[self.iCa],vm,method=self.mthd)
            #
            # self.cc_env[self.iH],self.cc_cells[self.iH],fH = \
            #     electrofuse(self.cc_env[self.iH],self.cc_cells[self.iH],self.Dm_cells[self.iH],self.tm,self.sacell,
            #         self.envV,self.volcell,self.zs[self.iH],vm,method=self.mthd)
            #
            # self.cc_env[self.iM],self.cc_cells[self.iM],fM = \
            #     electrofuse(self.cc_env[self.iM],self.cc_cells[self.iM],self.Dm_cells[self.iM],self.tm,self.sacell,
            #         self.envV,self.volcell,self.zs[self.iM],vm,method=self.mthd)

            # add the new concentration and voltage data to the time-storage matrices:
            self.cc_time.append(self.cc_cells)
            self.vm_time.append(vm)
            self.time.append(t)

def diffuse(cA,cB,Dc,d,sa,vola,volb,method=None):
    """
    Returns updated concentrations for diffusion between two
    connected volumes.

    Parameters
    ----------
    cA          Initial concentration of c in region A [mol/m3]
    cB          Initial concentration of c in region B [mol/m3]
    Dc          Diffusion constant of c  [m2/s]
    d           Distance between region A and region B [m]
    sa          Surface area separating region A and B [m2]
    vola        volume of region A [m3]
    volb        volume of region B [m3]
    dt          time step   [s]
    method      EULER or RK4 for Euler and Runge-Kutta 4
                integration methods, respectively.

    Returns
    --------
    cA2         Updated concentration of cA in region A [mol/m3]
    cB2         Updated concentration of cB in region B [mol/m3]
    flux        Chemical flux magnitude between region A and B [mol/s]

    """

    assert vola>0
    assert volb>0
    assert d >0

    volab = (vola + volb)/2
    qualityfactor = (Dc/d)*(sa/volab)*p.dt   # quality factor should be <1.0 for stable simulations

    flux = -sa*Dc*(cB - cA)/d

    if method == None or method == 0:

        dmol = sa*p.dt*Dc*(cB - cA)/d

        cA2 = cA + dmol/vola
        cB2 = cB - dmol/volb

        cA2 = check_c(cA2)
        cB2 = check_c(cB2)

    elif method == 1:

        k1 = sa*Dc*(cB - cA)/d

        k2 = sa*Dc*(cB - (cA + (1/2)*k1*p.dt))/d

        k3 = sa*Dc*(cB - (cA + (1/2)*k2*p.dt))/d

        k4 = sa*Dc*(cB - (cA + k3*p.dt))/d

        dmol = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

        cA2 = cA + dmol/vola
        cB2 = cB - dmol/volb

        cA2 = check_c(cA2)
        cB2 = check_c(cB2)

    return cA2, cB2, flux


def electrofuse(cA,cB,Dc,d,sa,vola,volb,zc,Vba,method=None):
    """
    Returns updated concentrations for electro-diffusion between two
    connected volumes. Note for cell work, 'b' is 'inside', 'a' is outside.
    and Vba is equivalent to Vmem.

    This function defaults to regular diffusion if Vba == 0.0

    Also, this function takes both matrix and singular values as input.
    However, if inputing matrix data, all inputs must be the same shape
    matrices. Likewise, if inputting singular values, all input must be
    singular values.

    Parameters
    ----------
    cA          Initial concentration of c in region A [mol/m3] (out)
    cB          Initial concentration of c in region B [mol/m3] (in)
    Dc          Diffusion constant of c  [m2/s]
    d           Distance between region A and region B [m]
    sa          Surface area separating region A and B [m2]
    vola        volume of region A [m3]
    volb        volume of region B [m3]
    zc          valence of ionic species c
    Vba         voltage between B and A as Vb - Va  [V]
    dt          time step   [s]
    method      EULER or RK4 for Euler and Runge-Kutta 4
                integration methods, respectively. 'RK4_AB' is an
                alternative implementation of RK4.

    Returns
    --------
    cA2         Updated concentration of cA in region A [mol/m3]
    cB2         Updated concentration of cB in region B [mol/m3]
    flux        Chemical flux magnitude between region A and B [mol/s]

    """

    alpha = (zc*Vba*p.F)/(p.R*p.T)

    volab = (vola + volb)/2
    qualityfactor = abs((Dc/d)*(sa/volab)*p.dt*alpha)   # quality factor should be <1.0 for stable simulations

    deno = 1 - np.exp(-alpha)   # calculate the denominator for the electrodiffusion equation,..

    if isinstance(deno,np.ndarray):   # if we have matrix data:

        izero = (deno==0).nonzero()     # get the indices of the zero and non-zero elements of the denominator
        inzero = (deno!=0).nonzero()

        # initialize data matrices to the same shape as input data
        dmol = np.zeros(deno.shape)
        cA2 = np.zeros(deno.shape)
        cB2 = np.zeros(deno.shape)
        flux = np.zeros(deno.shape)
        k1 = np.zeros(deno.shape)
        k2 = np.zeros(deno.shape)
        k3 = np.zeros(deno.shape)
        k4 = np.zeros(deno.shape)

        if len(deno[izero]):   # if there's anything in the izero array:
             # calculate the flux for those elements:
            flux[izero] = -sa[izero]*p.dt*Dc[izero]*(cB[izero] - cA[izero])/d[izero]


            if method == None or method == 0:

                dmol[izero] = sa[izero]*p.dt*Dc[izero]*(cB[izero] - cA[izero])/d[izero]

                cA2[izero] = cA[izero] + dmol[izero]/vola[izero]
                cB2[izero] = cB[izero] - dmol[izero]/volb[izero]

                cA2[izero] = check_c(cA2[izero])
                cB2[izero] = check_c(cB2[izero])

            elif method == 1:

                k1[izero] = sa[izero]*Dc[izero]*(cB[izero] - cA[izero])/d[izero]

                k2[izero] = sa[izero]*Dc[izero]*(cB[izero] - (cA[izero] + (1/2)*k1[izero]*p.dt))/d[izero]

                k3[izero] = sa[izero]*Dc[izero]*(cB[izero] - (cA[izero] + (1/2)*k2[izero]*p.dt))/d[izero]

                k4[izero] = sa[izero]*Dc[izero]*(cB[izero] - (cA[izero] + k3[izero]*p.dt))/d[izero]

                dmol[izero] = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

                cA2[izero] = cA[izero] + dmol[izero]/vola[izero]
                cB2[izero] = cB[izero] - dmol[izero]/volb[izero]

                cA2[izero] = check_c(cA2[izero])
                cB2[izero] = check_c(cB2[izero])

        if len(deno[inzero]):   # if there's any indices in the inzero array:

            # calculate the flux for those elements:
            flux[inzero] = -(sa[inzero]*Dc[inzero]/d[inzero])*alpha[inzero]*\
                           ((cB[inzero] - cA[inzero]*np.exp(-alpha[inzero]))/deno[inzero])


            if method == None or method == 0:

                dmol[inzero] = (sa[inzero]*p.dt*Dc[inzero]/d[inzero])*alpha[inzero]*\
                               ((cB[inzero] - cA[inzero]*np.exp(-alpha[inzero]))/deno[inzero])

                cA2[inzero] = cA[inzero] + dmol[inzero]/vola[inzero]
                cB2[inzero] = cB[inzero] - dmol[inzero]/volb[inzero]

                cA2[inzero] = check_c(cA2[inzero])
                cB2[inzero] = check_c(cB2[inzero])

            elif method == 1:

                k1[inzero] = (sa[inzero]*Dc[inzero]/d[inzero])*alpha[inzero]*\
                             (cB[inzero] - cA[inzero]*np.exp(-alpha[inzero]))/deno[inzero]

                k2[inzero] = (sa[inzero]*Dc[inzero]/d[inzero])*alpha[inzero]*\
                             (cB[inzero] - (cA[inzero] + (1/2)*k1[inzero]*p.dt)*np.exp(-alpha[inzero]))/deno[inzero]

                k3[inzero] = (sa[inzero]*Dc[inzero]/d[inzero])*alpha[inzero]*\
                             (cB[inzero] - (cA[inzero] + (1/2)*k2[inzero]*p.dt)*np.exp(-alpha[inzero]))/deno[inzero]

                k4[inzero] = (sa[inzero]*Dc[inzero]/d[inzero])*alpha[inzero]*\
                             (cB[inzero] - (cA[inzero] + k3[inzero]*p.dt)*np.exp(-alpha[inzero]))/deno[inzero]

                dmol[inzero] = (p.dt/6)*(k1[inzero] + 2*k2[inzero] + 2*k3[inzero] + k4[inzero])

                cA2[inzero] = cA[inzero] + dmol[inzero]/vola[inzero]
                cB2[inzero] = cB[inzero] - dmol[inzero]/volb[inzero]

                cA2[inzero] = check_c(cA2[inzero])
                cB2[inzero] = check_c(cB2[inzero])



    elif isinstance(deno,np.float64): # else if we have a single-valued data:

        if deno == 0.0:

            flux = -sa*p.dt*Dc*(cB - cA)/d

            if method == None or method == 0:

                dmol = sa*p.dt*Dc*(cB - cA)/d

                cA2 = cA + dmol/vola
                cB2 = cB - dmol/volb

                cA2 = check_c(cA2)
                cB2 = check_c(cB2)

            elif method == 1:

                k1 = sa*Dc*(cB - cA)/d

                k2 = sa*Dc*(cB - (cA + (1/2)*k1*p.dt))/d

                k3 = sa*Dc*(cB - (cA + (1/2)*k2*p.dt))/d

                k4 = sa*Dc*(cB - (cA + k3*p.dt))/d

                dmol = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

                cA2 = cA + dmol/vola
                cB2 = cB - dmol/volb

                cA2 = check_c(cA2)
                cB2 = check_c(cB2)

        else:

            flux = -(sa*Dc/d)*alpha*((cB - cA*np.exp(-alpha))/deno)

            if method == None or method == 0:

                dmol = (sa*p.dt*Dc/d)*alpha*((cB - cA*np.exp(-alpha))/deno)

                cA2 = cA + dmol/vola
                cB2 = cB - dmol/volb

                cA2 = check_c(cA2)
                cB2 = check_c(cB2)

            elif method == 1:

                k1 = (sa*Dc/d)*alpha*(cB - cA*np.exp(-alpha))/deno

                k2 = (sa*Dc/d)*alpha*(cB - (cA + (1/2)*k1*p.dt)*np.exp(-alpha))/deno

                k3 = (sa*Dc/d)*alpha*(cB - (cA + (1/2)*k2*p.dt)*np.exp(-alpha))/deno

                k4 = (sa*Dc/d)*alpha*(cB - (cA + k3*p.dt)*np.exp(-alpha))/deno

                dmol = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

                cA2 = cA + dmol/vola
                cB2 = cB - dmol/volb

                cA2 = check_c(cA2)
                cB2 = check_c(cB2)

    return cA2, cB2, flux


def pumpNaKATP(cNai,cNao,cKi,cKo,voli,volo,Vm,method=None):

    """
    Parameters
    ----------
    cNai            Concentration of Na+ inside the cell
    cNao            Concentration of Na+ outside the cell
    cKi             Concentration of K+ inside the cell
    cKo             Concentration of K+ outside the cell
    voli            Volume of the cell [m3]
    volo            Volume outside the cell [m3]
    Vm              Voltage across cell membrane [V]
    sa              Surface area of membrane
    method          EULER or RK4 solver

    Returns
    -------
    cNai2           Updated Na+ inside cell
    cNao2           Updated Na+ outside cell
    cKi2            Updated K+ inside cell
    cKo2            Updated K+ outside cell
    f_Na            Na+ flux (into cell +)
    f_K             K+ flux (into cell +)
    """

    delG_Na = p.R*p.T*np.log(cNao/cNai) - p.F*Vm
    delG_K = p.R*p.T*np.log(cKi/cKo) + p.F*Vm
    delG_NaKATP = p.deltaGATP - (3*delG_Na + 2*delG_K)
    delG = (delG_NaKATP/1000)

    alpha = p.alpha_NaK*step(delG,p.halfmax_NaK,p.slope_NaK)

    f_Na  = -alpha*cNai*cKo      #flux as [mol/s]
    f_K = -(2/3)*f_Na          # flux as [mol/s]

    if method == None or method == 0:

        dmol = -alpha*cNai*cKo*p.dt

        cNai2 = cNai + dmol/voli
        cNao2 = cNao - dmol/volo

        cKi2 = cKi - (2/3)*dmol/voli
        cKo2 = cKo + (2/3)*dmol/volo

        cNai2 = check_c(cNai2)
        cNao2 = check_c(cNao2)
        cKi2 = check_c(cKi2)
        cKo2 = check_c(cKo2)

    elif method == 1:

        k1 = alpha*cNai*cKo

        k2 = alpha*(cNai+(1/2)*k1*p.dt)*cKo

        k3 = alpha*(cNai+(1/2)*k2*p.dt)*cKo

        k4 = alpha*(cNai+ k3*p.dt)*cKo

        dmol = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

        cNai2 = cNai - dmol/voli
        cNao2 = cNao + dmol/volo

        cKi2 = cKi + (2/3)*dmol/voli
        cKo2 = cKo - (2/3)*dmol/volo

        cNai2 = check_c(cNai2)
        cNao2 = check_c(cNao2)
        cKi2 = check_c(cKi2)
        cKo2 = check_c(cKo2)

    return cNai2,cNao2,cKi2,cKo2, f_Na, f_K


def get_volt(q,sa):

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


def get_charge(concentrations,zs,vol):

    q = 0

    for conc,z in zip(concentrations,zs):
        q = q+ conc*z

    netcharge = p.F*q*vol

    return netcharge


def check_c(cA):
    """
    Does a quick check on two values (concentrations)
    and sets one to zero if it is below zero.

    """
    if isinstance(cA,np.float64):  # if we just have a singular value
        if cA < 0.0:
            cA2 = 0.0

    elif isinstance(cA,np.ndarray): # if we have matrix data
        isubzeros = (cA<0).nonzero()
        if isubzeros:  # if there's anything in the isubzeros matrix...
            cA[isubzeros] = 0.0

    return cA


def sigmoid(x,g,y_sat):
    """
    A sigmoidal function (logistic curve) allowing user
    to specify a saturation level (y_sat) and growth rate (g).

    Parameters
    ----------
    x            Input values, may be numpy array or float
    g            Growth rate
    y_sat        Level at which growth saturates

    Returns
    --------
    y            Numpy array or float of values

    """
    y = (y_sat*np.exp(g*x))/(y_sat + (np.exp(g*x)-1))
    return y

def hill(x,K,n):

    """
    The Hill equation (log-transformed sigmoid). Function ranges
    from y = 0 to +1.

    Parameters
    ----------
    x            Input values, may be numpy array or float. Note all x>0 !
    K            Value of x at which curve is 1/2 maximum (y=0.5)
    n            Hill co-efficient n<1 negative cooperativity, n>1 positive.

    Returns
    --------
    y            Numpy array or float of values

    """
    assert x.all() > 0

    y = x**n/((K**n)+(x**n))

    return y

def step(t,t_on,t_change):
    """
    A step function (bounded by 0 and 1) based on a logistic curve
    and allowing user to specify time for step to come on (t_on) and time for
    change from zero to one to happen.

    Parameters
    ----------
    t            Input values, may be numpy array or float
    t_on         Time step turns on
    t_change     Time for change from 0 to 1 (off to on)

    Returns
    --------
    y            Numpy array or float of values

    """
    g = (1/t_change)*10
    y = 1/(1 + (np.exp(-g*(t-t_on))))
    return y

def pulse(t,t_on,t_off,t_change):
    """
    A pulse function (bounded by 0 and 1) based on logistic curves
    and allowing user to specify time for step to come on (t_on) and time for
    change from zero to one to happen, and time for step to come off (t_change).

    Parameters
    ----------
    t            Input values, may be numpy array or float
    t_on         Time step turns on
    t_off        Time step turns off
    t_change     Time for change from 0 to 1 (off to on)

    Returns
    --------
    y            Numpy array or float of values

    """
    g = (1/t_change)*10
    y1 = 1/(1 + (np.exp(-g*(t-t_on))))
    y2 = 1/(1 + (np.exp(-g*(t-t_off))))
    y = y1 - y2
    return y



