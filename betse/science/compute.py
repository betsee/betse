#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

"""
About this module!



"""
#FIXME implement volt_calc
#FIXME check total concentration is conserved in diffuse and implement stability check to % change threshhold
# FIXME have pumpNaK update concentrations

import math
import numpy as np
from betse.science.parameters import params as p

def check_c(cA):
    """
    Does a quick check on two values (concentrations)
    and sets one to zero if it is below zero.

    """
    if cA < 0.0:
        cA2 = 0.0

    return cA

def pumpNaKATP(cNai,cNao,cKi,cKo,voli,volo,Vm,sa,method=None):

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
    method          'Euler' or 'RK4' solver

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

    f_Na  = alpha*cNai*cKo*sa      #flux as [mol/s]
    f_K = -(2/3)*f_Na*sa          # flux as [mol/s]

    if method == None or method == 'Euler':

        dmol = alpha*cNai*cKo*sa*p.dt

        cNai2 = cNai + dmol/voli
        cNao2 = cNao - dmol/volo

        cKi2 = cKi + (2/3)*dmol/voli
        cKo2 = cKo - (2/3)*dmol/volo

        cNai2 = check_c(cNai2)
        cNao2 = check_c(cNao2)
        cKi2 = check_c(cKi2)
        cKo2 = check_c(cKo2)

    elif method == 'RK4':

        k1 = alpha*cNai*cKo*sa

        k2 = alpha*(cNai+(1/2)*k1*p.dt)*cKo*sa

        k3 = alpha*(cNai+(1/2)*k2*p.dt)*cKo*sa

        k4 = alpha*(cNai+ k3*p.dt)*cKo*sa

        dmol = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

        cNai2 = cNai + dmol/voli
        cNao2 = cNao - dmol/volo

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
    method      'Euler' or 'RK4' for Euler and Runge-Kutta 4
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

    if method == None or method == 'Euler':

        dmol = sa*p.dt*Dc*(cB - cA)/d

        cA2 = cA + dmol/vola
        cB2 = cB - dmol/volb

        cA2 = check_c(cA2)
        cB2 = check_c(cB2)

    elif method == 'RK4':

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
    method      'Euler' or 'RK4' for Euler and Runge-Kutta 4
                integration methods, respectively. 'RK4_AB' is an
                alternative implementation of RK4.

    Returns
    --------
    cA2         Updated concentration of cA in region A [mol/m3]
    cB2         Updated concentration of cB in region B [mol/m3]
    flux        Chemical flux magnitude between region A and B [mol/s]

    """

    assert vola>0
    assert volb>0
    assert d >0

    alpha = (zc*Vba*p.F)/(p.R*p.T)

    volab = (vola + volb)/2
    qualityfactor = abs((Dc/d)*(sa/volab)*p.dt*alpha)   # quality factor should be <1.0 for stable simulations

    deno = 1 - math.exp(-alpha)

    if deno == 0.0:

        flux = -sa*p.dt*Dc*(cB - cA)/d

        if method == None or method == 'Euler':

            dmol = sa*p.dt*Dc*(cB - cA)/d

            cA2 = cA + dmol/vola
            cB2 = cB - dmol/volb

            cA2 = check_c(cA2)
            cB2 = check_c(cB2)

        elif method == 'RK4':

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

        flux = -(sa*Dc/d)*alpha*((cB - cA*math.exp(-alpha))/deno)

        if method == None or method == 'Euler':

            dmol = (sa*p.dt*Dc/d)*alpha*((cB - cA*math.exp(-alpha))/deno)

            cA2 = cA + dmol/vola
            cB2 = cB - dmol/volb

            cA2 = check_c(cA2)
            cB2 = check_c(cB2)

        elif method == 'RK4':

            k1 = (sa*Dc/d)*alpha*(cB - cA*math.exp(-alpha))/deno

            k2 = (sa*Dc/d)*alpha*(cB - (cA + (1/2)*k1*p.dt)*math.exp(-alpha))/deno

            k3 = (sa*Dc/d)*alpha*(cB - (cA + (1/2)*k2*p.dt)*math.exp(-alpha))/deno

            k4 = (sa*Dc/d)*alpha*(cB - (cA + k3*p.dt)*math.exp(-alpha))/deno

            dmol = (p.dt/6)*(k1 + 2*k2 + 2*k3 + k4)

            cA2 = cA + dmol/vola
            cB2 = cB - dmol/volb

            cA2 = check_c(cA2)
            cB2 = check_c(cB2)

    return cA2, cB2, flux

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

