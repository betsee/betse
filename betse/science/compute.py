#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

"""
About this module!



"""
#FIXME implement pumpNaK
#FIXME implement volt_calc
#FIXME check total concentration is conserved in diffuse and implement stability check to % change threshhold

import math
import numpy as np
from betse.science.parameters import params as p

def diffuse(cA,cB,Dc,d,sa,vola,volb,dt,method=None):
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
    qualityfactor = (Dc/d)*(sa/volab)*dt   # quality factor should be <1.0 for stable simulations

    flux = -sa*Dc*(cB - cA)/d

    if method == None or method == 'Euler':

        dmol = sa*dt*Dc*(cB - cA)/d

        cA2 = cA + dmol/vola
        cB2 = cB - dmol/volb

        if cA2 < 0.0 and cB2 < 0.0:
            cA2 =0.0
            cB2 =0.0
        elif cA2 < 0.0 and cB2 >=0.0:
            cA2 = 0.0
        elif cB2 < 0.0 and cA2 >=0.0:
            cB2 = 0.0

    elif method == 'RK4':

        k1 = sa*Dc*(cB - cA)/d

        k2 = sa*Dc*(cB - (cA + (1/2)*k1*dt))/d

        k3 = sa*Dc*(cB - (cA + (1/2)*k2*dt))/d

        k4 = sa*Dc*(cB - (cA + k3*dt))/d

        dmol = (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

        cA2 = cA + dmol/vola
        cB2 = cB - dmol/volb

        if cA2 < 0.0 and cB2 < 0.0:
            cA2 =0.0
            cB2 =0.0
        elif cA2 < 0.0 and cB2 >=0.0:
            cA2 = 0.0
        elif cB2 < 0.0 and cA2 >=0.0:
            cB2 = 0.0

    return cA2, cB2, flux


def electrofuse(cA,cB,Dc,d,sa,vola,volb,zc,Vba,T,dt,method=None):
    """
    Returns updated concentrations for electro-diffusion between two
    connected volumes.

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
    T           temperature [K]
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

    alpha = (zc*Vba*p.F)/(p.R*T)

    volab = (vola + volb)/2
    qualityfactor = abs((Dc/d)*(sa/volab)*dt*alpha)   # quality factor should be <1.0 for stable simulations

    deno = 1 - math.exp(-alpha)

    if deno == 0.0:

        flux = -sa*dt*Dc*(cB - cA)/d

        if method == None or method == 'Euler':

            dmol = sa*dt*Dc*(cB - cA)/d

            cA2 = cA + dmol/vola
            cB2 = cB - dmol/volb

            if cA2 < 0.0 and cB2 < 0.0:
                cA2 =0.0
                cB2 =0.0
            elif cA2 < 0.0 and cB2 >=0.0:
                cA2 = 0.0
            elif cB2 < 0.0 and cA2 >=0.0:
                cB2 = 0.0

        elif method == 'RK4':

            k1 = sa*Dc*(cB - cA)/d

            k2 = sa*Dc*(cB - (cA + (1/2)*k1*dt))/d

            k3 = sa*Dc*(cB - (cA + (1/2)*k2*dt))/d

            k4 = sa*Dc*(cB - (cA + k3*dt))/d

            dmol = (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

            cA2 = cA + dmol/vola
            cB2 = cB - dmol/volb

            if cA2 < 0.0 and cB2 < 0.0:
                cA2 =0.0
                cB2 =0.0
            elif cA2 < 0.0 and cB2 >=0.0:
                cA2 = 0.0
            elif cB2 < 0.0 and cA2 >=0.0:
                cB2 = 0.0

    else:

        flux = -(sa*Dc/d)*alpha*((cB - cA*math.exp(-alpha))/deno)

        if method == None or method == 'Euler':

            dmol = (sa*dt*Dc/d)*alpha*((cB - cA*math.exp(-alpha))/deno)

            cA2 = cA + dmol/vola
            cB2 = cB - dmol/volb

            if cA2 < 0.0 and cB2 < 0.0:
                cA2 =0.0
                cB2 =0.0
            elif cA2 < 0.0 and cB2 >=0.0:
                cA2 = 0.0
            elif cB2 < 0.0 and cA2 >=0.0:
                cB2 = 0.0

        elif method == 'RK4':

            k1 = (sa*Dc/d)*alpha*(cB - cA*math.exp(-alpha))/deno

            k2 = (sa*Dc/d)*alpha*(cB - (cA + (1/2)*k1*dt)*math.exp(-alpha))/deno

            k3 = (sa*Dc/d)*alpha*(cB - (cA + (1/2)*k2*dt)*math.exp(-alpha))/deno

            k4 = (sa*Dc/d)*alpha*(cB - (cA + k3*dt)*math.exp(-alpha))/deno

            dmol = (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

            cA2 = cA + dmol/vola
            cB2 = cB - dmol/volb

            if cA2 < 0.0 and cB2 < 0.0:
                cA2 =0.0
                cB2 =0.0
            elif cA2 < 0.0 and cB2 >=0.0:
                cA2 = 0.0
            elif cB2 < 0.0 and cA2 >=0.0:
                cB2 = 0.0

    return cA2, cB2, flux

