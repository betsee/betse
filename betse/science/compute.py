#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

"""
About this module!



"""

#FIXME implement electrofuse
#FIXME implement pumpNaK
#FIXME implement volt_calc

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
                integration methods, respectively

    Returns
    --------
    cA2         Updated concentration of cA in region A
    cB2         Updated concentration of cB in region B

    """

    assert vola>0
    assert volb>0
    assert d >0

    flux = -Dc*(cB - cA)/d

    if method == None or method == 'Euler':

        dcA = sa*dt*Dc*(cB - cA)/(d*vola)
        dcB = -sa*dt*Dc*(cB - cA)/(d*volb)
        cA2 = cA + dcA
        cB2 = cB + dcB

        if cA2 >= 0.0 and cB2 >= 0.0:
            pass
        elif cA2 >= 0.0:
            cB2 = 0.0
        elif cB2 >= 0.0:
            cA2 = 0.0

    elif method == 'RK4':

        k1A = sa*Dc*(cB - cA)/(d*vola)
        k1B = -sa*Dc*(cB - cA)/(d*volb)

        k2A = sa*Dc*(cB - (cA + (1/2)*k1A*dt))/(d*vola)
        k2B = -sa*Dc*((cB + (1/2)*k1B*dt) - cA)/(d*volb)

        k3A = sa*Dc*(cB - (cA + (1/2)*k2A*dt))/(d*vola)
        k3B = -sa*Dc*((cB + (1/2)*k2B*dt) - cA)/(d*volb)

        k4A = sa*Dc*(cB - (cA + k3A*dt))/(d*vola)
        k4B = -sa*Dc*((cB + k2B*dt) - cA)/(d*volb)

        cA2 = cA + (dt/6)*(k1A + 2*k2A + 2*k3A + k4A)
        cB2 = cB + (dt/6)*(k1B + 2*k2B + 2*k3B + k4B)

        if cA2 >= 0.0 and cB2 >= 0.0:
            pass
        elif cA2 >= 0.0:
            cB2 = 0.0
        elif cB2 >= 0.0:
            cA2 = 0.0

    return cA2, cB2, flux