#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

#FIXME: This submodule doesn't appear to be used by the codebase and should
#ideally be shifted elsewhere. Since this appears to be a collection of lambda
#functions derived from Channelpedia models, would the "doc/example" directory
#perhaps be a more suitable location for this?

import numpy as np

#-------------------------------------------------------------------

# Equivalent notations

#dm_dt = (mInf-m)/mTau
#dh_dt = (hInf-h)/hTau

# dm_dt = (alpha_m*(1- m) - beta_m*(m))

# m_o = mInf
# h_o = hInf

# mInf = mAlpha/(mAlpha + mBeta)
# mTau = 1/(mAlpha + mBeta)

#hInf = hAlpha/(hAlpha + hBeta)
# hTau = 1/(hAlpha + hBeta)


# NaV1.3------------------------------------------------------------

# Neonatal sodium channel

mpower = 3
hpower = 1

mpower = 3
hpower = 1

mAlpha = lambda V: (0.182 * ((V)- -26))/(1-(np.exp(-((V)- -26)/9)))
mBeta = lambda V: (0.124 * (-(V) -26))/(1-(np.exp(-(-(V) -26)/9)))

#
hInf = lambda V: 1 /(1 + np.exp((V-(-65.0))/8.1))
hTau = lambda V: 0.40 + (0.265*np.exp(-V/9.47))


# NaV1.6 -----------------------------------------------------------
# Persistent current sodium channel
mpower = 1
mInf = lambda V: 1.0000/(1+ np.exp(-0.03937*4.2*(V +17.000)))
mTau = 1


#NaV1.5--------------------------------------------------------------
# Cardiac, mutation cancer



# Na_Rat1-----------------------------------------------------------
mpower = 3
hpower = 1

mAlpha = lambda V: 0.091*(V+38)/(1-np.exp((-V-38)/5))
mBeta = lambda V: -0.062*(V+38)/(1-np.exp((V+38)/5))
hAlpha = lambda V: 0.016*np.exp((-55-V)/15)
hBeta = lambda V: 2.07/(np.exp((17-V)/21)+1)


# HH Na Channel------------------------------------------------------
mpower = 3
hpower = 1


mAlpha = lambda V: (0.1*(25-V))/(np.exp((25 - V)/10) -1.0)
mBeta = lambda V: 4.0 * (np.exp(-V/18))
hAlpha = lambda V: 0.07 * np.exp(-V/20)
hBeta = lambda V: 1/(np.exp((30- V)/10) + 1.0)


# Kv1.5  # cardiac, other tissues-----------------------------------

mpower = 1
hpower = 1

mInf = lambda V: 1.0000/(1+ np.exp(-(V + 6.0000)/6.4000))
mTau = lambda V: (-0.1163 * V) + 8.3300
hInf = lambda V: 1.0000/(1+ np.exp((V + 25.3000)/3.5000))
hTau = lambda V: (-15.5000 * V) + 1620.0000

# Kv1.2 # delayed rectifier -----------------------------------------
mpower = 1
hpower = 1

mInf = lambda V: 1.0000/(1+ np.exp(-(V +21.0000)/11.3943))
mTau = lambda V: 150.0000/(1+ np.exp((V + 67.5600)/34.1479))
hInf = lambda V: 1.0000/(1+ np.exp((V + 22.0000)/11.3943))
hTau = lambda V: 15000.0000/(1+ np.exp(-(V + 46.5600)/44.1479))

# Kv3.3 Fast K+ channel----------------------------------------------
# Rooma, Desai, et al 2008
mpower = 2
hpower = 1

mInf = lambda V: 1/(1+np.exp((V-35)/-7.3))
mTau = lambda V: 0.676808 +( 27.913114 / (1 + np.exp((V - 22.414149)/9.704638)))
hInf = lambda V: 0.25+( 0.75 /(1+ np.exp((V-(-28.293856))/29.385636)))
hTau = lambda V: 199.786728 + (2776.119438*np.exp(-V/7.309565))


# Kv Rat Fast ------------------------------------------------------------

mpower = 1
hpower = 1

mInf = lambda V: 1/(1 + np.exp(-(V+47)/29))
mTau = lambda V: (0.34+0.92*np.exp(-((V+71)/59)**2))
hInf = lambda V: 1/(1 + np.exp(-(V+56)/-10))
hTau = lambda V: (8+49*np.exp(-((V+73)/23)**2))

# HH K Channel-----------------------------------------------------

mpower = 4

mAlpha = lambda V: (0.01*(10 - V))/(np.exp((10 - V)/10) - 1.0)
mBeta = lambda V: 0.125 * (np.exp(-V/80))
# mInf = mAlpha/(mAlpha + mBeta)
# mTau = 1/(mAlpha + mBeta)


# Kir2.1-----------------------------------------------------------
# Inward rectifying K+ channel
mpower = 1
hpower = 2

mInf = lambda V: 1 /(1 + np.exp((V-(-96.48))/23.26))
mTau = lambda V: 3.7 + (-3.37 / (1 + np.exp((V - -32.9)/27.93)))
hInf = lambda V: 1 /(1 + np.exp((V-(-168.28))/-44.13))
hTau = lambda V: 0.85 + (306.3 / (1 + np.exp((V - -118.29)/-27.23)))

# HCN2----------------------------------------------------------------
mpower = 1

mInf = lambda V: 1.0000/(1+ np.exp((V + 99)/6.2))
mTau = 184.0000

# HCN4-----------------------------------------------------------------
mpower = 1

mInf = lambda V: 1.0000/(1+np.exp((V + 100)/9.6))
mTau = lambda V: 461.0000

# CaV3.1---------------------------------------------------------------
mpower = 1
hpower = 1

mInf = lambda V: 1 /(1+np.exp((V-(-42.921064))/-5.163208))
mTau = lambda V: -0.855809 + (1.493527 * np.exp(-V/27.414182))
hInf = lambda V: 1 /(1+ np.exp((V-(-72.907420))/4.575763))
hTau = lambda V: 9.987873 + (0.002883 * np.exp(-V/5.598574))


# CaV3.3--------------------------------------------------------------
mpower = 1
hpower = 1

mInf = lambda V: 1/(1 + np.exp((V- -45.454426)/-5.073015))
mTau = lambda V: 3.394938 +( 54.187616 / (1 + np.exp((V - -40.040397)/4.110392)))
hInf = lambda V: 1 /(1+ np.exp((V-(-74.031965))/8.416382))
hTau = lambda V: 109.701136 + (0.003816 * np.exp(-V/4.781719))

# L-Type Ca----------------------------------------------------------
mpower = 2
hpower = 1

mInf = lambda V: 1.0000/(1+ np.exp(((V-10)+30.000)/-6))
mTau = lambda V: 5.0000 + 20.0000/(1+np.exp(((V-10)+25.000)/5))
hInf = lambda V: 1.0000/(1+ np.exp(((V-10)+80.000)/6.4))
hTau = lambda V: 20.0000 + 50.0000/(1+np.exp(((V-10)+40.000)/7))
