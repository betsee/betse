#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
from betse.exceptions import BetseExceptionLambda
from betse.util.type import types



class Gap_Junction(object):
    """
    Defines functions controling gap junction voltage gating characteristics.

    From Harris et al. J of Neurosci.(1983) 3:79-100.
    For Ambystoma Mexicanum (Axolotl!) early embryo gap junction gating.

    """

    def __init__(self,p):

        # voltage-dependent rate constant for channel opening (1/ms, V in mV):
        self.alpha_gj = lambda V: 0.0013*np.exp(-0.077*(V*1e3 - p.gj_vthresh))

        # voltage-dependent rate constant for channel closing (1/ms, V in mV):
        self.beta_gj = lambda V: 0.0013*np.exp(0.14*(V*1e3 - p.gj_vthresh))

        self.gmin = 0.04

        self.beta_gj_p = lambda V: 0.0013*np.exp(0.14*(V*1e3 - p.gj_vthresh))/(1+50*0.0013*np.exp(0.14*(V*1e3 -
                                                                                                        p.gj_vthresh)))


#----------------------------------------------------------------------------------------------------------------------
#  Voltage Gated Sodium Channels
#----------------------------------------------------------------------------------------------------------------------

def vgNa_squid(self, sim, dyna, p):

    V = sim.vm[dyna.targets_vgNa]*1000 + 50.0

    alpha_m = (0.1*(25-V))/(np.exp((25-V)/10)-1)
    beta_m = 4.0*np.exp(-V/18)


    alpha_h = 0.07*np.exp(-V/20)
    beta_h = 1/(1 + np.exp((30-V)/10))

    # calculate m channels
    dyna.m_Na = (alpha_m*(1-dyna.m_Na) - beta_m*(dyna.m_Na))*p.dt*1e3 + dyna.m_Na

    dyna.h_Na = (alpha_h*(1-dyna.h_Na) - beta_h*(dyna.h_Na))*p.dt*1e3 + dyna.h_Na

    # as equations are sort of ill-behaved, threshhold to ensure 0 to 1 status
    inds_mNa_over = (dyna.m_Na > 1.0).nonzero()
    dyna.m_Na[inds_mNa_over] = 1.0

    inds_hNa_over = (dyna.h_Na > 1.0).nonzero()
    dyna.h_Na[inds_hNa_over] = 1.0

    inds_mNa_under = (dyna.m_Na < 0.0).nonzero()
    dyna.m_Na[inds_mNa_under] = 0.0

    inds_hNa_under = (dyna.h_Na < 0.0).nonzero()
    dyna.h_Na[inds_hNa_under] = 0.0

    gNa_max = 4.0e-14

    # Define ultimate activity of the vgNa channel:
    sim.Dm_vg[sim.iNa][dyna.targets_vgNa] = gNa_max*(dyna.m_Na**3)*(dyna.h_Na)

def vgNa_rat(self,sim,dyna,p):
    '''
    Handle all **targeted voltage-gated sodium channels** (i.e., only
    applicable to specific tissue profiles) specified by the passed
    user-specified parameters on the passed tissue simulation and cellular
    world for the passed time step.

    Channel model uses Hodgkin-Huxley model for voltage gated sodium channels.

    '''

    V = sim.vm[dyna.targets_vgNa]*1000

    alpha_m = (0.091*(V+38))/(1-np.exp((-V-38)/5))
    beta_m = (-0.062*(V+38))/(1-np.exp((V+38)/5))


    alpha_h = 0.016*np.exp((-55-V)/15)
    beta_h = 2.07/(1 + np.exp((17-V)/21))

    # calculate m channels
    dyna.m_Na = (alpha_m*(1-dyna.m_Na) - beta_m*(dyna.m_Na))*p.dt*1e3 + dyna.m_Na

    dyna.h_Na = (alpha_h*(1-dyna.h_Na) - beta_h*(dyna.h_Na))*p.dt*1e3 + dyna.h_Na

    # as equations are sort of ill-behaved, threshhold to ensure 0 to 1 status
    inds_mNa_over = (dyna.m_Na > 1.0).nonzero()
    dyna.m_Na[inds_mNa_over] = 1.0

    inds_hNa_over = (dyna.h_Na > 1.0).nonzero()
    dyna.h_Na[inds_hNa_over] = 1.0

    inds_mNa_under = (dyna.m_Na < 0.0).nonzero()
    dyna.m_Na[inds_mNa_under] = 0.0

    inds_hNa_under = (dyna.h_Na < 0.0).nonzero()
    dyna.h_Na[inds_hNa_under] = 0.0

    gNa_max = 5.0e-14

    # Define ultimate activity of the vgNa channel:
    sim.Dm_vg[sim.iNa][dyna.targets_vgNa] = gNa_max*(dyna.m_Na**3)*(dyna.h_Na)

def vgNa_1p6(self,sim,dyna,p):

    pass

#---------------------------------------------------------------------------------------------------------------------
# Voltage Gated Potassium Channels
#---------------------------------------------------------------------------------------------------------------------


def vgK_slow(self,sim,dyna,p):

    pass

def vgK_fast(self,sim,dyna,p):
    """
    A "fast" voltage gated potassium channel from:
    A Korngreen et. al; J. Physiol. (Lond.) 2000 Jun 15

    Has a very wide response making firing of action potentials
    difficult.

    """

    # detecting channels to turn on:
    V = sim.vm[dyna.targets_vgK]*1000

    m_Inf = 1/(1 + np.exp(-(V+47)/29))
    m_Tau = (0.34+0.92*np.exp(-((V+71)/59)**2))
    h_Inf = 1/(1 + np.exp(-(V+56)/-10))
    h_Tau = (8+49*np.exp(-((V+73)/23)**2))

    dyna.m_K = ((m_Inf - dyna.m_K)/m_Tau)*p.dt*1e3 + dyna.m_K
    dyna.h_K = ((h_Inf - dyna.h_K)/h_Tau)*p.dt*1e3 + dyna.h_K

    gK_max = 1.0e-14

    inds_mK_over = (dyna.m_K > 1.0).nonzero()
    dyna.m_K[inds_mK_over] = 1.0

    inds_mK_under = (dyna.m_K < 0.0).nonzero()
    dyna.m_K[inds_mK_under] = 0.0

    inds_hK_over = (dyna.h_K > 1.0).nonzero()
    dyna.h_K[inds_hK_over] = 1.0

    inds_hK_under = (dyna.h_K < 0.0).nonzero()
    dyna.h_K[inds_hK_under] = 0.0

    sim.Dm_vg[sim.iK][dyna.targets_vgK] = (dyna.m_K)*(dyna.h_K)*gK_max

def vgK_1p5(self,sim,dyna,p):

    pass

def vgK_squid(self, sim, dyna,p):

    pass

def vgK_Kir2p1(self,sim,dyna,p):

    pass

def vgPotassium_rat(dyna,sim,cells,p):
    '''
    Handle all **targeted voltage-gated potassium channels** (i.e., only
    applicable to specific tissue profiles) specified by the passed
    user-specified parameters on the passed tissue simulation and cellular
    world for the passed time step.
    '''
     # detecting channels to turn on:
    V = sim.vm[dyna.targets_vgK]*1000

    # m_Inf = 1/(1 + np.exp(-(V+47)/29))
    m_Inf_lambda = define_lambda("lambda V: 1/(1 + np.exp(-(V+47)/29))")
    m_Inf = m_Inf_lambda(V)
    m_Tau = (0.34+0.92*np.exp(-((V+71)/59)**2))
    h_Inf = 1/(1 + np.exp(-(V+56)/-10))
    h_Tau = (8+49*np.exp(-((V+73)/23)**2))

    dyna.m_K = ((m_Inf - dyna.m_K)/m_Tau)*p.dt*1e3 + dyna.m_K
    dyna.h_K = ((h_Inf - dyna.h_K)/h_Tau)*p.dt*1e3 + dyna.h_K

    gK_max = 1.0e-14

    inds_mK_over = (dyna.m_K > 1.0).nonzero()
    dyna.m_K[inds_mK_over] = 1.0

    inds_mK_under = (dyna.m_K < 0.0).nonzero()
    dyna.m_K[inds_mK_under] = 0.0

    inds_hK_over = (dyna.h_K > 1.0).nonzero()
    dyna.h_K[inds_hK_over] = 1.0

    inds_hK_under = (dyna.h_K < 0.0).nonzero()
    dyna.h_K[inds_hK_under] = 0.0

    sim.Dm_vg[sim.iK][dyna.targets_vgK] = (dyna.m_K)*(dyna.h_K)*gK_max

# calcium (and voltage) gated K+ channel-------------------------------------------------------------------------------
#
def cagPotassium(dyna,sim,cells,p):
    """
    Model of a high-conductance calcium activated potassium channel (SK channel), obtained from
    the allosteric model of Cox DH, Cui J, Aldrich RW. J Gen Physiology. 1997. 110: 257-281.

    """

    # get data on cytosolic calcium levels in target cells
    if p.sim_ECM is False:

        ca = sim.cc_cells[sim.iCa][dyna.targets_cagK]

    else:

        ca = sim.cc_cells[sim.iCa][cells.mem_to_cells][dyna.targets_cagK]

    # calculate different terms:
    t1 = 1 + ((ca*1e3)/10.22)
    t2 = 1 + ((ca*1e3)/0.89)

    # calculate probability of channel being open or closed:
    P = 1/(1+((t1/t2)**4)*6182*np.exp(-(1.64*p.F*sim.vm)/(p.R*sim.T)))

    # ensure proper probability behaviour:
    inds_P_over = (P > 1.0).nonzero()
    P[inds_P_over] = 1.0

    inds_P_under = (P < 0.0).nonzero()
    P[inds_P_under] = 0.0

    # calculate conductance of this potassium channel:
    sim.Dm_cag[sim.iK] = dyna.maxDmKcag*P


#----------------------------------------------------------------------------------------------------------------------
# Voltage gated Calcium Channels
#----------------------------------------------------------------------------------------------------------------------

def vgCalcium(dyna,sim,cells,p):

    v = 1e3*sim.vm  # FIXME must be at targets!


    # # Find areas where the differential equation is intrinsically ill-behaved:
    # truth_inds_ha = v < -69
    # truth_inds_hb = v > -71
    #
    # v_inds_h = (truth_inds_ha*truth_inds_hb).nonzero()
    #
    # truth_inds_ma = v < -19
    # truth_inds_mb = v > -21
    #
    # v_inds_m = (truth_inds_ma*truth_inds_mb).nonzero()
    #
    # # small correction constant on the voltage
    # corr_const = 1.0e-6
    #
    # v[v_inds_m] = v + corr_const
    # v[v_inds_h] = v + corr_const

    # model of: K P Carlin et. al; Eur. J. Neurosci. 2000 L-Type------------
    mInf = 1.0/(1+ np.exp(-((v - 10) + 30.0)/6))
    hInf = 1.0/(1+ np.exp(((v-10)+80.0)/6.4))
    mTau = 5.0 + 20.0/(1 + np.exp(((v-10)+25.0)/5))
    hTau = 20.0 + (50.0/(1 + np.exp(((v-10)+40.0)/7)))

    mpower = 2

    # model of: R B Avery et. al; J. Neurosci. 1996 Sep 15 L-Type--------

    # mInf = 1.0000/(1+ np.exp(-(v + 30.000)/6))
    # mTau = 10
    # hInf = 1.0000/(1+ np.exp((v + 80.000)/6.4))
    # hTau = 59

    # T-type calcium channel---------------------------------------
    # mInf = 1 /(1+np.exp((v-(-42.921064))/-5.163208))
    #
    # inds_vn10 = (v < -10).nonzero()
    # inds_vp10 = (v >= -10).nonzero()
    #
    # mTau = np.zeros(len(v))
    #
    # mTau[inds_vn10] = -0.855809 + (1.493527*np.exp(-v[inds_vn10]/27.414182))
    # mTau[inds_vp10] = 1.0
    #
    # # mTau = -0.855809 + (1.493527*np.exp(-v/27.414182))
    #
    # hInf = 1 /(1+np.exp((v-(-72.907420))/4.575763))
    # hTau = 9.987873 + (0.002883 * np.exp(-v/5.598574))
    #
    # mpower = 1

    #-------------------------------------------------------------------

    dyna.m_Ca = ((mInf - dyna.m_Ca)/mTau)*p.dt*1e3 + dyna.m_Ca
    dyna.h_Ca = ((hInf - dyna.h_Ca)/hTau)*p.dt*1e3 + dyna.h_Ca

    # print(dyna.m_Ca.max(),dyna.h_Ca.max())

    # dyna.m_Ca[v_inds_m] = dyna.m_Ca*(1 + corr_const)
    # dyna.m_Ca[v_inds_h] = dyna.m_Ca*(1 + corr_const)
    # dyna.h_Ca[v_inds_m] = dyna.h_Ca*(1 + corr_const)
    # dyna.h_Ca[v_inds_h] = dyna.h_Ca*(1 + corr_const)

    # inds_mCa_over = (dyna.m_Ca > 1.0).nonzero()
    # dyna.m_Ca[inds_mCa_over] = 1.0
    #
    # inds_mCa_under = (dyna.m_Ca < 0.0).nonzero()
    # dyna.m_Ca[inds_mCa_under] = 0.0
    #
    # inds_hCa_over = (dyna.h_Ca > 1.0).nonzero()
    # dyna.h_Ca[inds_hCa_over] = 1.0
    #
    # inds_hCa_under = (dyna.h_Ca < 0.0).nonzero()
    # dyna.h_Ca[inds_hCa_under] = 0.0

    P = (dyna.m_Ca**mpower)*(dyna.h_Ca)

    # print("P",P.max())
    #
    # print('-----')

    sim.Dm_vg[sim.iCa] = dyna.maxDmCa*P

def vgCalcium_init(dyna,sim,p):

    v = 1e3*sim.vm  # FIXME must be at targets!

    # model of: K P Carlin et. al; Eur. J. Neurosci. 2000--------------
    dyna.m_Ca = 1.0/(1+ np.exp(-((v - 10) + 30.0)/6))
    dyna.h_Ca = 1.0/(1+ np.exp(((v-10)+80.0)/6.4))

    # T-type channel:
    # dyna.m_Ca = 1 /(1+np.exp((v-(-42.921064))/-5.163208))
    # dyna.h_Ca = 1 /(1+np.exp((v-(-72.907420))/4.575763))

# defaults--------------------------------------------------------------------------------------------------------------

def vgSodium(dyna,sim,cells,p):
    '''
    Handle all **targeted voltage-gated sodium channels** (i.e., only
    applicable to specific tissue profiles) specified by the passed
    user-specified parameters on the passed tissue simulation and cellular
    world for the passed time step.

    Channel model uses Hodgkin-Huxley model for voltage gated sodium channels.

    '''

    V = sim.vm[dyna.targets_vgNa]*1000

        # add in random noise:
    # corr_const = 1.e-6*np.random.random(len(V))
    # corr_const2 = 1.e-6*np.random.random(len(V))

    #----HH squid model modified (sustains pulses)------------------------

    # V = V + 50
    #
    # mAlpha = (0.1*(25-V))/(np.exp((25-V)/10) -1.0)
    # mBeta = 4.0 * (np.exp(-V/18))
    # mInf = mAlpha/(mAlpha + mBeta)
    # mTau = 1/(mAlpha + mBeta)
    # hAlpha = 0.07 * np.exp(-V/20)
    # hBeta = 1/(np.exp((30-V)/10) + 1.0)
    # hInf = hAlpha/(hAlpha + hBeta)
    # hTau = 1/(hAlpha + hBeta)
    #
    # print(dyna.m_Na.min(),dyna.m_Na.max(),dyna.h_Na.min(),dyna.h_Na.max())


    #--------McCormack Model---------------(very persistant...)-----------------------

    # # Find areas where the differential equation is intrinsically ill-behaved:
    # v_inds_h = []
    #
    # truth_inds_ma = V < -37
    # truth_inds_mb = V > -39
    #
    # v_inds_m = (truth_inds_ma*truth_inds_mb).nonzero()
    #
    # # small correction constant on the voltage
    # corr_const = 1.0e-6
    #
    # V[v_inds_m] = V + corr_const
    # V[v_inds_h] = V + corr_const
    #
    # mAlpha = 0.091*(V+38)/(1-np.exp((-V-38)/5))
    # mBeta = -0.062*(V+38)/(1-np.exp((V+38)/5))
    # hAlpha = 0.016*np.exp((-55-V)/15)
    # hBeta = 2.07/(np.exp((17-V)/21)+1)
    #
    # mInf = mAlpha/(mAlpha + mBeta)
    # mTau = 1/(mAlpha + mBeta)
    #
    # hInf = hAlpha/(hAlpha + hBeta)
    # hTau = 1/(hAlpha + hBeta)

    #------NaS model Hammil 1991------(produces cool, well behaved AP with p.sim_ECM=True. Slow to activate. use!)
    #
    # Find areas where the differential equation is intrinsically ill-behaved:
    truth_inds_ha = V < -39
    truth_inds_hb = V > -41

    v_inds_h = (truth_inds_ha*truth_inds_hb).nonzero()

    truth_inds_ma = V < -24
    truth_inds_mb = V > -26

    v_inds_m = (truth_inds_ma*truth_inds_mb).nonzero()

    # small correction constant on the voltage
    corr_const = 1.0e-10

    V[v_inds_m] = V + corr_const
    V[v_inds_h] = V + corr_const
    #
    mAlpha = (0.182 * ((V-10.0)- -35.0))/(1-(np.exp(-((V-10.0)- -35.0)/9)))
    mBeta = (0.124 * (-(V-10.0) -35.0))/(1-(np.exp(-(-(V-10.0) -35.0)/9)))

    mInf = mAlpha/(mAlpha + mBeta)
    mTau = 1/(mAlpha + mBeta)
    hInf = 1.0/(1+np.exp((V- -65.0-10.0)/6.2))
    hTau = 1/((0.024 * ((V-10.0)- -50.0))/(1-(np.exp(-((V-10.0)- -50.0)/5))) +(0.0091 * (-(V-10.0) -
                                                                75.000123))/(1-(np.exp(-(-(V-10) - 75.000123)/5))))

    #--Neonatal NaV1.3---------------activates at very high depolarizations----
    # # Find areas where the differential equation is intrinsically ill-behaved:
    # v_inds_h = []
    #
    # truth_inds_ma = V < -25
    # truth_inds_mb = V > -27
    #
    # v_inds_m = (truth_inds_ma*truth_inds_mb).nonzero()
    #
    # # small correction constant on the voltage
    # corr_const = 1.0e-6
    #
    # V[v_inds_m] = V + corr_const
    # V[v_inds_h] = V + corr_const
    #
    # mAlpha = (0.182 * ((V)- -26))/(1-(np.exp(-((V)- -26)/9)))
    # mBeta = (0.124 * (-(V) -26))/(1-(np.exp(-(-(V) -26)/9)))
    # mInf = mAlpha/(mAlpha + mBeta)
    # mTau = 1/(mAlpha + mBeta)
    # hInf = 1 /(1+np.exp((V-(-65.0))/8.1))
    # hTau = 0.40 + (0.265 * np.exp(-V/9.47))

    #-----action-potential promising vgNa---------(slow, persistent, weak)
    # qt = 2.3**((34-21)/10)
    # mInf = 1.0/(1+np.exp((V - -52.6)/-4.6))
    # mAlpha = (0.182 * (V- -38))/(1-(np.exp(-(V- -38)/6)))
    # mBeta  = (0.124 * (-V -38))/(1-(np.exp(-(-V -38)/6)))
    # mTau = 6*(1/(mAlpha + mBeta))/qt
    # hInf = 1.0/(1+np.exp((V- -48.8)/10))
    # hAlpha = -2.88e-6 * (V + 17) / (1 - np.exp((V + 17)/4.63))
    # hBeta = 6.94e-6 * (V + 64.4) / (1 - np.exp(-(V + 64.4)/2.63))
    # hTau = (1/(hAlpha + hBeta))/qt

    # NaTa_t-------Costa 2002 -----(fast, strong, persistent)



    # Find areas where the differential equation is intrinsically in need of correction:
    # truth_inds_ma = V < -37
    # truth_inds_mb = V > -39
    #
    # v_inds_m = (truth_inds_ma*truth_inds_mb).nonzero()
    #
    # truth_inds_ha = V < -65
    # truth_inds_hb = V > -67
    #
    # v_inds_h = (truth_inds_ha*truth_inds_hb).nonzero()
    #
    # # small correction constant on the voltage
    # corr_const = 1.0e-6
    #
    # V[v_inds_m] = V + corr_const
    # V[v_inds_h] = V + corr_const

    # qt = 2.3**((34-21)/10)
    #
    # mAlpha = (0.182 * (V- -38))/(1-(np.exp(-(V- -38)/6)))
    # mBeta  = (0.124 * (-V -38))/(1-(np.exp(-(-V -38)/6)))
    # mTau = (1/(mAlpha + mBeta))/qt
    # mInf = mAlpha/(mAlpha + mBeta)
    # hAlpha = (-0.015 * (V- -66))/(1-(np.exp((V- -66)/6)))
    # hBeta  = (-0.015 * (-V -66))/(1-(np.exp((-V -66)/6)))
    # hTau = (1/(hAlpha + hBeta))/qt
    # hInf = hAlpha/(hAlpha + hBeta)

    #-----------------------------------------------------------------------------------

    # calculate m and h channel states:
    dyna.m_Na = ((mInf - dyna.m_Na)/mTau)*p.dt*1e3 + dyna.m_Na
    dyna.h_Na = ((hInf - dyna.h_Na)/hTau)*p.dt*1e3 + dyna.h_Na

    # correction strategy for voltages (to prevent stalling at nullclines)----------
    dyna.m_Na[v_inds_m] = dyna.m_Na*(1 + corr_const)
    dyna.m_Na[v_inds_h] = dyna.m_Na*(1 + corr_const)
    dyna.h_Na[v_inds_m] = dyna.h_Na*(1 + corr_const)
    dyna.h_Na[v_inds_h] = dyna.h_Na*(1 + corr_const)

    #--------------------------------------------------------------------------------

    # threshhold to ensure 0 to 1 status of gates
    # inds_mNa_over = (dyna.m_Na > 1.0).nonzero()
    # dyna.m_Na[inds_mNa_over] = 1.0
    #
    # inds_hNa_over = (dyna.h_Na > 1.0).nonzero()
    # dyna.h_Na[inds_hNa_over] = 1.0
    #
    # inds_mNa_under = (dyna.m_Na < 0.0).nonzero()
    # dyna.m_Na[inds_mNa_under] = 0.0
    #
    # inds_hNa_under = (dyna.h_Na < 0.0).nonzero()
    # dyna.h_Na[inds_hNa_under] = 0.0

    # final probability of open-closed state:

    P = (dyna.m_Na**3)*(dyna.h_Na)

    # Define ultimate activity of the vgNa channel:
    sim.Dm_vg[sim.iNa][dyna.targets_vgNa] = dyna.maxDmNa*P

def vgSodium_init(dyna,sim,p):

    V = sim.vm[dyna.targets_vgNa]*1000

     # Find areas where the differential equation is intrinsically ill-behaved:
    truth_inds_ha = V < -39
    truth_inds_hb = V > -41

    v_inds_h = (truth_inds_ha*truth_inds_hb).nonzero()

    truth_inds_ma = V < -24
    truth_inds_mb = V > -26

    v_inds_m = (truth_inds_ma*truth_inds_mb).nonzero()

    # small correction constant on the voltage
    corr_const = 1.0e-10

    V[v_inds_m] = V + corr_const
    V[v_inds_h] = V + corr_const
    #
    mAlpha = (0.182 * ((V-10.0)- -35.0))/(1-(np.exp(-((V-10.0)- -35.0)/9)))
    mBeta = (0.124 * (-(V-10.0) -35.0))/(1-(np.exp(-(-(V-10.0) -35.0)/9)))

    dyna.m_Na = mAlpha/(mAlpha + mBeta)
    dyna.h_Na = 1.0/(1+np.exp((V- -65.0-10.0)/6.2))


def vgPotassium(dyna,sim,cells,p):
    '''
    Handle all **targeted voltage-gated potassium channels** (i.e., only
    applicable to specific tissue profiles) specified by the passed
    user-specified parameters on the passed tissue simulation and cellular
    world for the passed time step.
    '''
     # detecting channels to turn on:
    V = sim.vm[dyna.targets_vgK]*1000

    # Hodgkin Huxley--------------------------

    # V= V + 15
    #
    # mAlpha = (0.01*(10-V))/(np.exp((10-V)/10) - 1.0)
    # mBeta = 0.125 * (np.exp(-V/80))
    # mInf = mAlpha/(mAlpha + mBeta)
    # mTau = 1/(mAlpha + mBeta)
    #
    # mexp = 4.0


    #dyna.n_K = (alpha_n*(1-dyna.n_K) - beta_n*dyna.n_K)*p.dt*1e3 + dyna.n_K

    # inds_nK_over = (dyna.n_K > 1.0).nonzero()
    # dyna.n_K[inds_nK_over] = 1.0
    #
    # inds_nK_under = (dyna.n_K < 0.0).nonzero()
    # dyna.n_K[inds_nK_under] = 0.0
    #
    # sim.Dm_vg[sim.iK][dyna.targets_vgK] = (dyna.n_K**4)*gK_max
    #
    # mexp = 1

    # Kv1.2-----------------------------------------

    mInf = 1.0000/(1+ np.exp(-(V +21.0000)/11.3943))
    mTau = 150.0000/(1+ np.exp((V + 67.5600)/34.1479))
    hInf = 1.0000/(1+ np.exp((V + 22.0000)/11.3943))
    hTau = 15000.0000/(1+ np.exp(-(V + 46.5600)/44.1479))

    mexp = 1

    # Kv1.5-------------------------------------------------

    # inds_lt50 = (V < 50).nonzero()
    # inds_gt50 = (V>=50).nonzero()
    # inds_lt100 = (V<100).nonzero()
    # inds_gt100 = (V > 100).nonzero()

    # mInf = 1.0000/(1+ np.exp(-(V + 6.0000)/6.4000))
    # mTau = (-0.1163 * V) + 8.3300
    # hInf = 1.0000/(1+ np.exp((V + 25.3000)/3.5000))
    # hTau = (-15.5000 * V) + 1620.0000
    #

    # need to redo this channel with the following:
    #		mInf = 1.0000/(1+ exp((v - -6.0000)/-6.4000))
		# if(v < 50){
		# 	mTau = (-0.1163 * v) + 8.3300
		# }
		# if(v >= 50){
		# 	mTau =  2
		# }
		# hInf = 1.0000/(1+ exp((v - -25.3000)/3.5000))
		# if(v < 100){
		# 	hTau = (-15.5000 * v) + 1620.0000
		# }
		# if(v >= 100){
		# 	hTau =  50
    #
    # mexp = 1


    #- Fast Kv3.3-------------------------------------------
    # mInf = 1/(1+np.exp((V-35)/-7.3))
    # mTau = 0.676808 +( 27.913114 / (1 + np.exp((V - 22.414149)/9.704638)))
    # hInf = 0.25+( 0.75 /(1+ np.exp((V-(-28.293856))/29.385636)))
    # hTau = 199.786728 + (2776.119438*np.exp(-V/7.309565))
    # #
    # mexp = 1

    #--K slow rat-------------------------------------------------------

    # mInf = (1/(1 + np.exp(-(V+14)/14.6)))
    # inds_lt50 = (V < -50).nonzero()
    # inds_gt50 = (V >= -50).nonzero()
    # mTau = np.zeros(len(V))
    # mTau[inds_lt50] = (1.25+175.03*np.exp(V*0.026))
    # mTau[inds_gt50] = (1.25+13*np.exp(-V*0.026))
    # hInf = 1/(1 + np.exp((V+54)/11))
    # hTau = 360+(1010+24*(V+55))*np.exp(-((V+75)/48)**2)
    #
    # mexp = 2

    #-------------------------------#
    dyna.m_K = ((mInf - dyna.m_K)/mTau)*p.dt*1e3 + dyna.m_K
    dyna.h_K = ((hInf - dyna.h_K)/hTau)*p.dt*1e3 + dyna.h_K
    # dyna.h_K[:] = 1

    # inds_mK_over = (dyna.m_K > 1.0).nonzero()
    # dyna.m_K[inds_mK_over] = 1.0
    #
    # inds_mK_under = (dyna.m_K < 0.0).nonzero()
    # dyna.m_K[inds_mK_under] = 0.0
    #
    # inds_hK_over = (dyna.h_K > 1.0).nonzero()
    # dyna.h_K[inds_hK_over] = 1.0
    #
    # inds_hK_under = (dyna.h_K < 0.0).nonzero()
    # dyna.h_K[inds_hK_under] = 0.0

    # open probability
    P =  (dyna.m_K**mexp)*(dyna.h_K)

    # correction factors to push out of stall points:
    # v_inds5 = (np.round(V,0) == -5).nonzero()
    #
    # P[v_inds5] = P - 1.0e-6
    #
    # inds_P_over = (P > 1.0).nonzero()
    # P[inds_P_over] = 1.0
    #
    # inds_P_under = (P < 0.0).nonzero()
    # P[inds_P_under] = 0.0

    sim.Dm_vg[sim.iK][dyna.targets_vgK] = P*dyna.maxDmK


def vgPotassium_init(dyna,sim,p):

     # detecting channels to turn on:
    V = sim.vm[dyna.targets_vgK]*1000

    dyna.m_K = 1.0000/(1+ np.exp(-(V +21.0000)/11.3943))
    dyna.h_K = 1.0000/(1+ np.exp((V + 22.0000)/11.3943))












#FIXME: Shift into a "betse.util" module. Hefty superstrings in the cleft!
def define_lambda(lambda_func: str, arg: object) -> 'lambda':
    '''
    Dynamically define and return a callable lambda function object defined by
    the passed string.

    The passed string _must_ conform to Python syntax for lambda functions. In
    particular, this string _must_ be prefixed by the identifier `lambda `
    (e.g., `lambda x: x**2`).

    Parameters
    ----------------------------
    lambda_func : str
        String defining the lambda function to be created and returned.

    Returns
    ----------------------------
    lambda
        Callable lambda function object defined by the passed string.
    '''
    assert types.is_str_nonempty(lambda_func), (
        types.assert_not_str_nonempty(lambda_func, 'Lambda function'))

    # If this is not a lambda function, raise an exception.
    if not lambda_func.startswith('lambda '):
        raise BetseExceptionLambda(
            'Lambda function not preceded by "lambda ":\n{}'.format(
                lambda_func))

    # Define and return this lambda function.
    exec('return {}'.format(lambda_func))


# Wastelands
#------------------------

    # m_Inf = call_lambda("lambda V: 1/(1 + np.exp(-(V+47)/29))", V)

# def call_lambda(lambda_func: str, arg: object) -> object:
#     '''
#     Dynamically call the lambda function defined by the passed string passed the
#     passed argument _and_ return the value returned by this call.
#
#     Parameters
#     ----------------------------
#     lambda_func : str
#         String defining the lambda function to be called.
#     arg : object
#         First and only argument to be passed to this lambda function.
#     '''
#     assert types.is_str_nonempty(lambda_func), (
#         types.assert_not_str_nonempty(lambda_func, 'Lambda function'))
#
#     # If this is not a lambda function, raise an exception.
#     if not lambda_func.startswith('lambda '):
#         raise BetseExceptionLambda(
#             'Lambda function not preceded by "lambda ":\n{}'.format(
#                 lambda_func))
#
#     # Call this lambda function and return the return value.
#     exec('return ({})(arg)'.format(lambda_func))

#def sigmoid(V,po,p1,p2,p3):
#     """
#      Template for typical channel voltage-dependent parameter function m_inf or tau_inf,
#      used in Hodgkin-Huxley channel models.
#
#      Parameters:
#      ------------
#      Input parameters describe the equation:
#
#      po + p1/(1 + exp((V - p2)/p3)
#
#      Returns
#      --------
#      f          The value of the sigmoid function
#     """
#
#     f = po + (p1/(1 + np.exp((V-p2)/p3)))
#
#     return f
#
# def alpha(V,po,p1,p2,p3):
#
#     f = (po*(p1-V))/(np.exp((p2-V)/p3)-1)
#
#     return f
#
# def beta(V, po, p1):
#
#     f = po*np.exp(-V/p1)





