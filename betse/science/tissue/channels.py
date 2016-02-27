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




def vgSodium(dyna,sim,cells,p):
    '''
    Handle all **targeted voltage-gated sodium channels** (i.e., only
    applicable to specific tissue profiles) specified by the passed
    user-specified parameters on the passed tissue simulation and cellular
    world for the passed time step.

    Channel model uses Hodgkin-Huxley model for voltage gated sodium channels.

    '''

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


def vgPotassium(dyna,sim,cells,p):
    '''
    Handle all **targeted voltage-gated potassium channels** (i.e., only
    applicable to specific tissue profiles) specified by the passed
    user-specified parameters on the passed tissue simulation and cellular
    world for the passed time step.
    '''
     # detecting channels to turn on:
    V = sim.vm[dyna.targets_vgK]*1000 + 20.0

    alpha_n = (0.01*(10 - V))/(np.exp((10-V)/10)-1)
    beta_n = 0.125*np.exp(-V/80)

    dyna.n_K = (alpha_n*(1-dyna.n_K) - beta_n*dyna.n_K)*p.dt*1e3 + dyna.n_K

    gK_max = 2.0e-14    # gK_max from 1/5 of gNa_max to 1/10th produces "bistable" state

    inds_nK_over = (dyna.n_K > 1.0).nonzero()
    dyna.n_K[inds_nK_over] = 1.0

    inds_nK_under = (dyna.n_K < 0.0).nonzero()
    dyna.n_K[inds_nK_under] = 0.0

    sim.Dm_vg[sim.iK][dyna.targets_vgK] = (dyna.n_K**4)*gK_max




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





