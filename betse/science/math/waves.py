#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.


import math

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation


# FIXME to integrate with BETSE, needs to be able to take in forcing function data (which comes from BETSE) and may
# need to be interpolated to this particular grid (along with 'c', by the way, which might be 0 outside of cell cluster)

# FIXME thinking about implementing the FitzHugh-Nagumo or van der Poll relaxation oscillator math as as simple model
# to demonstrate how phase oscillators might create Pilot Wave fields that go on to control the phase oscillators.
# Would want to compare to fMRI images of symmetric brain activity zones.

class WaveSolver(object):

    """
    WaveSolver provides a suite of tools to initialize and iterate a damped, forced wave function on a rectangular
    grid using the Finite Difference Time Domain method.

    Methods
    ---------



    Attributes
    ----------

    """

    def __init__(self,p,constant_c = True):

        # define parameters (these will eventually be linked to the parameters object p)
        self.c = 1     # wave speed [m/s] (can be an array of values defined over the rectangular grid)

        self.xmin = 0  # bounding dimensions of gridded space [m]
        self.xmax = 5
        self.ymin = 0
        self.ymax = 5

        self.extent_p = [self.xmin,self.xmax,self.ymin,self.ymax]  # list holding the world extents

        self.grid_xy_n = 150  # grid points in the x and y dimensions
        self.grid_t_n = 2500  # number of time points

        self.gamma = 1.0e-3  # damping factor [N s/m]

        self.delta_xy = (self.xmax/self.grid_xy_n)   # spacing in the x-y grid

        self.t_end = 20  # end time for simulation

        self.dt = (self.t_end/self.grid_t_n)  # spacing in time grid

        self.x_grid = np.linspace(self.xmin,self.xmax,self.grid_xy_n)
        self.y_grid = np.linspace(self.ymin,self.ymax,self.grid_xy_n)

        self.X, self.Y = np.meshgrid(self.x_grid,self.y_grid)  # gridded X and Y dimensions [m]

        self.fxy = np.exp(-10*((self.X-(self.xmax/2))**2 + (self.Y-(self.ymax/2))**2))  # forcing function amplitude

        f = 0.75 # forcing function frequency

        self.forcing_function = lambda t: np.cos(2*math.pi*f*t)

        self.K = (self.dt*self.c)/self.delta_xy  # Courant number (magnitude relates to stability)

        if self.K > 0.5:

            print('Warning! Your simulation is indicated to be unstable. '
                  'Increase number of time steps or decrease wave speed.')

        self.mShape = list(self.X.shape)
        self.mShape.append(self.grid_t_n)
        self.mShape = tuple(self.mShape)

        self.bounds = 'reflective'   # state what boundary condition applies to the box boundaries

        self.colormap = cm.coolwarm

    def init_matrices(self):  # FIXME next version of this must be able to take in an arbitrary F property to force U
        '''
        Initialize the main matrices of the simulation. Assumes that the initial intensity and velocity of the
        wave property U are zero, and that the substratum is influenced by a forcing function (or property) F.

        '''

        # initial conditions and first time step -- forced damped wave equation
        self.U = np.zeros(self.mShape)   # initialize the U matrix to hold all time steps

        self.F = np.zeros(self.mShape)   # initialize the forcing function F matrix to hold all time steps

        self.F[:,:,0] = self.fxy*self.forcing_function(0)  # create the first time-point of the F matrix

        self.U[:,:,0] = np.zeros(self.X.shape)    # the initial state of the U-property will be zero

        # the next state of the U-property will depend on the forcing function state

        self.U[:,:,1] = ((self.dt**2)/(2-self.dt*self.gamma))*self.F[:,:,0]

    def iter_wave(self,k,t):

        """
        Iterates the wave function under the forcing condition F for time t.

        Parameters
        -----------
        k           Index into the time vector
        t           Time value at index k

        """

        # Advance the forcing function (FIXME later versions will accept forcing *data*)

        self.F[:,:,k] = self.fxy*self.forcing_function(t)

        UU = self.U*self.c  # if c varies over space, this will allow the function to take account for that

        # Interior points:

        self.U[1:-1,1:-1,k+1] = (1/(1+self.gamma))*(self.gamma*self.U[1:-1,1:-1,k]+
                            ((self.dt/self.delta_xy)**2)*(UU[2:,1:-1,k] +
                            UU[0:-2,1:-1,k] + UU[1:-1,2:,k] + UU[1:-1,0:-2,k] - 4*UU[1:-1,1:-1,k]) +
                            2*self.U[1:-1,1:-1,k] - self.U[1:-1,1:-1,k-1] + (self.dt**2)*self.F[1:-1,1:-1,k])


        # Boundary conditions:

        # Dirchlet (fixed boundary, say 0)
        if self.bounds == 'fixed':

            self.U[:,0,k+1] = 0
            self.U[:,-1,k+1] = 0
            self.U[0,:,k+1] = 0
            self.U[-1,:,k+1] = 0


        # reflective boundary
        if self.bounds == 'reflective':

            aa = self.dt/self.delta_xy

            # left boundary, j = 0:
            self.U[1:-1,0,k+1] = (aa**2)*(self.UU[2:,0,k] + self.UU[0:-2,0,k] + 2*self.UU[1:-1,1,k] -
                                              4*self.UU[1:-1,0,k]) + 2*self.U[1:-1,0,k] - self.U[1:-1,0,k-1]

            # right boundary, j = jmax:
            self.U[1:-1,-1,k+1] = (aa**2)*(self.UU[2:,-1,k] + self.UU[0:-2,-1,k] + 2*self.UU[1:-1,-2,k] -
                                               4*self.UU[1:-1,-1,k]) + 2*self.U[1:-1,-1,k] - self.U[1:-1,-1,k-1]

            # bottom boundary, i = 0:
            self.U[0,1:-1,k+1] = (aa**2)*(2*self.UU[1,1:-1,k] + self.UU[0,2:,k] + self.UU[0,0:-2,k] -
                                              4*self.UU[0,1:-1,k]) + 2*self.U[0,1:-1,k] - self.U[0,1:-1,k-1]

            # top boundary, i = max:
            self.U[-1,1:-1,k+1] = (aa**2)*(2*self.UU[-2,1:-1,k] + self.UU[-1,2:,k] + self.UU[-1,0:-2,k] -
                                               4*self.UU[-1,1:-1,k]) + 2*self.U[-1,1:-1,k] - self.U[-1,1:-1,k-1]

        if self.bounds == 'open':

            # left boundary, j = 0:
            self.U[1:-1,0,k+1] = self.U[1:-1,1,k] + ((self.K-1)/(self.K+1))*(self.U[1:-1,1,k+1] - self.U[1:-1,0,k])

            # right boundary, j = jmax:
            self.U[1:-1,-1,k+1] = self.U[1:-1,-2,k] + ((1-self.K)/(1+self.K))*(self.U[1:-1,-1,k] - self.U[1:-1,-2,k+1])


            # bottom boundary, i = 0:
            self.U[0,1:-1,k+1] =  self.U[1,1:-1,k] + ((self.K-1)/(self.K+1))*(self.U[1,1:-1,k+1] - self.U[0,1:-1,k])

            # top boundary, i = max:
            self.U[-1,1:-1,k+1] = self.U[-2,1:-1,k] + ((1-self.K)/(1+self.K))*(self.U[-1,1:-1,k] - self.U[-2,1:-1,k+1])


class AnimateWaves(object):

    def __init__(self,sim, autoscale=True,min_val = -0.50,max_val=0.5):

        fig = plt.figure()

        data = sim.U

        if autoscale is True:

            data_max = np.max(sim.U)
            data_min = np.min(sim.U)

        else:
            data_max = max_val
            data_min = min_val

        l = plt.imshow(sim.U[:,:,0],origin='lower',extent=sim.extent_p,vmin=data_min,vmax=data_max,cmap=sim.colormap)
        plt.colorbar()

        plt.xlim(sim.xmin, sim.xmax)
        plt.ylim(sim.ymin, sim.ymax)
        plt.xlabel('x')
        plt.title('test')
        line_ani = animation.FuncAnimation(fig, self.update_line, sim.grid_t_n, fargs=(data, l),
            interval=50, blit=True)

    def update_line(num, data, line):
        line.set_data(data[:,:,num])
        return line

