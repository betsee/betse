#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from betse.science import filehandling as fh
from betse.science.math.modulate import gradient_x, gradient_xr, gradient_bitmap, double_x

class WernerSim(object):

    def __init__(self, fn):

        self.sim, self.cells, self.p = fh.loadSim(fn)

    def init_sim(self):
        # define parameters:
        h = 5

        # constants for activation/inhibition of substance A:
        alpha_A = 0.1
        k_A = 1.0e-3
        D_A = 0.5e-10

        # constants for activation/inhibition of substance B:
        alpha_B = 0.4
        k_B = 2.0e-3
        D_B = 1.5e-9

        # constants for substance E:
        k_E = 2 * k_A
        D_E = 5.0e-10
        alpha_E = 0.4 * alpha_A

        # end time sequence:
        endt = 3600
        dt = 1.0e-2

        tpoints = int(endt / dt)

        #
        ttt = np.linspace(0, endt, tpoints)

        tsample = ttt[0:-1:1000]

        gx, _ = gradient_bitmap(cells.cell_centres, cells, p)

        # initial conditions for cA and cB
        # cA = np.ones(sim.cdl)*(1 + 0.25*(cells.cell_centres[:,0]/cells.cell_centres[:,0].max()))
        cA = np.ones(sim.cdl) * gx + 0.1
        cB = np.zeros(sim.cdl)
        cE = 0.1 * np.ones(sim.cdl)

    def sim_loop(self):

        for tt in ttt:

            dAt = deltaA(cA, cB, cE, D_A, alpha_A, k_A, h)
            cA += dAt * dt

            dBt = deltaB(cA, cB, cE, D_B, alpha_B, k_B, h)
            cB += dBt * dt

            dEt = deltaE(cB, D_E, alpha_E, k_E)
            cE += dEt * dt

            if tt in tsample:
                ii += 1
                print(ii, ' of ', len(tsample))
                t_time.append(1 * tt)
                cA_time.append(1 * cA)
                cB_time.append(1 * cB)
                cE_time.append(1 * cE)




    def deltaA(cA, cB, cE, D_A, alpha_A, k_A, h):
        """
        Function calculating the change for substance A in the reaction diffusion system.
        """

        # calculate the activation-inhibition term:
        termAB = (cA ** h) / (cA ** h + cB ** h)

        # calculate the extender concentration decay term:
        beta_A = k_A * cE

        # calculate the gradient of concentration A:
        gcA = (cA[cells.cell_nn_i[:, 1]] - cA[cells.cell_nn_i[:, 0]]) / (cells.nn_len)

        # calculate the diffusive flux of A:
        fluxA = -D_A * gcA

        # enforce zero flux at outer boundary:
        fluxA[cells.bflags_mems] = 0.0

        # take the divergence of the flux to obtain the net change with time:
        div_fluxA = np.dot(cells.M_sum_mems, -fluxA * cells.mem_sa) / cells.cell_vol

        # calculate the change with time for the full reaction-diffusion expression:
        dAt = alpha_A * termAB - beta_A * cA + div_fluxA

        return dAt

    def deltaB(cA, cB, cE, D_B, alpha_B, k_B, h):
        """
        Function calculating the change for substance B in the reaction diffusion system.
        """

        # calculate the activation-inhibition term:
        termAB = (cA ** h) / (cA ** h + cB ** h)

        # calculate the extender concentration decay term:
        beta_B = k_B * cE

        # calculate the gradient of concentration B:
        gcB = (cB[cells.cell_nn_i[:, 1]] - cB[cells.cell_nn_i[:, 0]]) / (cells.nn_len)

        # calculate the diffusive flux of B:
        fluxB = -D_B * gcB

        # enforce zero flux at outer boundary:
        fluxB[cells.bflags_mems] = 0.0

        # take the divergence of the flux to obtain the net change with time:
        div_fluxB = np.dot(cells.M_sum_mems, -fluxB * cells.mem_sa) / cells.cell_vol

        # calculate the change with time for the full reaction-diffusion expression:
        dBt = alpha_B * termAB - beta_B * cB + div_fluxB

        return dBt

    def deltaE(cB, D_E, alpha_E, k_E):
        """
        Function calculating the change for substance E in the reaction diffusion system.
        """

        # calculate the gradient of concentration E:
        gcE = (cE[cells.cell_nn_i[:, 1]] - cE[cells.cell_nn_i[:, 0]]) / (cells.nn_len)

        # calculate the diffusive flux of E:
        fluxE = -D_E * gcE

        # enforce zero flux at outer boundary:
        fluxE[cells.bflags_mems] = 0.0

        # take the divergence of the flux to obtain the net change with time:
        div_fluxE = np.dot(cells.M_sum_mems, -fluxE * cells.mem_sa) / cells.cell_vol

        # calculate the change with time for the full reaction-diffusion expression:
        dEt = alpha_E - k_E * cB * cE + div_fluxE

        return dEt


    def init_plots(self):

        xc = p.um * cells.cell_centres[:, 0]
        yc = p.um * cells.cell_centres[:, 1]

        xm = p.um * cells.mem_mids_flat[:, 0]
        ym = p.um * cells.mem_mids_flat[:, 1]

        xec = p.um * cells.ecm_mids[:, 0]
        yec = p.um * cells.ecm_mids[:, 1]

        xenv = p.um * cells.xypts[:, 0]
        yenv = p.um * cells.xypts[:, 1]

        xyaxis = [p.um * cells.xmin, p.um * cells.xmax, p.um * cells.ymin, p.um * cells.ymax]

        verts = p.um * cells.cell_verts

        fig = plt.figure()
        ax = plt.subplot(111)
        col = PolyCollection(verts, cmap='RdBu_r')
        col.set_array(cA_time[-1])
        ax.add_collection(col)
        fig.colorbar(col)
        plt.axis('equal')
        plt.axis(xyaxis)
        plt.show()

    def save(self):

        MMi = np.asarray([cA_time2, cB_time2, cE_time2])
        np.save('Werner_Nov/Werner_1H_2EC.npy', MMi)

