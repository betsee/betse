#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level facilities for displaying and/or saving all enabled plots _and_
animations.
'''

#FIXME: For safety, most "== 1"-style tests in this module should be converted
#to "is True"-style tests instead. Into the trackless reaches of ice and snow!
#FIXME: I believe I've finally tracked down the issue relating to the following
#runtime "pyplot" warning:
#
#    pyplot.py:424: RuntimeWarning: More than 20 figures have been opened. Figures created through the pyplot interface (`matplotlib.pyplot.figure`) are retained until explicitly closed and may consume too much memory. (To control this warning, see the rcParam `figure.max_open_warning`).
#
#The issue is that whenever we call a viz.plot* function below (e.g.,
#viz.plotSingleCellCData()), we localize the figure returned by that function.
#That figure will then be garbage collected on whichever of the following
#occurs last: (A) the corresponding local variable goes out of scope and (B)
#the corresponding plot window is closed by the user. Normally, neither would
#be a problem. Except this function is 1,400 lines long, which means that each
#figure's local variable effectively *NEVER* goes out of scope for the duration
#of plotting. Thus, figures will only be garbage collected *AFTER* this
#function terminates -- which is pretty much unacceptable.
#
#There are a couple solutions, thankfully. The simplest would simply be to stop
#localizing figures returned by viz.plot*() functions for all unused figure
#locals. The harder but probably more ideal solution would be to refactor all
#viz.plot*() functions to stop returning figures altogether. Since the current
#figure is *ALWAYS* accessible via the "matplotlib.pyplot.gcf()" getter (i.e.,
#[g]et[c]urrent[f]igure), there's no nead to explicitly return figures at all.
#
#In the unlikely event that we actually use a figure local, we need to either:
#
#* Extract all plotting logic related to that figure into a new helper function
#  of this module, ensuring that local will go out of scope immediately after
#  plotting that figure in that function.
#* Manually destroy that local with either "figure = None" or "del figure".
#
#Undomesticated unicorns running into the carefree sunset!

# ....................{ IMPORTS                            }....................
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib.collections import LineCollection
from betse.exceptions import BetseSimConfigException
from betse.science.visual.plot import plotutil as viz
from betse.science.visual.anim.animpipe import pipeline_anims
from betse.util.io.log import logs
from betse.util.type import types
from scipy.ndimage.filters import gaussian_filter

# ....................{ PIPELINES                          }....................
#FIXME: Refactor the "plot_type" parameter to be somewhat less crazy. This
#parameter is passed by the "simrunner" submodule. Ideally, this parameter
#should be a typesafe enum rather than a non-typesafe string or, preferably,
#simply go away entirely. Related commentary follows on how to best achieve the
#latter goal.
#FIXME: I don't quite grok our usage of "sim.run_sim". This undocumented
#attribute appears to be internally set by the Simulator.run_phase_sans_ecm()
#method. That makes sense; however, what's the parallel "p.run_sim" attribute
#for, then?  Interestingly, the "SimRunner" class sets "p.run_sim" as follows:
#
#* To "False" if an initialization is being performed.
#* To "True" if a simulation is being performed.
#
#This doesn't seem quite ideal, however. Ideally, there would exist one and only
#one attribute whose value is an instance of a multi-state "SimPhaseType" class
#rather than two binary boolean attributes. Possible enum values might include:
#
#* "SimPhaseType.seed" when seeding a new cluster.
#* "SimPhaseType.init" when initializing a seeded cluster.
#* "SimPhaseType.sim" when simulating an initialized cluster.
#
#This attribute would probably exist in the "Simulator" class -- say, as
#"sim.phase". In light of that, consider the following refactoring:
#
#* Define a new "SimPhaseType" class in the "sim" module with the above
#  attributes.
#* Define a new "Simulator.phase" attribute initialized to None.
#* Replace all existing uses of the "p.run_sim" and "sim.run_sim" booleans with
#  "sim.phase" instead. Note that only the:
#  * "SimRunner" class sets "p.run_sim".
#  * "Simulator" class sets "sim.run_sim".
#
#Note also the:
#
#* "plot_type" parameter passed to the pipeline_plots() function by the
#  "SimRunner" class.
#* The seemingly duplicate "p.plot_type" attribute internally set by the
#  pipeline_plots() function, which is frankly crazy.
#
#Both parameters should probably receive similar treatment and be replaced
#entirely by use of the new "sim.phase" attribute.
#
#Wonder temptress at the speed of light and the sound of love!

#FIXME: Shift this function into a new "betse.science.visual.visualpipe"
#submodule.

def pipeline_results(
    sim: 'Simulator',
    cells: 'Cells',
    p: 'Parameters',
    plot_type: str = 'init',
) -> None:
    '''
    Serially (i.e., in series) display and/or save all enabled plots and
    animations for the passed simulation phase (e.g., `init`, `sim`).

    Parameters
    ----------------------------
    sim : Simulator
        Current simulation.
    cells : Cells
        Current cell cluster.
    p : Parameters
        Current simulation configuration.
    plot_type : str
        String constant corresponding to the current simulation phase. Valid
        values include:
        * `init`, for plotting simulation initialization results.
        * `sim`, for plotting simulation run results.
    '''
    assert types.is_simulator(sim), types.assert_not_simulator(sim)
    assert types.is_cells(cells), types.assert_not_parameters(cells)
    assert types.is_parameters(p), types.assert_not_parameters(p)

    #FIXME: This is terrible. I don't even.
    p.plot_type = plot_type

    # Display and/or save all plots.
    pipeline_plots(sim, cells, p)

    # Display and/or save all animations.
    pipeline_anims(sim, cells, p)

    #FIXME: What is this? What requires showing? Are we finalizing some
    #previously displayed visual artifact? We suspect this to be safely
    #jettisoned deadweight, but... let's verify that, please.

    # If displaying plots and animations, display... something? I guess?
    if p.turn_all_plots_off is False:
        plt.show()
    # Else, log the directory to which results were exported.
    else:
        #FIXME: This is terrible. I don't even. For one, this logic is
        #duplicated below by the pipeline_plots() function. For another, the
        #"p.sim_results" and "p.sim_results" parameters should probably simply
        #be aggregated into a single "sim.export_dirname" parameter
        #corresponding to the export directory for the current phase.
        if p.plot_type == 'sim':
            export_dirname = p.sim_results
        elif p.plot_type == 'init':
            export_dirname = p.init_results

        logs.log_info('Results exported to: %s', export_dirname)

# ....................{ PIPELINES ~ plots                  }....................
def pipeline_plots(
    sim: 'Simulator',
    cells: 'Cells',
    p: 'Parameters',
) -> None:
    '''
    Serially (i.e., in series) display and/or save all enabled plots for the
    current simulation phase if animations are enabled _or_ noop otherwise.

    Parameters
    ----------------------------
    sim : Simulator
        Current simulation.
    cells : Cells
        Current cell cluster.
    p : Parameters
        Current simulation configuration.
    '''
    assert types.is_simulator(sim), types.assert_not_simulator(sim)
    assert types.is_cells(cells), types.assert_not_parameters(cells)
    assert types.is_parameters(p), types.assert_not_parameters(p)

    # If post-simulation plots are disabled, noop.
    if not p.plot.is_after_sim:
       return

    if p.autosave is True:
        if p.plot_type == 'sim':
            images_path = p.sim_results
        elif p.plot_type == 'init':
            images_path = p.init_results

        image_cache_dir = os.path.expanduser(images_path)
        os.makedirs(image_cache_dir, exist_ok=True)
        savedImg = os.path.join(image_cache_dir, 'fig_')

    # check that the plot cell is in range of the available cell indices:
    if p.plot_cell not in cells.cell_i:
        raise BetseSimConfigException(
            'The "plot cell" defined in the "results" section of your '
            'configuration file does not exist in your cluster. '
            'Choose a plot cell number smaller than the maximum cell number.')

    if p.sim_ECM is True:
        plot_cell_ecm = cells.cell_to_mems[p.plot_cell][0]  # convert from cell to mem index
    else:
        plot_cell_ecm = p.plot_cell

    if p.exportData is True:
        viz.exportData(cells, sim, p)

    if p.exportData2D is True:
        for i, t in enumerate(sim.time):
            simdata = 1.0e3*sim.vm_ave_time[i]
            viz.export2dData(i, simdata, cells, p)

        # for i, t in enumerate(sim.time):
        #     simdata_x = 1.0e3*sim.pol_x_time[i]
        #     viz.export2dData(i, simdata_x, cells, p, foldername = 'Polarization_x', filebit = 'Pol_x')
        #
        #     simdata_y = 1.0e3 * sim.pol_y_time[i]
        #     viz.export2dData(i, simdata_y, cells, p, foldername='Polarization_y', filebit='Pol_y')

    #-------------------------------------------------------------------------------------------------------------------
    #               SINGLE CELL DATA GRAPHS
    #-------------------------------------------------------------------------------------------------------------------

    if p.plot_single_cell_graphs is True:
        # Plot cell sodium concentration versus time.

        mem_i = cells.cell_to_mems[p.plot_cell][0]

        figConcsNa, axConcsNa = viz.plotSingleCellCData(
            sim.cc_time, sim.time, sim.iNa, p.plot_cell, fig=None,
            ax=None, lncolor='g', ionname='Na+')

        titNa = 'Sodium concentration in cell ' + str(p.plot_cell)
        axConcsNa.set_title(titNa)

        if p.autosave is True:
            savename1 = savedImg + 'concNa_time' + '.png'
            plt.savefig(savename1, dpi=300, format='png', transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # Plot cell potassium concentration versus time.
        figConcsK, axConcsK = viz.plotSingleCellCData(
            sim.cc_time, sim.time, sim.iK, p.plot_cell, fig=None,
            ax=None, lncolor='b', ionname='K+')

        titK = 'Potassium concentration in cell ' + str(p.plot_cell)
        axConcsK.set_title(titK)

        if p.autosave is True:
            savename1 = savedImg + 'concK_time' + '.png'
            plt.savefig(savename1,dpi=300,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # plot-cell anion (bicarbonate) concentration vs time:

        figConcsM, axConcsM = viz.plotSingleCellCData(
            sim.cc_time, sim.time, sim.iM, p.plot_cell,
            fig=None,
            ax=None,
            lncolor='r',
            ionname='M-',
        )

        titM = 'M Anion concentration in cell ' + str(p.plot_cell)
        axConcsM.set_title(titM)

        if p.autosave is True:
            savename1 = savedImg + 'concM_time' + '.png'
            plt.savefig(savename1,dpi=300,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # Plot single cell Vmem vs time.
        figVt, axVt = viz.plotSingleCellVData(
            sim, mem_i, p, fig=None, ax=None, lncolor='k')
        titV = 'Voltage (Vmem) in cell ' + str(p.plot_cell)
        axVt.set_title(titV)

        if p.autosave is True:
            savename2 = savedImg + 'Vmem_time' + '.png'
            plt.savefig(savename2,dpi=300,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # Plot fast-Fourier-transform (fft) of Vmem.
        figFFT, axFFT = viz.plotFFT(
            sim.time, sim.vm_time, mem_i, lab="Power")
        titFFT = 'Fourier transform of Vmem in cell ' + str(p.plot_cell)
        axFFT.set_title(titFFT)

        if p.autosave is True:
            savename = savedImg + 'FFT_time' + '.png'
            plt.savefig(savename,dpi=300,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # plot rate of Na-K-ATPase pump vs time:
        figNaK = plt.figure()
        axNaK = plt.subplot()

        if p.sim_ECM is False:
            pump_rate = [pump_array[p.plot_cell] for pump_array in sim.rate_NaKATP_time]

        else:
            pump_rate = [pump_array[plot_cell_ecm] for pump_array in sim.rate_NaKATP_time]

        axNaK.plot(sim.time, pump_rate)
        axNaK.set_xlabel('Time [s]')
        axNaK.set_ylabel('Pumping Rate of Na+ Out of Cell [mol/(m2 s)]')
        axNaK.set_title('Rate of NaK-ATPase pump in cell: ' + str(p.plot_cell) )

        if p.autosave is True:
            savename = savedImg + 'NaKATPaseRaTE_' + '.png'
            plt.savefig(savename,dpi=300,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        #--------------------------------------------------------

        # Plot cell trans-membrane current vs time.
        figI = plt.figure()
        axI = plt.subplot(111)

        if p.sim_ECM is False:
            Imem = [100*memArray[p.plot_cell] for memArray in sim.I_mem_time]

        else:
            Imem = []  # initialize a total cell current storage vector
            mems_for_plotcell = cells.cell_to_mems[p.plot_cell]  # get membranes for the plot cell

            for t in range(len(sim.time)):
                memArray = sim.I_mem_time[t]
                # get the current components at each membrane (net current, not density)
                Ixo = memArray[mems_for_plotcell]*cells.mem_vects_flat[mems_for_plotcell,2]*cells.mem_sa[mems_for_plotcell]
                Iyo = memArray[mems_for_plotcell]*cells.mem_vects_flat[mems_for_plotcell,3]*cells.mem_sa[mems_for_plotcell]
                # add components of current at each membrane (this takes account for in-out directionality)
                Ix = np.sum(Ixo)
                Iy = np.sum(Iyo)

                # get the total magnitude of net current and divide by cell surface area to return to density:
                Io = np.sqrt(Ix**2 + Iy**2)/cells.cell_sa[p.plot_cell]
                Imem.append(100*Io)

        axI.plot(sim.time, Imem)
        axI.set_title('Transmembrane current density for cell ' + str(p.plot_cell) )
        axI.set_xlabel('Time [s]')
        axI.set_ylabel('Current density [uA/cm2]')

        if p.autosave is True:
            savename = savedImg + 'Imem_time' + '.png'
            plt.savefig(savename,dpi=300,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        # optional 1D plots--------------------------------------------------------------------------------------------

        # hydrostatic pressure in cells:

        if np.mean(sim.P_cells_time) != 0.0:
            p_hydro = [arr[p.plot_cell] for arr in sim.P_cells_time]
            figOP = plt.figure()
            axOP = plt.subplot(111)
            axOP.plot(sim.time, p_hydro)
            axOP.set_xlabel('Time [s]')
            axOP.set_ylabel('Hydrostatic Pressure [Pa]')
            axOP.set_title('Hydrostatic pressure in cell ' + str(p.plot_cell) )

            if p.autosave is True:
                savename = savedImg + 'HydrostaticP_' + '.png'
                plt.savefig(savename,dpi=300,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

        # Plot cell calcium vs time (if Ca enabled in ion profiles).
        if p.ions_dict['Ca'] ==1:
            figA, axA = viz.plotSingleCellCData(
                sim.cc_time, sim.time, sim.iCa, p.plot_cell,
                fig=None,
                ax=None,
                lncolor='g',
                ionname='Ca2+ cell',
            )
            titCa = 'Cytosolic Ca2+ in cell index ' + str(p.plot_cell)
            axA.set_title(titCa)

            if p.autosave is True:
                savename3 = savedImg + 'cytosol_Ca_time' + '.png'
                plt.savefig(savename3,dpi=300,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)


        if p.deform_osmo is True:
            p_osmo = [arr[p.plot_cell] for arr in sim.osmo_P_delta_time]
            figOP = plt.figure()
            axOP = plt.subplot(111)
            axOP.plot(sim.time, p_osmo)
            axOP.set_xlabel('Time [s]')
            axOP.set_ylabel('Osmotic Pressure [Pa]')
            axOP.set_title('Osmotic pressure in cell ' + str(p.plot_cell) )

            if p.autosave is True:
                savename = savedImg + 'OsmoticP_' + '.png'
                plt.savefig(savename,dpi=300,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

        # Total displacement in cell.
        if p.deformation is True and sim.run_sim is True:
            # Extract time-series deformation data for the plot cell.
            dx = np.asarray([arr[p.plot_cell] for arr in sim.dx_cell_time])
            dy = np.asarray([arr[p.plot_cell] for arr in sim.dy_cell_time])

            # Get the total magnitude.
            disp = np.sqrt(dx**2 + dy**2)

            figD = plt.figure()
            axD = plt.subplot(111)
            axD.plot(sim.time, p.um*disp)
            axD.set_xlabel('Time [s]')
            axD.set_ylabel('Displacement [um]')
            axD.set_title('Displacement of cell ' + str(p.plot_cell) )

            if p.autosave is True:
                savename = savedImg + 'Displacement_' + '.png'
                plt.savefig(savename,dpi=300,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

    #-------------------------------------------------------------------------------------------------------------------
    #                       2D Data Map Plotting
    #-------------------------------------------------------------------------------------------------------------------

    if p.plot_vm2d is True:

        figV, axV, cbV = viz.plotPrettyPolyData(1000*sim.vm_time[-1],
            sim, cells, p,
            clrAutoscale=p.autoscale_Vmem,
            clrMin=p.Vmem_min_clr,
            clrMax=p.Vmem_max_clr,
            number_cells=p.enumerate_cells,
            clrmap=p.default_cm,
            current_overlay=False,
            plotIecm=p.IecmPlot,
        )

        figV.suptitle('Final Vmem',fontsize=14, fontweight='bold')
        axV.set_xlabel('Spatial distance [um]')
        axV.set_ylabel('Spatial distance [um]')
        cbV.set_label('Voltage mV')

        if p.autosave is True:
            savename5 = savedImg + 'final_Vmem_2D' + '.png'
            plt.savefig(savename5,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        figVa, axVa, cbVa = viz.plotPolyData(
            sim, cells, p,
            zdata=1000*sim.vm_ave,
            clrAutoscale=p.autoscale_Vmem,
            clrMin=p.Vmem_min_clr,
            clrMax=p.Vmem_max_clr,
            number_cells=p.enumerate_cells,
            clrmap=p.default_cm,
            current_overlay=False,
            plotIecm=p.IecmPlot,
        )

        # axVa.quiver(p.um*cells.cell_centres[:,0], p.um*cells.cell_centres[:,1], sim.pol_cell_x, sim.pol_cell_y)

        figVa.suptitle('Final Average Vmem', fontsize=14, fontweight='bold')
        axVa.set_xlabel('Spatial distance [um]')
        axVa.set_ylabel('Spatial distance [um]')
        cbVa.set_label('Voltage [mV]')

        if p.autosave is True:
            savename5 = savedImg + 'final_AverageVmem_2D' + '.png'
            plt.savefig(savename5, format='png', transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)


        if p.sim_ECM is True:

            plt.figure()
            plt.imshow(
                1e3*sim.v_env.reshape(cells.X.shape),
                origin='lower',
                extent=[p.um * cells.xmin, p.um * cells.xmax, p.um * cells.ymin, p.um * cells.ymax],
                cmap=p.default_cm,
            )
            plt.colorbar()
            plt.title('Environmental Voltage [mV]')

            if p.autosave is True:
                savename10 = savedImg + 'Final_environmental_V' + '.png'
                plt.savefig(savename10, format='png', transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

    if p.GHK_calc is True:
        figV_ghk, axV_ghk, cbV_ghk = viz.plotPolyData(
            sim, cells, p,
            zdata=1000*sim.vm_GHK_time[-1],
            clrAutoscale=p.autoscale_Vmem,
            clrMin=p.Vmem_min_clr,
            clrMax=p.Vmem_max_clr,
            number_cells=p.enumerate_cells,
            clrmap=p.default_cm,
            current_overlay = False,
            plotIecm=p.IecmPlot,
        )

        figV_ghk.suptitle('Final Vmem using Goldman Equation',fontsize=14, fontweight='bold')
        axV_ghk.set_xlabel('Spatial distance [um]')
        axV_ghk.set_ylabel('Spatial distance [um]')
        cbV_ghk.set_label('Voltage [mV]')

        if p.autosave is True:
            savename5 = savedImg + 'final_Vmem_GHK_2D' + '.png'
            plt.savefig(savename5,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    #-------------------------------------------------------------------------------------------------------------------

    if p.plot_ca2d is True and p.ions_dict['Ca'] == 1:


        figCa, axCa, cbCa = viz.plotPolyData(sim, cells, p, zdata=sim.cc_time[-1][sim.iCa]*1e6, number_cells=p.enumerate_cells,
                         clrAutoscale=p.autoscale_Ca, clrMin=p.Ca_min_clr, clrMax=p.Ca_max_clr,
                         clrmap=p.default_cm)


        axCa.set_title('Final cytosolic Ca2+')
        axCa.set_xlabel('Spatial distance [um]')
        axCa.set_ylabel('Spatial distance [um]')
        cbCa.set_label('Concentration nmol/L')

        if p.autosave is True:
            savename8 = savedImg + 'final_Ca_2D' + '.png'
            plt.savefig(savename8,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        if p.sim_ECM is True:


            if p.smooth_level == 0.0:

                cc_Ca = gaussian_filter(sim.cc_env[sim.iCa].reshape(cells.X.shape), 1.0)

            else:

                cc_Ca = sim.cc_env[sim.iCa].reshape(cells.X.shape)


            plt.figure()
            plt.imshow(
                cc_Ca,
                origin='lower',
                extent=[p.um * cells.xmin, p.um * cells.xmax, p.um * cells.ymin, p.um * cells.ymax],
                cmap=p.default_cm,
            )
            plt.colorbar()
            plt.title('Environmental Calcium [mmol/L]')

            if p.autosave is True:
                savename10 = savedImg + 'Final_environmental_calcium' + '.png'
                plt.savefig(savename10, format='png', transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)


    if p.plot_pH2d is True and p.ions_dict['H'] == 1:
        pHdata = -np.log10(1e-3*sim.cc_time[-1][sim.iH])

        figH, axH, cbH = viz.plotPolyData(sim, cells, p, zdata=pHdata, number_cells=p.enumerate_cells,
                         clrAutoscale=p.autoscale_pH, clrMin=p.pH_min_clr, clrMax=p.pH_max_clr,
                         clrmap=p.default_cm)

        # figH, axH, cbH = viz.plotPrettyPolyData(pHdata, sim,cells,p,
        #     number_cells= p.enumerate_cells, clrAutoscale = p.autoscale_pH,
        #     clrMin = p.pH_min_clr, clrMax = p.pH_max_clr, clrmap = p.default_cm)

        axH.set_title('Final cytosolic pH')
        axH.set_xlabel('Spatial distance [um]')
        axH.set_ylabel('Spatial distance [um]')
        cbH.set_label('pH')

        if p.autosave is True:
            savename8 = savedImg + 'final_pH_2D' + '.png'
            plt.savefig(savename8,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    #----plot 2D pump data--------------------------------------------------------------------------------
    pumpData = sim.rate_NaKATP*1e9

    figPump, axPump, cbPump = viz.plotPrettyPolyData(pumpData, sim, cells, p,
        number_cells=p.enumerate_cells, clrmap=p.default_cm)

    axPump.set_title('Final Na/K-ATPase Pump Rate')
    axPump.set_xlabel('Spatial distance [um]')
    axPump.set_ylabel('Spatial distance [um]')
    cbPump.set_label('Pump Na+ Flux (nmol/m2*s)')

    if p.autosave is True:
        savename8 = savedImg + 'final_NaKPump_2D' + '.png'
        plt.savefig(savename8, format='png', transparent=True)

    if p.turn_all_plots_off is False:
        plt.show(block=False)

    #------------------------------------------------------------------------------------------------------------------

    if p.plot_I2d is True and p.calc_J is True:

        figI, axI, cbI = viz.plotStreamField(
            100*sim.J_cell_x, 100*sim.J_cell_y,
            cells,
            p,
            plot_ecm = False,
            title = 'Intracellular Current Density',
            cb_title = 'Current Density [uA/cm2]',
            show_cells = False,
            colorAutoscale = p.autoscale_I2d,
            minColor=p.I_min_clr,
            maxColor=p.I_max_clr,
        )


        axI.set_xlabel('Spatial distance [um]')
        axI.set_ylabel('Spatial distance [um]')
        cbI.set_label('Current Density [uA/cm2]')

        if p.autosave is True:
            savename10 = savedImg + 'Final_Current_gj' + '.png'
            plt.savefig(savename10,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        if p.sim_ECM is True:

            figI2, axI2, cbI2 = viz.plotStreamField(
                100 * sim.J_env_x, 100 * sim.J_env_y,
                cells,
                p,
                plot_ecm=True,
                title='Extracellular Current Density',
                cb_title='Current Density [uA/cm2]',
                show_cells=False,
                colorAutoscale=p.autoscale_I2d,
                minColor=p.I_min_clr,
                maxColor=p.I_max_clr,
            )


            axI2.set_xlabel('Spatial distance [um]')
            axI2.set_ylabel('Spatial distance [um]')
            cbI2.set_label('Extracellular Current Density [uA/cm2]')

            if p.autosave is True:
                savename11 = savedImg + 'Final_Current_extracellular' + '.png'
                plt.savefig(savename11,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

    #-------------------------------------------------------------------------------------------------------------------

    if p.plot_Efield is True:

        if p.sim_ECM is True:

            viz.plotVectField(sim.E_env_x,sim.E_env_y,cells,p,plot_ecm = True,
                title='Final Electric Field', cb_title = 'Electric Field [V/m]',
                colorAutoscale = p.autoscale_Efield, minColor = p.Efield_min_clr,
                maxColor = p.Efield_max_clr)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

            if p.autosave is True:
                if p.sim_ECM is True:
                    savename = savedImg + 'Final_Electric_Field_ECM' + '.png'
                    plt.savefig(savename,format='png',transparent=True)

        viz.plotVectField(
            sim.E_gj_x,sim.E_gj_y,cells,p,plot_ecm = False,
            title='Final Electric Field',
            cb_title='Electric Field [V/m]',
            colorAutoscale=p.autoscale_Efield,
            minColor=p.Efield_min_clr,
            maxColor=p.Efield_max_clr,
        )

        if p.autosave is True:
            savename = savedImg + 'Final_Electric_Field_GJ' + '.png'
            plt.savefig(savename,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    #------------------------------------------------------------------------------------------------------------------
    if p.plot_P is True and np.mean(sim.P_cells_time) != 0.0:

        figP, axP, cbP = viz.plotPolyData(sim, cells,p,zdata=sim.P_cells,number_cells=p.enumerate_cells,
        clrAutoscale = p.autoscale_P, clrMin = p.P_min_clr, clrMax = p.P_max_clr, clrmap = p.default_cm)

        axP.set_title('Final Pressure in Cell Network')
        axP.set_xlabel('Spatial distance [um]')
        axP.set_ylabel('Spatial distance [um]')
        cbP.set_label('Pressure [Pa]')

        if p.autosave is True:
            savename13 = savedImg + 'final_P_2D_gj' + '.png'
            plt.savefig(savename13,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

    #------------------------------------------------------------------------------------------------------------------

    if p.deformation is True and sim.run_sim is True:
        viz.plotStreamField(
            p.um*sim.dx_cell_time[-1],
            p.um*sim.dy_cell_time[-1],
            cells, p,
            plot_ecm=False,
            title='Final Displacement of Cell Collective',
            cb_title='Displacement [um]',
            show_cells=p.showCells,
            colorAutoscale=p.autoscale_Deformation_ani,
            minColor=p.Deformation_ani_min_clr,
            maxColor=p.Deformation_ani_max_clr,
        )

        if p.autosave is True:
            savename13 = savedImg + 'final_displacement_2D' + '.png'
            plt.savefig(savename13,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)


    if (p.plot_Vel is True and p.fluid_flow is True):
        viz.plotStreamField(
            (1e9)*sim.u_cells_x,
            (1e9)*sim.u_cells_y,
            cells, p,
            plot_ecm=False,
            title='Final Fluid Velocity in Cell Collective',
            cb_title='Velocity [nm/s]',
            colorAutoscale=p.autoscale_Vel,
            minColor=p.Vel_min_clr,
            maxColor=p.Vel_max_clr,
        )

        if p.autosave is True:
            savename13 = savedImg + 'final_vel_2D_gj' + '.png'
            plt.savefig(savename13,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        if p.sim_ECM is True:
            viz.plotStreamField(
                (1e6)*sim.u_env_x,
                (1e6)*sim.u_env_y,
                cells, p, plot_ecm=True,
                title='Final Fluid Velocity in Cell Collective',
                cb_title='Velocity [um/s]',
                colorAutoscale=p.autoscale_Vel,
                minColor=p.Vel_min_clr,
                maxColor=p.Vel_max_clr,
            )

            if p.autosave is True:
                savename13 = savedImg + 'final_vel_2D_env' + '.png'
                plt.savefig(savename13,format='png',transparent=True)

            if p.turn_all_plots_off is False:
                plt.show(block=False)

    # if p.gj_flux_sensitive is True or p.v_sensitive_gj is True:
    # viz.plotMemData(cells,p,zdata=sim.rho_gj,clrmap=p.default_cm)
    fig_x = plt.figure()
    ax_x = plt.subplot(111)
    con_segs = cells.nn_edges
    connects = p.um*np.asarray(con_segs)
    collection = LineCollection(connects, array=sim.gjopen, cmap= p.background_cm, linewidths=2.0)
    ax_x.add_collection(collection)
    # collection.set_clim(0, 1)
    cb = fig_x.colorbar(collection)
    plt.axis('equal')
    plt.axis([cells.xmin*p.um,cells.xmax*p.um,cells.ymin*p.um,cells.ymax*p.um])


    cb.set_label('Relative Permeability')
    ax_x.set_xlabel('Spatial x [um]')
    ax_x.set_ylabel('Spatial y [um')
    ax_x.set_title('Final Gap Junction Relative Permeability')

    if p.autosave is True:
        savename = savedImg + 'final_gjState' + '.png'
        plt.savefig(savename,format='png',transparent=True)

    if p.turn_all_plots_off is False:
        plt.show(block=False)

    if p.sim_eosmosis is True and sim.run_sim is True:
        viz.plotMemData(cells,p,zdata=sim.rho_pump,clrmap=p.default_cm)
        plt.xlabel('Spatial Dimension [um]')
        plt.ylabel('Spatial Dimension [um]')
        plt.title('Membrane ion pump density factor')

        if p.autosave is True:
            savename = savedImg + 'final_pumps_2D' + '.png'
            plt.savefig(savename,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)

        viz.plotMemData(cells,p,zdata=sim.rho_channel,clrmap=p.default_cm)
        plt.xlabel('Spatial Dimension [um]')
        plt.ylabel('Spatial Dimension [um]')
        plt.title('Membrane ion channel density factor')

        if p.autosave is True:
            savename = savedImg + 'final_channels_2D' + '.png'
            plt.savefig(savename,format='png',transparent=True)

        if p.turn_all_plots_off is False:
            plt.show(block=False)
