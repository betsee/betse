#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level facilities for **pipelining** (i.e., iteratively displaying and/or
exporting) post-simulation plots.
'''

#FIXME: I believe I've finally tracked down the issue relating to the following
#runtime "pyplot" warning:
#
#    pyplot.py:424: RuntimeWarning: More than 20 figures have been opened. Figures created through the pyplot interface (`matplotlib.pyplot.figure`) are retained until explicitly closed and may consume too much memory. (To control this warning, see the rcParam `figure.max_open_warning`).
#
#The issue is that whenever we call a plotutil.plot* function below (e.g.,
#plotutil.plotSingleCellCData()), we localize the figure returned by that
#function.  That figure will then be garbage collected on whichever of the
#following occurs last: (A) the corresponding local variable goes out of scope
#and (B) the corresponding plot window is closed by the user. Normally, neither
#would be a problem. Except this function is 1,400 lines long, which means that
#each figure's local variable effectively *NEVER* goes out of scope for the
#duration of plotting. Thus, figures will only be garbage collected *AFTER* this
#function terminates -- which is pretty much unacceptable.
#
#There are a couple solutions, thankfully. The simplest would simply be to stop
#localizing figures returned by plotutil.plot*() functions for all unused figure
#locals. The harder but probably more ideal solution would be to refactor all
#plotutil.plot*() functions to stop returning figures altogether. Since the current
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
#FIXME: The above analysis for the following runtime "pyplot" warning is, sadly,
#completely wrong.
#
#    pyplot.py:424: RuntimeWarning: More than 20 figures have been opened. Figures created through the pyplot interface (`matplotlib.pyplot.figure`) are retained until explicitly closed and may consume too much memory. (To control this warning, see the rcParam `figure.max_open_warning`).
#
#The underlying cause is simple. The pyplot.figure() function internally caches
#created figures into a global cache of such figures. There are many reasons why
#this is a bad idea, including the aforementioned warning as well as the
#probably non-thread-safety of this approach. The solution, of course, is to
#simply manually instantiate Figure() instances directly rather than call the
#pyplot.figure() function: e.g.,
#
#    # Instead of this...
#    pyplot.figure()
#    subplot = pyplot.subplot()
#
#    # ...do this.
#    fig = Figure()
#    subplot = fig.subplot()  # does this work? no idea.
#
#The obvious downside of the above approach, of course, is proper creation of
#non-blocking plots -- which will probably require creation and use of some sort
#of thread-safe cache. Assuming a one-to-one mapping is preserved between each
#non-blocking plot and each thread, the simplest mechanism would be to simply
#cache that plot's figure as an attribute of that thread. Sweet, no?

# ....................{ IMPORTS                            }....................
import matplotlib
import numpy as np
from betse.exceptions import BetseSimConfigException
from betse.lib.matplotlib import mplutil
from betse.lib.matplotlib.matplotlibs import mpl_config
from betse.science.config.export.confvisabc import SimConfVisualListable
from betse.science.simulate.pipe import piperunreq
from betse.science.simulate.pipe.pipeabc import SimPipelinerExportABC
from betse.science.simulate.pipe.piperun import exporter_metadata
from betse.science.simulate.simphase import SimPhaseABC, SimPhaseKind
from betse.science.visual.plot import plotutil
from betse.util.path import dirs, paths
from betse.util.type.types import type_check, IterableTypes
from matplotlib import pyplot as pyplot
from matplotlib.collections import LineCollection
from scipy.ndimage.filters import gaussian_filter

# ....................{ SUBCLASSES                         }....................
class PlotCellsPipeliner(SimPipelinerExportABC):
    '''
    **Post-simulation plot pipeline** (i.e., class iteratively creating all
    post-simulation plots requested by the current simulation configuration).
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, label_singular='plot', **kwargs)

        # If this index is not that of an actual cell, raise an exception.
        if self._phase.p.plot_cell not in self._phase.cells.cell_i:
            raise BetseSimConfigException(
                'Plot cell index {} invalid '
                '(i.e., not in range [{}, {}]).'.format(
                    self._phase.p.plot_cell,
                    self._phase.cells.cell_i[0],
                    self._phase.cells.cell_i[-1],
                ))

        # If saving post-simulation plots...
        if self._phase.p.plot.is_after_sim_save:
            # Create the top-level directory containing these plots if needed.
            dirs.make_unless_dir(self._phase.save_dirname)

    # ..................{ SUPERCLASS                         }..................
    @property
    def is_enabled(self) -> bool:
        return self._phase.p.plot.is_after_sim


    @property
    def _runners_conf(self) -> IterableTypes:
        return self._phase.p.plot.after_sim_pipeline

    # ..................{ EXPORTERS ~ cell : current         }..................
    #FIXME: Force every currently optional "conf" parameter to be mandatory.
    #Specifically, excise "= None" everywhere below.
    @exporter_metadata(categories=(
        'Single Cell', 'Transmembrane Current Density',))
    def export_cell_current_membrane(self, conf: SimConfVisualListable = None) -> None:
        '''
        Plot all transmembrane current densities for all time steps for the
        single cell indexed by the current simulation configuration.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        pyplot.figure()
        axI = pyplot.subplot(111)

        if self._phase.p.sim_ECM:
            # Total cell current storage vector.
            Imem = []

            # Indices of all membranes for this cell.
            mems_for_plotcell = self._phase.cells.cell_to_mems[
                self._phase.p.plot_cell]

            for t in range(len(self._phase.sim.time)):
                memArray = self._phase.sim.I_mem_time[t]

                # X and Y components of the net current at each membrane.
                Ixo = (
                    memArray[mems_for_plotcell] *
                    self._phase.cells.mem_vects_flat[mems_for_plotcell,2] *
                    self._phase.cells.mem_sa[mems_for_plotcell]
                )
                Iyo = (
                    memArray[mems_for_plotcell] *
                    self._phase.cells.mem_vects_flat[mems_for_plotcell,3] *
                    self._phase.cells.mem_sa[mems_for_plotcell]
                )

                # X and Y components of the net current over all membranes,
                # implicitly accounting for in-out directionality.
                Ix = np.sum(Ixo)
                Iy = np.sum(Iyo)

                # Current density, obtained by dividing the total magnitude of
                # the net current by this cell's surface area.
                Io = np.sqrt(Ix**2 + Iy**2) / self._phase.cells.cell_sa[
                    self._phase.p.plot_cell]
                Imem.append(100*Io)
        else:
            Imem = [
                100*memArray[self._phase.p.plot_cell]
                for memArray in self._phase.sim.I_mem_time
            ]

        axI.plot(self._phase.sim.time, Imem)
        axI.set_xlabel('Time [s]')
        axI.set_ylabel('Current density [uA/cm2]')
        axI.set_title(
            'Transmembrane current density for cell {}'.format(
                self._phase.p.plot_cell))

        # Export this plot to disk and/or display.
        self._export(basename='Imem_time')

    # ..................{ EXPORTERS ~ cell : ion             }..................
    @exporter_metadata(
        categories=('Single Cell', 'Ion Concentration', 'M anion'),
        requirements={piperunreq.ION_M,},
    )
    def export_cell_ion_m_anion(self, conf: SimConfVisualListable = None) -> None:
        '''
        Plot all M anion (i.e., M-) ion concentrations for all time steps for
        the single cell indexed by the current simulation configuration.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figConcsM, axConcsM = plotutil.plotSingleCellCData(
            self._phase.sim.cc_time,
            self._phase.sim.time,
            self._phase.sim.iM,
            self._phase.p.plot_cell,
            fig=None, ax=None,
            lncolor='r',
            ionname='M-',
        )
        axConcsM.set_title(
            'M Anion concentration in cell {}'.format(self._phase.p.plot_cell))

        # Export this plot to disk and/or display.
        self._export(basename='concM_time')


    @exporter_metadata(
        categories=('Single Cell', 'Ion Concentration', 'Potassium'),
        requirements={piperunreq.ION_K,},
    )
    def export_cell_ion_potassium(self, conf: SimConfVisualListable = None) -> None:
        '''
        Plot all potassium (i.e., K+) ion concentrations for all time steps for
        the single cell indexed by the current simulation configuration.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figConcsK, axConcsK = plotutil.plotSingleCellCData(
            self._phase.sim.cc_time,
            self._phase.sim.time,
            self._phase.sim.iK,
            self._phase.p.plot_cell,
            fig=None, ax=None,
            lncolor='b',
            ionname='K+',
        )
        axConcsK.set_title(
            'Potassium concentration in cell {}'.format(
                self._phase.p.plot_cell))

        # Export this plot to disk and/or display.
        self._export(basename='concK_time')


    @exporter_metadata(
        categories=('Single Cell', 'Ion Concentration', 'Sodium'),
        requirements={piperunreq.ION_NA,},
    )
    def export_cell_ion_sodium(self, conf: SimConfVisualListable = None) -> None:
        '''
        Plot all sodium (i.e., Na+) ion concentrations for all time steps for
        the single cell indexed by the current simulation configuration.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        # Plot cell sodium concentration versus time.
        figConcsNa, axConcsNa = plotutil.plotSingleCellCData(
            self._phase.sim.cc_time,
            self._phase.sim.time,
            self._phase.sim.iNa,
            self._phase.p.plot_cell,
            fig=None,
            ax=None,
            lncolor='g',
            ionname='Na+',
        )
        axConcsNa.set_title(
            'Sodium concentration in cell {}'.format(self._phase.p.plot_cell))

        # Export this plot to disk and/or display.
        self._export(basename='concNa_time')

    # ..................{ EXPORTERS ~ pressure               }..................
    @exporter_metadata(
        categories=('Single Cell', 'Pressure Total',),
        requirements={piperunreq.PRESSURE_TOTAL,},
    )
    def export_cell_pressure_total(self, conf: SimConfVisualListable = None) -> None:
        '''
        Plot the **total cellular pressure** (i.e., summation of the cellular
        mechanical and osmotic pressure) for all time steps for the single cell
        indexed by the current simulation configuration.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        p_hydro = tuple(
            arr[self._phase.p.plot_cell]
            for arr in self._phase.sim.P_cells_time)

        pyplot.figure()
        axOP = pyplot.subplot(111)
        axOP.plot(self._phase.sim.time, p_hydro)
        axOP.set_xlabel('Time [s]')
        axOP.set_ylabel('Total Pressure [Pa]')
        axOP.set_title(
            'Total pressure in cell {}'.format(self._phase.p.plot_cell))

        # Export this plot to disk and/or display.
        self._export(basename='P_time')

    # ..................{ EXPORTERS ~ cell : pump            }..................
    #FIXME: Actually plot the average of these rates. Currently, this method
    #only plots rates for a single arbitrarily selected membrane of this cell.
    @exporter_metadata(categories=(
        'Single Cell', 'Na-K-ATPase Pump Rate',))
    def export_cell_pump_nakatpase(self, conf: SimConfVisualListable = None) -> None:
        '''
        Plot the averages of the Na-K-ATPase membrane pump rates for all time
        steps for the single cell indexed by the current simulation
        configuration.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        pyplot.figure()
        axNaK = pyplot.subplot()

        if self._phase.p.sim_ECM:
            #FIXME: Generalize to average across all membranes of this cell.
            pump_array_index = self._phase.cells.cell_to_mems[
                self._phase.p.plot_cell][0]
        else:
            pump_array_index = self._phase.p.plot_cell

        pump_rate = tuple(
            pump_array[pump_array_index]
            for pump_array in self._phase.sim.rate_NaKATP_time)

        axNaK.plot(self._phase.sim.time, pump_rate)
        axNaK.set_xlabel('Time [s]')
        axNaK.set_ylabel('Pumping Rate of Na+ Out of Cell [mol/(m2 s)]')
        axNaK.set_title(
            'Rate of NaK-ATPase pump in cell: {}'.format(
                self._phase.p.plot_cell))

        # Export this plot to disk and/or display.
        self._export(basename='NaKATPase_rate')

    # ..................{ EXPORTERS ~ cell : voltage         }..................
    #FIXME: Actually plot the average of these Vmems. Currently, this method
    #only plots Vmems for a single arbitrarily selected membrane of this cell.
    @exporter_metadata(categories=(
        'Single Cell', 'Transmembrane Voltages', 'Average',))
    def export_cell_voltage_membrane(self, conf: SimConfVisualListable = None) -> None:
        '''
        Plot the averages of all transmembrane voltages for all time steps for
        the single cell indexed by the current simulation configuration.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        #FIXME: Generalize to average across all membranes of this cell.
        mem_i = self._phase.cells.cell_to_mems[self._phase.p.plot_cell][0]

        # Plot single cell Vmem vs time.
        figVt, axVt = plotutil.plotSingleCellVData(
            self._phase.sim,
            mem_i,
            self._phase.p,
            fig=None, ax=None,
            lncolor='k',
        )
        axVt.set_title(
            'Voltage (Vmem) in cell {}'.format(self._phase.p.plot_cell))

        # Export this plot to disk and/or display.
        self._export(basename='Vmem_time')


    #FIXME: Actually plot the average of these Vmems. Currently, this method
    #only plots Vmems for a single arbitrarily selected membrane of this cell.
    @exporter_metadata(categories=(
        'Single Cell',
        'Transmembrane Voltages',
        'Fast Fourier Transform (FFT)',
    ))
    def export_cell_voltage_membrane_fft(
        self, conf: SimConfVisualListable = None) -> None:
        '''
        Plot the fast Fourier transform (FFT) of the averages of all
        transmembrane voltages for all time steps for the single cell indexed by
        the current simulation configuration.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        #FIXME: Generalize to average across all membranes of this cell.
        mem_i = self._phase.cells.cell_to_mems[self._phase.p.plot_cell][0]

        figFFT, axFFT = plotutil.plotFFT(
            self._phase.sim.time, self._phase.sim.vm_time, mem_i, lab="Power")
        axFFT.set_title(
            'Fourier transform of Vmem in cell {}'.format(
                self._phase.p.plot_cell))

        # Export this plot to disk and/or display.
        self._export(basename='Vmem_FFT_time')

    # ..................{ PRIVATE ~ preparers                }..................
    def _export_prep(self) -> None:
        '''
        Prepare to export the current plot.
        '''

        #FIXME: DRY. This functionality perfectly duplicates the
        #AnimCellsWhileSolving.__enter__() method, which is bad. To resolve
        #this:
        #
        #* Shift that method into the "VisualCellsABC" superclass.
        #* Refactor all plots to subclass that superclass.

        # Id displaying this plot, do so in a non-blocking manner.
        if self._phase.p.plot.is_after_sim_show:
            # If the current matplotlib backend supports "true" non-blocking
            # behavior, prefer this non-deprecated approach.
            if mpl_config.is_backend_current_nonblockable():
                pyplot.show(block=False)
            # Else, fallback to the deprecated approach guaranteed to apply to
            # all matplotlib backends.
            else:
                # pass
                matplotlib.interactive(True)
                # pyplot.show()


    @type_check
    def _export(self, basename: str) -> None:
        '''
        Export the current plot to the current screen if displaying plots and/or
        to a file with the passed basename if saving plots.

        Parameters
        -----------
        basename : str
            Basename excluding filetype of the plot to be exported.
        '''

        #FIXME: DRY. This functionality perfectly duplicates the
        #AnimCellsWhileSolving.__exit__() method, which is bad.

        # Id displaying this plot *AND* the current matplotlib backend
        # fails to support "true" non-blocking behavior...
        if (self._phase.p.plot.is_after_sim_show and
            not mpl_config.is_backend_current_nonblockable()):
            # Update all artists displayed by this plot.
            pyplot.draw()

            #FIXME: DRY. This functionality perfectly duplicates the
            #VisualCellsABC._show_frame(() method, which is also bad.

            # Temporarily yield the time slice for the smallest amount of time
            # required by the current matplotlib backend to handle queued events in
            # the GUI-specific event loop of the current process.
            with mplutil.deprecations_ignored():
                pyplot.pause(0.0001)

            # Disable the "fake" non-blocking behavior enabled by the prior
            # _export_prep() call.
            matplotlib.interactive(False)

        # If saving this plot...
        if self._phase.p.plot.is_after_sim_save:
            # Filetype of the file to be saved.
            filetype = self._phase.p.plot.image_filetype

            # Absolute path of the file to be saved.
            filename = paths.join(
                self._phase.save_dirname,
                'fig_{}.{}'.format(basename, filetype))

            # Save this plot to this file.
            pyplot.savefig(
                filename,
                dpi=self._phase.p.plot.image_dpi,
                format=filetype,
                transparent=True,
            )

# ....................{ OBSOLETE                           }....................
#FIXME: Replace *ALL* functionality defined below with the "PlotCellsPipeliner"
#class defined above.
@type_check
def pipeline(phase: SimPhaseABC) -> None:
    '''
    Display and/or save all currently enabled plots for the passed simulation
    phase.

    Parameters
    -----------
    phase: SimPhaseABC
        Current simulation phase.
    '''

    # Post-simulation animation pipeline producing all such animations.
    pipeliner = PlotCellsPipeliner(phase)

    #FIXME: Replace *ALL* logic below with the following single call:
    #    pipeliner.run()
    #FIXME: When doing so, note that *ALL* uses of hardcoded plot-specific
    #parameter options (e.g., "self._phase.p.I_ani_min_clr") will need to be
    #refactored to use the general-purpose settings for the current plot.
    #FIXME: Likewise, refactor tests to exercise the new dynamic pipeline schema
    #rather than the obsolete hardcoded schema.

    # If post-simulation plots are disabled, noop.
    if not phase.p.plot.is_after_sim:
       return

    #FIXME: Excise this once no longer required.
    # Substring prefixing the absolute path of each plot created below.
    savedImg = paths.join(phase.save_dirname, 'fig_')

    #-------------------------------------------------------------------------------------------------------------------
    #               SINGLE CELL DATA GRAPHS
    #-------------------------------------------------------------------------------------------------------------------

    #FIXME: Replace these local variable placeholders with the equivalent
    #"phase.sim", "phase.cells", and "phase.p". To do so, continue iteratively
    #pushing these declarations further and further down this function until
    #they are no longer required at all.
    sim   = phase.sim
    cells = phase.cells
    p     = phase.p

    if p.plot_single_cell_graphs:
        # Plot all cell transmembrane voltages.
        pipeliner.export_cell_voltage_membrane()
        pipeliner.export_cell_voltage_membrane_fft()

        # Plot all cell transmembrane current densities.
        pipeliner.export_cell_current_membrane()

        # Plot all Na-K-ATPase pump rates.
        pipeliner.export_cell_pump_nakatpase()

        # If M anions are enabled, plot all cell M anion concentrations.
        if phase.p.ions_dict['M'] == 1:
            pipeliner.export_cell_ion_m_anion()

        # If potassium is enabled, plot all cell potassium concentrations.
        if phase.p.ions_dict['K'] == 1:
            pipeliner.export_cell_ion_potassium()

        # If sodium is enabled, plot all cell sodium concentrations.
        if phase.p.ions_dict['Na'] == 1:
            pipeliner.export_cell_ion_sodium()

        # If pressure is enabled, plot all cell pressure totals.
        if piperunreq.PRESSURE_TOTAL.is_satisfied(phase):
            pipeliner.export_cell_pressure_total()

        # optional 1D plots--------------------------------------------------------------------------------------------


        # Plot cell calcium vs time (if Ca enabled in ion profiles).
        if p.ions_dict['Ca'] ==1:
            figA, axA = plotutil.plotSingleCellCData(
                sim.cc_time, sim.time, sim.iCa, p.plot_cell,
                fig=None,
                ax=None,
                lncolor='g',
                ionname='Ca2+ cell',
            )
            titCa = 'Cytosolic Ca2+ in cell index ' + str(p.plot_cell)
            axA.set_title(titCa)

            if p.plot.is_after_sim_save is True:
                savename3 = savedImg + 'cytosol_Ca_time' + '.png'
                pyplot.savefig(savename3,dpi=300,format='png',transparent=True)

            if phase.p.plot.is_after_sim_show:
                pyplot.show(block=False)


        if p.deform_osmo is True:
            p_osmo = [arr[p.plot_cell] for arr in sim.osmo_P_delta_time]
            pyplot.figure()
            axOP = pyplot.subplot(111)
            axOP.plot(sim.time, p_osmo)
            axOP.set_xlabel('Time [s]')
            axOP.set_ylabel('Osmotic Pressure [Pa]')
            axOP.set_title('Osmotic pressure in cell ' + str(p.plot_cell) )

            if p.plot.is_after_sim_save is True:
                savename = savedImg + 'OsmoticP_' + '.png'
                pyplot.savefig(savename,dpi=300,format='png',transparent=True)

            if phase.p.plot.is_after_sim_show:
                pyplot.show(block=False)

        # Total displacement in cell.
        if phase.p.deformation is True and phase.kind is SimPhaseKind.SIM:
            # Extract time-series deformation data for the plot cell.
            dx = np.asarray([arr[p.plot_cell] for arr in sim.dx_cell_time])
            dy = np.asarray([arr[p.plot_cell] for arr in sim.dy_cell_time])

            # Get the total magnitude.
            disp = np.sqrt(dx**2 + dy**2)

            pyplot.figure()
            axD = pyplot.subplot(111)
            axD.plot(sim.time, p.um*disp)
            axD.set_xlabel('Time [s]')
            axD.set_ylabel('Displacement [um]')
            axD.set_title('Displacement of cell ' + str(p.plot_cell) )

            if p.plot.is_after_sim_save is True:
                savename = savedImg + 'Displacement_' + '.png'
                pyplot.savefig(savename,dpi=300,format='png',transparent=True)

            if phase.p.plot.is_after_sim_show:
                pyplot.show(block=False)

    #-------------------------------------------------------------------------------------------------------------------
    #                       2D Data Map Plotting
    #-------------------------------------------------------------------------------------------------------------------

    if p.plot_vm2d is True:

        figV, axV, cbV = plotutil.plotPrettyPolyData(1000*sim.vm_time[-1],
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

        if p.plot.is_after_sim_save is True:
            savename5 = savedImg + 'final_Vmem_2D' + '.png'
            pyplot.savefig(savename5,format='png',transparent=True)

        if phase.p.plot.is_after_sim_show:
            pyplot.show(block=False)

        figVa, axVa, cbVa = plotutil.plotPolyData(
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

        if p.plot.is_after_sim_save is True:
            savename5 = savedImg + 'final_AverageVmem_2D' + '.png'
            pyplot.savefig(savename5, format='png', transparent=True)

        if phase.p.plot.is_after_sim_show:
            pyplot.show(block=False)


        if p.sim_ECM is True:

            vv = sim.v_env.reshape(cells.X.shape)
            vv = gaussian_filter(vv, 1, mode = 'constant')

            pyplot.figure()
            pyplot.imshow(
                1e3*vv,
                origin='lower',
                extent=[p.um * cells.xmin, p.um * cells.xmax, p.um * cells.ymin, p.um * cells.ymax],
                cmap=p.default_cm,
            )
            pyplot.colorbar()
            pyplot.title('Environmental Voltage [mV]')

            if p.plot.is_after_sim_save is True:
                savename10 = savedImg + 'Final_environmental_V' + '.png'
                pyplot.savefig(savename10, format='png', transparent=True)

            if phase.p.plot.is_after_sim_show:
                pyplot.show(block=False)

    if p.GHK_calc is True:
        figV_ghk, axV_ghk, cbV_ghk = plotutil.plotPolyData(
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

        if p.plot.is_after_sim_save is True:
            savename5 = savedImg + 'final_Vmem_GHK_2D' + '.png'
            pyplot.savefig(savename5,format='png',transparent=True)

        if phase.p.plot.is_after_sim_show:
            pyplot.show(block=False)

    #-------------------------------------------------------------------------------------------------------------------

    if p.plot_ca2d is True and p.ions_dict['Ca'] == 1:


        figCa, axCa, cbCa = plotutil.plotPolyData(sim, cells, p, zdata=sim.cc_time[-1][sim.iCa]*1e6, number_cells=p.enumerate_cells,
                         clrAutoscale=p.autoscale_Ca, clrMin=p.Ca_min_clr, clrMax=p.Ca_max_clr,
                         clrmap=p.default_cm)


        axCa.set_title('Final cytosolic Ca2+')
        axCa.set_xlabel('Spatial distance [um]')
        axCa.set_ylabel('Spatial distance [um]')
        cbCa.set_label('Concentration nmol/L')

        if p.plot.is_after_sim_save is True:
            savename8 = savedImg + 'final_Ca_2D' + '.png'
            pyplot.savefig(savename8,format='png',transparent=True)

        if phase.p.plot.is_after_sim_show:
            pyplot.show(block=False)

        if p.sim_ECM is True:


            if p.smooth_level == 0.0:

                cc_Ca = gaussian_filter(sim.cc_env[sim.iCa].reshape(cells.X.shape), 1.0)

            else:

                cc_Ca = sim.cc_env[sim.iCa].reshape(cells.X.shape)


            pyplot.figure()
            pyplot.imshow(
                cc_Ca,
                origin='lower',
                extent=[p.um * cells.xmin, p.um * cells.xmax, p.um * cells.ymin, p.um * cells.ymax],
                cmap=p.default_cm,
            )
            pyplot.colorbar()
            pyplot.title('Environmental Calcium [mmol/L]')

            if p.plot.is_after_sim_save is True:
                savename10 = savedImg + 'Final_environmental_calcium' + '.png'
                pyplot.savefig(savename10, format='png', transparent=True)

            if phase.p.plot.is_after_sim_show:
                pyplot.show(block=False)


    if p.plot_pH2d is True and p.ions_dict['H'] == 1:
        pHdata = -np.log10(1e-3*sim.cc_time[-1][sim.iH])

        figH, axH, cbH = plotutil.plotPolyData(sim, cells, p, zdata=pHdata, number_cells=p.enumerate_cells,
                         clrAutoscale=p.autoscale_pH, clrMin=p.pH_min_clr, clrMax=p.pH_max_clr,
                         clrmap=p.default_cm)

        # figH, axH, cbH = plotutil.plotPrettyPolyData(pHdata, sim,cells,p,
        #     number_cells= p.enumerate_cells, clrAutoscale = p.autoscale_pH,
        #     clrMin = p.pH_min_clr, clrMax = p.pH_max_clr, clrmap = p.default_cm)

        axH.set_title('Final cytosolic pH')
        axH.set_xlabel('Spatial distance [um]')
        axH.set_ylabel('Spatial distance [um]')
        cbH.set_label('pH')

        if p.plot.is_after_sim_save is True:
            savename8 = savedImg + 'final_pH_2D' + '.png'
            pyplot.savefig(savename8,format='png',transparent=True)

        if phase.p.plot.is_after_sim_show:
            pyplot.show(block=False)

    #----plot 2D pump data--------------------------------------------------------------------------------
    pumpData = sim.rate_NaKATP*1e9

    figPump, axPump, cbPump = plotutil.plotPrettyPolyData(pumpData, sim, cells, p,
        number_cells=p.enumerate_cells, clrmap=p.default_cm)

    axPump.set_title('Final Na/K-ATPase Pump Rate')
    axPump.set_xlabel('Spatial distance [um]')
    axPump.set_ylabel('Spatial distance [um]')
    cbPump.set_label('Pump Na+ Flux (nmol/m2*s)')

    if p.plot.is_after_sim_save is True:
        savename8 = savedImg + 'final_NaKPump_2D' + '.png'
        pyplot.savefig(savename8, format='png', transparent=True)

    if phase.p.plot.is_after_sim_show:
        pyplot.show(block=False)

    #------------------------------------------------------------------------------------------------------------------

    if p.plot_I2d is True and p.calc_J is True:

        figI, axI, cbI = plotutil.plotStreamField(
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

        if p.plot.is_after_sim_save is True:
            savename10 = savedImg + 'Final_Current_gj' + '.png'
            pyplot.savefig(savename10,format='png',transparent=True)

        if phase.p.plot.is_after_sim_show:
            pyplot.show(block=False)

        if p.sim_ECM is True:

            figI2, axI2, cbI2 = plotutil.plotStreamField(
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

            if p.plot.is_after_sim_save is True:
                savename11 = savedImg + 'Final_Current_extracellular' + '.png'
                pyplot.savefig(savename11,format='png',transparent=True)

            if phase.p.plot.is_after_sim_show:
                pyplot.show(block=False)

    #-------------------------------------------------------------------------------------------------------------------

    if p.plot_Efield is True:

        if p.sim_ECM is True:

            plotutil.plotVectField(sim.E_env_x,sim.E_env_y,cells,p,plot_ecm = True,
                title='Final Electric Field', cb_title = 'Electric Field [V/m]',
                colorAutoscale = p.autoscale_Efield, minColor = p.Efield_min_clr,
                maxColor = p.Efield_max_clr)

            if phase.p.plot.is_after_sim_show:
                pyplot.show(block=False)

            if p.plot.is_after_sim_save is True:
                if p.sim_ECM is True:
                    savename = savedImg + 'Final_Electric_Field_ECM' + '.png'
                    pyplot.savefig(savename,format='png',transparent=True)

        plotutil.plotVectField(
            sim.E_gj_x,sim.E_gj_y,cells,p,plot_ecm = False,
            title='Final Electric Field',
            cb_title='Electric Field [V/m]',
            colorAutoscale=p.autoscale_Efield,
            minColor=p.Efield_min_clr,
            maxColor=p.Efield_max_clr,
        )

        if p.plot.is_after_sim_save is True:
            savename = savedImg + 'Final_Electric_Field_GJ' + '.png'
            pyplot.savefig(savename,format='png',transparent=True)

        if phase.p.plot.is_after_sim_show:
            pyplot.show(block=False)

    #------------------------------------------------------------------------------------------------------------------
    if p.plot_P and np.mean(sim.P_cells_time) != 0.0:

        figP, axP, cbP = plotutil.plotPolyData(
            sim, cells, p,
            zdata=sim.P_cells,
            number_cells=p.enumerate_cells,
            clrAutoscale=p.autoscale_P,
            clrMin=p.P_min_clr,
            clrMax=p.P_max_clr,
            clrmap=p.default_cm,
        )

        axP.set_title('Final Pressure in Cell Network')
        axP.set_xlabel('Spatial distance [um]')
        axP.set_ylabel('Spatial distance [um]')
        cbP.set_label('Pressure [Pa]')

        if p.plot.is_after_sim_save is True:
            savename13 = savedImg + 'final_P_2D_gj' + '.png'
            pyplot.savefig(savename13,format='png',transparent=True)

        if phase.p.plot.is_after_sim_show:
            pyplot.show(block=False)

    #------------------------------------------------------------------------------------------------------------------

    if (p.plot_Deformation and
        p.deformation and
        phase.kind is SimPhaseKind.SIM):
        plotutil.plotStreamField(
            p.um*sim.dx_cell_time[-1],
            p.um*sim.dy_cell_time[-1],
            cells, p,
            plot_ecm=False,
            title='Final Displacement of Cell Collective',
            cb_title='Displacement [um]',
            show_cells=p.showCells,
            colorAutoscale=p.autoscale_Deformation,
            minColor=p.Deformation_min_clr,
            maxColor=p.Deformation_max_clr,
        )

        if p.plot.is_after_sim_save is True:
            savename13 = savedImg + 'final_displacement_2D' + '.png'
            pyplot.savefig(savename13,format='png',transparent=True)

        if phase.p.plot.is_after_sim_show:
            pyplot.show(block=False)


    if p.plot_Vel and p.fluid_flow:
        plotutil.plotStreamField(
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

        if p.plot.is_after_sim_save is True:
            savename13 = savedImg + 'final_vel_2D_gj' + '.png'
            pyplot.savefig(savename13,format='png',transparent=True)

        if phase.p.plot.is_after_sim_show:
            pyplot.show(block=False)

        if p.sim_ECM is True:
            plotutil.plotStreamField(
                (1e6)*sim.u_env_x,
                (1e6)*sim.u_env_y,
                cells, p, plot_ecm=True,
                title='Final Fluid Velocity in Cell Collective',
                cb_title='Velocity [um/s]',
                colorAutoscale=p.autoscale_Vel,
                minColor=p.Vel_min_clr,
                maxColor=p.Vel_max_clr,
            )

            if p.plot.is_after_sim_save is True:
                savename13 = savedImg + 'final_vel_2D_env' + '.png'
                pyplot.savefig(savename13,format='png',transparent=True)

            if phase.p.plot.is_after_sim_show:
                pyplot.show(block=False)

    # if p.gj_flux_sensitive is True or p.v_sensitive_gj is True:
    # plotutil.plotMemData(cells,p,zdata=sim.rho_gj,clrmap=p.default_cm)
    fig_x = pyplot.figure()
    ax_x = pyplot.subplot(111)
    con_segs = cells.nn_edges
    connects = p.um*np.asarray(con_segs)
    collection = LineCollection(connects, array=sim.gjopen, cmap= p.background_cm, linewidths=2.0)
    ax_x.add_collection(collection)
    # collection.set_clim(0, 1)
    cb = fig_x.colorbar(collection)
    pyplot.axis('equal')
    pyplot.axis([cells.xmin*p.um,cells.xmax*p.um,cells.ymin*p.um,cells.ymax*p.um])

    cb.set_label('Relative Permeability')
    ax_x.set_xlabel('Spatial x [um]')
    ax_x.set_ylabel('Spatial y [um')
    ax_x.set_title('Final Gap Junction Relative Permeability')

    if p.plot.is_after_sim_save is True:
        savename = savedImg + 'final_gjState' + '.png'
        pyplot.savefig(savename,format='png',transparent=True)

    if phase.p.plot.is_after_sim_show:
        pyplot.show(block=False)

    if p.sim_eosmosis is True and phase.kind is SimPhaseKind.SIM:
        plotutil.plotMemData(cells,p,zdata=sim.rho_pump,clrmap=p.default_cm)
        pyplot.xlabel('Spatial Dimension [um]')
        pyplot.ylabel('Spatial Dimension [um]')
        pyplot.title('Membrane ion pump density factor')

        if p.plot.is_after_sim_save:
            savename = savedImg + 'final_pumps_2D' + '.png'
            pyplot.savefig(savename,format='png',transparent=True)

        if phase.p.plot.is_after_sim_show:
            pyplot.show(block=False)

        plotutil.plotMemData(cells,p,zdata=sim.rho_channel,clrmap=p.default_cm)
        pyplot.xlabel('Spatial Dimension [um]')
        pyplot.ylabel('Spatial Dimension [um]')
        pyplot.title('Membrane ion channel density factor')

        if p.plot.is_after_sim_save:
            savename = savedImg + 'final_channels_2D' + '.png'
            pyplot.savefig(savename,format='png',transparent=True)

        if phase.p.plot.is_after_sim_show:
            pyplot.show(block=False)
