#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Post-simulation cell cluster plot pipelines** (i.e., pipelines plotting
simulated data of all cells in the cell cluster).
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.science.config.visual.confvisabc import SimConfVisualCellsListItem
from betse.science.simulate.pipe import piperunreq
from betse.science.simulate.pipe.piperun import piperunner
from betse.science.visual.plot import plotutil
from betse.science.visual.plot.pipe.plotpipeabc import PlotPipeABC
from betse.util.type.types import IterableTypes
from matplotlib import pyplot as pyplot
from matplotlib.collections import LineCollection
from scipy.ndimage.filters import gaussian_filter

# ....................{ SUBCLASSES                         }....................
class PlotCellsPipe(PlotPipeABC):
    '''
    **Post-simulation cell cluster plot pipeline** (i.e., object iteratively
    displaying and/or saving all plots depicting all cells of the cell cluster,
    produced after initialization and simulation as specified by the current
    simulation configuration).
    '''

    # ..................{ INITIALIZERS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, label_singular='cell cluster plot', **kwargs)

    # ..................{ SUPERCLASS                         }..................
    @property
    def _runners_conf(self) -> IterableTypes:
        return self._phase.p.plot.plots_cells_after_sim

    # ..................{ EXPORTERS ~ channel                }..................
    # @piperunner(
    #     categories=('Ion Channel', 'Density Factor',),
    #     requirements={piperunreq.ELECTROOSMOSIS,},
    # )
    # def export_channel_density(self, conf: SimConfVisualCellsListItem) -> None:
    #     '''
    #     Plot all cell membrane ion channel density factors for the cell cluster
    #     at the last time step.
    #     '''
    #
    #     # Prepare to export the current plot.
    #     self._export_prep()
    #
    #     plotutil.plotMemData(
    #         self._phase.cells,
    #         self._phase.p,
    #         zdata=self._phase.sim.rho_channel,
    #         clrmap=self._phase.p.default_cm,
    #     )
    #     pyplot.xlabel('Spatial Dimension [um]')
    #     pyplot.ylabel('Spatial Dimension [um]')
    #     pyplot.title('Membrane ion channel density factor')
    #
    #     # Export this plot to disk and/or display.
    #     self._export(basename='final_channels_2D')

    # ..................{ EXPORTERS ~ current                }..................
    @piperunner(
        categories=('Current Density', 'Intracellular',),
        # requirements={piperunreq.VOLTAGE_POLARITY,},
    )
    def export_currents_intra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all intracellular current densities for the cell cluster at the
        last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figI, axI, cbI = plotutil.plotStreamField(
            100*self._phase.sim.J_cell_x,
            100*self._phase.sim.J_cell_y,
            self._phase.cells,
            self._phase.p,
            plot_ecm=False,
            title='Intracellular Current Density',
            cb_title='Current Density [uA/cm2]',
            show_cells=False,
            colorAutoscale=conf.is_color_autoscaled,
            minColor=conf.color_min,
            maxColor=conf.color_max,
        )

        axI.set_xlabel('Spatial distance [um]')
        axI.set_ylabel('Spatial distance [um]')
        cbI.set_label('Current Density [uA/cm2]')

        # Export this plot to disk and/or display.
        self._export(basename='Final_Current_gj')


    @piperunner(categories=('Current Density', 'Extracellular',))
    def export_currents_extra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all extracellular current densities for the cell cluster
        environment at the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figI2, axI2, cbI2 = plotutil.plotStreamField(
            100*self._phase.sim.J_env_x,
            100*self._phase.sim.J_env_y,
            self._phase.cells,
            self._phase.p,
            plot_ecm=True,
            title='Extracellular Current Density',
            cb_title='Current Density [uA/cm2]',
            show_cells=False,
            colorAutoscale=conf.is_color_autoscaled,
            minColor=conf.color_min,
            maxColor=conf.color_max,
        )

        axI2.set_xlabel('Spatial distance [um]')
        axI2.set_ylabel('Spatial distance [um]')
        cbI2.set_label('Extracellular Current Density [uA/cm2]')

        # Export this plot to disk and/or display.
        self._export(basename='Final_Current_extracellular')

    # ..................{ EXPORTERS ~ deform                 }..................
    @piperunner(
        categories=('Deformation', 'Total',),
        requirements={piperunreq.DEFORM,},
    )
    def export_deform_total(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all **total cellular displacements** (i.e., summations of all
        cellular deformations due to galvanotropic and osmotic pressure body
        forces) for the cell cluster at the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        plotutil.plotStreamField(
            self._phase.p.um*self._phase.sim.dx_cell_time[-1],
            self._phase.p.um*self._phase.sim.dy_cell_time[-1],
            self._phase.cells, self._phase.p,
            plot_ecm=False,
            title='Final Displacement of Cell Collective',
            cb_title='Displacement [um]',
            show_cells=self._phase.p.showCells,
            colorAutoscale=conf.is_color_autoscaled,
            minColor=conf.color_min,
            maxColor=conf.color_max,
        )

        # Export this plot to disk and/or display.
        self._export(basename='final_displacement_2D')

    # ..................{ EXPORTERS ~ electric               }..................
    @piperunner(
        categories=('Electric Field', 'Intracellular',))
    def export_electric_intra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all intracellular electric field lines for the cell cluster at the
        last time step.
        '''

        # Prepare to export the electric plot.
        self._export_prep()

        plotutil.plotVectField(
            self._phase.sim.E_cell_x,
            self._phase.sim.E_cell_y,
            self._phase.cells,
            self._phase.p,
            plot_ecm=False,
            title='Final Electric Field',
            cb_title='Electric Field [V/m]',
            colorAutoscale=conf.is_color_autoscaled,
            minColor=conf.color_min,
            maxColor=conf.color_max,
        )

        # Export this plot to disk and/or display.
        self._export(basename='Final_Electric_Field_GJ')


    @piperunner(
        categories=('Electric Field', 'Extracellular',),
        requirements={piperunreq.ECM,},
    )
    def export_electric_extra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all extracellular electric field lines for the cell cluster
        environment at the last time step.
        '''

        # Prepare to export the electric plot.
        self._export_prep()

        plotutil.plotVectField(
            self._phase.sim.E_env_x,
            self._phase.sim.E_env_y,
            self._phase.cells,
            self._phase.p,
            plot_ecm=True,
            title='Final Electric Field',
            cb_title='Electric Field [V/m]',
            colorAutoscale=conf.is_color_autoscaled,
            minColor=conf.color_min,
            maxColor=conf.color_max,
        )

        # Export this plot to disk and/or display.
        self._export(basename='Final_Electric_Field_ECM')

    # ..................{ EXPORTERS ~ fluid                  }..................
    @piperunner(
        categories=('Fluid Flow', 'Intracellular',),
        requirements={piperunreq.FLUID,},
    )
    def export_fluid_intra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all intracellular fluid flow field lines for the cell cluster at
        the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        plotutil.plotStreamField(
            1e6*self._phase.sim.u_cells_x,
            1e6*self._phase.sim.u_cells_y,
            self._phase.cells, self._phase.p,
            plot_ecm=False,
            title='Final Fluid Velocity in Cell Collective',
            cb_title='Velocity [um/s]',
            colorAutoscale=conf.is_color_autoscaled,
            minColor=conf.color_min,
            maxColor=conf.color_max,
        )

        # Export this plot to disk and/or display.
        self._export(basename='final_vel_2D_gj')


    @piperunner(
        categories=('Fluid Flow', 'Extracellular',),
        requirements={piperunreq.FLUID, piperunreq.ECM,},
    )
    def export_fluid_extra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all extracellular fluid flow field lines for the cell cluster
        environment at the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        plotutil.plotStreamField(
            1e9*self._phase.sim.u_env_x,
            1e9*self._phase.sim.u_env_y,
            self._phase.cells, self._phase.p,
            plot_ecm=True,
            title='Final Fluid Velocity in Environment',
            cb_title='Velocity [nm/s]',
            colorAutoscale=conf.is_color_autoscaled,
            minColor=conf.color_min,
            maxColor=conf.color_max,
        )

        # Export this plot to disk and/or display.
        self._export(basename='final_vel_2D_env')

    # ..................{ EXPORTERS ~ ion : calcium          }..................
    @piperunner(
        categories=('Ion Concentration', 'Calcium', 'Intracellular',),
        requirements={piperunreq.ION_CALCIUM,},
    )
    def export_ion_calcium_intra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all intracellular calcium (i.e., Ca2+) ion concentrations for the
        cell cluster at the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figCa, axCa, cbCa = plotutil.plotPolyData(
            self._phase.sim, self._phase.cells, self._phase.p,
            zdata=self._phase.sim.cc_time[-1][self._phase.sim.iCa]*1e6,
            number_cells=self._phase.p.enumerate_cells,
            clrmap=self._phase.p.default_cm,
            clrAutoscale=conf.is_color_autoscaled,
            clrMin=conf.color_min,
            clrMax=conf.color_max,
        )

        axCa.set_title('Final cytosolic Ca2+')
        axCa.set_xlabel('Spatial distance [um]')
        axCa.set_ylabel('Spatial distance [um]')
        cbCa.set_label('Concentration nmol/L')

        # Export this plot to disk and/or display.
        self._export(basename='final_Ca_2D')


    @piperunner(
        categories=('Ion Concentration', 'Calcium', 'Extracellular',),
        requirements={piperunreq.ION_CALCIUM, piperunreq.ECM,},
    )
    def export_ion_calcium_extra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all extracellular calcium (i.e., Ca2+) ion concentrations for the
        cell cluster environment at the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        cc_Ca = self._phase.sim.cc_env[self._phase.sim.iCa].reshape(
                self._phase.cells.X.shape)

        pyplot.figure()
        pyplot.imshow(
            cc_Ca,
            origin='lower',
            extent=self._cells_extent,
            cmap=self._phase.p.default_cm,
        )
        pyplot.colorbar()
        pyplot.title('Environmental Calcium [mmol/L]')

        # Export this plot to disk and/or display.
        self._export(basename='Final_environmental_calcium')

    # ..................{ EXPORTERS ~ ion : hydrogen         }..................
    # @piperunner(
    #     categories=('Ion Concentration', 'Hydrogen', 'Intracellular',),
    #     requirements={piperunreq.ION_HYDROGEN,},
    # )
    # def export_ion_hydrogen_intra(self, conf: SimConfVisualCellsListItem) -> None:
    #     '''
    #     Plot all intracellular hydrogen (i.e., H+) ion concentrations for the
    #     cell cluster at the last time step.
    #     '''
    #
    #     # Prepare to export the current plot.
    #     self._export_prep()
    #
    #     pHdata = -np.log10(1e-3*self._phase.sim.cc_time[-1][self._phase.sim.iH])
    #
    #     figH, axH, cbH = plotutil.plotPolyData(
    #         self._phase.sim, self._phase.cells, self._phase.p,
    #         zdata=pHdata,
    #         number_cells=self._phase.p.enumerate_cells,
    #         clrmap=self._phase.p.default_cm,
    #         clrAutoscale=conf.is_color_autoscaled,
    #         clrMin=conf.color_min,
    #         clrMax=conf.color_max,
    #     )
    #
    #     # figH, axH, cbH = plotutil.plotPrettyPolyData(pHdata, sim,cells,p,
    #     #     number_cells= p.enumerate_cells, clrAutoscale = p.autoscale_pH,
    #     #     clrMin = p.pH_min_clr, clrMax = p.pH_max_clr, clrmap = p.default_cm)
    #
    #     axH.set_title('Final cytosolic pH')
    #     axH.set_xlabel('Spatial distance [um]')
    #     axH.set_ylabel('Spatial distance [um]')
    #     cbH.set_label('pH')
    #
    #     # Export this plot to disk and/or display.
    #     self._export(basename='final_pH_2D')

    # ..................{ EXPORTERS ~ junction               }..................
    @piperunner(
        categories=('Gap Junction', 'Connectivity State',))
    def export_junction_state(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all **gap junction connectivity states** (i.e., relative
        permeabilities of the gap junctions connecting all cell membranes) for
        the cell cluster at the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        fig_x = pyplot.figure()
        ax_x = pyplot.subplot(111)

        con_segs = self._phase.cells.nn_edges
        connects = self._phase.p.um*np.asarray(con_segs)
        collection = LineCollection(
            connects,
            array=self._phase.sim.gjopen,
            cmap=self._phase.p.background_cm,
            linewidths=2.0,
        )
        # collection.set_clim(0, 1)

        ax_x.add_collection(collection)
        cb = fig_x.colorbar(collection)
        pyplot.axis('equal')
        pyplot.axis(self._cells_extent)

        cb.set_label('Relative Permeability')
        ax_x.set_xlabel('Spatial x [um]')
        ax_x.set_ylabel('Spatial y [um')
        ax_x.set_title('Final Gap Junction Relative Permeability')

        # Export this plot to disk and/or display.
        self._export(basename='final_gjState')

    # ..................{ EXPORTERS ~ microtubule            }..................
    @piperunner(categories=('Microtubule', 'Coherence',))
    def export_microtubule(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot the coherence of all cellular microtubules for the cell cluster at
        the last time step.
        '''

        # Prepare to export the microtubules plot.
        self._export_prep()

        pyplot.figure()
        ax = pyplot.subplot(111)

        plotutil.mem_quiver(
            self._phase.sim.mtubes.mtubes_x,
            self._phase.sim.mtubes.mtubes_y,
            ax,
            self._phase.cells,
            self._phase.p,
        )

        ax.set_xlabel('X-Distance [um]')
        ax.set_ylabel('Y-Distance [um]')
        ax.set_title('Microtubule arrangement in cells')

        # Export this plot to disk and/or display.
        self._export(basename='Final_Microtubules')

    # ..................{ EXPORTERS ~ pump                   }..................
    # @piperunner(
    #     categories=('Ion Pump', 'Density Factor',),
    #     requirements={piperunreq.ELECTROOSMOSIS,},
    # )
    # def export_pump_density(self, conf: SimConfVisualCellsListItem) -> None:
    #     '''
    #     Plot all cell membrane ion pump density factors for the cell cluster at
    #     the last time step.
    #     '''
    #
    #     # Prepare to export the current plot.
    #     self._export_prep()
    #
    #     plotutil.plotMemData(
    #         self._phase.cells, self._phase.p,
    #         zdata=self._phase.sim.rho_pump,
    #         clrmap=self._phase.p.default_cm,
    #     )
    #     pyplot.xlabel('Spatial Dimension [um]')
    #     pyplot.ylabel('Spatial Dimension [um]')
    #     pyplot.title('Membrane ion pump density factor')
    #
    #     # Export this plot to disk and/or display.
    #     self._export(basename='final_pumps_2D')


    @piperunner(categories=('Ion Pump', 'Pump Rate', 'Na-K-ATPase',))
    def export_pump_nakatpase(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all cell membrane Na-K-ATPase pump rates for the cell cluster at
        the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        pumpData = self._phase.sim.rate_NaKATP*1e9

        figPump, axPump, cbPump = plotutil.plotPrettyPolyData(
            pumpData,
            self._phase.sim, self._phase.cells, self._phase.p,
            number_cells=self._phase.p.enumerate_cells,
            clrmap=self._phase.p.default_cm,
        )

        axPump.set_title('Final Na/K-ATPase Pump Rate')
        axPump.set_xlabel('Spatial distance [um]')
        axPump.set_ylabel('Spatial distance [um]')
        cbPump.set_label('Pump Na+ Flux (nmol/m2*s)')

        # Export this plot to disk and/or display.
        self._export(basename='final_NaKPump_2D')

    # ..................{ EXPORTERS ~ pressure               }..................
    @piperunner(
        categories=('Pressure', 'Total',),
        requirements={piperunreq.PRESSURE_TOTAL,},
    )
    def export_pressure_total(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all **cellular pressure totals** (i.e., summations of all cellular
        mechanical and osmotic pressures) for the cell cluster at the last time
        step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figP, axP, cbP = plotutil.plotPolyData(
            self._phase.sim, self._phase.cells, self._phase.p,
            zdata=self._phase.sim.P_cells,
            number_cells=self._phase.p.enumerate_cells,
            clrmap=self._phase.p.default_cm,
            clrAutoscale=conf.is_color_autoscaled,
            clrMin=conf.color_min,
            clrMax=conf.color_max,
        )

        axP.set_title('Final Pressure in Cell Network')
        axP.set_xlabel('Spatial distance [um]')
        axP.set_ylabel('Spatial distance [um]')
        cbP.set_label('Pressure [Pa]')

        # Export this plot to disk and/or display.
        self._export(basename='final_P_2D_gj')

    # ..................{ EXPORTERS ~ voltage                }..................
    @piperunner(
        categories=('Voltage', 'Extracellular',),
        requirements={piperunreq.ECM,},
    )
    def export_voltage_extra(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all extracellular voltages for the cell cluster environment at the
        last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        vv = self._phase.sim.v_env.reshape(self._phase.cells.X.shape)

        pyplot.figure()
        pyplot.imshow(
            1e3*vv,
            origin='lower',
            extent=self._cells_extent,
            cmap=self._phase.p.default_cm,
        )
        pyplot.colorbar()
        pyplot.title('Environmental Voltage [mV]')

        # Export this plot to disk and/or display.
        self._export(basename='Final_environmental_V')


    @piperunner(
        categories=('Voltage', 'Polarity',),
        requirements={piperunreq.VOLTAGE_POLARITY,},
    )
    def export_voltage_polarity(
        self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all cellular voltage polarities for the cell cluster at the last
        time step.
        '''

        # Prepare to export the polarization plot.
        self._export_prep()

        pyplot.figure()

        # Plot a background Vmem mesh.
        fig, ax, cb = plotutil.plotPrettyPolyData(
            1000*self._phase.sim.vm_time[-1],
            self._phase.sim, self._phase.cells, self._phase.p,
            number_cells=self._phase.p.enumerate_cells,
            current_overlay=False,
            plotIecm=self._phase.p.IecmPlot,
            clrmap=self._phase.p.default_cm,
            clrAutoscale=conf.is_color_autoscaled,
            clrMin=conf.color_min,
            clrMax=conf.color_max,
        )

        # Calculate the Vmem polarity vectors.
        polm = self._phase.sim.vm - (
            self._phase.sim.vm_ave_time[-1][self._phase.cells.mem_to_cells])
        polx = polm*self._phase.cells.mem_vects_flat[:,2]
        poly = polm*self._phase.cells.mem_vects_flat[:,3]

        pcx = np.dot(
            self._phase.cells.M_sum_mems,
            polx*self._phase.cells.mem_sa) / self._phase.cells.cell_sa
        pcy = np.dot(
            self._phase.cells.M_sum_mems,
            poly*self._phase.cells.mem_sa) / self._phase.cells.cell_sa

        plotutil.cell_quiver(pcx, pcy, ax, self._phase.cells, self._phase.p)

        ax.set_xlabel('X-Distance [um]')
        ax.set_ylabel('Y-Distance [um]')
        ax.set_title('Cell Vmem polarity')

        # Export this plot to disk and/or display.
        self._export(basename='Final_Polarity')

    # ..................{ EXPORTERS ~ voltage : vmem         }..................
    @piperunner(
        categories=('Voltage', 'Transmembrane', 'Actual',))
    def export_voltage_membrane(self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all transmembrane voltages (Vmem) for the cell cluster at the last
        time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figV, axV, cbV = plotutil.plotPrettyPolyData(
            1000*self._phase.sim.vm_time[-1],
            self._phase.sim, self._phase.cells, self._phase.p,
            number_cells=self._phase.p.enumerate_cells,
            current_overlay=False,
            plotIecm=self._phase.p.IecmPlot,
            clrmap=self._phase.p.default_cm,
            clrAutoscale=conf.is_color_autoscaled,
            clrMin=conf.color_min,
            clrMax=conf.color_max,
        )

        figV.suptitle('Final Vmem', fontsize=14, fontweight='bold')
        axV.set_xlabel('Spatial distance [um]')
        axV.set_ylabel('Spatial distance [um]')
        cbV.set_label('Voltage mV')

        # Export this plot to disk and/or display.
        self._export(basename='final_Vmem_2D')


    @piperunner(
        categories=('Voltage', 'Transmembrane', 'Average',))
    def export_voltage_membrane_average(
        self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot the averages of all transmembrane voltages (Vmem) for the cell
        cluster at the last time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figVa, axVa, cbVa = plotutil.plotPolyData(
            self._phase.sim, self._phase.cells, self._phase.p,
            zdata=1000*self._phase.sim.vm_ave,
            number_cells=self._phase.p.enumerate_cells,
            current_overlay=False,
            plotIecm=self._phase.p.IecmPlot,
            clrmap=self._phase.p.default_cm,
            clrAutoscale=conf.is_color_autoscaled,
            clrMin=conf.color_min,
            clrMax=conf.color_max,
        )

        # axVa.quiver(
        #     p.um*cells.cell_centres[:,0],
        #     p.um*cells.cell_centres[:,1], sim.pol_cell_x, sim.pol_cell_y)

        figVa.suptitle('Final Average Vmem', fontsize=14, fontweight='bold')
        axVa.set_xlabel('Spatial distance [um]')
        axVa.set_ylabel('Spatial distance [um]')
        cbVa.set_label('Voltage [mV]')

        # Export this plot to disk and/or display.
        self._export(basename='final_AverageVmem_2D')


    @piperunner(
        categories=('Voltage', 'Transmembrane', 'GHK',),
        requirements={piperunreq.VOLTAGE_MEMBRANE_GHK,},
    )
    def export_voltage_membrane_ghk(
        self, conf: SimConfVisualCellsListItem) -> None:
        '''
        Plot all transmembrane voltages (Vmem) calculated by the
        Goldman-Hodgkin-Katz (GHK) equation for the cell cluster at the last
        time step.
        '''

        # Prepare to export the current plot.
        self._export_prep()

        figV_ghk, axV_ghk, cbV_ghk = plotutil.plotPolyData(
            self._phase.sim, self._phase.cells, self._phase.p,
            zdata=1000*self._phase.sim.vm_GHK_time[-1],
            number_cells=self._phase.p.enumerate_cells,
            current_overlay=False,
            plotIecm=self._phase.p.IecmPlot,
            clrmap=self._phase.p.default_cm,
            clrAutoscale=conf.is_color_autoscaled,
            clrMin=conf.color_min,
            clrMax=conf.color_max,
        )

        figV_ghk.suptitle(
            'Final Vmem using Goldman Equation', fontsize=14, fontweight='bold')
        axV_ghk.set_xlabel('Spatial distance [um]')
        axV_ghk.set_ylabel('Spatial distance [um]')
        cbV_ghk.set_label('Voltage [mV]')

        # Export this plot to disk and/or display.
        self._export(basename='final_Vmem_GHK_2D')
