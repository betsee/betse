#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**Post-simulation single-cell plot pipelines** (i.e., pipelines plotting
simulated data of a single cell of the cell cluster).
'''

# ....................{ IMPORTS                           }....................
import numpy as np
from betse.exceptions import BetseSimConfException
from betse.science.config.export.visual.confexpvisplot import (
    SimConfExportPlotCell)
from betse.science.phase.phasecls import SimPhase
from betse.science.phase.require import phasereqs
from betse.science.pipe.export.plot.pipeexpplotabc import SimPipeExportPlotABC
from betse.science.pipe.piperun import piperunner
from betse.science.visual.plot import plotutil
from betse.util.type.descriptor.descs import classproperty_readonly
from betse.util.type.types import type_check, SequenceTypes
from matplotlib import pyplot

# ....................{ SUBCLASSES                        }....................
class SimPipeExportPlotCell(SimPipeExportPlotABC):
    '''
    **Post-simulation single-cell plot pipeline** (i.e., object iteratively
    displaying and/or saving all plots specific to a single cell of the cell
    cluster, produced after initialization and simulation as specified by the
    current simulation configuration).
    '''

    # ..................{ SUPERCLASS ~ properties           }..................
    @classproperty_readonly
    def _NOUN_SINGULAR(cls) -> str:
        return 'single-cell plot'

    # ..................{ INITIALIZERS                      }..................
    @type_check
    def init(self, phase: SimPhase) -> None:

        # Initialize our superclass for the current call to the _run() method.
        super().init(phase)

        # If this index is not that of an actual cell, raise an exception.
        if phase.p.visual.single_cell_index not in phase.cells.cell_i:
            raise BetseSimConfException(
                'Plot cell index {} invalid '
                '(i.e., not in range [{}, {}]).'.format(
                    phase.p.visual.single_cell_index,
                    phase.cells.cell_i[0],
                    phase.cells.cell_i[-1]))

    # ..................{ SUPERCLASS                        }..................
    @type_check
    def iter_runners_conf(self, phase: SimPhase) -> SequenceTypes:
        return phase.p.plot.plots_cell_after_sim

    # ..................{ EXPORTERS ~ cell : current        }..................
    @piperunner(
        categories=('Current Density', 'Transmembrane',),
        requirements=phasereqs.ELECTRIC_CURRENT_MEMBRANE,
    )
    def export_currents_membrane(
        self, phase: SimPhase, conf: SimConfExportPlotCell) -> None:
        '''
        Plot all transmembrane current densities for the single cell indexed by
        the current simulation configuration over all sampled time steps.
        '''

        # Prepare to export the current plot.
        self._export_prep(phase)

        pyplot.figure()
        axI = pyplot.subplot(111)

        if phase.p.is_ecm:
            # Total cell current storage vector.
            Imem = []

            # Indices of all membranes for this cell.
            mems_for_plotcell = phase.cells.cell_to_mems[
                phase.p.visual.single_cell_index]

            for t in range(len(phase.sim.time)):
                memArray = phase.sim.I_mem_time[t]

                # X and Y components of the net current at each membrane.
                Ixo = (
                    memArray[mems_for_plotcell] *
                    phase.cells.mem_vects_flat[mems_for_plotcell,2] *
                    phase.cells.mem_sa[mems_for_plotcell]
                )
                Iyo = (
                    memArray[mems_for_plotcell] *
                    phase.cells.mem_vects_flat[mems_for_plotcell,3] *
                    phase.cells.mem_sa[mems_for_plotcell]
                )

                # X and Y components of the net current over all membranes,
                # implicitly accounting for in-out directionality.
                Ix = np.sum(Ixo)
                Iy = np.sum(Iyo)

                # Current density, obtained by dividing the total magnitude of
                # the net current by this cell's surface area.
                Io = np.sqrt(Ix**2 + Iy**2) / phase.cells.cell_sa[
                    phase.p.visual.single_cell_index]
                Imem.append(100*Io)
        else:
            Imem = [
                100*memArray[phase.p.visual.single_cell_index]
                for memArray in phase.sim.I_mem_time
            ]

        axI.plot(phase.sim.time, Imem)
        axI.set_xlabel('Time [s]')
        axI.set_ylabel('Current density [uA/cm2]')
        axI.set_title(
            'Transmembrane current density for cell {}'.format(
                phase.p.visual.single_cell_index))

        # Export this plot to disk and/or display.
        self._export(phase=phase, basename='Imem_time')

    # ..................{ EXPORTERS ~ cell : deform         }..................
    @piperunner(
        categories=('Deformation', 'Total',),
        requirements=phasereqs.DEFORM,
    )
    def export_deform_total(
        self, phase: SimPhase, conf: SimConfExportPlotCell) -> None:
        '''
        Plot all **total cellular displacements** (i.e., summations of all
        cellular deformations due to galvanotropic and osmotic pressure body
        forces) for the single cell indexed by the current simulation
        configuration over all sampled time steps.
        '''

        # Prepare to export the current plot.
        self._export_prep(phase)

        # Extract time-series deformation data for the plot cell.
        dx = np.asarray([
            arr[phase.p.visual.single_cell_index]
            for arr in phase.sim.dx_cell_time])
        dy = np.asarray([
            arr[phase.p.visual.single_cell_index]
            for arr in phase.sim.dy_cell_time])

        # Get the total magnitude.
        disp = np.sqrt(dx**2 + dy**2)

        pyplot.figure()
        axD = pyplot.subplot(111)
        axD.plot(phase.sim.time, phase.p.um*disp)
        axD.set_xlabel('Time [s]')
        axD.set_ylabel('Displacement [um]')
        axD.set_title('Displacement of cell {}'.format(
            phase.p.visual.single_cell_index))

        # Export this plot to disk and/or display.
        self._export(phase=phase, basename='Displacement_time')

    # ..................{ EXPORTERS ~ cell : ion            }..................
    @piperunner(
        categories=('Ion Concentration', 'M anion'),
        requirements=phasereqs.ION_M_ANION,
    )
    def export_ion_m_anion(
        self, phase: SimPhase, conf: SimConfExportPlotCell) -> None:
        '''
        Plot all M anion (i.e., M-) ion concentrations for the single cell
        indexed by the current simulation configuration over all sampled time
        steps.
        '''

        # Prepare to export the current plot.
        self._export_prep(phase)

        figConcsM, axConcsM = plotutil.plotSingleCellCData(
            phase.sim.cc_time,
            phase.sim.time,
            phase.sim.iM,
            phase.p.visual.single_cell_index,
            fig=None, ax=None,
            lncolor='r',
            ionname='M-',
        )
        axConcsM.set_title(
            'M Anion concentration in cell {}'.format(phase.p.visual.single_cell_index))

        # Export this plot to disk and/or display.
        self._export(phase=phase, basename='concM_time')


    @piperunner(
        categories=('Ion Concentration', 'Potassium'),
        requirements=phasereqs.ION_POTASSIUM,
    )
    def export_ion_potassium(
        self, phase: SimPhase, conf: SimConfExportPlotCell) -> None:
        '''
        Plot all potassium (i.e., K+) ion concentrations for the single cell
        indexed by the current simulation configuration over all sampled time
        steps.
        '''

        # Prepare to export the current plot.
        self._export_prep(phase)

        figConcsK, axConcsK = plotutil.plotSingleCellCData(
            phase.sim.cc_time,
            phase.sim.time,
            phase.sim.iK,
            phase.p.visual.single_cell_index,
            fig=None, ax=None,
            lncolor='b',
            ionname='K+',
        )
        axConcsK.set_title(
            'Potassium concentration in cell {}'.format(
                phase.p.visual.single_cell_index))

        # Export this plot to disk and/or display.
        self._export(phase=phase, basename='concK_time')


    @piperunner(
        categories=('Ion Concentration', 'Sodium'),
        requirements=phasereqs.ION_SODIUM,
    )
    def export_ion_sodium(
        self, phase: SimPhase, conf: SimConfExportPlotCell) -> None:
        '''
        Plot all sodium (i.e., Na+) ion concentrations for the single cell
        indexed by the current simulation configuration over all sampled time
        steps.
        '''

        # Prepare to export the current plot.
        self._export_prep(phase)

        # Plot cell sodium concentration versus time.
        figConcsNa, axConcsNa = plotutil.plotSingleCellCData(
            phase.sim.cc_time,
            phase.sim.time,
            phase.sim.iNa,
            phase.p.visual.single_cell_index,
            fig=None,
            ax=None,
            lncolor='g',
            ionname='Na+',
        )
        axConcsNa.set_title(
            'Sodium concentration in cell {}'.format(phase.p.visual.single_cell_index))

        # Export this plot to disk and/or display.
        self._export(phase=phase, basename='concNa_time')

    # ..................{ EXPORTERS ~ cell : ion : calcium  }..................
    @piperunner(
        categories=('Ion Concentration', 'Calcium', 'Cellular',),
        requirements=phasereqs.ION_CALCIUM,
    )
    def export_ion_calcium(
        self, phase: SimPhase, conf: SimConfExportPlotCell) -> None:
        '''
        Plot all calcium (i.e., Ca2+) ion concentrations for the single cell
        indexed by the current simulation configuration over all sampled time
        steps.
        '''

        # Prepare to export the current plot.
        self._export_prep(phase)

        figA, axA = plotutil.plotSingleCellCData(
            phase.sim.cc_time,
            phase.sim.time,
            phase.sim.iCa,
            phase.p.visual.single_cell_index,
            fig=None,
            ax=None,
            lncolor='g',
            ionname='Ca2+ cell',
        )
        axA.set_title('Cytosolic Ca2+ in cell {}'.format(
            phase.p.visual.single_cell_index))

        # Export this plot to disk and/or display.
        self._export(phase=phase, basename='cytosol_Ca_time')


    @piperunner(
        categories=('Ion Concentration', 'Calcium', 'Endoplasmic Reticulum',),
        requirements=phasereqs.ION_CALCIUM_DYNAMICS,
    )
    def export_ion_calcium_er(
        self, phase: SimPhase, conf: SimConfExportPlotCell) -> None:
        '''
        Plot all calcium (i.e., Ca2+) ion concentrations for the endoplasmic
        reticulum of the single cell indexed by the current simulation
        configuration over all sampled time steps.
        '''

        # Prepare to export the current plot.
        self._export_prep(phase)

        # One-dimensional Numpy array of all ion concentrations for the
        # endoplasmic reticulum of this cell over all sampled time steps.
        times_cell_ion_calcium_er = [
            cells_ion_calcium_er[phase.p.visual.single_cell_index]
            for cells_ion_calcium_er in phase.sim.endo_retic.Ca_er_time
        ]

        # Plot this array.
        pyplot.figure()
        pyplot.plot(phase.sim.time, times_cell_ion_calcium_er)

        # Export this plot to disk and/or display.
        self._export(phase=phase, basename='CaER')


    @piperunner(
        categories=('Voltage', 'Transmembrane', 'Endoplasmic Reticulum',),
        requirements=phasereqs.ION_CALCIUM_DYNAMICS,
    )
    def export_voltage_membrane_er(
        self, phase: SimPhase, conf: SimConfExportPlotCell) -> None:
        '''
        Plot all transmembrane voltages across the endoplasmic reticulum of the
        single cell indexed by the current simulation configuration over all
        sampled time steps.
        '''

        # Prepare to export the current plot.
        self._export_prep(phase)

        # One-dimensional Numpy array of all ion concentrations for the
        # endoplasmic reticulum of this cell over all sampled time steps.
        times_cell_vmem_er = [
            cells_vmem_er[phase.p.visual.single_cell_index]
            for cells_vmem_er in phase.sim.endo_retic.ver_time
        ]

        # Plot this array.
        pyplot.figure()
        pyplot.plot(phase.sim.time, times_cell_vmem_er)

        # Export this plot to disk and/or display.
        self._export(phase=phase, basename='VER')

    # ..................{ EXPORTERS ~ pressure              }..................
    @piperunner(
        categories=('Pressure', 'Osmotic',),
        requirements=phasereqs.PRESSURE_OSMOTIC,
    )
    def export_pressure_osmotic(
        self, phase: SimPhase, conf: SimConfExportPlotCell) -> None:
        '''
        Plot the osmotic cellular pressure for the single cell indexed by the
        current simulation configuration over all sampled time steps.
        '''

        # Prepare to export the current plot.
        self._export_prep(phase)

        p_osmo = tuple(
            arr[phase.p.visual.single_cell_index]
            for arr in phase.sim.osmo_P_delta_time)

        pyplot.figure()
        axOP = pyplot.subplot(111)
        axOP.plot(phase.sim.time, p_osmo)
        axOP.set_xlabel('Time [s]')
        axOP.set_ylabel('Osmotic Pressure [Pa]')
        axOP.set_title(
            'Osmotic pressure in cell {}'.format(phase.p.visual.single_cell_index))

        # Export this plot to disk and/or display.
        self._export(phase=phase, basename='OsmoticP_time')


    @piperunner(
        categories=('Pressure', 'Total',),
        requirements=phasereqs.PRESSURE_TOTAL,
    )
    def export_pressure_total(
        self, phase: SimPhase, conf: SimConfExportPlotCell) -> None:
        '''
        Plot the **total cellular pressure** (i.e., summation of the cellular
        mechanical and osmotic pressure) for the single cell indexed by the
        current simulation configuration over all sampled time steps.
        '''

        # Prepare to export the current plot.
        self._export_prep(phase)

        p_hydro = tuple(
            arr[phase.p.visual.single_cell_index]
            for arr in phase.sim.P_cells_time)

        pyplot.figure()
        axOP = pyplot.subplot(111)
        axOP.plot(phase.sim.time, p_hydro)
        axOP.set_xlabel('Time [s]')
        axOP.set_ylabel('Total Pressure [Pa]')
        axOP.set_title(
            'Total pressure in cell {}'.format(phase.p.visual.single_cell_index))

        # Export this plot to disk and/or display.
        self._export(phase=phase, basename='P_time')

    # ..................{ EXPORTERS ~ pump                  }..................
    #FIXME: Actually plot the average of these rates. Currently, this method
    #only plots rates for a single arbitrarily selected membrane of this cell.
    @piperunner(
        categories=('Pump Rate', 'Na-K-ATPase',),
        requirements=phasereqs.PUMP_NAKATPASE,
    )
    def export_pump_nakatpase(
        self, phase: SimPhase, conf: SimConfExportPlotCell) -> None:
        '''
        Plot the averages of the Na-K-ATPase membrane pump rates for the single
        cell indexed by the current simulation configuration over all time
        steps.
        '''

        # Prepare to export the current plot.
        self._export_prep(phase)

        pyplot.figure()
        axNaK = pyplot.subplot()

        if phase.p.is_ecm:
            #FIXME: Generalize to average across all membranes of this cell.
            pump_array_index = phase.cells.cell_to_mems[
                phase.p.visual.single_cell_index][0]
        else:
            pump_array_index = phase.p.visual.single_cell_index

        pump_rate = tuple(
            pump_array[pump_array_index]
            for pump_array in phase.sim.rate_NaKATP_time)

        axNaK.plot(phase.sim.time, pump_rate)
        axNaK.set_xlabel('Time [s]')
        axNaK.set_ylabel('Pumping Rate of Na+ Out of Cell [mol/(m2 s)]')
        axNaK.set_title(
            'Rate of NaK-ATPase pump in cell: {}'.format(
                phase.p.visual.single_cell_index))

        # Export this plot to disk and/or display.
        self._export(phase=phase, basename='NaKATPase_rate')

    # ..................{ EXPORTERS ~ voltage               }..................
    #FIXME: Actually plot the average of these Vmems. Currently, this method
    #only plots Vmems for a single arbitrarily selected membrane of this cell.
    @piperunner(categories=('Voltage', 'Transmembrane', 'Average',))
    def export_voltage_membrane(
        self, phase: SimPhase, conf: SimConfExportPlotCell) -> None:
        '''
        Plot the averages of all transmembrane voltages for the single cell
        indexed by the current simulation configuration over all sampled time
        steps.
        '''

        # Prepare to export the current plot.
        self._export_prep(phase)

        #FIXME: Generalize to average across all membranes of this cell.
        mem_i = phase.cells.cell_to_mems[phase.p.visual.single_cell_index][0]

        # Plot single cell Vmem vs time.
        figVt, axVt = plotutil.plotSingleCellVData(
            phase.sim,
            mem_i,
            phase.p,
            fig=None, ax=None,
            lncolor='k',
        )
        axVt.set_title(
            'Voltage (Vmem) in cell {}'.format(phase.p.visual.single_cell_index))

        # Export this plot to disk and/or display.
        self._export(phase=phase, basename='Vmem_time')


    #FIXME: Actually plot the average of these Vmems. Currently, this method
    #only plots Vmems for a single arbitrarily selected membrane of this cell.
    @piperunner(categories=(
        'Voltage', 'Transmembrane', 'Fast Fourier Transform (FFT)',))
    def export_voltage_membrane_fft(
        self, phase: SimPhase, conf: SimConfExportPlotCell) -> None:
        '''
        Plot the fast Fourier transform (FFT) of the averages of all
        transmembrane voltages for the single cell indexed by the current
        simulation configuration over all sampled time steps.
        '''

        # Prepare to export the current plot.
        self._export_prep(phase)

        #FIXME: Generalize to average across all membranes of this cell.
        mem_i = phase.cells.cell_to_mems[phase.p.visual.single_cell_index][0]

        figFFT, axFFT = plotutil.plotFFT(
            phase.sim.time, phase.sim.vm_time, mem_i, lab="Power")
        axFFT.set_title('Fourier transform of Vmem in cell {}'.format(
            phase.p.visual.single_cell_index))

        # Export this plot to disk and/or display.
        self._export(phase=phase, basename='Vmem_FFT_time')
