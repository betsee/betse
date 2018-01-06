#!/usr/bin/env python3
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level facilities for **pipelining** (i.e., iteratively exporting)
post-simulation files containing raw simulation data, typically as plaintext
spreadsheets in comma-separated value (CSV) format.
'''

#FIXME: Log each attempt to create a CSV file in this submodule.

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.exceptions import BetseMethodUnimplementedException
from betse.lib.numpy import nparray
from betse.science.export import expmath
from betse.science.simulate.pipe.pipeabc import SimPipeExportABC
from betse.science.simulate.pipe.piperun import piperunner
from betse.science.simulate.simphase import SimPhase, SimPhaseKind
from betse.science.visual.plot.plotutil import cell_ave
from betse.util.path import dirs, pathnames
from betse.util.type.call.memoizers import property_cached
from betse.util.type.mapping.mapcls import OrderedArgsDict
from betse.util.type.types import type_check, IterableTypes, SequenceTypes
from numpy import ndarray

# ....................{ SUBCLASSES                         }....................
class SimPipelinerExportCSV(SimPipeExportABC):
    '''
    **Post-simulation CSV export pipeline** (i.e., class iteratively creating
    all post-simulation plaintext spreadsheets in comma-separated value (CSV)
    format requested by the current simulation configuration).
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, label_singular='CSV', **kwargs)

    # ..................{ SUPERCLASS                         }..................
    #FIXME: Implement this property *AFTER* generalizing the configuration file
    #format to support CSV pipelines.
    @property
    def _runners_conf(self) -> IterableTypes:

        raise BetseMethodUnimplementedException()

    # ..................{ EXPORTERS ~ cell                   }..................
    @piperunner(categories=('Single Cell', 'Raw Data'))
    def export_cell_raw(self) -> None:
        '''
        Save a plaintext file in comma-separated value (CSV) format containing
        several cell-specific time series (e.g., ion concentrations, membrane
        voltages, voltage-gated ion channel pump rates) for the single cell
        whose index indexed by the ``plot cell index`` entry for the current
        simulation configuration.
        '''

        # 0-based index of the cell to serialize time data for.
        cell_index = self._phase.p.plot_cell

        # Sequence of key-value pairs containing all simulation data to be
        # exported for this cell, suitable for passing to the
        # OrderedArgsDict.__init__() method calleb below.
        csv_column_name_values = []

        # One-dimensional Numpy array of null data of the required length,
        # suitable for use as CSV column data for columns whose corresponding
        # simulation feature (e.g., deformations) is disabled.
        column_data_empty = np.zeros(len(self._phase.sim.time))

        # ................{ TIME STEPS                       }..................
        csv_column_name_values.extend(('time_s', self._phase.sim.time))

        # ................{ VMEM                             }..................
        csv_column_name_values.extend(('Vmem_mV', self._cell_times_vmems))

        # ................{ VMEM ~ goldman                   }..................
        if self._phase.p.GHK_calc:
            vm_goldman = expmath.upscale_units_milli([
                vm_GHK_time_cells[cell_index]
                for vm_GHK_time_cells in self._phase.sim.vm_GHK_time])
        else:
            vm_goldman = column_data_empty

        csv_column_name_values.extend(('Goldman_Vmem_mV', vm_goldman))

        # ................{ Na K PUMP RATE                   }..................
        if self._phase.p.is_ecm:
            pump_rate = [
                pump_array[self._phase.cells.cell_to_mems[cell_index][0]]
                for pump_array in self._phase.sim.rate_NaKATP_time]
        else:
            pump_rate = [
                pump_array[cell_index]
                for pump_array in self._phase.sim.rate_NaKATP_time]

        csv_column_name_values.extend((
            'NaK-ATPase_Rate_mol/m2s', pump_rate))

        # ................{ ION CONCENTRATIONS               }..................
        # Create the header starting with cell concentrations.
        for i in range(len(self._phase.sim.ionlabel)):
            csv_column_name = 'cell_{}_mmol/L'.format(
                self._phase.sim.ionlabel[i])
            cc_m = [arr[i][cell_index] for arr in self._phase.sim.cc_time]
            csv_column_name_values.extend((csv_column_name, cc_m))

        # ................{ MEMBRANE PERMEABILITIES          }..................
        # create the header starting with membrane permeabilities
        for i in range(len(self._phase.sim.ionlabel)):
            if self._phase.p.is_ecm:
                dd_m = [
                    arr[i][self._phase.cells.cell_to_mems[cell_index][0]]
                    for arr in self._phase.sim.dd_time
                ]
            else:
                dd_m = [arr[i][cell_index] for arr in self._phase.sim.dd_time]

            csv_column_name = 'Dm_{}_m2/s'.format(self._phase.sim.ionlabel[i])
            csv_column_name_values.extend((csv_column_name, dd_m))

        # ................{ TRANSMEMBRANE CURRENTS           }..................
        if self._phase.p.is_ecm:
            Imem = [
                memArray[self._phase.cells.cell_to_mems[cell_index][0]]
                for memArray in self._phase.sim.I_mem_time]
        else:
            Imem = [
                memArray[cell_index] for memArray in self._phase.sim.I_mem_time]

        csv_column_name_values.extend(('I_A/m2', Imem))

        # ................{ HYDROSTATIC PRESSURE             }..................
        p_hydro = [arr[cell_index] for arr in self._phase.sim.P_cells_time]
        csv_column_name_values.extend(('HydroP_Pa', p_hydro))

        # ................{ OSMOTIC PRESSURE                 }..................
        if self._phase.p.deform_osmo:
            p_osmo = [
                arr[cell_index] for arr in self._phase.sim.osmo_P_delta_time]
        else:
            p_osmo = column_data_empty

        csv_column_name_values.extend(('OsmoP_Pa', p_osmo))

        # ................{ DEFORMATION                      }..................
        if (
            self._phase.p.deformation and
            self._phase.kind is SimPhaseKind.SIM
        ):
            # Extract time-series deformation data for the plot cell:
            dx = nparray.from_iterable([
                arr[cell_index] for arr in self._phase.sim.dx_cell_time])
            dy = nparray.from_iterable([
                arr[cell_index] for arr in self._phase.sim.dy_cell_time])

            # Get the total magnitude.
            disp = expmath.upscale_coordinates(np.sqrt(dx ** 2 + dy ** 2))
        else:
            disp = column_data_empty

        csv_column_name_values.extend(('Displacement_um', disp))

        # ................{ CSV EXPORT                       }..................
        # Ordered dictionary mapping from CSV column names to data arrays.
        csv_column_name_to_values = OrderedArgsDict(*csv_column_name_values)

        # Absolute path of the CSV file to export.
        csv_filename = pathnames.join(
            self._phase.save_dirname, 'ExportedData.csv')

        # Export this data to this CSV file.
        nparray.write_csv(
            filename=csv_filename,
            column_name_to_values=csv_column_name_to_values)

    # ..................{ EXPORTERS ~ cell : vmem            }..................
    #FIXME: This exporter appears to currently be broken for the non-ECM
    #case, despite working for the ECM case. The exception being raised is:
    #
    #    BetseSequenceException: Column "FFT_Vmem" length 9 differs from length
    #    5 of prior columns.

    @piperunner(categories=('Single Cell', 'Vmem FFT'))
    def export_cell_vmem_fft(self) -> None:
        '''
        Save a plaintext file in comma-separated value (CSV) format containing
        the finite Fourier transform (FFT) of all transmembrane voltages for all
        sampled time steps spatially situated at the centre of the cell indexed
        by the ``plot cell index`` entry for the current simulation
        configuration.
        '''

        #FIXME: Remove this after correcting this exporter.
        raise BetseMethodUnimplementedException()

        # Number of sampled time steps.
        sample_size = len(self._phase.sim.time)

        # Time in seconds between each sampled time step.
        sample_spacing = self._phase.sim.time[1] - self._phase.sim.time[0]

        cell_data = (1/sample_size) * (
            self._cell_times_vmems - np.mean(self._cell_times_vmems))

        # FFT of voltage.
        f_axis = np.fft.rfftfreq(sample_size, d=sample_spacing)
        fft_data_o = np.fft.rfft(cell_data)

        #FIXME: Numpy already provides a function for obtaining the magnitude of
        #a complex array: np.absolute(). The following math reduces to simply:
        #
        #    fft_data = np.absolute(fft_data_o)
        fft_data = np.sqrt(np.real(fft_data_o)**2 + np.imag(fft_data_o)**2)
        print('f_axis: {}'.format(f_axis))
        print('fft_data: {}'.format(fft_data))

        # Ordered dictionary mapping from CSV column names to data arrays.
        csv_column_name_to_values = OrderedArgsDict(
            'frequency_Hz', f_axis,
            'FFT_Vmem', fft_data,
        )

        # Absolute path of the CSV file to export.
        csv_filename = pathnames.join(
            self._phase.save_dirname, 'ExportedData_FFT.csv')

        # Export this data to this CSV file.
        nparray.write_csv(
            filename=csv_filename,
            column_name_to_values=csv_column_name_to_values)

    # ..................{ EXPORTERS ~ cells : vmem           }..................
    @piperunner(categories=('Cell Cluster', 'Transmembrane Voltages'))
    def export_cells_vmems(self) -> None:
        '''
        Save one plaintext file in comma-separated value (CSV) format for each
        sequence of transmembrane voltages spatially situated at cell centres
        for each sampled time step of the current simulation phase.
        '''

        # One-dimensional Numpy array of all upscaled cell voltages.
        cells_times_vmems = expmath.upscale_units_milli(
            self._phase.sim.vm_ave_time)

        # Save all membrane voltages for this time step to a unique CSV file.
        self._export_cells_times_data(
            cells_times_data=cells_times_vmems,
            csv_column_name='Vmem [mV]',
            csv_dir_basename='Vmem2D_TextExport',
            csv_basename_prefix='Vmem2D_',
        )

    # ..................{ PRIVATE ~ exporters                }..................
    @type_check
    def _export_cells_times_data(
        self,
        cells_times_data: SequenceTypes,
        csv_column_name: str,
        csv_dir_basename: str,
        csv_basename_prefix: str,
    ) -> None:
        '''
        Save a plaintext file in comma-separated value (CSV) format containing
        the passed two-dimensional Numpy array of arbitrary simulation data
        spatially situated at cell centres for all sampled time steps of the
        current simulation phase.

        Parameters
        ----------------------------
        cells_times_data : ndarray
            Two-dimensional Numpy array of arbitrary simulation data spatially
            situated at the centres of all cells for all sampled time steps.
        csv_column_name : str
            Name of the column containing this data in all CSV-formatted files
            exported by this method.
        csv_dir_basename : str
            Basename of the directory containing all CSV-formatted files
            exported by this method.
        csv_basename_prefix : str
            Substring prefixing the basenames of all CSV-formatted files
            exported by this method.
        '''

        # Absolute path of the directory containing all CSV-formatted files
        # exported by this method and creating this directory if needed.
        csv_dirname = pathnames.join(self._phase.save_dirname, csv_dir_basename)
        dirs.make_unless_dir(csv_dirname)

        # One-dimensional Numpy arrays of the X and Y coordinates (respectively)
        # of the centres of all cells.
        cell_centres_x = expmath.upscale_coordinates(
            self._phase.cells.cell_centres[:,0])
        cell_centres_y = expmath.upscale_coordinates(
            self._phase.cells.cell_centres[:,1])

        # For the 0-based index of each sampled time step...
        for time_step in range(len(self._phase.sim.time)):
            # Basename of the CSV-formatted file exported for this time step.
            csv_basename = '{}{}.csv'.format(csv_basename_prefix, time_step)

            # Absolute path of this file.
            csv_filename = pathnames.join(csv_dirname, csv_basename)

            # Ordered dictionary mapping from CSV column names to data arrays.
            csv_column_name_to_values = OrderedArgsDict(
                'x [um]', cell_centres_x,
                'y [um]', cell_centres_y,
                csv_column_name, cells_times_data[time_step],
            )

            # Export this data to this CSV file.
            nparray.write_csv(
                filename=csv_filename,
                column_name_to_values=csv_column_name_to_values)

    # ..................{ PRIVATE ~ properties               }..................
    @property_cached
    def _cell_times_vmems(self) -> ndarray:
        '''
        One-dimensional Numpy array of all transmembrane voltages for each
        sampled time step spatially situated at the centre of the single cell
        indexed by the ``plot cell index`` entry for the current simulation
        configuration.
        '''

        # 0-based index of the cell to serialize time data for.
        cell_index = self._phase.p.plot_cell

        if self._phase.p.is_ecm:
            cell_times_vmems = []
            for vm_at_mem in self._phase.sim.vm_time:
                vm_t = expmath.upscale_units_milli(
                    cell_ave(self._phase.cells,vm_at_mem)[cell_index])
                cell_times_vmems.append(vm_t)
        else:
            cell_times_vmems = expmath.upscale_units_milli(
                self._phase.sim.vm_time)

        return nparray.from_iterable(cell_times_vmems)

# ....................{ PIPELINES                          }....................
@type_check
def pipeline(phase: SimPhase) -> None:
    '''
    Save all currently enabled data exports (e.g., spreadsheets in
    comma-separated value form) for the passed simulation phase.

    Parameters
    ----------------------------
    phase: SimPhase
        Current simulation phase.
    '''

    # Post-simulation CSV pipeline producing all such CSVs.
    pipeliner = SimPipelinerExportCSV(phase)

    if phase.p.exportData:
        pipeliner.export_cell_raw()

        #FIXME: Reenable after correcting for non-ECM usage.
        # pipeliner.export_cell_vmem_fft()

    if phase.p.exportData2D:
        pipeliner.export_cells_vmems()
