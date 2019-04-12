#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2019 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level facilities for **pipelining** (i.e., iteratively exporting)
post-simulation files containing raw simulation data, typically as plaintext
spreadsheets in comma-separated value (CSV) format.
'''

# ....................{ IMPORTS                           }....................
import numpy as np
from betse.lib.numpy import nparray, npcsv
from betse.science.config.export.confexpcsv import SimConfExportCSV
from betse.science.math import mathunit
from betse.science.phase.phasecls import SimPhase
from betse.science.enum.enumphase import SimPhaseKind
from betse.science.phase.require import phasereqs
from betse.science.pipe.export.pipeexpabc import SimPipeExportABC
from betse.science.pipe.piperun import piperunner
from betse.science.visual.plot.plotutil import cell_ave
# from betse.util.io.log import logs
from betse.util.path import dirs, pathnames
from betse.util.type.descriptor.descs import classproperty_readonly
from betse.util.type.iterable.mapping.mapcls import OrderedArgsDict
from betse.util.type.types import (
    type_check, NumpyArrayType, SequenceTypes, StrOrNoneTypes,)

# ....................{ SUBCLASSES                        }....................
class SimPipeExportCSVs(SimPipeExportABC):
    '''
    **Post-simulation CSV export pipeline** (i.e., class iteratively creating
    all post-simulation plaintext spreadsheets in comma-separated value (CSV)
    format requested by the current simulation configuration).
    '''

    # ..................{ SUPERCLASS ~ properties           }..................
    @classproperty_readonly
    def _NOUN_SINGULAR(cls) -> str:
        return 'CSV file'

    # ..................{ SUPERCLASS                        }..................
    @type_check
    def iter_runners_conf(self, phase: SimPhase) -> SequenceTypes:
        return phase.p.csv.csvs_after_sim

    @type_check
    def _is_enabled(self, phase: SimPhase) -> bool:
        return phase.p.csv.is_after_sim_save

    # ..................{ EXPORTERS ~ cell                  }..................
    #FIXME: The requirements list should (arguably) be refined from the
    #coarse-grained "SOLVER_FULL" requirement to the exact list of fine-grained
    #requirements required by this exporter.
    @piperunner(
        categories=('Single Cell', 'Time Series'),
        requirements=phasereqs.SOLVER_FULL,
    )
    def export_cell_series(
        self, phase: SimPhase, conf: SimConfExportCSV) -> None:
        '''
        Save a plaintext file in comma-separated value (CSV) format containing
        several cell-specific time series (e.g., ion concentrations, membrane
        voltages, voltage-gated ion channel pump rates) for the single cell
        indexed ``plot cell index`` in the current simulation configuration.
        '''

        # 0-based index of the cell to serialize time data for.
        cell_index = phase.p.visual.single_cell_index

        # Sequence of key-value pairs containing all simulation data to be
        # exported for this cell, suitable for passing to the
        # OrderedArgsDict.__init__() method calleb below.
        csv_column_name_values = []

        # One-dimensional Numpy array of null data of the required length,
        # suitable for use as CSV column data for columns whose corresponding
        # simulation feature (e.g., deformations) is disabled.
        column_data_empty = np.zeros(len(phase.sim.time))

        # ................{ TIME STEPS                      }..................
        csv_column_name_values.extend(('time_s', phase.sim.time))

        # ................{ VMEM                            }..................
        csv_column_name_values.extend(
            ('Vmem_mV', self._get_cell_times_vmems(phase)))

        # ................{ VMEM ~ goldman                  }..................
        if phase.p.GHK_calc:
            vm_goldman = mathunit.upscale_units_milli([
                vm_GHK_time_cells[cell_index]
                for vm_GHK_time_cells in phase.sim.vm_GHK_time])
        else:
            vm_goldman = column_data_empty

        csv_column_name_values.extend(('Goldman_Vmem_mV', vm_goldman))

        # ................{ Na K PUMP RATE                  }..................
        if phase.p.is_ecm:
            pump_rate = [
                pump_array[phase.cells.cell_to_mems[cell_index][0]]
                for pump_array in phase.sim.rate_NaKATP_time]
        else:
            pump_rate = [
                pump_array[cell_index]
                for pump_array in phase.sim.rate_NaKATP_time]

        csv_column_name_values.extend((
            'NaK-ATPase_Rate_mol/m2s', pump_rate))

        # ................{ ION CONCENTRATIONS              }..................
        # Create the header starting with cell concentrations.
        for i in range(len(phase.sim.ionlabel)):
            csv_column_name = 'cell_{}_mmol/L'.format(
                phase.sim.ionlabel[i])
            cc_m = [arr[i][cell_index] for arr in phase.sim.cc_time]
            csv_column_name_values.extend((csv_column_name, cc_m))

        # ................{ MEMBRANE PERMEABILITIES         }..................
        # Create the header starting with membrane permeabilities.
        for i in range(len(phase.sim.ionlabel)):
            if phase.p.is_ecm:
                dd_m = [
                    arr[i][phase.cells.cell_to_mems[cell_index][0]]
                    for arr in phase.sim.dd_time
                ]
            else:
                dd_m = [arr[i][cell_index] for arr in phase.sim.dd_time]

            csv_column_name = 'Dm_{}_m2/s'.format(phase.sim.ionlabel[i])
            csv_column_name_values.extend((csv_column_name, dd_m))

        # ................{ TRANSMEMBRANE CURRENTS          }..................
        if phase.p.is_ecm:
            Imem = [
                memArray[phase.cells.cell_to_mems[cell_index][0]]
                for memArray in phase.sim.I_mem_time]
        else:
            Imem = [
                memArray[cell_index]
                for memArray in phase.sim.I_mem_time
            ]

        csv_column_name_values.extend(('I_A/m2', Imem))

        # ................{ HYDROSTATIC PRESSURE            }..................
        p_hydro = [arr[cell_index] for arr in phase.sim.P_cells_time]
        csv_column_name_values.extend(('HydroP_Pa', p_hydro))

        # ................{ OSMOTIC PRESSURE                }..................
        if phase.p.deform_osmo:
            p_osmo = [
                arr[cell_index] for arr in phase.sim.osmo_P_delta_time]
        else:
            p_osmo = column_data_empty

        csv_column_name_values.extend(('OsmoP_Pa', p_osmo))

        # ................{ DEFORMATION                     }..................
        if (
            phase.p.deformation and
            phase.kind is SimPhaseKind.SIM
        ):
            # Extract time-series deformation data for the plot cell:
            dx = nparray.from_iterable([
                arr[cell_index] for arr in phase.sim.dx_cell_time])
            dy = nparray.from_iterable([
                arr[cell_index] for arr in phase.sim.dy_cell_time])

            # Get the total magnitude.
            disp = mathunit.upscale_coordinates(np.sqrt(dx ** 2 + dy ** 2))
        else:
            disp = column_data_empty

        csv_column_name_values.extend(('Displacement_um', disp))

        # ................{ CSV EXPORT                      }..................
        # Ordered dictionary mapping from CSV column names to data arrays.
        csv_column_name_to_values = OrderedArgsDict(*csv_column_name_values)

        # Export this data to this CSV file.
        npcsv.write_csv(
            filename=self._get_csv_filename(
                phase=phase, basename_sans_filetype='ExportedData'),
            column_name_to_values=csv_column_name_to_values,
        )

    # ..................{ EXPORTERS ~ cell : vmem           }..................
    #FIXME: This exporter appears to currently be broken for the non-ECM
    #case, despite working for the ECM case. The exception being raised is:
    #
    #    BetseSequenceException: Column "FFT_Vmem" length 9 differs from length
    #    5 of prior columns.
    @piperunner(
        categories=('Single Cell', 'Voltage', 'FFT'),

        #FIXME: Eliminate this requirement after resolving the above issue.
        requirements=phasereqs.VOLTAGE_EXTRA,
    )
    def export_cell_vmem_fft(
        self, phase: SimPhase, conf: SimConfExportCSV) -> None:
        '''
        Save a plaintext file in comma-separated value (CSV) format containing
        the finite Fourier transform (FFT) of all transmembrane voltages for
        all sampled time steps spatially situated at the centre of the single
        cell indexed ``plot cell index`` in the current simulation
        configuration.
        '''

        # Number of sampled time steps.
        sample_size = len(phase.sim.time)

        # Time in seconds between each sampled time step.
        sample_spacing = phase.sim.time[1] - phase.sim.time[0]

        # One-dimensional Numpy array of all transmembrane voltages (Vmems) for
        # each sampled time step spatially situated at the centre of the single
        # cell indexed by the "plot cell index" entry for the passed phase.
        cell_times_vmems = self._get_cell_times_vmems(phase)

        cell_data = (1/sample_size) * (
            cell_times_vmems - np.mean(cell_times_vmems))

        # FFT of voltage.
        f_axis = np.fft.rfftfreq(sample_size, d=sample_spacing)
        fft_data_o = np.fft.rfft(cell_data)

        #FIXME: Numpy already provides a function for obtaining the magnitude of
        #a complex array: np.absolute(). The following math reduces to simply:
        #    fft_data = np.absolute(fft_data_o)
        fft_data = np.sqrt(np.real(fft_data_o)**2 + np.imag(fft_data_o)**2)
        # print('f_axis: {}'.format(f_axis))
        # print('fft_data: {}'.format(fft_data))

        # Ordered dictionary mapping from CSV column names to data arrays.
        csv_column_name_to_values = OrderedArgsDict(
            'frequency_Hz', f_axis,
            'FFT_Vmem', fft_data,
        )

        # Export this data to this CSV file.
        npcsv.write_csv(
            filename=self._get_csv_filename(
                phase=phase, basename_sans_filetype='ExportedData_FFT'),
            column_name_to_values=csv_column_name_to_values,
        )

    # ..................{ EXPORTERS ~ cells : vmem          }..................
    @piperunner(categories=('Cell Cluster', 'Voltage', 'Transmembrane',))
    def export_cells_vmem(
        self, phase: SimPhase, conf: SimConfExportCSV) -> None:
        '''
        Save one plaintext file in comma-separated value (CSV) format
        containing all transmembrane voltages (Vmem) upscaled and averaged from
        all cell membranes onto cell centres for each sampled time step of the
        current simulation phase.
        '''

        # Two-dimensional Numpy array of all transmembrane voltages.
        cells_times_vmems = mathunit.upscale_units_milli(
            phase.sim.vm_ave_time)

        # Export this data to this CSV file.
        self._export_cells_times_data(
            phase=phase,
            cells_times_data=cells_times_vmems,
            csv_column_name='Vmem [mV]',
            csv_dir_basename='Vmem2D_TextExport',
            csv_basename_prefix='Vmem2D_',
        )

    # ..................{ PRIVATE ~ getters                 }..................
    @type_check
    def _get_csv_filename(
        self,

        # Mandatory parameters.
        phase: SimPhase,
        basename_sans_filetype: str,

        # Optional parameters.
        dirname: StrOrNoneTypes = None,
    ) -> str:
        '''
        Absolute filename of a CSV file to be exported with the passed basename
        excluding suffixing ``.``-prefixed filetype (which this method appends
        to this basename as configured by the current simulation configuration)
        residing in the directory with the passed dirname.

        Parameters
        ----------
        phase : SimPhase
            Current simulation phase.
        basename_sans_filetype : str
            Basename (excluding suffixing ``.``-prefixed filetype) of this
            file.
        dirname : StrOrNoneTypes
            Absolute pathname of the directory containing this file. If this
            directory does *not* already exist, this method creates this
            directory. Defaults to ``None``, in which case this pathname
            defaults to the top-level directory containing all CSV files
            exported for the current simulation phase.

        Returns
        ----------
        str
            Absolute filename of a CSV file to be exported with this basename.
        '''

        # If unpassed, default this directory to the top-level directory
        # containing all CSV files exported for this simulation phase.
        if dirname is None:
            dirname = phase.export_dirname

        # Create this directory if needed.
        dirs.make_unless_dir(dirname)

        # Basename of this file.
        basename = '{}.{}'.format(
            basename_sans_filetype, phase.p.csv.filetype)

        # Create and return the absolute filename of this file.
        return pathnames.join(dirname, basename)

    # ..................{ PRIVATE ~ exporters               }..................
    @type_check
    def _export_cells_times_data(
        self,
        phase: SimPhase,
        cells_times_data: SequenceTypes,
        csv_column_name: str,
        csv_dir_basename: str,
        csv_basename_prefix: str,
    ) -> None:
        '''
        Save one plaintext file in comma-separated value (CSV) format
        containing arbitrary simulation data spatially situated at cell centres
        for each sampled time step of the current simulation phase.

        Parameters
        ----------
        phase : SimPhase
            Current simulation phase.
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

        # Absolute pathname of the directory containing all CSV files
        # specifically exported by this method.
        csv_dirname = pathnames.join(
            phase.export_dirname, csv_dir_basename)

        # One-dimensional Numpy arrays of the X and Y coordinates
        # (respectively) of the centres of all cells.
        cell_centres_x = mathunit.upscale_coordinates(
            phase.cells.cell_centres[:,0])
        cell_centres_y = mathunit.upscale_coordinates(
            phase.cells.cell_centres[:,1])

        # For the 0-based index of each sampled time step...
        for time_step in range(len(phase.sim.time)):
            # Basename of the CSV-formatted file exported for this time step,
            # excluding suffixing "."-prefixed filetype.
            csv_basename_sans_filetype = '{}{}'.format(
                csv_basename_prefix, time_step)

            # Absolute filename of this CSV file.
            csv_filename = self._get_csv_filename(
                phase=phase,
                basename_sans_filetype=csv_basename_sans_filetype,
                dirname=csv_dirname,
            )

            # Ordered dictionary mapping from CSV column names to data arrays.
            csv_column_name_to_values = OrderedArgsDict(
                'x [um]', cell_centres_x,
                'y [um]', cell_centres_y,
                csv_column_name, cells_times_data[time_step],
            )

            # Export this data to this CSV file.
            npcsv.write_csv(
                filename=csv_filename,
                column_name_to_values=csv_column_name_to_values)

    # ..................{ PRIVATE ~ properties              }..................
    @type_check
    def _get_cell_times_vmems(self, phase: SimPhase) -> NumpyArrayType:
        '''
        One-dimensional Numpy array of all transmembrane voltages for each
        sampled time step spatially situated at the centre of the single cell
        indexed by the ``plot cell index`` entry specified by the passed
        simulation phase.

        Parameters
        ----------
        phase : SimPhase
            Current simulation phase.
        '''

        # 0-based index of the cell to serialize time data for.
        cell_index = phase.p.visual.single_cell_index

        if phase.p.is_ecm:
            cell_times_vmems = []
            for vm_at_mem in phase.sim.vm_time:
                vm_t = mathunit.upscale_units_milli(
                    cell_ave(phase.cells,vm_at_mem)[cell_index])
                cell_times_vmems.append(vm_t)
        else:
            cell_times_vmems = mathunit.upscale_units_milli(
                phase.sim.vm_time)

        return nparray.from_iterable(cell_times_vmems)
