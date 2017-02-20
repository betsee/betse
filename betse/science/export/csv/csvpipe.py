#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
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
from betse.lib.numpy import arrays
from betse.science.simulate.simphase import SimPhaseABC
from betse.science.simulate.simpipeabc import SimPipelayerABC
# from betse.util.io.log import logs
from betse.science.visual.plot.plotutil import cell_ave
from betse.util.path import dirs, paths
from betse.util.type.call.memoizers import property_cached
from betse.util.type.types import type_check, SequenceTypes
from numpy import ndarray

# ....................{ SUBCLASSES                         }....................
class SimPipelayerCSV(SimPipelayerABC):
    '''
    **Post-simulation CSV pipeline** (i.e., class iteratively creating all
    post-simulation plaintext spreadsheets in comma-separated value (CSV) format
    requested by the current simulation configuration).
    '''

    # ..................{ INITIALIZERS                       }..................
    @type_check
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(
            *args,
            label_verb='Saving',
            label_singular='CSV',
            label_plural='CSVs',
            **kwargs)

    # ..................{ SUPERCLASS                         }..................
    #FIXME: Implement this property *AFTER* generalizing the configuration file
    #format to support CSV pipelines..
    @property
    def runner_names(self) -> SequenceTypes:
        '''
        Sequence of the names of all post-simulation animations enabled by this
        simulation configuration.
        '''

        raise BetseMethodUnimplementedException()

    # ..................{ RUNNERS ~ cell                     }..................
    def run_cell_data_all(self) -> None:
        '''
        Save a plaintext file in comma-separated value (CSV) format containing
        several cell-specific time series (e.g., ion concentrations, membrane
        voltages, voltage-gated ion channel pump rates) for the single cell
        whose index is given by the ``plot cell index`` entry for the current
        simulation configuration.
        '''

        # Prepare to serialize the current CSV file.
        self._prep_run()

        cc_cell = []
        dd_cell = []

        # 0-based index of the cell to serialize time data for.
        ci = self._phase.p.plot_cell

        # create the header, first entry will be time:
        headr = 'time_s'
        t = np.asarray(self._phase.sim.time)

        #-----------Vmem----------------------------------------------

        #FIXME: Replace these local variable placeholders with the equivalent
        #"phase.sim", "phase.cells", and "phase.p". To do so, continue iteratively
        #pushing these declarations further and further down this function until
        #they are no longer required at all.
        sim   = self._phase.sim
        cells = self._phase.cells
        p     = self._phase.p

        # next entry will be Vm:
        headr = headr + ',' + 'Vmem_mV'

        # golman Vmem------------------------------------------------------

        if p.GHK_calc is True:
            vm_goldman = [x[p.plot_cell]*1000 for x in sim.vm_GHK_time]
        else:
            vm_goldman = np.zeros(len(sim.time))

        vm_goldman = np.asarray(vm_goldman)
        # next entry will be Vm:
        headr = headr + ',' + 'Goldman_Vmem_mV'

        # ------Na K pump rate-----------------------------------------------

        # Na-K-pump rate:
        if p.sim_ECM is False:
            pump_rate = [pump_array[p.plot_cell] for pump_array in sim.rate_NaKATP_time]

        else:
            pump_rate = [pump_array[cells.cell_to_mems[p.plot_cell][0]] for pump_array in sim.rate_NaKATP_time]

        pump_rate = np.asarray(pump_rate)
        headr = headr + ',' + 'NaK-ATPase_Rate_mol/m2s'

        #----------cell concentrations---------------------------------------

        # create the header starting with cell concentrations
        for i in range(0,len(sim.ionlabel)):
            label = sim.ionlabel[i]
            headr = headr + ',' + 'cell_' + label + '_mmol/L'
            cc_m = [arr[i][ci] for arr in sim.cc_time]
            cc_m = np.asarray(cc_m)
            cc_cell.append(cc_m)

        cc_cell = np.asarray(cc_cell)

        #----------membrane permeabilities---------------------------------------

        # create the header starting with membrane permeabilities
        if p.sim_ECM is False:
            for i in range(0,len(sim.ionlabel)):
                label = sim.ionlabel[i]
                headr = headr + ',' + 'Dm_' + label + '_m2/s'
                dd_m = [arr[i][ci] for arr in sim.dd_time]
                dd_m = np.asarray(dd_m)
                dd_cell.append(dd_m)

        else:
            for i in range(0,len(sim.ionlabel)):
                label = sim.ionlabel[i]
                headr = headr + ',' + 'Dm_' + label + '_m2/s'
                dd_m = [arr[i][cells.cell_to_mems[ci][0]] for arr in sim.dd_time]
                dd_m = np.asarray(dd_m)
                dd_cell.append(dd_m)

        dd_cell = np.asarray(dd_cell)

        #FIXME: The local "Imem" list defined below is never used anywhere. This
        #list should either be saved to this CSV or removed entirely. Indeed, this
        #doesn't seem quite right at all; we append the currents-specific header
        #"I_A/m2" here but fail to save this currents-specific list to the CSV at
        #the end of this function. Doesn't this imply that the CSV will erroneously
        #contain more headers than actual columns?

        #----Transmembrane currents--------------------------
        if p.sim_ECM is False:
            Imem = [memArray[p.plot_cell] for memArray in sim.I_mem_time]
        else:
            Imem = [
                memArray[cells.cell_to_mems[p.plot_cell][0]]
                for memArray in sim.I_mem_time]

        headr = headr + ',' + 'I_A/m2'

        #----Hydrostatic pressure--------
        p_hydro = [arr[p.plot_cell] for arr in sim.P_cells_time]

        headr = headr + ',' + 'HydroP_Pa'

        # ---Osmotic pressure-----------
        if p.deform_osmo is True:

            p_osmo = [arr[p.plot_cell] for arr in sim.osmo_P_delta_time]

        else:
            p_osmo = np.zeros(len(sim.time))

        headr = headr + ',' + 'OsmoP_Pa'

        # total deformation ---------------------------------------
        if p.deformation is True and sim.run_sim is True:

            # extract time-series deformation data for the plot cell:
            dx = np.asarray([arr[p.plot_cell] for arr in sim.dx_cell_time])
            dy = np.asarray([arr[p.plot_cell] for arr in sim.dy_cell_time])

            # get the total magnitude:
            disp = p.um*np.sqrt(dx**2 + dy**2)

        else:
            disp = np.zeros(len(sim.time))

        headr = headr + ',' + 'Displacement_um'

        # # cluster polarization vector------------------------------
        # polar_x = sim.Pol_tot_x_time
        # headr = headr + ',' + 'Polarization x_C/m'
        # polar_y = sim.Pol_tot_y_time
        # headr = headr + ',' + 'Polarization y_C/m'

        dataM = np.column_stack((
            t,
            self._cell_times_vmems,
            vm_goldman,
            pump_rate,
            cc_cell.T,
            dd_cell.T,
            Imem,
            p_hydro,
            p_osmo,
            disp,
        ))

        savedData = paths.join(self._phase.save_dirname, 'ExportedData.csv')

        #FIXME: Wrap with a new betse.lib.numpy.numpys.write_csv() utility function.
        #To avoid desynchronization issues between the number of header columns and
        #the number of all following columns, this utility function should:
        #
        #* Accept a sequence of header names rather than a ","-delimited string
        #  concatenating such names.
        #* Validate that the length of this sequence is equal to the length of each
        #  sequence in the second dimension of the passed column data sequence.
        np.savetxt(savedData, dataM, delimiter = ',', header=headr)

    # ..................{ RUNNERS ~ cell : vmem              }..................
    def run_cell_vmem_fft(self) -> None:
        '''
        Save a plaintext file in comma-separated value (CSV) format containing
        the finite Fourier transform (FFT) of all transmembrane voltages for all
        sampled time steps spatially situated at the centre of the cell whose
        index is given by the ``plot cell index`` entry for the current
        simulation configuration.
        '''

        # Prepare to serialize the current CSV file.
        self._prep_run()

        # FFT of voltage :
        sample_size = len(self._phase.sim.time)
        sample_spacing = self._phase.sim.time[1] - self._phase.sim.time[0]

        cell_data = (1/sample_size) * (
            self._cell_times_vmems - np.mean(self._cell_times_vmems))

        f_axis = np.fft.rfftfreq(sample_size, d=sample_spacing)
        fft_data_o = np.fft.rfft(cell_data)
        fft_data = np.sqrt(np.real(fft_data_o)**2 + np.imag(fft_data_o)**2)

        headr2 = 'frequency_Hz'
        headr2 = headr2 + ',' + 'FFT_Vmem'

        dataFFT = np.column_stack((f_axis, fft_data))

        savedData_FFT = paths.join(
            self._phase.save_dirname, 'ExportedData_FFT.csv')

        np.savetxt(savedData_FFT, dataFFT, delimiter = ',', header=headr2)

    # ..................{ RUNNERS ~ cells : vmem             }..................
    def run_cells_vmems(self) -> None:
        '''
        Save one plaintext file in comma-separated value (CSV) format for each
        sequence of transmembrane voltages spatially situated at cell centres
        for each sampled time step of the current simulation phase.
        '''

        # Prepare to serialize the current CSV file.
        self._prep_run()

        # For the 0-based index of each sampled time step...
        for time_step in range(len(self._phase.sim.time)):
            # One-dimensional Numpy array of all cell voltages.
            cells_vmem = 1.0e3*self._phase.sim.vm_ave_time[time_step]

            # Save all membrane voltages for this time step to a unique CSV file.
            self._save_cells_vmem(
                time_step=time_step,
                cells_vmem=cells_vmem,
            )

            #     simdata_x = 1.0e3*sim.pol_x_time[i]
            #     save_cells_vmem(i, simdata_x, cells, p, foldername = 'Polarization_x', filebit = 'Pol_x')
            #
            #     simdata_y = 1.0e3 * sim.pol_y_time[i]
            #     save_cells_vmem(i, simdata_y, cells, p, foldername='Polarization_y', filebit='Pol_y')


    @type_check
    def _save_cells_vmem(
        self,
        time_step: int,
        cells_vmem: SequenceTypes,

        #FIXME: Document the following parameters.
        foldername: str = 'Vmem2D_TextExport',
        filebit: str = 'Vmem2D_',
    ) -> None:
        '''
        Save a plaintext file in comma-separated value (CSV) format containing
        all transmembrane voltages spatially situated at cell centres for the
        passed sampled time step of the current simulation phase.

        Parameters
        ----------------------------
        time_step : int
            Current sampled time step.
        cells_vmem : ndarray
            One-dimensional Numpy array of all transmembrane voltages spatially
            situated at the centres of all cells for this time step.
        '''

        filepath = paths.join(self._phase.save_dirname, foldername)
        dirs.make_unless_dir(filepath)

        basename = '{}{}.csv'.format(filebit, time_step)
        savedData_2d = paths.join(filepath, basename)

        dataM = np.column_stack((
            self._phase.p.um * self._phase.cells.cell_centres[:,0],
            self._phase.p.um * self._phase.cells.cell_centres[:,1],
            cells_vmem,
        ))
        hdr = 'x [um], y [um], Vmem [mV]'

        #FIXME: Wrap with a new betse.lib.numpy.numpys.write_csv() utility function.
        #See above for details.
        np.savetxt(savedData_2d, dataM, delimiter=',', header=hdr)

    # ....................{ PREPARERS                          }....................
    def _prep_run(self) -> None:
        '''
        Prepare to serialize the current CSV file.

        Specifically, this method:

        * Logs the current attempt to run the calling runner.
        * Creates the directory containing the CSV file subsequently created by
          this runner if needed.
        '''

        # Log the current attempt to run the calling runner.
        self._log_run()

        # Create the directory containing the current CSV file if needed.
        dirs.make_unless_dir(self._phase.save_dirname)

    # ....................{ PROPERTIES                         }....................
    @property_cached
    def _cell_times_vmems(self) -> ndarray:
        '''
        One-dimensional Numpy array of all transmembrane voltages for each
        sampled time step spatially situated at the centre of the cell
        whose index is given by the ``plot cell index`` entry for the current
        simulation configuration.
        '''

        # 0-based index of the cell to serialize time data for.
        cell_index = self._phase.p.plot_cell

        if self._phase.p.sim_ECM:
            cell_times_vmems = []
            for vm_at_mem in self._phase.sim.vm_time:
                vm_t = 1000*cell_ave(self._phase.cells,vm_at_mem)[cell_index]
                cell_times_vmems.append(vm_t)
        else:
            cell_times_vmems = [
                arr[cell_index]*1000 for arr in self._phase.sim.vm_time]

        return arrays.from_sequence(cell_times_vmems)

# ....................{ PIPELINES                          }....................
@type_check
def pipeline(phase: SimPhaseABC) -> None:
    '''
    Save all currently enabled data exports (e.g., spreadsheets in
    comma-separated value form) for the passed simulation phase.

    Parameters
    ----------------------------
    phase: SimPhaseABC
        Current simulation phase.
    '''

    # Post-simulation CSV pipeline producing all such CSVs.
    pipelayer = SimPipelayerCSV(phase)

    if phase.p.exportData:
        pipelayer.run_cell_data_all()
        pipelayer.run_cell_vmem_fft()

    if phase.p.exportData2D:
        pipelayer.run_cells_vmems()
