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
from betse.science.simulate.simphase import SimPhaseABC
# from betse.util.io.log import logs
from betse.science.visual.plot.plotutil import cell_ave
from betse.util.path import dirs, paths
from betse.util.type.types import type_check, SequenceTypes

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

    if phase.p.exportData:
        save_cell_times_data(phase)

    if phase.p.exportData2D:
        save_cells_times_vmem(phase)

# ....................{ SAVERS                             }....................
#FIXME: Split all logic saving the FFT-specific CSV file saved by this function
#into a new save_cell_vmem_fft() function.
@type_check
def save_cell_times_data(phase: SimPhaseABC) -> None:
    '''
    Save several cell-specific time series (e.g., ion concentrations, membrane
    voltages, voltage-gated ion channel pump rates) for the single cell whose
    index is given by the ``plot cell index`` entry for the passed simulation
    phase's configuration, as plaintext files in comma-separated value (CSV)
    format.

    Parameters
    ----------------------------
    phase: SimPhaseABC
        Current simulation phase.
    '''

    # Create the top-level directory containing these exports if needed.
    dirs.make_unless_dir(phase.save_dirname)

    savedData     = paths.join(phase.save_dirname, 'ExportedData.csv')
    savedData_FFT = paths.join(phase.save_dirname, 'ExportedData_FFT.csv')

    cc_cell = []
    dd_cell = []

    ci = phase.p.plot_cell  # index of cell to get time data for

    # create the header, first entry will be time:
    headr = 'time_s'
    t = np.asarray(phase.sim.time)

    #-----------Vmem----------------------------------------------

    #FIXME: Replace these local variable placeholders with the equivalent
    #"phase.sim", "phase.cells", and "phase.p". To do so, continue iteratively
    #pushing these declarations further and further down this function until
    #they are no longer required at all.
    sim   = phase.sim
    cells = phase.cells
    p     = phase.p

    # next entry will be Vm:
    headr = headr + ',' + 'Vmem_mV'

    if p.sim_ECM is False:
        vm = [arr[ci]*1000 for arr in sim.vm_time]

    else:
        vm = []
        for vm_at_mem in sim.vm_time:
            vm_t = 1000*cell_ave(cells,vm_at_mem)[ci]
            vm.append(vm_t)

    vm = np.asarray(vm)

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
        Imem = [memArray[cells.cell_to_mems[p.plot_cell][0]] for memArray in sim.I_mem_time]

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

    # FFT of voltage :
    sample_size = len(sim.time)
    sample_spacing = sim.time[1] - sim.time[0]

    cell_data = (1/sample_size)*(vm - np.mean(vm))

    f_axis = np.fft.rfftfreq(sample_size, d=sample_spacing)
    fft_data_o = np.fft.rfft(cell_data)
    fft_data = np.sqrt(np.real(fft_data_o)**2 + np.imag(fft_data_o)**2)

    dataM = np.column_stack((t,vm,vm_goldman,pump_rate,cc_cell.T, dd_cell.T,
                             p_hydro,p_osmo,disp))

    headr2 = 'frequency_Hz'
    headr2 = headr2 + ',' + 'FFT_Vmem'

    dataFFT = np.column_stack((f_axis,fft_data))

    #FIXME: Wrap with a new betse.lib.numpy.numpys.write_csv() utility function.
    #To avoid desynchronization issues between the number of header columns and
    #the number of all following columns, this utility function should:
    #
    #* Accept a sequence of header names rather than a ","-delimited string
    #  concatenating such names.
    #* Validate that the length of this sequence is equal to the length of each
    #  sequence in the second dimension of the passed column data sequence.
    np.savetxt(savedData, dataM, delimiter = ',', header=headr)
    np.savetxt(savedData_FFT, dataFFT, delimiter = ',', header=headr2)

# ....................{ SAVERS ~ cells : vmem              }....................
@type_check
def save_cells_times_vmem(phase: SimPhaseABC) -> None:
    '''
    Iteratively save each one-dimensional Numpy array of all transmembrane
    voltages spatially situated at the centres of all cells for the passed
    simulation phase and each sampled time step as a plaintext file in
    comma-separated value (CSV) format.

    Parameters
    ----------------------------
    phase: SimPhaseABC
        Current simulation phase.
    '''

    # For the 0-based index of each sampled time step...
    for time_step in range(len(phase.sim.time)):
        # One-dimensional Numpy array of all cell voltages.
        voltages_cell_centre = 1.0e3 * phase.sim.vm_ave_time[time_step]

        # Save all membrane voltages for this time step to a unique CSV file.
        _save_cells_vmem(
            phase=phase,
            time_step=time_step,
            voltages_cell_centre=voltages_cell_centre,
        )

        #     simdata_x = 1.0e3*sim.pol_x_time[i]
        #     save_cells_vmem(i, simdata_x, cells, p, foldername = 'Polarization_x', filebit = 'Pol_x')
        #
        #     simdata_y = 1.0e3 * sim.pol_y_time[i]
        #     save_cells_vmem(i, simdata_y, cells, p, foldername='Polarization_y', filebit='Pol_y')


@type_check
def _save_cells_vmem(
    phase: SimPhaseABC,
    time_step: int,
    voltages_cell_centre: SequenceTypes,

    #FIXME: Document the following parameters.
    foldername: str = 'Vmem2D_TextExport',
    filebit: str = 'Vmem2D_',
) -> None:
    '''
    Save the one-dimensional Numpy array of all transmembrane voltages spatially
    situated at the centres of all cells for the passed simulation phase and
    sampled time step as a plaintext file in comma-separated value (CSV) format.

    Parameters
    ----------------------------
    phase: SimPhaseABC
        Current simulation phase.
    time_step : int
        Current sampled time step.
    voltages_cell_centre : ndarray
        One-dimensional Numpy array of all transmembrane voltages spatially
        situated at the centres of all cells for this time step.
    '''

    filename = filebit + str(time_step) + '.csv'
    filepath = paths.join(phase.save_dirname, foldername)

    dirs.make_unless_dir(filepath)

    savedData_2d = paths.join(filepath, filename)

    dataM = np.column_stack((
        phase.p.um*phase.cells.cell_centres[:,0],
        phase.p.um*phase.cells.cell_centres[:,1],
        voltages_cell_centre,
    ))
    hdr = 'x [um], y [um], Vmem [mV]'

    #FIXME: Wrap with a new betse.lib.numpy.numpys.write_csv() utility function.
    #See above for details.
    np.savetxt(savedData_2d, dataM, delimiter=',', header=hdr)
