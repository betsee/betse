#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level vector field subcache functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.lib.numpy import nparray
from betse.science.export import expmath
from betse.science.math.cache.cacheabc import SimPhaseCacheABC
from betse.science.math.vector.veccls import VectorCellsCache
from betse.science.math.vector.vecfldcls import VectorFieldCellsCache
from betse.util.type.decorator.decmemo import property_cached

# FIXME: Sess, why the heck would you put these plotting functions into math?!?!?

# ....................{ SUBCLASSES                         }....................
class SimPhaseCacheVectorFieldCells(SimPhaseCacheABC):
    '''
    Simulation phase-specific vector field subcache, persisting all previously
    constructed vector fields for a single simulation phase.
    '''

    # ..................{ PROPERTIES ~ currents              }..................
    @property_cached
    def currents_intra(self) -> VectorFieldCellsCache:
        '''
        Vector field cache of all intracellular current densities over all time
        steps of the current simulation phase, originally spatially situated at
        cell centres.

        For readability of units in exported visuals (e.g., plots), this cache
        additionally upscales these densities to units of uA/cm^2.
        '''

        return VectorFieldCellsCache(
            x=VectorCellsCache(
                phase=self._phase,
                times_cells_centre=expmath.upscale_units_centi(
                    self._phase.sim.I_cell_x_time)),
            y=VectorCellsCache(
                phase=self._phase,
                times_cells_centre=expmath.upscale_units_centi(
                    self._phase.sim.I_cell_y_time)),
        )


    @property_cached
    def currents_extra(self) -> VectorFieldCellsCache:
        '''
        Vector field cache of all exracellular current densities over all time
        steps of the current simulation phase, originally spatially situated at
        environmental grid space centres.

        For readability of units in exported visuals (e.g., plots), this cache
        additionally upscales these densities to units of uA/cm^2.
        '''

        # Create, return, and cache this field.
        return VectorFieldCellsCache(
            x=VectorCellsCache(
                phase=self._phase,
                times_grids_centre=expmath.upscale_units_centi(
                    self._phase.sim.I_tot_x_time)),
            y=VectorCellsCache(
                phase=self._phase,
                times_grids_centre=expmath.upscale_units_centi(
                    self._phase.sim.I_tot_y_time)),
        )

    # ..................{ PROPERTIES ~ deform                }..................
    #FIXME: Raise an exception unless deformations are enabled. To do so sanely,
    #we'll want to define a new @phase_property_cached decorator accepting an
    #optional "requirements" parameter, much like the existing @piperunner
    #decorator. (For now, simply ignore this for simplicity.)
    @property_cached
    def deform_total(self) -> VectorFieldCellsCache:
        '''
        Vector field cache of all upscaled **total cellular displacements**
        (i.e., summations of all cellular deformations due to galvanotropic and
        osmotic pressure body forces) over all sampled time steps of the current
        simulation phase, originally spatially situated at cell centres.

        For readability of units in exported visuals (e.g., plots), this cache
        additionally upscales these displacements from meters to micrometers.
        '''

        return VectorFieldCellsCache(
            x=VectorCellsCache(
                phase=self._phase,
                times_cells_centre=expmath.upscale_coordinates(
                    self._phase.sim.dx_cell_time)),
            y=VectorCellsCache(
                phase=self._phase,
                times_cells_centre=expmath.upscale_coordinates(
                    self._phase.sim.dy_cell_time)),
        )

    # ..................{ PROPERTIES ~ electric              }..................
    @property_cached
    def electric_intra(self) -> VectorFieldCellsCache:
        '''
        Vector field cache of the intracellular electric field over all sampled
        time steps of the current simulation phase, originally spatially
        situated at cell membrane midpoints.
        '''

        return VectorFieldCellsCache(
            x=VectorCellsCache(
                phase=self._phase,
                times_membranes_midpoint=self._phase.sim.efield_gj_x_time),
            y=VectorCellsCache(
                phase=self._phase,
                times_membranes_midpoint=self._phase.sim.efield_gj_y_time),
        )


    @property_cached
    def electric_extra(self) -> VectorFieldCellsCache:
        '''
        Vector field cache of the extracellular electric field over all sampled
        time steps of the current simulation phase, originally spatially
        situated at environmental grid space centres.

        Raises
        ----------
        BetseSimConfException
            If this simulation has disabled extracellular spaces.
        '''

        # If extracellular spaces are disabled, raise an exception.
        self._phase.p.die_unless_ecm()

        # Create, return, and cache this field.
        return VectorFieldCellsCache(
            x=VectorCellsCache(
                phase=self._phase,
                times_grids_centre=self._phase.sim.efield_ecm_x_time),
            y=VectorCellsCache(
                phase=self._phase,
                times_grids_centre=self._phase.sim.efield_ecm_y_time),
        )

    # ..................{ PROPERTIES ~ microtubule           }..................
    @property_cached
    def microtubule(self) -> VectorFieldCellsCache:
        '''
        Vector field cache of all cellular microtubules over all time steps of
        the current simulation phase, originally spatially situated at cell
        membrane midpoints.
        '''

        return VectorFieldCellsCache(
            x=VectorCellsCache(
                phase=self._phase,
                times_membranes_midpoint=self._phase.sim.mtubes_x_time),
            y=VectorCellsCache(
                phase=self._phase,
                times_membranes_midpoint=self._phase.sim.mtubes_y_time),
        )

    # ..................{ PROPERTIES ~ voltage               }..................
    @property_cached
    def voltage_polarity(self) -> VectorFieldCellsCache:
        '''
        Vector field cache of all cellular voltage polarities over all time
        steps of the current simulation phase, originally spatially situated at
        cell membrane midpoints.
        '''

        # Two-dimensional Numpy arrays of all transmembrane voltages (Vmem) and
        # Vmem averages across all cell membranes over all time steps.
        vm_time     = nparray.from_iterable(self._phase.sim.vm_time)
        vm_ave_time = nparray.from_iterable(self._phase.sim.vm_ave_time)

        # Two-dimensional Numpy array of all transmembrane voltage polarity
        # vector magnitudes whose:
        #
        # * First dimension indexes each time step.
        # * Second dimension indexes each cell membrane such that each element
        #   is the magnitude of the polarity of the transmembrane voltage across
        #   that membrane for this time step, where this magnitude is defined as
        #   the difference between:
        #   * The transmembrane voltage across that membrane for this time step.
        #   * The average transmembrane voltage across all membranes of the cell
        #     containing that membrane for this time step.
        polarity_membranes_midpoint_magnitudes = (
            vm_time - vm_ave_time[:,self._phase.cells.mem_to_cells])

        # Two-dimensional Numpy arrays of the X and Y components of all Vmem
        # polarity vectors, spatially situated at cell membrane midpoints.
        polarity_membranes_midpoint_x = (
            polarity_membranes_midpoint_magnitudes *
            self._phase.cells.membranes_normal_unit_x)
        polarity_membranes_midpoint_y = (
            polarity_membranes_midpoint_magnitudes *
            self._phase.cells.membranes_normal_unit_y)

        # Two-dimensional Numpy arrays of the X and Y components of all Vmem
        # polarity vectors, spatially situated at cell centres.
        polarity_cells_centre_x = (
            self._phase.cells.map_membranes_midpoint_to_cells_centre(
                polarity_membranes_midpoint_x * self._phase.cells.mem_sa) /
            self._phase.cells.cell_sa)
        polarity_cells_centre_y = (
            self._phase.cells.map_membranes_midpoint_to_cells_centre(
                polarity_membranes_midpoint_y * self._phase.cells.mem_sa) /
            self._phase.cells.cell_sa)

        # Create, return, and cache this vector field.
        return VectorFieldCellsCache(
            x=VectorCellsCache(
                phase=self._phase, times_cells_centre=polarity_cells_centre_x),
            y=VectorCellsCache(
                phase=self._phase, times_cells_centre=polarity_cells_centre_y))
