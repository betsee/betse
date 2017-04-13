#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
High-level vector field subcache functionality.
'''

# ....................{ IMPORTS                            }....................
from betse.science.math.vector.vecfldcls import VectorFieldCellsCache
from betse.science.math.vector.veccls import VectorCellsCache
from betse.science.simulate.cache.cacheabc import SimPhaseCacheABC
from betse.util.type.call.memoizers import property_cached
# from betse.util.type.types import type_check

# ....................{ CONSTANTS                          }....................
_MAGNITUDE_FACTOR_CURRENTS = 100
'''
Factor by which to multiply each magnitude of each vector in each vector field
of intra- and/or extracellular current densities, producing magnitude in units
of uA/cm^2.
'''

# ....................{ SUBCLASSES                         }....................
class SimPhaseCacheVectorField(SimPhaseCacheABC):
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
        '''

        return VectorFieldCellsCache(
            x=VectorCellsCache(
                phase=self._phase,
                times_cells_centre=self._phase.sim.I_cell_x_time),
            y=VectorCellsCache(
                phase=self._phase,
                times_cells_centre=self._phase.sim.I_cell_y_time),
            magnitude_factor=_MAGNITUDE_FACTOR_CURRENTS,
        )


    @property_cached
    def currents_extra(self) -> VectorFieldCellsCache:
        '''
        Vector field cache of all exracellular current densities over all time
        steps of the current simulation phase, originally spatially situated at
        environmental grid space centres.

        Raises
        ----------
        BetseSimConfigException
            If this simulation has disabled extracellular spaces.
        '''

        # If extracellular spaces are disabled, raise an exception.
        self._phase.p.die_unless_ecm()

        # Create and return this field.
        return VectorFieldCellsCache(
            x=VectorCellsCache(
                phase=self._phase,
                times_grids_centre=self._phase.sim.I_tot_x_time),
            y=VectorCellsCache(
                phase=self._phase,
                times_grids_centre=self._phase.sim.I_tot_y_time),
            magnitude_factor=_MAGNITUDE_FACTOR_CURRENTS,
        )

    # ..................{ PROPERTIES ~ electric              }..................
    @property_cached
    def electric_intra(self) -> VectorFieldCellsCache:
        '''
        Vector field cache of the electric field across all intracellular spaces
        over all time steps of the current simulation phase, originally
        spatially situated at cell membrane midpoints.
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
        Vector field cache of the electric field across all extracellular spaces
        over all time steps of the current simulation phase, originally
        spatially situated at environmental grid space centres.

        Raises
        ----------
        BetseSimConfigException
            If this simulation has disabled extracellular spaces.
        '''

        # If extracellular spaces are disabled, raise an exception.
        self._phase.p.die_unless_ecm()

        # Create and return this field.
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

    # ..................{ PROPERTIES ~ vmem                  }..................
    @property_cached
    def voltage_membrane_polarity(self) -> VectorFieldCellsCache:
        '''
        Vector field cache of all transmembrane voltage (Vmem) polarities over
        all time steps of the current simulation phase, originally spatially
        situated at cell membrane midpoints.
        '''

        #FIXME: Is "sim.vm_ave_time" a list? If so, we'll need to Numpify here.

        # Two-dimensional Numpy array of all Vmem polarity magnitudes whose:
        #
        # * First dimension indexes each time step.
        # * Second dimension indexes each cell membrane such that each element
        #   is the magnitude of the polarity of the transmembrane voltage across
        #   that membrane for this time step, where this magnitude is defined as
        #   the difference between:
        #   * The transmembrane voltage across that membrane for this time step.
        #   * The average transmembrane voltage across all membranes of the cell
        #     containing that membrane for this time step.
        polm = (
            self._phase.sim.vm -
            self._phase.sim.vm_ave_time[:,self._phase.cells.mem_to_cells])

        #FIXME: What does multiplication by "polm" do here?

        # X and Y coordinates of all normal unit vectors orthogonal to all
        # tangent unit vectors of all cell membranes.
        polx = polm * self._phase.cells.mem_vects_flat[:,2]
        poly = polm * self._phase.cells.mem_vects_flat[:,3]

        #FIXME: Document us up.
        pcx = self._phase.map_membranes_midpoint_to_cells_centre(
            polx*self._phase.cells.mem_sa) / self._phase.cells.cell_sa
        pcy = self._phase.map_membranes_midpoint_to_cells_centre(
            poly*self._phase.cells.mem_sa) / self._phase.cells.cell_sa

        # Create, return, and cache this vector field.
        return VectorFieldCellsCache(
            x=VectorCellsCache(phase=self._phase, times_cells_centre=pcx),
            y=VectorCellsCache(phase=self._phase, times_cells_centre=pcy))
