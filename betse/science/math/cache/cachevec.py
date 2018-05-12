#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2018 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
High-level vector subcache functionality.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.exceptions import BetseSimVectorException
from betse.science.export import expmath
from betse.science.math.cache.cacheabc import SimPhaseCacheABC
from betse.science.math.vector.veccls import VectorCellsCache
from betse.util.type.decorator.decmemo import property_cached

# ....................{ SUBCLASSES                         }....................
class SimPhaseCacheVectorCells(SimPhaseCacheABC):
    '''
    Simulation phase-specific vector subcache, persisting all previously
    constructed vectors for a single simulation phase.

    Attributes
    ----------
    layer : SimPhaseCacheLayerCellsVector
        Subcache of all vector-based layers constructed for this phase.
    '''

    # ..................{ INITIALIZORS                       }..................
    def __init__(self, *args, **kwargs) -> None:

        # Initialize our superclass with all passed parameters.
        super().__init__(*args, **kwargs)

        # Avoid circular import dependencies.
        from betse.science.math.cache.cachelyrvec import (
            SimPhaseCacheLayerCellsVector)

        # Classify all subcaches imported above.
        self.layer = SimPhaseCacheLayerCellsVector(self._phase)

    # ..................{ PROPERTIES ~ ions                  }..................
    #FIXME: Decorate to require that calcium ions be enabled.
    @property_cached
    def ion_calcium_intra(self) -> VectorCellsCache:
        '''
        Vector cache of all upscaled cellular calcium ion (Ca2+) concentrations
        over all sampled time steps of the current simulation phase, originally
        spatially situated at cell centres.
        '''

        return self._make_ion_intra_cache(ion_index=self._phase.sim.iCa)


    #FIXME: Decorate to require that hydrogen ions be enabled.
    @property_cached
    def ion_hydrogen_intra(self) -> VectorCellsCache:
        '''
        Vector cache of all logarithmically scaled cellular hydrogen ion (H+)
        concentrations over all sampled time steps of the current simulation
        phase, originally spatially situated at cell centres.

        For generality, these concentrations are logarathmically scaled to
        correspond exactly to cellular pH.
        '''

        # Create and return this cache, which logarithmically rather than
        # multiplicatively scales and hence *CANNOT* defer to the standardized
        # _make_ion_cache() method.
        return VectorCellsCache(
            phase=self._phase,
            times_cells_centre=[
                -np.log10(1.0e-3 * ions_concentration[self._phase.sim.iH])
                for ions_concentration in self._phase.sim.cc_time])

    # ..................{ PROPERTIES ~ voltage               }..................
    @property_cached
    def voltage_extra(self) -> VectorCellsCache:
        '''
        Vector cache of all upscaled **extracellular voltages** (i.e., voltages
        across all environmental grid spaces) over all sampled time steps of the
        current simulation phase, originally spatially situated at grid space
        centres.

        For readability of units in exported visuals (e.g., plots), voltages are
        cached upscaled from volts (V) to millivolts (mV).
        '''

        # Shape of the extracellular voltages array to be cached. While the
        # original "venv_time" is a two-dimensional array indexed first by
        # sampled time steps and then by flattened grid spaces, this is a
        # three-dimensional array indexed first by sampled time steps and then
        # by nonflattened grid space rows and columns as unique dimensions. Why?
        # Because layers plotting extracellular data assume the latter.
        voltage_extra_shape = (
            len(self._phase.sim.venv_time),) + self._phase.cells.X.shape

        # Create, return, and create this cache, both upscaled and reshaped as
        # detailed above.
        return VectorCellsCache(
            phase=self._phase,
            times_grids_centre=expmath.upscale_units_milli(
                self._phase.sim.venv_time).reshape(voltage_extra_shape))


    @property_cached
    def voltage_membrane(self) -> VectorCellsCache:
        '''
        Vector cache of all upscaled **transmembrane voltages** (i.e., voltages
        across all gap junctions connecting intracellular membranes) over all
        sampled time steps of the current simulation phase, originally spatially
        situated at cell membrane midpoints.

        For readability of units in exported visuals (e.g., plots), voltages are
        cached upscaled from volts (V) to millivolts (mV).
        '''

        return VectorCellsCache(
            phase=self._phase,
            times_membranes_midpoint=expmath.upscale_units_milli(
                self._phase.sim.vm_time))

    # ..................{ PRIVATE ~ ions                     }..................
    def _make_ion_intra_cache(self, ion_index: int) -> VectorCellsCache:
        '''
        Create and return a vector cache of all upscaled concentrations of the
        ion with the passed index over all sampled time steps of the current
        simulation phase, originally spatially situated at cell centres.

        Parameters
        ----------
        ion_index : int
            0-based index of the ion to cache concentrations for.

        Returns
        ----------
        VectorCellsCache
            Cache of all upscaled concentrations of this ion.

        Raises
        ----------
        BetseSimVectorException
            If no ion with this index exists (i.e., if this index is *not* in
            the range ``[0, len(self._phase.sim.cc_time[0]))``).
        '''

        # If no ion with this index exists, raise an exception.
        if not (0 <= ion_index < len(self._phase.sim.cc_time[0])):
            raise BetseSimVectorException(
                'Ion with index {} not found '
                '(i.e., not in range [0, {})).'.format(
                    ion_index, len(self._phase.sim.cc_time[0])))

        # Create and return this cache.
        return VectorCellsCache(
            phase=self._phase,
            times_cells_centre=expmath.upscale_units_micro(
                ions_concentration[ion_index]
                for ions_concentration in self._phase.sim.cc_time))
