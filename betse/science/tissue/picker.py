#!/usr/bin/env python3
# Copyright 2015 by Alexis Pietak & Cecil Curry
# See "LICENSE" for further details.

'''
Class hierarchy collectively implementing various methods for assigning a subset
of the total cell population to the corresponding tissue profile.
'''

# ....................{ IMPORTS                            }....................
import random
from abc import ABCMeta, abstractmethod
from betse.exceptions import (
    BetseExceptionMethodUnimplemented, BetseExceptionParameters)
from betse.science import toolbox
from betse.util.path import files, paths
from betse.util.type import types

# ....................{ BASE                               }....................
class TissuePicker(object, metaclass = ABCMeta):
    '''
    Abstract base class of all tissue matching classes.

    Instances of this class assign a subset of all cells matching
    subclass-specific criteria (e.g., explicit indexing, randomized selection,
    spatial location) to the corresponding tissue profile.
    '''

    @abstractmethod
    def get_cell_indices(
        self, cells, p, ignoreECM: bool = False):
        '''
        Get a Numpy array of the indices of all cells selected by this picker.

        Parameters
        ---------------------------------
        cells : Cells
            Instance of the `Cells` class.
        p : Parameters
            Instance of the `Parameters` class.
        ignoreECM : bool
            If `True`, electromagnetism (e.g., `p.sim_ECM`) will be ignored;
            else, electromagnetism will be simulated. Defaults to `False`.

        Returns
        ---------------------------------
        ndarray
            See method synopsis above.
        '''
        pass


    @abstractmethod
    def get_removal_cell_indices(self, cells):
        '''
        Get a Numpy array of the indices of all cells selected by this picker,
        suitable for permanently removing these cells from the cluster.

        Parameters
        ---------------------------------
        cells : Cells
            Instance of the `Cells` class.

        Returns
        ---------------------------------
        ndarray
            See method synopsis above.
        '''
        pass


    #FIXME: No idea what a "cluster mask" specifically is. Improve docstring.
    #You'll thank me in the pancake-laden morning.
    @abstractmethod
    def get_removal_cluster_mask(self, cells):
        '''
        Get a ??? excluding all cells selected by this picker, intended to
        replace the cluster mask when permanently removing these cells from the
        cluster.

        Parameters
        ---------------------------------
        cells : Cells
            Instance of the `Cells` class.

        Returns
        ---------------------------------
        ndarray
            See method synopsis above.
        '''
        pass


class TissuePickerAll(TissuePicker):
    '''
    All-inclusive tissue picker.

    This matcher unconditionally matches _all_ cells.
    '''

    def get_cell_indices(
        self, cells, p, ignoreECM: bool = False):
        assert types.is_cells(cells),  types.assert_not_cells(cells)
        assert types.is_parameters(p), types.assert_not_parameters(p)
        assert types.is_bool(ignoreECM), types.assert_not_bool(ignoreECM)

        # If either not simulating *OR* ignoring electromagnetism, do so.
        if p.sim_ECM is False or ignoreECM is True:
            target_inds = cells.cell_i

        # Else, simulate electromagnetism.
        else:
            target_inds = cells.cell_to_mems[cells.cell_i]
            target_inds, _, _ = toolbox.flatten(target_inds)

        return target_inds


    def get_removal_cell_indices(self, cells):
        return cells.cell_i

    #FIXME: Since I have no idea how to return a cluster mask encompassing the
    #entire cell cluster, this class does *NOT* currently support removal.
    def get_removal_cluster_mask(self, cells):
        raise BetseExceptionMethodUnimplemented()

# ....................{ BITMAP                             }....................
class TissuePickerBitmap(TissuePicker):
    '''
    Bitmap-specific tissue picker.

    This matcher matches all cells residing inside the colored pixel area
    defined by an associated bitmap file.

    Attributes
    ----------------------------
    filename : str
        Absolute path of this bitmap.
    '''

    def __init__(self, filename, dirname):
        '''
        Initialize this matcher.

        Parameters
        ----------------------------
        filename : str
            Absolute or relative path of the desired bitmap. If relative (i.e.,
            _not_ prefixed by a directory separator), this path will be
            canonicalized into an absolute path relative to the directory
            containing the current simulation's configuration file.
        dirname : str
            Absolute path of the directory containing the path of the bitmap to
            be loaded (i.e., `filename`). If that path is relative, that path
            will be prefixed by this path to convert that path into an absolute
            path; otherwise, this path will be ignored.
        '''
        assert types.is_str(filename), types.assert_nonstr(filename)
        assert types.is_str( dirname), types.assert_nonstr( dirname)

        # If this is a relative path, convert this into an absolute path
        # relative to the directory containing the source configuration file.
        if paths.is_relative(filename):
            filename = paths.join(dirname, filename)

        # If this absolute path is *NOT* an existing file, raise an exception.
        files.die_unless_file(filename)

        # Persist this path.
        self.filename = filename


    def get_cell_indices(
        self, cells, p, ignoreECM: bool = False):
        assert types.is_cells(cells),  types.assert_not_cells(cells)
        assert types.is_parameters(p), types.assert_not_parameters(p)
        assert types.is_bool(ignoreECM), types.assert_not_bool(ignoreECM)

        # Calculate the indices of all cells residing inside this bitmap.
        bitmask = self._get_bitmapper(cells)
        target_inds = bitmask.good_inds

        # If simulating electromagnetism and at least one cell matches...
        if p.sim_ECM is True and ignoreECM is False and len(target_inds):
            target_inds = cells.cell_to_mems[target_inds]
            target_inds,_,_ = toolbox.flatten(target_inds)

        return target_inds


    def get_removal_cell_indices(self, cells):
        return self._get_bitmapper(cells).good_inds

    def get_removal_cluster_mask(self, cells):
        return cells.cluster_mask - self._get_bitmapper(cells).clipping_matrix

    # ..................{ BITMAP ~ private                   }..................
    def _get_bitmapper(self, cells):
        '''
        Get an instance of the `BitMapper` object providing the indices of all
        cells residing inside this bitmap.

        Parameters
        ---------------------------------
        cells : Cells
            Instance of the `Cells` class.
        '''
        assert types.is_cells(cells), types.assert_not_cells(cells)

        # Avoid circular import dependencies.
        from betse.science.tissue.bitmapper import BitMapper

        # Return the desired bitmap object.
        bitmapper = BitMapper(
            self, cells.xmin, cells.xmax, cells.ymin, cells.ymax)
        bitmapper.clipPoints(cells.cell_centres[:,0], cells.cell_centres[:,1])
        return bitmapper

# ....................{ INDICES                            }....................
class TissuePickerIndices(TissuePicker):
    '''
    Indices-specific tissue picker.

    This matcher matches all cells with the listed indices.

    Attributes
    ----------------------------
    indices : collections.Sequence
        Sequence (e.g., list, tuple) of the indices of all cells to be matched.
    '''

    def __init__(self, indices):
        '''
        Initialize this matcher.

        Parameters
        ----------------------------
        indices : collections.Sequence
            See the class docstring.
        '''
        assert types.is_sequence_nonstr(indices),\
            types.assert_not_sequence_nonstr(indices)
        self.indices = indices


    def get_cell_indices(
        self, cells, p, ignoreECM: bool = False):
        assert types.is_cells(cells),  types.assert_not_cells(cells)
        assert types.is_parameters(p), types.assert_not_parameters(p)
        assert types.is_bool(ignoreECM), types.assert_not_bool(ignoreECM)

        # If either not simulating *OR* ignoring electromagnetism, do so.
        if p.sim_ECM is False or ignoreECM is True:
            target_inds = self.indices

        # Else, simulate electromagnetism.
        else:
            target_inds = cells.cell_to_mems[self.indices]
            target_inds,_,_ = toolbox.flatten(target_inds)

        return target_inds


    def get_removal_cell_indices(self, cells):
        return self.indices

    #FIXME: Unconvinced this is in any way right. Shouldn't we return a
    #clipping mask excluding all cells with the listed indices? *shrugs alot*
    def get_removal_cluster_mask(self, cells):
        return None

# ....................{ RANDOM                             }....................
class TissuePickerRandom(TissuePicker):
    '''
    Randomized cell matcher.

    This matcher randomly matches a percentage of cells.

    Attributes
    ----------------------------
    percentage : {int, float}
        Percentage of the total cell population to be randomly matched as an
        integer or float in the range `[0,0, 100.0]`.
    '''

    def __init__(self, percentage):
        '''
        Initialize this matcher.

        Parameters
        ----------------------------
        percentage : {int, float}
            See the class docstring.
        '''
        assert types.is_numeric(percentage),\
            types.assert_not_numeric(percentage)

        # If this is not a valid percentage, raise an exception. This is
        # important enough to always test rather than defer to assertions.
        if not 0.0 <= percentage <= 100.0:
            raise BetseExceptionParameters(
                '{} not in the range [0.0, 100.0].'.format(percentage))

        self.percentage = percentage


    def get_cell_indices(
        self, cells, p, ignoreECM: bool = False):
        assert types.is_cells(cells),  types.assert_not_cells(cells)
        assert types.is_parameters(p), types.assert_not_parameters(p)
        assert types.is_bool(ignoreECM), types.assert_not_bool(ignoreECM)

        data_length = len(cells.cell_i)
        data_fraction = int((self.percentage/100)*data_length)
        random.shuffle(cells.cell_i)
        target_inds_cell = [cells.cell_i[x] for x in range(0,data_fraction)]

        # If either not simulating *OR* ignoring electromagnetism, do so.
        if p.sim_ECM is False or ignoreECM is True:
            target_inds = target_inds_cell

        # Else, simulate electromagnetism.
        else:
            target_inds = cells.cell_to_mems[target_inds_cell]
            target_inds,_,_ = toolbox.flatten(target_inds)

        return target_inds


    # This class does *NOT* currently support removal.
    def get_removal_cell_indices(self, cells):
        raise BetseExceptionMethodUnimplemented()

    def get_removal_cluster_mask(self, cells):
        raise BetseExceptionMethodUnimplemented()
