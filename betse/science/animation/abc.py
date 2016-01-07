#!/usr/bin/env python3
# Copyright 2014-2016 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Abstract base classes of all Matplotlib-based animation classes.
'''

# ....................{ IMPORTS                            }....................
from abc import ABCMeta #, abstractmethod, abstractstaticmethod
from betse.exceptions import BetseExceptionParameters
# from betse.lib.matplotlib import mpl
from betse.util.path import dirs, paths
from betse.util.type import types

# ....................{ BASE                               }....................
#FIXME: Privatize all public attributes declared below. Raving river madness!
#FIXME: Document the "clrAutoscale", "clrMin", and "clrMax" attributes. Sizzle!
#FIXME: Refactor the "savedAni" attribute into a template containing exactly one
#"{"- and "}"-delimited string, to satisfy our "FrameFileWriter" API. Currently,
#this attribute is merely the prefix for all frame image files to be saved.

class Animation(object, metaclass=ABCMeta):
    '''
    Abstract base class of all animation classes.

    Instances of this class animate the spatial distribution of modelled
    variables (e.g., Vmem) over all time steps of the simulation.

    Attributes
    ----------
    sim : Simulation
        Current simulation.
    cells : Cells
        Current cell cluster.
    p : Parameters
        Current simulation configuration.
    ani_repeat : bool
        `True` if this animation is to be indefinitely repeated (and hence
        require manual closing) _or_ `False` if this animation is to be
        performed only once.
    clrAutoscale : bool
        ???.
    clrMin : float
        ???.
    clrMax : float
        ???.
    colormap : Colormap
        Matplotlib colormap with which to create this animation's colorbar.
    savedAni : str
        Path prefix for all frame image files to be saved.
    saveFile : str
        Basename prefix of all frame image files to be saved.
    saveFolder : str
        Basename of the subdirectory in the phase-specific results directory
        to which all animation results will be saved.
    '''
    pass

    # ..................{ CONCRETE ~ private                 }..................
    def __init__(
        self,

        # Mandatory parameters.
        sim: 'Simulator',
        cells: 'Cells',
        p: 'Parameters',
        clrAutoscale: bool,
        clrMin: float,
        clrMax: float,
        saveFolder: str,
        saveFile: str,

        # Optional parameters.
        colormap: 'Colormap' = None,
    ) -> None:
        '''
        Initialize this animation.

        Parameters
        ----------
        sim : Simulator
            Current simulation.
        cells : Cells
            Current cell cluster.
        p : Parameters
            Current simulation configuration.
        clrAutoscale : bool
            ???.
        clrMin : float
            ???.
        clrMax : float
            ???.
        saveFolder : str
            Basename of the subdirectory in the phase-specific results directory
            to which all animation results will be saved.
        saveFile : str
            Basename prefix of all frame image files to be saved.
        colormap : Colormap
            Matplotlib colormap to be used in this animation's colorbar.
        '''
        # Validate core parameters.
        assert types.is_simulator(sim), types.assert_not_simulator(sim)
        assert types.is_cells(cells), types.assert_not_parameters(cells)
        assert types.is_parameters(p), types.assert_not_parameters(p)

        # Default unpassed parameters.
        if colormap is None:
            colormap = p.default_cm

        # Validate all remaining parameters *AFTER* defaulting parameters.
        assert types.is_bool(clrAutoscale), types.assert_not_bool(clrAutoscale)
        assert types.is_numeric(clrMin), types.assert_not_numeric(clrMin)
        assert types.is_numeric(clrMax), types.assert_not_numeric(clrMax)
        assert types.is_matplotlib_colormap(colormap), (
            types.assert_not_matplotlib_colormap(colormap))
        assert types.is_str_nonempty(saveFolder), (
            types.assert_not_str_nonempty(saveFolder, 'Save directory'))
        assert types.is_str_nonempty(saveFile), (
            types.assert_not_str_nonempty(saveFile, 'Save frame file prefix'))

        # Classify *AFTER* validating parameters.
        self.sim = sim
        self.cells = cells
        self.p = p
        self.clrAutoscale = clrAutoscale
        self.clrMin = clrMin
        self.clrMax = clrMax
        self.colormap = colormap
        self.saveFolder = saveFolder
        self.saveFile = saveFile

        # Define all remaining attributes.
        self.ani_repeat = True

        # Initialize animation saving *AFTER* defining all attribute defaults.
        self._init_saving()


    def _init_saving(self) -> None:
        '''
        Initialize this animation for platform-compatible file saving if enabled
        by the current simulation configuration or noop otherwise.
        '''

        # If animation saving is disabled, noop.
        if self.p.saveAnimations is False:
            return

        # Path of the phase-specific parent directory of the subdirectory to
        # which these files will be saved.
        phase_dirname = None
        if self.p.plot_type == 'sim':
            phase_dirname = self.p.sim_results
        elif self.p.plot_type == 'init':
            phase_dirname = self.p.init_results
        else:
            raise BetseExceptionParameters(
                'Animation saving unsupported during the "{}" phase.'.format(
                    self.p.plot_type))

        #FIXME: Refactor all calls to os.makedirs() everywhere similarly.
        # Path of the subdirectory to which these files will be saved, creating
        # this subdirectory and all parents thereof if needed.
        images_dirname = paths.join(phase_dirname, self.saveFolder)
        images_dirname = dirs.canonicalize_and_make_unless_dir(images_dirname)

        # Path of the file to be saved.
        self.savedAni = paths.join(images_dirname, self.saveFile)

        # Force animations to *NOT* repeat.
        self.ani_repeat = False
