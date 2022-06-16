#!/usr/bin/env python3
# --------------------( LICENSE                           )--------------------
# Copyright 2014-2022 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Classes to easily run BETSE simulations from an external script.
'''

# ....................{ IMPORTS                            }....................
import numpy as np
from betse.science.simrunner import SimRunner
from betse.science import filehandling as fh
from betse.util.path import files
from betse.science.phase.phasecls import SimPhase
from betse.science.enum.enumphase import SimPhaseKind
from betse.science.parameters import Parameters as p
from betse.util.io.log import logs
from betse.util.io.log.conf import logconf
from betse.util.io.log.logenum import LogLevel
from scipy import interpolate
from betse.lib.pil import pilnumpy
from betse.lib.pil.pilnumpy import ImageModeType
from betse.science.math import finitediff as fd

# ....................{ Main                            }....................

class BetseWrapper(object):
    """
    Object allowing for simple creation of BETSE cell cluster and simulation object
    that can be easily worked with in an external script using BETSE as a dependency.

    This class creates (or optionally loads) a BETSE cells object, runs (or optionally loads)
    an init phase, and optionally runs a sim phase.

    """

    def __init__(self, *args, **kwargs):

        self.model_init(*args, **kwargs)  # Initializes the model when called

    def model_init(self, config_filename, new_mesh=True, verbose=True,
                   run_init=True, run_sim=False):
        """
        Initializes the BETSE modelling object, which includes creating or loading
        a cell cluster, running an init phase simulation and running a sim phase simulation.

        Parameters:
        --------------
        config_filename : str
            Full path to the configuration file for the simulation.
        new_mesh : bool
            Whether to generate a whole new mesh (True) or try to load a saved one (False).
        verbose: bool
            Spit out comments (True) or stay silent (False).
        run_init : bool
            Whether to run through the BETSE initialization (True) or not (False).
        run_sim : bool
            whether to run a BETSE simulation phase (True) or not (False).

        """

        # Make an instance of the BETSE 'parameters' object based on
        # settings in the configuration file supplied:

        self.p = p.make(config_filename)

        self.verbose = verbose  # save verbosity setting

        log_config = logconf.get_log_conf()

        if verbose:
            log_config.handler_stdout.setLevel(LogLevel.INFO)

        else:
            log_config = logconf.get_log_conf()

            # Reduce logging verbosity to improve readability.
            log_config.handler_stdout.setLevel(LogLevel.WARNING)

        self._make_mesh(new_mesh=new_mesh)  # make or load a BETSE cell cluster

        if run_init:
            self._init_runner(runsim=run_sim)

        if self.verbose is True:
            logs.log_info("Successfully generated betse model!")

    def _make_mesh(self, new_mesh=False):
        """
        Generates or loads a saved 2D Voronoi mesh for the BETSE simulation, based
        on settings in the supplied config file.

        Parameters
        ----------
        new_mesh : bool
            Whether to generate a whole new mesh (True) or try to load a saved one (False).


        """

        if new_mesh is True:  # If 'new mesh' is requested, make a whole new cell cluster

            if self.verbose is True:
                logs.log_info("Creating a new 2D Grid.")

            # Create a new grid:

            self.simrun = SimRunner(self.p)  # Call an instance of BETSE's "SimRunner'
            phase = self.simrun.seed()  # Go through the making of a whole cell cluster
            self.phase = phase  # save the phase object (which contains the cell cluster
            # to the DemoSim object)


        else:  # otherwise, if a new mesh isn't needed, go ahead and load a saved cluster,
            # if it exists!

            if not files.is_file(self.p.seed_pickle_filename):  # If it doesn't exist...
                if self.verbose is True:
                    logs.log_warning("File not found; Creating a new 2D Grid")

                # Make a new mesh

                self.simrun = SimRunner(self.p)
                phase = self.simrun.seed()  # Go ahead and make a new cluster
                self.phase = phase

            else:  # Otherwise, load the saved cell cluster:
                if self.verbose is True:
                    logs.log_info("Loading a 2D Grid from file.")

                # Load from previous creation:
                cells, _ = fh.loadWorld(self.p.seed_pickle_filename)

                self.simrun = SimRunner(self.p)

                # Simulation phase, created after unpickling these objects above

                self.phase = SimPhase(
                    kind=SimPhaseKind.SEED,
                    p=self.p,
                    cells=cells,
                )

                # Initialize core simulation data structures.

                self.phase.sim.init_core(self.phase)
                self.phase.dyna.init_profiles(self.phase)

        # assign commonly used dimensional parameters to the DemoSim object:
        self._assign_shorts(self.phase.cells)

    def _init_runner(self, runsim=False):
        """
        Run the BETSE cell cluster through the initialization specified in
        the loaded configuration file.

        Parameters
        ----------
        runsim : bool

        """

        self.phase = self.simrun.init()

        if runsim:
            self.phase = self.simrun.sim()

        # short forms of commonly used data structures from BETSE sim object, based on
        # ion settings for 'basic' profile:

        # Data defined on cell centrepoints (self.xc, self.yc):
        self.vm_ave = self.phase.sim.vm_ave  # Vmem state averaged over a whole cell [V]
        self.cc_cells = self.phase.sim.cc_cells  # Array of cytosolic ion conc arrays [mol/m^3]
        self.Ex = self.phase.sim.E_cell_x  # x-component of electric field [V/m]
        self.Ey = self.phase.sim.E_cell_y  # y-component of electric field [V/m]
        self.rho_cells = self.phase.sim.rho_cells  # charge density [C/m^3]

        # Data defined on membrane midpoints (self.xm, self.ym):

        self.vm = self.phase.sim.vm  # Vmem defined at cell edges (membranes) [V]
        self.Dm_cells = self.phase.sim.Dm_cells  # Array of mem diffusion constant arrays [m^2/s]

        # Data defined on environmental points (self.xenv, self.yenv)

        self.cc_env = self.phase.sim.cc_env  # Array of extracellular ion conc arrays [mol/m^3]

        # Single-value data (one value for each ion included in the simulation)

        self.D_free = self.phase.sim.D_free  # Array of media diffusion constant arrays [m^2/s]

        # Numerical indices to access specific ions in cc_cells, Dm_cells, cc_env, etc

        self.iNa = self.phase.sim.iNa  # Index of Na ion in array-of-arrays
        self.iK = self.phase.sim.iK  # Index of K ion in array-of-arrays
        self.iM = self.phase.sim.iM  # Index of M ion in array-of-arrays
        self.iP = self.phase.sim.iP  # Index of P ion in array-of-arrays

        # re-assign commonly used dimensional parameters:
        self._assign_shorts(self.phase.cells)

    def _assign_shorts(self, cells):
        """
        Assigns short forms to commonly used mesh components.

        Parameters
        -----------
        cells: a Betse cell cluster object

        """

        # points of cell centres:

        self.xc = cells.cell_centres[:, 0]
        self.yc = cells.cell_centres[:, 1]

        # membrane midpoints for edges of each cell:

        self.xmem = cells.mem_mids_flat[:, 0]
        self.ymem = cells.mem_mids_flat[:, 1]

        # points of the square environmental grid:

        self.xenv = cells.xypts[:, 0]
        self.yenv = cells.xypts[:, 1]

        # global size limits of the whole environment and cell cluster:

        self.xyaxis = np.asarray([cells.xmin, cells.xmax, cells.ymin, cells.ymax])

        # vertices of each polygon representing each cell in the cluster:
        self.verts = cells.cell_verts

        # Normal vectors to membrane edges
        self.nx = cells.mem_vects_flat[:, 2]
        self.ny = cells.mem_vects_flat[:, 3]

        # true extracellular matrix points (point shared between two membrane mids):
        self.xec = cells.ecm_mids[:, 0]
        self.yec = cells.ecm_mids[:, 1]

        # assign length of cell (cdl), mems (mdl) and env (edl) of mesh:

        self.cdl = len(cells.cell_i)
        self.mdl = len(cells.mem_i)
        self.edl = len(cells.ecm_mids)
        self.envdl = len(cells.xypts)

    def interp_bitmap_to_cells(self, bitmap_filename, to_mems=True, smooth=False):
        '''
        Interpolate a greyscale bitmap (supplied in the bitmap_filename string) to
        the cell cluster. The resulting interpolation is an array of floats, where
        0.0 is black and 255.0 is white.

        Parameters
        -------------
        bitmap_filename : str
            Path and filename to the bitmap to read and interpolate as a greyscale image.
        to_mems : bool
            If True, interpolates the bitmap to membrane domains of the cell cluster. If
            False, interpolates the bitmap to cell centres of the cell cluster.
        smooth : bool
            Apply smoothing to the interpolated bitmap (True) or keep raw interpolation (False).

        Returns
        -------
        ndarray
            An array (shaped to membrane domains or for to_mems == True and cell centres for
            to_mems=False) of bitmap greyscale values interpolated over the cell cluster. Values of
            0.0 are pure black while 255.0 are pure white.

        '''

        cells = self.phase.cells

        xmi = cells.xmin
        ymi = cells.ymin
        xma = cells.xmax
        yma = cells.ymax

        if to_mems:
            xmap = cells.map_mem2ecm

        else:
            xmap = cells.map_cell2ecm

        xx = np.linspace(cells.xmin, cells.xmax, cells.X.shape[1])
        yy = np.linspace(cells.ymin, cells.ymax, cells.X.shape[0])

        # Three-dimensional Numpy array of the RGB-ordered integer components of all
        # pixels loaded from this image.
        a1o = pilnumpy.load_image(filename=bitmap_filename, mode=ImageModeType.GRAYSCALE_FLOAT)

        a1 = np.asarray(a1o, dtype=np.float64)

        a1 = np.flipud(a1)

        xa = np.linspace(xmi, xma, a1.shape[1])
        ya = np.linspace(ymi, yma, a1.shape[0])

        spline_F = interpolate.interp2d(xa, ya, a1, kind='linear', fill_value=0.0)
        fe = spline_F(xx, yy)

        if smooth:
            fe = fd.integrator(fe, sharp=0.5)  # smooth a little to avoid bizarre visual effects

        f = fe.ravel()[xmap]

        return f