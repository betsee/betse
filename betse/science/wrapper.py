#!/usr/bin/env python3
# --------------------( LICENSE                            )--------------------
# Copyright 2014-2022 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
Classes to easily run BETSE simulations from an external script.
'''

# ....................{ IMPORTS                            }....................
import os
import numpy as np
from beartype import beartype
from beartype.typing import Optional
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
from betse.lib.numpy.npcsv import write_csv
from betse.lib.pil import pilnumpy
from betse.lib.pil.pilnumpy import ImageModeType
from betse.science.math import finitediff as fd
from betse.lib import libs
from betse.science.chemistry.netplot import plot_master_network

# ....................{ CLASSES                            }....................
class BetseWrapper(object):
    '''
    Object allowing for simple creation of BETSE cell cluster and simulation
    object that can be easily worked with in an external script using BETSE as a
    dependency.

    This class creates (or optionally loads) a BETSE cells object, runs (or
    optionally loads) an init phase, and optionally runs a sim phase.
    '''

    # ..................{ INITIALIZERS                       }..................
    @beartype
    def __init__(
        self,

        # Mandatory parameters.
        config_filename: str,

        # Optional parameters.
        log_filename: Optional[str] = None,

        #FIXME: Validate with something resembling:
        #    log_level: Optional[Union[Literal['ALL'], ..., Literal['None']]] = None,
        log_level: Optional[str] = None,
    ) -> None:
        '''
        Initialization routines for the BETSE Wrapper.

        Parameters
        ----------
        config_filename : str
            Path to the config filename to run the betse simulation.
        log_filename : Optional[str]
            Absolute or relative filename to log simulation messages into.
            Defaults to ``None``, in which case the default filename is
            used.
        log_level : Optional[str]
            Valid options are:

            * 'ALL'.
            * 'DEBUG'.
            * 'INFO'.
            * 'WARNING'.
            * 'ERROR'.
            * 'CRITICAL'.
            * 'NONE'.

            String applying custom settings to the logging level of the
            BetseWrapper. If the log_level is specified by one of the above
            strings, then it overrides the verbose key in BetseWrapper methods
            with the specified logging function.
        '''

        # Classify all passed parameters.
        self._config_filename = config_filename
        self._log_filename = log_filename

        if log_level is not None:
            self._log_level = getattr(logs.LogLevel, log_level)
        else:
            self._log_level = None

    # ..................{ RUNNERS                            }..................
    @beartype
    def run_pipeline(
        self,
        new_mesh: bool = True,
        verbose: bool = True,
        run_init: bool = True,
        run_sim: bool = False,
    ) -> None:
        '''
        Runs an entire BETSE modelling pipeline -- including creating or loading
        a cell cluster, running an init phase simulation, and running a sim
        phase simulation.

        Parameters:
        --------------
        config_filename : str
            Full path to the configuration file for the simulation.
        new_mesh : bool
            Whether to generate a whole new mesh (True) or try to load a saved
            one (False).
        verbose: bool
            Spit out comments (True) or stay silent (False).
        run_init : bool
            Whether to run through the BETSE initialization (True) or not
            (False).
        run_sim : bool
            whether to run a BETSE simulation phase (True) or not (False).
        '''

        # Make an instance of the BETSE 'parameters' object based on settings in
        # the configuration file supplied.
        self.p = p.make(self._config_filename)

        self._set_logging(verbose=verbose)

        # Make or load a BETSE cell cluster.
        self._make_mesh(new_mesh=new_mesh)

        if run_init:
            self._init_runner(runsim=run_sim)

        if self.verbose is True:
            logs.log_info("Successfully run betse pipeline!")


    @beartype
    def run_seed(self, verbose: bool = False):
        '''
        Initializes the BETSE modelling object, which includes creating
        a cell cluster.

        Parameters:
        --------------
        '''

        self.p = p.make(self._config_filename)

        self.verbose = verbose  # save verbosity setting

        self._set_logging(verbose=verbose)

        self._make_mesh(new_mesh=True)  # make a BETSE cell cluster

        if self.verbose is True:
            logs.log_info(
                "Successfully created a new BETSE cell cluster object!")


    @beartype
    def run_init(self, new_mesh: bool = True, verbose: bool = False):
        '''
        Initializes the BETSE modelling object, which includes creating or loading
        a cell cluster, and running an init phase simulation.

        Parameters:
        --------------
        new_mesh : bool
            Whether to generate a whole new mesh (True) or try to load a saved one (False).
        verbose: bool
            Spit out comments (True) or stay silent (False).
        '''

        # Make an instance of the BETSE 'parameters' object based on
        # settings in the configuration file supplied:

        self.p = p.make(self._config_filename)

        self.verbose = verbose  # save verbosity setting

        self._set_logging(verbose=verbose)

        self._make_mesh(new_mesh=new_mesh)  # make or load a BETSE cell cluster

        self._init_runner(runsim=False)

        if self.verbose is True:
            logs.log_info("Successfully run initialization on BETSE model!")


    @beartype
    def run_sim(self, verbose: bool = False):
        '''
        Loads a previously-made BETSE cell cluster and init phase simulation to run
        the sim phase simulation.

        Parameters:
        --------------
        verbose: bool
            Spit out comments (True) or stay silent (False).
        '''

        # Make an instance of the BETSE 'parameters' object based on
        # settings in the configuration file supplied:

        self.p = p.make(self._config_filename)

        self.verbose = verbose  # save verbosity setting

        self._set_logging(verbose=verbose)

        self._make_mesh(new_mesh=False)  # load a BETSE cell cluster

        self._sim_runner() # Run the BETSE simulation

        if self.verbose is True:
            logs.log_info("Successfully run simulation on BETSE model!")

    # ..................{ LOADERS                            }..................
    #FIXME: Docstring us up, please. Flying churros at midnight!
    @beartype
    def load_seed(self, verbose: bool = False):
        '''
        '''

        self.p = p.make(self._config_filename)

        self.verbose = verbose  # save verbosity setting

        self._set_logging(verbose=verbose)

        if not files.is_file(self.p.seed_pickle_filename):  # If file doesn't exist...
            if self.verbose is True:
                logs.log_warning("File not found. Run seed to create a cell cluster.")

        else:  # Otherwise, load the saved cell cluster:
            if self.verbose is True:
                logs.log_info("Loading 2D Grid from file.")

            # Load from previous creation:
            cells, _ = fh.loadWorld(self.p.seed_pickle_filename)

            self.simrun = SimRunner(self.p)

            # Simulation phase, created after unpickling these objects above
            self.phase = SimPhase(
                kind=SimPhaseKind.SEED,
                p=self.p,
                cells=cells,
            )


    @beartype
    def load_init(self, verbose: bool = False):
        '''
        '''

        self.p = p.make(self._config_filename)

        self.verbose = verbose  # save verbosity setting

        self._set_logging(verbose=verbose)

        if not files.is_file(self.p.init_pickle_filename):  # If file doesn't exist...
            if self.verbose is True:
                logs.log_warning("File not found. Run init to create an initialization.")

        else:  # Otherwise, load the saved cell cluster:
            if self.verbose is True:
                logs.log_info("Loading BETSE init from file.")

            sim, cells, _ = fh.loadSim(self.p.init_pickle_filename)

            self.simrun = SimRunner(self.p)

            # Simulation phase, created after unpickling these objects above
            self.phase = SimPhase(
                kind=SimPhaseKind.INIT,
                p=self.p,
                cells=cells,
                sim=sim,
            )

            if self.verbose is True:
                logs.log_info("Successfully loaded init of BETSE model!")


    @beartype
    def load_sim(self, verbose: bool = False):
        '''
        '''

        self.p = p.make(self._config_filename)

        self.verbose = verbose  # save verbosity setting

        self._set_logging(verbose=verbose)

        if not files.is_file(self.p.sim_pickle_filename):  # If file doesn't exist...
            if self.verbose is True:
                logs.log_warning("File not found. Run sim to create a simulation.")

        else:  # Otherwise, load the saved cell cluster:
            if self.verbose is True:
                logs.log_info("Loading BETSE sim from file.")

            sim, cells, _ = fh.loadSim(self.p.sim_pickle_filename)

            self.simrun = SimRunner(self.p)

            # Simulation phase, created after unpickling these objects above
            self.phase = SimPhase(
                kind=SimPhaseKind.SIM,
                p=self.p,
                cells=cells,
                sim=sim,
            )

            if self.verbose is True:
                logs.log_info("Successfully loaded sim of BETSE model!")


    @beartype
    def load_simgrn(self, verbose: bool = False):
        '''
        '''

        self.p = p.make(self._config_filename)

        self.verbose = verbose  # save verbosity setting

        self._set_logging(verbose=verbose)

        if not files.is_file(self.p.sim_pickle_filename):  # If file doesn't exist...
            if self.verbose is True:
                logs.log_warning("File not found. Run sim_grn to create a new grn simulation.")

        else:  # Otherwise, load the saved cell cluster:
            if self.verbose is True:
                logs.log_info("Loading simulated grn from file.")

            init, _, _ = fh.loadSim(self.p.init_pickle_filename)
            grn, cells, _ = fh.loadSim(self.p.grn_pickle_filename)

            # Simulation phase.
            self.phase = SimPhase(
                kind=SimPhaseKind.INIT,
                cells=cells,
                sim=init,
                p=self.p)

            self.phase.sim.grn.core = grn

            if self.verbose is True:
                logs.log_info(
                    "Successfully loaded simulated grn of BETSE model!")


    @beartype
    def run_sim_grn(self, new_mesh: bool = True, verbose: bool = False):
        '''
        Run only the BETSE GRN of the model (no bioelectricity).
        '''

        self.p = p.make(self._config_filename)

        self.verbose = verbose  # save verbosity setting

        self._set_logging(verbose=verbose)

        self._make_mesh(new_mesh=new_mesh)  # load a BETSE cell cluster

        self.phase = self.simrun.sim_grn()  # Run the BETSE simulation as GRN-only

        if self.verbose is True:
            logs.log_info(
                "Successfully run GRN-only simulation on BETSE model!")

    # ..................{ GETTERS                            }..................
    @beartype
    def get_network(self, verbose: bool = True, as_networkx: bool = False):
        '''
        Returns BETSE's GRN (if simulated) as a pydot (as_networkx=False) or
        networkx (as_networkx = True) Digraph object.
        '''

        if self.phase.sim.grn is not None:
            # Working with BETSE's networks:
            # Access the gene regulatory network core:
            grn = self.phase.sim.grn.core

            graph_net = plot_master_network(grn, self.p)

            if as_networkx is True:
                networkx = libs.import_runtime_optional('networkx')
                # Convert the pydot graph to a networkx file:
                graph_net = networkx.nx_pydot.from_pydot(graph_net)

        else:
            if verbose is True:
                logs.log_info("No GRN present in this BETSE model.")
            graph_net = None

        return graph_net


    @beartype
    def get_betse_grn(self, verbose: bool = True):
        '''
        Returns BETSE's GRN modelling object (if a GRN is simulated in the BETSE
        model).
        '''

        if self.phase.sim.grn is not None:
            # Working with BETSE's networks:
            # Access the gene regulatory network core:
            grn = self.phase.sim.grn.core

        else:
            if verbose is True:
                logs.log_info("No GRN present in this BETSE model.")
            grn = None

        return grn


    @beartype
    def get_connected_grn_elements(self, verbose: bool = False):
        '''
        Returns a list of sets of connected nodes of BETSE's GRN (if simulated).
        '''

        graph_network = self.get_network(verbose=verbose, as_networkx=True)

        if graph_network is not None:
            graph_network_o = graph_network.to_undirected()
            networkx = libs.import_runtime_optional('networkx')
            connected_elements = sorted(networkx.connected_components(graph_network_o))

            if verbose is True:
                print(connected_elements)

        else:
            if verbose is True:
                logs.log_info("No GRN present in this BETSE model.")
            connected_elements = None

        return connected_elements

    # ..................{ PLOTTERS                           }..................
    @beartype
    def plot_network(self, verbose: bool = True):
        '''
        Exports an svg of BETSE's GRN (if simulated).
        '''

        if self.phase.sim.grn is not None:
            # Working with BETSE's networks:
            # Access the gene regulatory network core:
            grn = self.phase.sim.grn.core

            graph_pydot = plot_master_network(grn, self.p)

            # Save the pydot graph to an svg file:
            # Initialize saving:
            grn.init_saving(self.phase.cells, self.p, plot_type='init', nested_folder_name='GRN')

            # Optionally print the location of the image path using: print(grn.imagePath)
            savename = os.path.join(grn.imagePath[0:-4], 'OptimizedNetworkGraph.svg')
            graph_pydot.write_svg(savename, prog='dot')

            if verbose is True:
                logs.log_info(f"Model GRN network image saved to {savename}.")

        else:
            if verbose is True:
                logs.log_info("No GRN present in this BETSE model.")

    # ..................{ OTHERS                             }..................
    @beartype
    def analyze_network(
        self,
        verbose: bool = True,
        plot_network: bool = False,
        save_csv: bool = True,
    ):
        '''
        Analyzes the node connectivity (node degree) of a BETSE GRN (if
        simulated) and optionally exports results to csv (save_csv = True). The
        graph of the GRN can also be optionally exported (plot_network=True).
        '''

        if self.phase.sim.grn is not None:
            # Import optional dependencies for working with networks:
            # Defer importation of optional runtime dependencies until necessary.
            pydot, networkx = libs.import_runtime_optional('pydot', 'networkx')

            # Working with BETSE's networks:
            # Access the gene regulatory network core:
            grn = self.phase.sim.grn.core

            graph_pydot = plot_master_network(grn, self.p)

            if plot_network is True:
                # Save the pydot graph to an svg file:
                # Initialize saving:
                grn.init_saving(self.phase.cells, self.p, plot_type='init', nested_folder_name='GRN')

                # Optionally print the location of the image path using: print(grn.imagePath)
                savename = os.path.join(grn.imagePath[0:-4], 'OptimizedNetworkGraph.svg')
                graph_pydot.write_svg(savename, prog='dot')

                if verbose is True:
                    logs.log_info(f"Model GRN network image saved to {savename}.")

            # Convert the pydot graph to a networkx file:
            graph_network = networkx.nx_pydot.from_pydot(graph_pydot)

            # Compute the connectivity degree of graph nodes:
            degree_list = np.asarray(sorted(networkx.degree(graph_network)))

            # determine the connectivity degree values as an array of floats:
            degree_vals = np.asarray(degree_list[:, 1], dtype=float)

            # get the indices into the degree list that will sort the connectivity from high to low:
            inds_sort = np.flip(np.argsort(degree_vals))

            # get the node names and degree values sorted from highest to lowest connectivity degree:
            sorted_nodes = degree_list[inds_sort, 0]
            sorted_vals = degree_list[inds_sort, 1]

            if verbose is True:
                for node, val in zip(sorted_nodes, sorted_vals):
                    logs.log_info(f'{repr(node)}, {repr(val)}')


            if save_csv is True:
                savename = os.path.join(grn.imagePath[0:-4], 'NetworkConnectivity.csv')
                column_name_to_values = {
                    'Node name': sorted_nodes,
                    'Connectivity': sorted_vals,
                }
                column_name_to_format = {
                    'Node name': '%s',
                    'Connectivity': '%s',
                }

                write_csv(
                    filename=savename,
                    column_name_to_values=column_name_to_values,
                    column_name_to_format=column_name_to_format,
                )

                if verbose is True:
                    logs.log_info(f"Analyzed BETSE model GRN network and exported csv to {savename}.")

        else:
            if verbose is True:
                logs.log_info("No GRN present in this BETSE model.")

        return sorted_nodes, sorted_vals


    @beartype
    def interp_bitmap_to_cells(
        self,
        bitmap_filename: str,
        to_mems: bool = True,
        smooth: bool = False,
    ):
        '''
        Interpolate a greyscale bitmap (supplied in the bitmap_filename string)
        to the cell cluster. The resulting interpolation is an array of floats,
        where
        0.0 is black and 255.0 is white.

        Parameters
        -------------
        bitmap_filename : str
            Path and filename to the bitmap to read and interpolate as a
            greyscale image.
        to_mems : bool
            If True, interpolates the bitmap to membrane domains of the cell
            cluster. If False, interpolates the bitmap to cell centres of the
            cell cluster.
        smooth : bool
            Apply smoothing to the interpolated bitmap (True) or keep raw
            interpolation (False).

        Returns
        -------
        ndarray
            An array (shaped to membrane domains or for to_mems == True and cell
            centres for to_mems=False) of bitmap greyscale values interpolated
            over the cell cluster. Values of
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
    # ..................{ PRIVATE                            }..................
    @beartype
    def _set_logging(self, verbose: bool = False) -> None:
        '''
        Set the logging properties of the BetseWrapper.
        '''

        self.verbose = verbose  # save verbosity setting

        # Logging configuration singleton.
        log_config = logconf.get_log_conf()

        # If the user passed a log filename, reconfigure our logging
        # configuration to log to this file.
        if self._log_filename is not None:
            log_config.filename = self._log_filename
        # Else, the user passed *NO* log filename. In this case, accept the
        # current default log filename.

        # If the user passed *NO* log level, default to a log level
        # corresponding to the passed verbosity.
        if self._log_level is None:
            if verbose:
                log_config.handler_stdout.setLevel(LogLevel.INFO)
            else:
                log_config.handler_stdout.setLevel(LogLevel.WARNING)
        # Else, the user has set a log level. In this case, apply it.
        else:
            log_config.handler_stdout.setLevel(self._log_level)


    @beartype
    def _make_mesh(self, new_mesh: bool = False):
        '''
        Generates or loads a saved 2D Voronoi mesh for the BETSE simulation,
        based on settings in the supplied config file.

        Parameters
        ----------
        new_mesh : bool
            Whether to generate a whole new mesh (True) or try to load a saved
            one (False).
        '''

        if new_mesh is True:  # If 'new mesh' is requested, make a whole new cell cluster

            if self.verbose is True:
                logs.log_info("Creating a new cell cluster.")

            # Create a new grid:

            self.simrun = SimRunner(self.p)  # Call an instance of BETSE's "SimRunner'
            phase = self.simrun.seed()  # Go through the making of a whole cell cluster
            self.phase = phase  # save the phase object (which contains the cell cluster
            # to the DemoSim object)

        else:  # otherwise, if a new mesh isn't needed, go ahead and load a saved cluster,
            # if it exists!

            if not files.is_file(self.p.seed_pickle_filename):  # If it doesn't exist...
                if self.verbose is True:
                    logs.log_warning("File not found; Creating a new cell cluster...")

                # Make a new mesh
                self.simrun = SimRunner(self.p)
                phase = self.simrun.seed()  # Go ahead and make a new cluster
                self.phase = phase

            else:  # Otherwise, load the saved cell cluster:
                if self.verbose is True:
                    logs.log_info("Loading a cell cluster from file.")

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


    @beartype
    def _init_runner(self, runsim: bool = False):
        '''
        Run the BETSE cell cluster through the initialization specified in
        the loaded configuration file.

        Parameters
        ----------
        runsim : bool
        '''

        self.phase = self.simrun.init()

        if runsim:
            self.phase = self.simrun.sim()

        # short forms of commonly used data structures from BETSE sim object, based on
        # ion settings for 'basic' profile:

        # re-assign commonly used dimensional parameters:
        self._assign_shorts(self.phase.cells)


    def _sim_runner(self):
        '''
        Run the BETSE cell cluster through the initialization specified in
        the loaded configuration file.

        Parameters
        ----------
        runsim : bool
        '''

        self.phase = self.simrun.sim()

        # short forms of commonly used data structures from BETSE sim object, based on
        # ion settings for 'basic' profile:

        # re-assign commonly used dimensional parameters:
        self._assign_shorts(self.phase.cells)


    def _assign_shorts(self, cells):
        '''
        Assigns short forms to commonly used mesh components.

        Parameters
        -----------
        cells : ???
            BETSE cell cluster object.
        '''

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
