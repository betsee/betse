#!/usr/bin/env python3
# Copyright 2014-2017 by Alexis Pietak & Cecil Curry.
# See "LICENSE" for further details.

'''
**In-simulation animation** (i.e., animation produced *while* rather than
*after* solving a simulation) subclasses.
'''

# ....................{ IMPORTS                            }....................
import matplotlib
import numpy as np
from betse.science.visual import visualutil
from betse.science.visual.anim.animabc import AnimCellsABC
# from betse.util.io.log import logs
from betse.util.type.types import type_check, SequenceTypes
from matplotlib import pyplot

# ....................{ SUBCLASSES                         }....................
#FIXME: Rename "_cell_data_plot" to "_cell_body_plot".
#FIXME: Rename "_cell_edges_plot" to "_cell_edge_plot".

class AnimCellsWhileSolving(AnimCellsABC):
    '''
    In-simulation animation.

    This class animates arbitrary cell data as a time series plotted over the
    cell cluster (e.g., cell membrane voltage as a function of time) *while*
    rather than *after* solving a simulation.

    Attributes
    ----------
    _cell_edges_plot : LineCollection
        Plot of the current or prior frame's cell edges.
    _cell_data_plot : Collection
        Plot of the current or prior frame's cell contents.
    _cell_verts_id : int
        Unique identifier for the array of cell vertices (i.e.,
        `cells.cell_verts`) when plotting the current or prior frame.
        Retaining this identifier permits the `_plot_frame_figure()` method to
        efficiently detect and respond to physical changes (e.g., deformation
        forces, cutting events) in the fundamental structure of the previously
        plotted cell cluster.
    _is_colorbar_autoscaling_telescoped : bool
        `True` if colorbar autoscaling is permitted to increase but _not_
        decrease the colorbar range _or_ `False` otherwise (i.e., if
        colorbar autoscaling is permitted to both increase and decrease the
        colorbar range). Such telescoping assists in emphasizing stable
        long-term patterns in cell data at a cost of deemphasizing unstable
        short-term patterns. If colorbar autoscaling is disabled (i.e.,
        `is_color_autoscaled` is `False`), this will be ignored.
    _is_time_step_first : bool
        `True` only if the current frame being animated is the first.
    '''


    @type_check
    def __init__(
        self,

        # Mandatory parameters.
        p: 'betse.science.parameters.Parameters',

        # Optional parameters.

        #FIXME: Permit this option to be configured via a new configuration
        #option in the YAML file.
        is_colorbar_autoscaling_telescoped: bool = False,
        *args, **kwargs
    ) -> None:
        '''
        Initialize this in-simulation animation.

        Parameters
        ----------
        p : Parameters
            Current simulation configuration.
        is_colorbar_autoscaling_telescoped : optional[bool]
            `True` if colorbar autoscaling is permitted to increase but _not_
            decrease the colorbar range _or_ `False` otherwise (i.e., if
            colorbar autoscaling is permitted to both increase and decrease the
            colorbar range). Such telescoping assists in emphasizing stable
            long-term patterns in cell data at a cost of deemphasizing unstable
            short-term patterns. If colorbar autoscaling is disabled (i.e.,
            `is_color_autoscaled` is `False`), this will be ignored. Defaults
            to `False`.

        See the superclass `__init__()` method for all remaining parameters.
        '''

        # Initialize the superclass.
        super().__init__(
            p=p,

            # Prevent the superclass from overlaying electric current or
            # concentration flux. Although this class does *NOT* animate a
            # streamplot, the time series required to plot this overlay is
            # unavailable until after the simulation ends.
            is_current_overlayable=False,

            # Save and show this mid-simulation animation only if this
            # configuration has enabled doing so.
            is_save=p.anim.is_while_sim_save,
            is_show=p.anim.is_while_sim_show,

            # Save this mid-simulation animation to a different parent
            # directory than that to which the corresponding post-simulation
            # animation is saved.
            save_dir_parent_basename='anim_while_solving',

            # Pass all remaining arguments as is to our superclass.
            *args, **kwargs
        )
        # logs.log_debug('Showing mid-simulation animation: {}'.format(self._is_show))

        # Classify all remaining parameters.
        self._is_colorbar_autoscaling_telescoped = (
            is_colorbar_autoscaling_telescoped)

        # "True" only if the current frame being animated is the first.
        self._is_time_step_first = True

        # Unique identifier for the array of cell vertices. (See docstring.)
        self._cell_verts_id = id(self._cells.cell_verts)

        #FIXME: This is a temp change until we get this right.
        #FIXME: Refactor to call the new
        #Cells.map_membranes_midpoint_to_cells_centre() method instead.

        # average the voltage to the cell centre
        vm_o = np.dot(
            self._cells.M_sum_mems, self._sim.vm) / self._cells.num_mems

        # self._cell_time_series = self.sim.vm_time
        self._cell_time_series = self._sim.vm_ave_time

        # cell_data_current = self.sim.vm
        cell_data_current = vm_o

        # Upscaled cell data for the first frame.
        cell_data = visualutil.upscale_cell_data(cell_data_current)

        # Collection of cell polygons with animated voltage data.
        #
        # If *NOT* simulating extracellular spaces, only animate intracellular
        # spaces.
        self._cell_data_plot = self._plot_cells_sans_ecm(
            cell_data=cell_data)

        # Perform all superclass plotting preparation immediately *BEFORE*
        # plotting this animation's first frame.
        self._prep_figure(
            color_mappables=self._cell_data_plot,
            color_data=cell_data,
        )

        # Id displaying this animation, do so in a non-blocking manner.
        if self._is_show:
            pyplot.show(block=False)

            #FIXME: Unfortunately, the "block=False" combination leveraged above
            #is backend-specific and seemingly unsupported by at least the
            #"Qt4Agg" backend. Under this backend, an alternative approach is
            #required. Unfortunately, all attempts to implement such an approach
            #have thus far failed. Clearly, something is not quite right.
            #FIXME: At the very least, we should:
            #
            #* Add a new tester to the
            #  "betse.lib.matplotlib.matplotlibs.MatplotlibConfig" class
            #  returning True only if the current backend is non-blockable. This
            #  is at least the case for the "TkAgg" and "Agg" backends. Any
            #  others?
            #* Call this tester here and, if False, log a non-fatal warning
            #  informing the user that the mid-simulation animation will
            #  probably not be displayed due to unresolved matplotlib issues.
            #FIXME: O.K.; the implementation below actually works now. It's
            #non-ideal, however, and results in matplotlib deprecation warnings
            #resembling:
            #
            #    [betse] /usr/lib64/python3.4/site-packages/matplotlib/backend_bases.py:2437: MatplotlibDeprecationWarning: Using default event loop until function specific to this GUI is implemented
            #
            #Is that bad? Not particularly, but it *IS* annoying. Can we avoid
            #this? Yes, but only for backends known to support the experimental
            #"block=False" parameter (e.g., "TkAgg"). Specifically:
            #
            #* If the current backend supports this parameter, preferably call:
            #
            #     pyplot.show(block=False)
            #
            #* Else, fallback to:
            #
            #     matplotlib.interactive(True)
            #     pyplot.show()
            #FIXME: Unfortunately, the fact that we need to temporarily enable
            #"matplotlib.interactive(True)" *ONLY* for this duration of this
            #animation implies that we need to generalize this animation with
            #context manager functionality: specifically, we need to implement
            #the following special methods:
            #
            #* __enter__(), performing the logic detailed in the prior FIXME.
            #* __exit__(), calling "matplotlib.interactive(False)" *ONLY* if the
            #  prior __enter__() call called "matplotlib.interactive(True)".
            #
            #Properly implementing these methods is somewhat subtle, of course.
            #FIXME: After implementing non-blocking correctly here, generalize
            #this logic to the superclass by adding a new "is_blocking=True"
            #parameter to the AnimCellsABC.__init__() method. Ideally, *ALL*
            #animations should default to non-blocking behaviour; for safety,
            #however, blocking should remain the default until well-tested.
            #
            #This implies that all animations will become context managers.
            #While that's absolutely the correct long-term design, we currently
            #use post-simulation animations everywhere as if they are *NOT*
            #context managers. Since refactoring such usage is highly
            #non-trivial, we'll want to wait on this refactoring for a bit.

            # logs.log_debug('Enabling non-blocking animation behaviour...')
            # matplotlib.interactive(True)
            # pyplot.show()

    # ..................{ PLOTTERS                           }..................
    def _plot_frame_figure(self) -> None:

        # Upscaled cell data for the current time step.
        cell_data = visualutil.upscale_cell_data(
            self._cell_time_series[self._time_step])

        #FIXME: Duplicated from above. What we probably want to do is define a
        #new _get_cell_data() method returning this array in a centralized
        #manner callable both here and above. Rays of deluded beaming sunspray!

        # If the unique identifier for the array of cell vertices has *NOT*
        # changed, the cell cluster has *NOT* fundamentally changed and need
        # only be updated with this time step's cell data.
        if self._cell_verts_id == id(self._cells.cell_verts):
            # loggers.log_info(
            #     'Updating animation "{}" cell plots...'.format(self._type))
            self._update_cell_plots(cell_data)
        # Else, the cell cluster has fundamentally changed (e.g., due to
        # physical deformations or cutting events) and must be recreated.
        else:
            # loggers.log_info(
            #     'Reviving animation "{}" cell plots...'.format(self._type))

            # Prevent subsequent calls to this method from erroneously
            # recreating the cell cluster again.
            self._cell_verts_id = id(self._cells.cell_verts)

            # Recreate the cell cluster.
            self._revive_cell_plots(cell_data)

        # Update the color bar with the content of the cell body plot *AFTER*
        # possibly recreating this plot above.
        if self._is_color_autoscaled:
            cell_data_vm = cell_data

            # If autoscaling this colorbar in a telescoping manner and this is
            # *NOT* the first time step, do so.
            #
            # If this is the first time step, the previously calculated minimum
            # and maximum colors are garbage and thus *MUST* be ignored.
            if (self._is_colorbar_autoscaling_telescoped and
                not self._is_time_step_first):
                self._color_min = min(self._color_min, np.min(cell_data_vm))
                self._color_max = max(self._color_max, np.max(cell_data_vm))
            # Else, autoscale this colorbar in an unrestricted manner.
            else:
                self._color_min = np.min(cell_data_vm)
                self._color_max = np.max(cell_data_vm)

                # If autoscaling this colorbar in a telescoping manner, do so
                # for subsequent time steps.
                self._is_time_step_first = False

            # Autoscale the colorbar to these colors.
            self._rescale_color_mappables()

        # If displaying this frame...
        if self._is_show:
            # If the current event loop is idle, draw this frame; else, noop.
            # This is the OO-style equivalent to calling pyplot.draw().
            self._figure.canvas.draw_idle()
            # self._figure.canvas.draw()

            #FIXME: Duplicate call in the superclass plot_frame() method.
            #Ideally, this method's implementation should be entirely integrated
            #into the superclass plot_frame() method and then excised.

            # Temporarily yield the time slice for the minimum amount of time
            # required by the current matplotlib backend and operating system
            # for responding to queued events in the GUI eventloop of the
            # current process. Failing to do so reliably results in unresponsive
            # animations under poorly implemented backends (e.g., "Qt4Agg") and
            # the POSIX-incompatible Windows process model. For further details,
            # see also:
            #
            #     https://gitlab.com/betse/betse/issues/9
            #     https://github.com/matplotlib/matplotlib/issues/2134/
            pyplot.pause(0.0001)


    @type_check
    def _update_cell_plots(self, cell_data: SequenceTypes) -> None:
        '''
        Update _without_ recreating all cell plots for this time step with the
        passed array of arbitrary cell data.

        This method is intended to be called _unless_ physical changes
        (e.g., deformation forces, cutting events) in the underlying structure
        of the cell cluster have occurred for this simulation time step.

        Parameters
        -----------
        cell_data : SequenceTypes
            Arbitrary cell data defined on an environmental grid to be plotted.
        '''

        self._update_cell_plot_sans_ecm(
            cell_plot=self._cell_data_plot,
            cell_data=cell_data)


    @type_check
    def _revive_cell_plots(self, cell_data: SequenceTypes) -> None:
        '''
        Recreate all cell plots for this time step with the passed array of
        arbitrary cell data.

        This method is intended to be called in response to physical changes
        (e.g., deformation forces, cutting events) in the underlying structure
        of the cell cluster for this simulation time step. This method is both
        inefficient and destructive, and should be called only when needed.

        Parameters
        -----------
        cell_data : SequenceTypes
            Arbitrary cell data defined on an environmental grid to be plotted.
        '''

        self._cell_data_plot = self._revive_cell_plots_sans_ecm(
            cell_plot=self._cell_data_plot,
            cell_data=cell_data)
